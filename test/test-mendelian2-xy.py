#!/usr/bin/env python3
"""CLI tests of the compiled native mendelian2 evaluator; no output rewriting."""
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / '.local' / 'validation' / 'mendelian2-xy'


class NativeXY(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        ignored = subprocess.run(['git', '-C', str(ROOT), 'check-ignore', '.local/'], capture_output=True)
        tracked = subprocess.run(['git', '-C', str(ROOT), 'ls-files', '.local'], capture_output=True, check=True)
        if ignored.returncode or tracked.stdout:
            raise RuntimeError('.local must be ignored and contain no tracked files')
        WORK.mkdir(parents=True, exist_ok=True)
        cls.env = dict(os.environ, BCFTOOLS_PLUGINS=str(ROOT / 'plugins'))
        cls.beds = WORK / 'beds'
        cls.beds.mkdir(exist_ok=True)
        for chrom in ('X', 'Y'):
            for kind, text in [('PAR', '0\t10\n20\t30\n'), ('nonPAR', '10\t15\n16\t20\n30\t1000\n')]:
                with gzip.open(cls.beds / f'assembly.chr{chrom}_{kind}.bed.gz', 'wt') as out:
                    out.write(''.join(f'chr{chrom}\t{line}\n' for line in text.strip().splitlines()))
        (cls.beds / 'par_regions.tsv').write_text('chrom\tstart\tend\tregion\nX\t0\t10\tPAR1\nX\t20\t30\tPAR2\nY\t0\t10\tPAR1\nY\t20\t30\tPAR2\n')

    def call(self, rows: list[str], sex: str = 'male', masked: bool = True,
             mode: str = 'a', extra: list[str] | None = None,
             resources: bool = True, pfm: bool = True,
             header: str = '') -> subprocess.CompletedProcess[str]:
        path = WORK / 'input.vcf'
        path.write_text('##fileformat=VCFv4.2\n' + ''.join(f'##contig=<ID={c},length=1000>\n' for c in ('X', 'Y', 'chrX', 'chrY', '1')) +
                        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' + header +
                        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tChild\tFather\tMother\n' + '\n'.join(rows) + '\n')
        cmd = [str(ROOT / 'bcftools'), '+mendelian2', str(path), '--child-sex', sex, '-m', mode]
        if pfm:
            cmd += ['-p', 'Child,Father,Mother']
        if resources:
            cmd += ['--par-dir', str(self.beds)]
        if masked:
            cmd += ['--ypar-masked']
        return subprocess.run(cmd + (extra or []), env=self.env, text=True, capture_output=True)

    @staticmethod
    def row(chrom: str = 'X', pos: int = 12, child: str = '1', dad: str = '.',
            mom: str = '0/1', ref: str = 'A', alt: str = 'C') -> str:
        return f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{child}\t{dad}\t{mom}'

    def check(self, row: str, reason: str = 'pass', status: str = 'consistent',
              region: str | None = None, **kwargs: object) -> dict[str, str]:
        result = self.call([row], **kwargs)
        self.assertEqual(result.returncode, 0, result.stderr)
        out = [line.split('\t') for line in result.stdout.splitlines() if not line.startswith('#')]
        self.assertEqual(len(out), 1, result.stdout)
        info = dict(token.split('=', 1) for token in out[0][7].split(';'))
        self.assertEqual(info['XY_REASON'], reason)
        self.assertEqual(info['XY_STATUS'], status)
        if region:
            self.assertEqual(info['XY_REGION'], region)
        self.assertEqual(out[0][9:], row.split('\t')[9:], 'Original GT must be preserved')
        return info

    def test_inheritance_and_required_members(self) -> None:
        for chrom in ('X', 'chrX'):
            self.check(self.row(chrom=chrom), region='X_nonPAR')
            self.check(self.row(chrom=chrom, dad='1/2'), region='X_nonPAR')
            self.check(self.row(chrom=chrom, child='0', mom='1/1'), 'mendelian_violation', 'violation')
            self.check(self.row(chrom=chrom, child='1/1'), region='X_nonPAR')
            self.check(self.row(chrom=chrom, child='1|1'), region='X_nonPAR')
        self.check(self.row(chrom='Y', dad='1', mom='.'), region='Y_nonPAR')
        self.check(self.row(chrom='Y', dad='1/1', mom='0/1'), region='Y_nonPAR')
        self.check(self.row(chrom='Y', child='0', dad='1'), 'mendelian_violation', 'violation')
        self.check(self.row(child='0/1', dad='1', mom='0/0'), sex='female')
        self.check(self.row(child='1/1', dad='1', mom='0/0'), 'mendelian_violation', 'violation', sex='female')
        self.check(self.row(child='0/1', dad='0', mom='0/0'), 'mendelian_violation', 'violation', sex='female')
        self.check(self.row(child='0/1', dad='1/1', mom='0/0'), sex='female')

    def test_exclusions(self) -> None:
        for gt in ('.', './.', '0/.', './1'):
            self.check(self.row(child=gt), 'required_gt_missing', 'not_evaluable')
        self.check(self.row(child='0/1'), 'haploid_heterozygous', 'not_evaluable')
        self.check(self.row(chrom='Y', dad='0/1'), 'haploid_heterozygous', 'not_evaluable')
        self.check(self.row(mom='1'), 'unexpected_ploidy', 'not_evaluable')
        self.check(self.row(child='1/1/1'), 'unexpected_ploidy', 'not_evaluable')
        self.check(self.row(child='0', mom='0/0', dad='1'), 'all_required_reference', 'not_evaluable')
        self.check(self.row(chrom='Y', child='0', dad='0', mom='1/1'), 'all_required_reference', 'not_evaluable')
        self.check(self.row(chrom='Y'), 'female_y_nonpar', 'unsupported', sex='female')
        self.check(self.row(alt='<DEL>'), 'unsupported_variant', 'unsupported')
        self.check(self.row(), 'unknown_child_sex', 'not_evaluable', sex='unknown')
        self.check(self.row(), 'missing_region_resources', 'not_evaluable', resources=False, masked=False)
        self.check(self.row(pos=16), 'unclassified_region', 'not_evaluable')
        self.check(self.row(pos=15, ref='AAA', alt='A'), 'boundary_crossing', 'not_evaluable')

    def test_par_matrix_and_boundaries(self) -> None:
        for sex in ('male', 'female', 'unknown'):
            for masked in (False, True):
                for chrom in ('X', 'Y'):
                    for pos, par in ((5, 1), (25, 2)):
                        reason = 'unknown_child_sex' if sex == 'unknown' else 'y_par_excluded' if chrom == 'Y' else 'par_requires_ypar_masked' if not masked else 'pass'
                        self.check(self.row(chrom=chrom, pos=pos, child='0/1', dad='0/1'), reason,
                                   'consistent' if reason == 'pass' else 'not_evaluable', f'{chrom}_PAR{par}', sex=sex, masked=masked)
        self.check(self.row(pos=10, child='0/1', dad='0/1', alt='AC'), region='X_PAR1')
        self.check(self.row(pos=10, ref='AA', alt='A'), 'boundary_crossing', 'not_evaluable')
        self.check(self.row(pos=10, ref='AA', alt='AAC', child='0/1', dad='0/1'), region='X_PAR1')
        manifest = self.beds / 'par_regions.tsv'
        text = manifest.read_text()
        manifest.unlink()
        try:
            self.check(self.row(pos=25, child='0/1', dad='0/1'), region='X_PAR')
        finally:
            manifest.write_text(text)

    def test_autosomes_and_multiallelic(self) -> None:
        self.check(self.row(chrom='1', child='0/1', dad='1/1', mom='0/0'), sex='unknown', region='autosome')
        self.check(self.row(chrom='1', child='0/0', dad='0/0', mom='0/0'), 'all_required_reference', 'not_evaluable')
        self.check(self.row(child='2', mom='1/2', alt='C,G'))
        self.check(self.row(child='2', mom='0/1', alt='C,G'), 'mendelian_violation', 'violation')

    def test_counts_match_annotations(self) -> None:
        rows = [self.row(), self.row(child='0', mom='1/1'), self.row(child='0/1'), self.row(child='0', mom='0/0')]
        result = self.call(rows, mode='c')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('XY_REGION\tX_nonPAR\t2\t1\t0.5', result.stdout)
        self.assertIn('XY_REGION\tY_nonPAR\t0\t0\tNA', result.stdout)
        for reason in ('pass', 'mendelian_violation', 'haploid_heterozygous', 'all_required_reference'):
            self.assertIn(f'XY_REASON\tX_nonPAR\t{reason}\t1', result.stdout)
        for mode, expected in [('e', 1), ('g', 1), ('m', 2), ('E', 3), ('M', 2)]:
            result = self.call(rows, mode=mode)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(len([x for x in result.stdout.splitlines() if not x.startswith('#')]), expected)

    def test_existing_info_headers(self) -> None:
        for tag in ('XY_STATUS', 'XY_REASON', 'XY_REGION'):
            for number, datatype in (('0', 'Flag'), ('2', 'String'), ('.', 'String'), ('1', 'Integer')):
                header = f'##INFO=<ID={tag},Number={number},Type={datatype},Description="Existing">\n'
                result = self.call([self.row()], header=header)
                self.assertNotEqual(result.returncode, 0, header)
                self.assertIn(f'Incompatible INFO/{tag}', result.stderr)
        header = ''.join(f'##INFO=<ID={tag},Number=1,Type=String,Description="Existing">\n'
                         for tag in ('XY_STATUS', 'XY_REASON', 'XY_REGION'))
        row = self.row().replace('\t.\tGT\t', '\tXY_STATUS=violation;XY_REASON=stale;XY_REGION=autosome\tGT\t')
        self.check(row, header=header, region='X_nonPAR')

    def test_cli_and_bad_resources(self) -> None:
        for extra in (['--ypar-masked=false'], ['--child-sex', 'other'], ['--rules', 'GRCh38'], ['--rules-file', 'rules'], ['-m', 'd']):
            self.assertNotEqual(self.call([self.row()], extra=extra).returncode, 0, str(extra))
        self.assertNotEqual(self.call([self.row()], resources=False).returncode, 0)
        ped = WORK / 'input.ped'
        ped.write_text('F Child Father Mother 2\n')
        self.assertNotEqual(self.call([self.row()], pfm=False, extra=['-P', str(ped)]).returncode, 0)
        bed = self.beds / 'assembly.chrX_nonPAR.bed.gz'
        original = bed.read_bytes()
        for content in ('X\t-1\t10\n', 'X\t20\t10\n', 'Y\t0\t20\n', 'X\t8\t12\n', 'X\t10\t20\nX\t19\t30\n', ''):
            with gzip.open(bed, 'wt') as out:
                out.write(content)
            self.assertNotEqual(self.call([self.row()]).returncode, 0, content)
        bed.write_bytes(original)
        duplicate = self.beds / 'duplicate.chrX_nonPAR.bed'
        duplicate.write_text('X\t10\t20\n')
        try:
            self.assertNotEqual(self.call([self.row()]).returncode, 0)
        finally:
            duplicate.unlink()
        bed.unlink()
        try:
            self.assertNotEqual(self.call([self.row()]).returncode, 0)
        finally:
            bed.write_bytes(original)


if __name__ == '__main__':
    unittest.main(verbosity=2)
