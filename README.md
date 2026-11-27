[![Build Status](https://api.cirrus-ci.com/github/samtools/bcftools.svg?branch=develop)](https://api.cirrus-ci.com/github/samtools/bcftools)
[![Build status](https://ci.appveyor.com/api/projects/status/yanx5wnsdsqm4pay?svg=true)](https://ci.appveyor.com/project/samtools/bcftools)
[![Github All Releases](https://img.shields.io/github/downloads/samtools/bcftools/total.svg)](https://github.com/samtools/bcftools/releases/latest)

This is the official development repository for BCFtools. It contains all the vcf* commands which previously lived in the htslib repository (such as vcfcheck, vcfmerge, vcfisec, etc.) and the samtools BCF calling from bcftools subdirectory of samtools.

For a full documentation, see [bcftools GitHub page](http://samtools.github.io/bcftools/).

Other useful links:
------------------

File format specifications live on [HTS-spec GitHub page](http://samtools.github.io/hts-specs/) [htslib](https://github.com/samtools/htslib) [samtools](https://github.com/samtools/samtools) [tabix](https://github.com/samtools/tabix)

### Citing

Please cite this paper when using BCFtools for your publications. http://samtools.github.io/bcftools/howtos/publications.html

> Twelve years of SAMtools and BCFtools </br>
> Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li </br>
> _GigaScience_, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008

```
@article{10.1093/gigascience/giab008,
    author = {Danecek, Petr and Bonfield, James K and Liddle, Jennifer and Marshall, John and Ohan, Valeriu and Pollard, Martin O and Whitwham, Andrew and Keane, Thomas and McCarthy, Shane A and Davies, Robert M and Li, Heng},
    title = "{Twelve years of SAMtools and BCFtools}",
    journal = {GigaScience},
    volume = {10},
    number = {2},
    year = {2021},
    month = {02},
    abstract = "{SAMtools and BCFtools are widely used programs for processing and analysing high-throughput sequencing data. They include tools for file format conversion and manipulation, sorting, querying, statistics, variant calling, and effect analysis amongst other methods.The first version appeared online 12 years ago and has been maintained and further developed ever since, with many new features and improvements added over the years. The SAMtools and BCFtools packages represent a unique collection of tools that have been used in numerous other software projects and countless genomic pipelines.Both SAMtools and BCFtools are freely available on GitHub under the permissive MIT licence, free for both non-commercial and commercial use. Both packages have been installed \\&gt;1 million times via Bioconda. The source code and documentation are available from https://www.htslib.org.}",
    issn = {2047-217X},
    doi = {10.1093/gigascience/giab008},
    url = {https://doi.org/10.1093/gigascience/giab008},
    note = {giab008},
    eprint = {https://academic.oup.com/gigascience/article-pdf/10/2/giab008/36332246/giab008.pdf},
}
```

### Support

If you have found a bug or would like a new feature, please report the same in the GitHub [BCFtools](https://github.com/samtools/bcftools/issues) issue tracker.

For any security related issue, please send a mail to [samtools@sanger.ac.uk](mailto:samtools@sanger.ac.uk) instead of reporting in the GitHub issue tracker.

## Native XY Mendelian evaluation (`XY` branch)

### 1. XY inheritance principles

This small-variant policy assumes an XX mother and an XY father. In X-PAR, the child inherits one allele from each parent and all three required GTs are diploid. Evaluation in X-PAR requires a declaration that the calling pipeline masked Y-PAR. Y-PAR records are excluded to avoid counting the same PAR inheritance twice. `--ypar-masked` records that declaration; it does not verify the reference, alignment or original genotype calls. Without it, X-PAR is excluded while nonPAR can still be evaluated.

| Region | Child sex | Required members and ploidy | Transmission |
| --- | --- | --- | --- |
| X-PAR, Y-PAR masking declared | male or female | mother 2, father 2, child 2 | One allele from each parent |
| X-nonPAR | male | mother 2, child 1 | Maternal only; father is irrelevant |
| X-nonPAR | female | mother 2, father 1, child 2 | Maternal and paternal contributions occupy distinct child allele slots |
| Y-nonPAR | male | father 1, child 1 | Paternal only; mother is irrelevant |
| Y-nonPAR | female | none | Unsupported, not a Mendelian violation |
| Y-PAR | any | none | Excluded |

The evaluator checks only the child and the parents listed for that region. Their GTs must be complete, have the expected number of alleles and include at least one nonreference allele across those people. Missing required GTs are excluded, never converted to reference. A homozygous pseudo-diploid call such as `1/1` may be interpreted internally as haploid `1`; the original GT remains in the output. A heterozygous GT in a required haploid member is excluded, not counted as an inheritance error. Unknown child sex excludes XY evaluation without changing the original autosome evaluator.

BED intervals for the chosen assembly determine which inheritance rule applies. The complete reference span must fit intervals that use the same rule. For an insertion whose alternate allele starts with the full reference allele, only the anchor position is used. Records crossing a model boundary or a covered interval into a gap are excluded. The evaluator does not assign a rule to gaps, and GRCh38 coordinates must not be reused for another assembly.

### 2. What the upstream baseline does

The fixed baseline is [`e9ad5ec2`](https://github.com/YuanfengZhang/bcftools/tree/e9ad5ec2fa91fecc060aec645c8eeb4b96c8e579), whose `version.sh` declares source version 1.25. This fork has no 1.24 tag. The original implementation is [`plugins/mendelian2.c`](https://github.com/YuanfengZhang/bcftools/blob/e9ad5ec2fa91fecc060aec645c8eeb4b96c8e579/plugins/mendelian2.c). Its GRCh38 rules include these exact entries:

```c
            "   1X  X:1-9999                M\n"
            "   1X  X:2781480-155701381     M\n"
            "   1X  Y:1-57227415            F\n"
            "   2X  Y:1-57227415            .\n"
```

Thus the predefined male Y rule covers the full listed Y interval, including Y-PAR. Unlisted regions default to diploid inheritance from both parents. In `collect_stats()`, rules are selected by overlap with the record, without requiring its entire span to fit the selected model:

```c
    int i,itr_set = regidx_overlap(args->rules, bcf_seqname(args->hdr,rec),rec->pos,rec->pos+rec->rlen-1, args->itr);
```

In `parse_ped()`, a PED sex other than 1 takes the `2X` branch:

```c
        if ( khash_str2int_get(args->str2sex_id, sex==1 ? "1X" : "2X", &trio->sex_id)<0 )
            error("Missing the sex \"%s\", it's not in the rules :-/\n",sex==1 ? "1X" : "2X");
```

The native algorithm encodes alleles as bits. A haploid child is compared with the selected parent's allele set using `parent&kid1`. For diploid transmission it tests both assignments of the two child alleles:

```c
        int nal_mom = parse_gt(&args->gt_arr[ngt*trio->idx[iMOM]],ngt,&mom,&mom,als,&nals);
        int nal_dad = parse_gt(&args->gt_arr[ngt*trio->idx[iDAD]],ngt,&dad,&dad,als,&nals);
        if ( (kid1&dad && kid2&mom) || (kid1&mom && kid2&dad) )
```

That transmission test remains useful. The XY branch adds checks for whether PAR can be evaluated, whether child sex is known, whether the whole record fits one inheritance rule, and whether each participating parent has the expected number of GT alleles. The baseline folds a homozygous diploid child in a haploid model, but counts a heterozygous child there as a Mendelian error. The baseline can also use partial information when a parent is missing; the XY policy requires complete GTs for every participating person. These checks decide which records reach the existing transmission test and which parent supplies an allele.

### 3. What the XY branch changes and how to run it

The implementation at [`3b9e0e5b`](https://github.com/YuanfengZhang/bcftools/tree/3b9e0e5b692e6f81cfe3d38dba52b79bc54a02ec) adds resource loading in `plugins/mendelian2-xy.h` and checks that select the participating parents and exclude unsuitable XY records in `plugins/mendelian2-xy-eval.h`, called from `plugins/mendelian2.c`. Eligible XY records reuse the existing `parse_gt()` and `allele_bit()` functions and the original bitmask comparisons for haploid or diploid transmission. Every record outside X and Y still calls the original `collect_stats()` even when the new options are present. Its native annotations and counters retain their upstream semantics, including partial-missing and reference-only behavior. No new option means the original upstream mode throughout. The plugin does no normalization or haplotype replay and does not change GTs.

The implementation files changed from the fixed upstream baseline have these responsibilities:

| File | Change and responsibility |
| --- | --- |
| `plugins/mendelian2.c` | Extends the existing plugin with reading CLI options, checking child sex against PED or PFM input, validating INFO definitions, choosing the XY checks when requested and printing counts by region; retains the original evaluator. |
| `plugins/mendelian2-xy.h` | Defines the XY settings, recognizes chromosome aliases, reads and validates the four BED files, checks coverage across adjacent intervals and reads optional PAR names. |
| `plugins/mendelian2-xy-eval.h` | Selects the participating parents and expected GT allele counts, excludes records that fail those checks, then uses the original bitmask comparison. Sends records outside X and Y to the original evaluator and records results for annotation and counts. |
| `plugins/mendelian2.mk` | New dependencies on both policy headers so incremental builds rebuild the plugin after policy changes. |
| `test/test-mendelian2-xy.py` | New CLI regression driver for the compiled plugin, including XY cases, invalid BED files and INFO definitions, preserved GTs and 729 autosome genotype combinations; performs no replacement evaluation. |
| `doc/bcftools.txt` | Extends the existing manual with the behavior of the new options, commands and limitations. |
| `.gitignore` | Excludes local build dependencies, fixtures and logs under `/.local/`. |

The evaluator first finds the region, then selects which parents participate and how many GT alleles each person must have. A zero in the expected-ploidy array means that parent is ignored. It checks region and sex exclusions before testing inheritance. Source: [`plugins/mendelian2-xy-eval.h:24-35`](https://github.com/YuanfengZhang/bcftools/blob/3b9e0e5b692e6f81cfe3d38dba52b79bc54a02ec/plugins/mendelian2-xy-eval.h#L24-L35) at `3b9e0e5b`.

```c
    int region=xy_region(p,chr,rec->pos,end,&boundary), reason=0;
    if (chr && !p->sex) reason=2;
    else if (chr && !p->dir) reason=3;
    else if (boundary) reason=5;
    else if (region==8) reason=4;
    else if (chr==2 && region!=7) reason=7;
    else if (chr==1 && region!=3 && !p->masked) reason=6;
    else if (region==7 && p->sex==2) reason=8;
    int expected[3]={2,2,2};
    if (region==3 && p->sex==1) { expected[iDAD]=0; expected[iKID]=1; }
    else if (region==3 && p->sex==2) expected[iDAD]=1;
    else if (region==7) { expected[iMOM]=0; expected[iDAD]=expected[iKID]=1; }
```

Once all participating people have complete GTs with the expected allele counts, and at least one has a nonreference allele, the evaluator calls the original parser that encodes alleles as bits. It then uses the original comparison for haploid or diploid transmission. For a diploid child, one allele must match the mother and the other must match the father; this also applies to female X-nonPAR. Source: [`plugins/mendelian2-xy-eval.h:78-91`](https://github.com/YuanfengZhang/bcftools/blob/3b9e0e5b692e6f81cfe3d38dba52b79bc54a02ec/plugins/mendelian2-xy-eval.h#L78-L91) at `3b9e0e5b`.

```c
                uint8_t kid1,kid2,mom=0,dad=0,parent=0;
                int als[7]={0}, nals=1;
                parse_gt(&args->gt_arr[stride*trio->idx[iKID]],stride,&kid1,&kid2,als,&nals);
                int consistent;
                if (expected[iKID]==1) {
                    int ipar=expected[iMOM]?iMOM:iDAD;
                    parse_gt(&args->gt_arr[stride*trio->idx[ipar]],stride,&parent,&parent,als,&nals);
                    consistent=parent&kid1;
                } else {
                    parse_gt(&args->gt_arr[stride*trio->idx[iMOM]],stride,&mom,&mom,als,&nals);
                    parse_gt(&args->gt_arr[stride*trio->idx[iDAD]],stride,&dad,&dad,als,&nals);
                    consistent=(kid1&dad && kid2&mom) || (kid1&mom && kid2&dad);
                }
                r=consistent?0:1;
```

Records outside X and Y immediately call the original evaluator. Its result flags and counters are unchanged. The added XY fields describe those results without imposing the new rules for which XY records can be evaluated. Source: [`plugins/mendelian2-xy-eval.h:8-18`](https://github.com/YuanfengZhang/bcftools/blob/3b9e0e5b692e6f81cfe3d38dba52b79bc54a02ec/plugins/mendelian2-xy-eval.h#L8-L18) at `3b9e0e5b`.

```c
    if (!chr) {
        int before_alt=args->trio[0].stats.ngood_alt;
        int ret=collect_stats(args,rec);
        int reason=(ret&HAS_MERR)?1:(ret&HAS_MISS)?9:(ret&HAS_GOOD)?0:10;
        if (reason==0 && args->trio[0].stats.ngood_alt==before_alt) reason=12;
        p->counts[9][reason]++;
        bcf_update_info_string(args->hdr_out,rec,"XY_STATUS",reason==0?"consistent":reason==1?"violation":"not_evaluable");
        bcf_update_info_string(args->hdr_out,rec,"XY_REASON",xy_reasons[reason]);
        bcf_update_info_string(args->hdr_out,rec,"XY_REGION","autosome");
        return ret;
    }
```

The count command computes N and E from the reason counters also updated during annotation, separately for each region. Only pass and mendelian_violation contribute to N; zero N prints NA. These are record counts, not haplotype replay events. Source: [`plugins/mendelian2.c:842-850`](https://github.com/YuanfengZhang/bcftools/blob/3b9e0e5b692e6f81cfe3d38dba52b79bc54a02ec/plugins/mendelian2.c#L842-L850) at `3b9e0e5b`.

```c
    if (args->xy.enabled) {
        fprintf(log_fh,"# XY policy v1; counts are required-member nonreference records, not haplotype replay events\n");
        fprintf(log_fh,"# XY_REGION\tregion\tN_evaluable\tE_violation\terror_rate\n");
        for (int r=0;r<10;r++) {
            uint64_t n=args->xy.counts[r][0]+args->xy.counts[r][1], e=args->xy.counts[r][1];
            fprintf(log_fh,"XY_REGION\t%s\t%"PRIu64"\t%"PRIu64"\t",xy_regions[r],n,e);
            if (n) fprintf(log_fh,"%.8g\n",(double)e/n); else fprintf(log_fh,"NA\n");
            for (int j=0;j<15;j++) if (args->xy.counts[r][j]) fprintf(log_fh,"XY_REASON\t%s\t%s\t%"PRIu64"\n",xy_regions[r],xy_reasons[j],args->xy.counts[r][j]);
        }
```

The explicit CLI is:

- `--child-sex male|female|unknown`; omission in explicit mode means unknown.
- `--par-dir DIR`; four unique, nonempty files ending in
  `chrX_PAR.bed[.gz]`, `chrX_nonPAR.bed[.gz]`, `chrY_PAR.bed[.gz]` and
  `chrY_nonPAR.bed[.gz]`. Both plain BED and gzip BED work without an index.
- `--ypar-masked`; optional pure flag, default off.
  `--ypar-masked=false` is rejected. The flag requires valid `--par-dir`.

BED coordinates are zero-based, half-open and nonnegative. X and chrX are aliases, as are Y and chrY. Intervals cannot overlap, including between PAR and nonPAR files. Optional `par_regions.tsv` has columns `chrom start end region`, with an optional tab-separated header of those exact names and region names PAR1 or PAR2. It assigns names only to PAR records already covered by the BED files. Without a containing manifest interval, the label is X_PAR or Y_PAR; fragmented BED rows never become PAR numbers by row order.

Use the standalone `+mendelian2` run interface shown below: its options go directly after the input, without the extra `--` separator used by some other bcftools plugins. `Child,Father,Mother` is the `-p` order. For a source build, set `BCFTOOLS_PLUGINS` to the built plugin directory first.

```sh
export BCFTOOLS_PLUGINS="$PWD/plugins"

# Annotate calls made with declared Y-PAR masking; allow X-PAR evaluation.
./bcftools +mendelian2 input.vcf.gz -p Child,Father,Mother \
    --child-sex male --par-dir beds --ypar-masked -m a -o annotated.vcf

# Count with the same native policy.
./bcftools +mendelian2 input.vcf.gz -p Child,Father,Mother \
    --child-sex male --par-dir beds --ypar-masked -m c -o counts.tsv

# Without the masking declaration, evaluate nonPAR and exclude PAR.
./bcftools +mendelian2 input.vcf.gz -p Child,Father,Mother \
    --child-sex female --par-dir beds -m a -o female.annotated.vcf
```

#### Compatibility and limits

Use one trio at a time; split multi-family PED input. Known PED sex or a `1X:`/`2X:` prefix must agree with `--child-sex`. The new mode cannot be combined with `--rules`, `--rules-file` or genotype-deleting `-m d`.

An invalid BED directory causes an error. Without a directory, XY records are skipped. Symbolic alleles, breakends and spanning-deletion `*` are unsupported. The plugin compares alleles within each record; it cannot match equivalent changes split across different records. It needs no FASTA and does not mask references, realign reads, call genotypes or add quality/depth thresholds.

#### Results and counts

- `XY_STATUS` says whether a record is `consistent`, `violation`, `not_evaluable` or `unsupported`. Only a checked inheritance failure is a violation.
- `XY_REASON` explains the result, such as `pass`, `mendelian_violation`, missing calls, invalid ploidy or an excluded region. The reserved reason `unsupported_replay` is not used by this plugin.
- `XY_REGION` identifies X/Y PAR or nonPAR, `unclassified` or `autosome`. PAR1/PAR2 names are used when supplied by the manifest.

These are single-value String INFO fields. Compatible existing fields are replaced; incompatible definitions cause an error.

In `-m c` output, `XY_REGION region N E error_rate` gives evaluated records (N), violations (E) and E/N. A record needs a nonreference allele in at least one participating person to count toward XY N. When N is zero, the rate is NA. `XY_REASON region reason count` lists the reasons by region.

Autosomal MERR, MGOOD, MMISS and summary counts retain their original meaning. Added XY fields describe that result without changing it. An all-reference autosomal record can therefore be good in the original counter while carrying `all_required_reference`. Always inspect `XY_STATUS` for XY eligibility; absence of MERR does not prove a record was checked.

The regression driver `python3 -I test/test-mendelian2-xy.py` tests the compiled plugin without rewriting its decisions.

## Usage

Build this fork's XY branch with a C compiler, make and the development libraries required by htslib (see `INSTALL`). The tested implementation is `3b9e0e5b`, with htslib `36ef0a19` and its htscodecs submodule `b9fc194f`. Use the pinned sources when reproducing that build; the XY branch may receive later updates. The commands below assume the bcftools source directory is the current directory and the pinned htslib source, including htscodecs, is at `../htslib`.

```sh
make HTSDIR=../htslib -j8
export BCFTOOLS_PLUGINS="$PWD/plugins"
./bcftools +mendelian2 -h
```

Prepare one VCF containing the child, father and mother, with FORMAT/GT calls from the same reference assembly and compatible allele representations. This plugin does not merge separate sample VCFs or normalize their alleles. Obtain the joint or merged VCF from your calling workflow and inspect its sample names with `./bcftools query -l trio.vcf.gz`. Supply those exact names in child,father,mother order with `-p`. Input may be plain VCF, gzip VCF or BCF; sequential reading does not require an index unless you also request indexed region access.

You can instead supply `-P trio.ped`, using family information from your sample metadata. The first five whitespace-separated columns are family ID (ignored by the plugin), child sample ID, father sample ID, mother sample ID and child sex (1 male, 2 female, 0 unknown). A sixth phenotype column is optional. For example, this line describes samples named Child, Father and Mother:

```text
FAM1 Child Father Mother 1 0
```

Use either `-p Child,Father,Mother` or `-P trio.ped`, not both. Explicit XY evaluation requires one unambiguous trio. Choose `--child-sex` from the child's known sex; it must agree with a known PED sex. Use `unknown` when sex is unknown rather than assuming female.

Obtain the four XY BED files (X PAR, X nonPAR, Y PAR and Y nonPAR) from the assembly-specific XY directory in the official [GIAB genome stratifications v3.6 release](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/). Use the directory matching the VCF and reference assembly, and keep assemblies in separate directories.

Use `bedtools maskfasta` with the assembly-matched Y PAR BED to generate a YPARmasked FASTA, then index it with `samtools faidx`. Keep the reference assembly and contig names consistent with the VCF and BED files.

Each directory must contain exactly one matching file for each suffix `chrX_PAR.bed`, `chrX_nonPAR.bed`, `chrY_PAR.bed` and `chrY_nonPAR.bed`, optionally followed by `.gz`. The original GIAB filenames match these suffixes. Plain and gzip files are both accepted, but keeping compressed and uncompressed copies of the same class in one directory makes the match ambiguous. The intervals must be nonempty, valid zero-based half-open BED intervals and must not overlap within a chromosome, even across the PAR and nonPAR files.

The optional `par_regions.tsv` in the same directory supplies PAR1 and PAR2 names. Omit it if you do not need separate names: eligible PAR records will then be labeled X_PAR or Y_PAR. If you supply it, use literal tabs and assembly-specific canonical intervals, not one PAR number per BED fragment. For CHM13, the following writes the corresponding manifest:

```sh
printf 'chrom\tstart\tend\tregion\nchrX\t0\t2394410\tPAR1\nchrX\t153925834\t154259566\tPAR2\nchrY\t0\t2458320\tPAR1\nchrY\t62122809\t62460029\tPAR2\n' \
    > beds/CHM13/par_regions.tsv
```

This manifest names covered PAR records but never adds eligible bases. Do not copy the CHM13 manifest into the GRCh38 directory. A GRCh38 manifest must use canonical GRCh38 intervals from its assembly resources; the plugin works without it.

For a male child whose VCF was called with Y-PAR masking, annotate and count with the same options:

```sh
./bcftools +mendelian2 trio.vcf.gz -p Child,Father,Mother \
    --child-sex male --par-dir beds/GRCh38 --ypar-masked \
    -m a -o trio.xy.vcf

./bcftools +mendelian2 trio.vcf.gz -p Child,Father,Mother \
    --child-sex male --par-dir beds/GRCh38 --ypar-masked \
    -m c -o trio.xy.counts.tsv
```

For CHM13 calls, replace `beds/GRCh38` with `beds/CHM13`. For a female child, pass `--child-sex female` and matching pedigree information. To use a PED file, replace only `-p Child,Father,Mother` with `-P trio.ped`. Do not add an extra `--` separator to this plugin's run interface or combine these options with `--rules` or `--rules-file`.

If the VCF was not called with Y-PAR masking, omit the flag. The following still evaluates X-nonPAR and Y-nonPAR for a male child but excludes X-PAR and Y-PAR:

```sh
./bcftools +mendelian2 trio.vcf.gz -p Child,Father,Mother \
    --child-sex male --par-dir beds/GRCh38 \
    -m a -o trio.nonpar.vcf
```

Do not pass `--ypar-masked=false`; the option takes no value. This plugin needs no reference FASTA for XY evaluation. The VCF must already have been called against the intended assembly. Supplying BED files or setting the flag does not mask a FASTA, realign reads, recalculate GTs or prove how the VCF was generated. Keep the calling workflow's reference and masking provenance separately.

Inspect `XY_STATUS`, `XY_REASON` and `XY_REGION` in the annotated VCF. `consistent` and `violation` identify evaluated inheritance results; `not_evaluable` and `unsupported` records are excluded from the XY N denominator. The reason distinguishes missing GTs, unexpected allele counts, unknown sex, region exclusions and unsupported alleles. The plugin preserves the original GTs. In the count file, `XY_REGION` rows give N evaluated records, E violations and E/N, with NA when N is zero; `XY_REASON` rows count each reason by region. These additions do not change upstream autosome annotations or summary counters. Symbolic alleles, breakends and spanning-deletion stars are unsupported on the explicit XY route; the plugin does not compare haplotypes across records.
