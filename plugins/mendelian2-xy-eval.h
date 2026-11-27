/* One evaluator feeds both annotations and native counters. */
static int collect_xy_stats(args_t *args, bcf1_t *rec)
{
    xy_policy_t *p=&args->xy;
    int chr=xy_chr(bcf_seqname(args->hdr,rec)), boundary=0;
    /* Non-sex chromosomes keep the upstream evaluator and native counters.
     * XY metadata is observational only on this route. */
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
    bcf_unpack(rec,BCF_UN_STR);
    int64_t end=rec->pos+rec->rlen;
    int pure_insertion=rec->n_allele>1;
    for (int a=1;a<rec->n_allele;a++) if (strlen(rec->d.allele[a])<=strlen(rec->d.allele[0]) || strncmp(rec->d.allele[a],rec->d.allele[0],strlen(rec->d.allele[0]))) pure_insertion=0;
    if (pure_insertion) end=rec->pos+1;
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
    if (!reason) {
        for (int a=0;a<rec->n_allele;a++) {
            const char *al=rec->d.allele[a];
            if (!*al) reason=13;
            for (;*al;al++) if (!strchr("ACGTNacgtn",*al)) reason=13;
        }
    }
    int nvalues=bcf_get_genotypes(args->hdr,rec,&args->gt_arr,&args->ngt_arr);
    int nsamples=bcf_hdr_nsamples(args->hdr), stride=nvalues>0?nvalues/nsamples:0;
    args->nmerr=args->nmiss=args->ngood=args->nrule=0;
    int ret=0, final_reason=reason;
    for (int t=0;t<args->ntrio;t++) {
        trio_t *trio=&args->trio[t]; int r=reason;
        int lengths[3]={0}, allele[3][2]={{0}}, missing=0, hetero=0, ploidy=0, nonref=0;
        trio->has_merr=0;
        if (trio->fail) { trio->stats.nfail++; continue; }
        if (!r) {
            for (int m=0;m<3;m++) if (expected[m]) {
                if (!stride) { missing=1; continue; }
                int32_t *gt=&args->gt_arr[stride*trio->idx[m]];
                for (int j=0;j<stride && gt[j]!=bcf_int32_vector_end;j++) {
                    lengths[m]++;
                    if (bcf_gt_is_missing(gt[j])) missing=1;
                    else {
                        int a=bcf_gt_allele(gt[j]);
                        if (a<0 || a>=rec->n_allele) missing=1;
                        if (a>0) nonref=1;
                        if (j<2) allele[m][j]=a;
                    }
                }
                if (!lengths[m]) missing=1;
                if (expected[m]==1 && lengths[m]==2 && allele[m][0]!=allele[m][1]) hetero=1;
                if (lengths[m]!=expected[m] && !(expected[m]==1 && lengths[m]==2 && allele[m][0]==allele[m][1])) ploidy=1;
            }
            if (missing) r=9;
            else if (hetero) r=11;
            else if (ploidy) r=10;
            else if (!nonref) r=12;
            else {
                /* Reuse upstream allele_bit/parse_gt and transmission tests.
                 * Eligibility above guarantees complete GT and supported ploidy;
                 * pseudo2 collapse is confined to the child/source bit view. */
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
            }
        }
        p->counts[region][r]++;
        final_reason=r;
        if (!r) { args->ngood++; trio->stats.ngood++; trio->stats.ngood_alt++; ret|=HAS_GOOD; }
        else if (r==1) { args->nmerr++; trio->stats.nmerr++; trio->has_merr=1; ret|=HAS_MERR; }
        else { args->nmiss++; trio->stats.nmiss++; ret|=HAS_MISS; }
    }
    const char *status=final_reason==0?"consistent":final_reason==1?"violation":(final_reason==8 || final_reason==13 || final_reason==14)?"unsupported":"not_evaluable";
    bcf_update_info_string(args->hdr_out,rec,"XY_STATUS",status);
    bcf_update_info_string(args->hdr_out,rec,"XY_REASON",xy_reasons[final_reason]);
    bcf_update_info_string(args->hdr_out,rec,"XY_REGION",xy_regions[region]);
    return ret;
}
