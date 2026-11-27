/* Explicit assembly-independent XY resources for mendelian2. */
#include <dirent.h>
#include <limits.h>
typedef struct { int chr, kind; int64_t beg, end; int par; } xy_interval_t;
typedef struct {
    int enabled, sex, sex_set, masked;
    char *dir;
    xy_interval_t *iv, *manifest;
    int niv, nmanifest;
    uint64_t counts[10][15];
} xy_policy_t;
static const char *xy_regions[] = {"X_PAR1","X_PAR2","X_PAR","X_nonPAR","Y_PAR1","Y_PAR2","Y_PAR","Y_nonPAR","unclassified","autosome"};
static const char *xy_reasons[] = {"pass","mendelian_violation","unknown_child_sex","missing_region_resources","unclassified_region","boundary_crossing","par_requires_ypar_masked","y_par_excluded","female_y_nonpar","required_gt_missing","unexpected_ploidy","haploid_heterozygous","all_required_reference","unsupported_variant","unsupported_replay"};
static int xy_chr(const char *s) { if (!strncmp(s,"chr",3)) s+=3; return !strcmp(s,"X") ? 1 : !strcmp(s,"Y") ? 2 : 0; }
static int xy_suffix(const char *s,const char *end) { size_t a=strlen(s),b=strlen(end); return a>=b && !strcmp(s+a-b,end); }
static void xy_read_bed(xy_policy_t *p,const char *path,int chr,int kind,int manifest)
{
    htsFile *fp=hts_open(path,"r");
    if (!fp) error("Cannot open XY resource: %s\n",path);
    kstring_t line={0,0,0}; int n=0;
    while (hts_getline(fp,KS_SEP_LINE,&line)>=0) {
        char seq[256],name[32]="",start[64],stop[64],tail; long long b,e;
        if (!line.l || line.s[0]=='#') continue;
        if (manifest && !strcmp(line.s,"chrom\tstart\tend\tregion")) continue;
        int nc=sscanf(line.s,"%255s %63s %63s %31s %c",seq,start,stop,name,&tail);
        char *bend=NULL,*eend=NULL;
        errno=0;
        b=nc>=3?strtoll(start,&bend,10):-1;
        int invalid=errno || !bend || bend==start || *bend;
        errno=0;
        e=nc>=3?strtoll(stop,&eend,10):-1;
        invalid |= errno || !eend || eend==stop || *eend;
        if (nc<3 || invalid || b<0 || e<=b || !xy_chr(seq) || (!manifest && xy_chr(seq)!=chr)) error("Invalid XY BED row in %s: %s\n",path,line.s);
        if (manifest && (nc!=4 || (strcmp(name,"PAR1") && strcmp(name,"PAR2") && strcmp(name,"X_PAR1") && strcmp(name,"X_PAR2") && strcmp(name,"Y_PAR1") && strcmp(name,"Y_PAR2")))) error("Invalid PAR manifest row in %s\n",path);
        xy_interval_t iv={.chr=xy_chr(seq),.kind=kind,.beg=b,.end=e,.par=manifest ? name[strlen(name)-1]-'0' : 0};
        xy_interval_t **arr=manifest?&p->manifest:&p->iv; int *num=manifest?&p->nmanifest:&p->niv;
        for (int i=0;i<*num;i++) if ((*arr)[i].chr==iv.chr && b<(*arr)[i].end && e>(*arr)[i].beg) error("Overlapping XY resource intervals in %s\n",path);
        *arr=realloc(*arr,(*num+1)*sizeof(**arr)); if (!*arr) error("Out of memory\n"); (*arr)[(*num)++]=iv; n++;
    }
    free(line.s); if (hts_close(fp)) error("Cannot close %s\n",path);
    if (!n) error("Empty XY resource: %s\n",path);
}
static void xy_load(xy_policy_t *p)
{
    if (!p->dir) return;
    const char *suffix[]={"chrX_PAR.bed","chrX_nonPAR.bed","chrY_PAR.bed","chrY_nonPAR.bed"};
    DIR *d=opendir(p->dir); if (!d) error("Cannot open --par-dir %s\n",p->dir);
    char *paths[4]={0}; struct dirent *ent;
    while ((ent=readdir(d))) for (int i=0;i<4;i++) {
        char gz[64]; snprintf(gz,sizeof(gz),"%s.gz",suffix[i]);
        if (!xy_suffix(ent->d_name,suffix[i]) && !xy_suffix(ent->d_name,gz)) continue;
        if (paths[i]) error("Ambiguous XY resource for %s\n",suffix[i]);
        size_t len=strlen(p->dir)+strlen(ent->d_name)+2; paths[i]=malloc(len); snprintf(paths[i],len,"%s/%s",p->dir,ent->d_name);
    }
    closedir(d);
    for (int i=0;i<4;i++) { if (!paths[i]) error("Missing XY resource: %s\n",suffix[i]); xy_read_bed(p,paths[i],i<2?1:2,i%2?1:0,0); free(paths[i]); }
    char path[PATH_MAX]; snprintf(path,sizeof(path),"%s/par_regions.tsv",p->dir);
    if (!access(path,F_OK)) xy_read_bed(p,path,0,0,1);
}
/* Coverage uses the union of adjacent fragments, never a complement. */
static int xy_cover(xy_interval_t *iv,int n,int chr,int kind,int64_t b,int64_t e)
{
    int64_t pos=b;
    while (pos<e) { int64_t next=pos; for (int i=0;i<n;i++) if (iv[i].chr==chr && iv[i].kind==kind && iv[i].beg<=pos && iv[i].end>next) next=iv[i].end; if (next==pos) return 0; pos=next; }
    return 1;
}
static int xy_region(xy_policy_t *p,int chr,int64_t b,int64_t e,int *boundary)
{
    *boundary=0; if (!chr) return 9;
    for (int kind=0;kind<2;kind++) if (xy_cover(p->iv,p->niv,chr,kind,b,e)) {
        if (kind) return chr==1?3:7;
        for (int i=0;i<p->nmanifest;i++) if (p->manifest[i].chr==chr && b>=p->manifest[i].beg && e<=p->manifest[i].end) return (chr==1?0:4)+p->manifest[i].par-1;
        return chr==1?2:6;
    }
    for (int i=0;i<p->niv;i++) if (p->iv[i].chr==chr && b<p->iv[i].end && e>p->iv[i].beg) *boundary=1;
    return 8;
}
