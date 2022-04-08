/*    ###################################################################
     ##                                                                   ##
     ##       #####################  ICeCAP ##########################    ##
     ##       ##                   Hi-Resolution-Hi-C @NML          ##    ##
     ##       ##                                                    ##    ##
     ##       ##          Algorithm by Gabriele Migliorini,         ##    ##
     ##       ########################################################    ##
     ##                                                                   ##
     ################################################################# */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h> 

#define PACKAGE "ICECap"
#define VERSION "0.01"
#define CHROMOSOMES 22
#define WIN_SHIFT1 1
#define RESOLUTION 1000

int rf[CHROMOSOMES],gp[CHROMOSOMES],bp[CHROMOSOMES];
int size,WIN_SIZE,CUT,GRID,pchic;
int chrom,chrom2,chrom_start,chrom_end,allchrom;
int *map_count,*chromgrid,*chrombait;
int **enzyme,**enriched,**gridpoints,**baitpoints;
int **hic;

double *weight,*weight2;
double *avg,*avg2,*dilute,*map;

char *lane,*folderN,*folder,*folderref,*folderrefmap,*folderfrag,*folderrefC,*folderwrite,*BAITS,*BAITFLAG,*PATH,*PATH2;
char *ALLELE,*PWEIGHT,*PATHBOWTIE,*blunt,*stats,*outprefix,*path;
char *outdir,*outdir2,*outdir3,*indir,*REFERENCE,*folder1,*TRANS;
char *ENZYME; 
char JUNCTION[20];
char stamp[200],stamp2[200];
char stringarray[15000000][100];
FILE *fmap;
FILE *fout;
FILE *fout2,*fout2b,*fout3,*fout4,*fout4b,*fout5,*fout6,*fchic;

char* concat(char *s1, char *s2)
{
    size_t len1 = strlen(s1);
    size_t len2 = strlen(s2);

    char *result = malloc(len1+len2+1);
        memcpy(result, s1, len1);
            memcpy(result+len1, s2, len2+1);
	    return result;          }
   
void append(char* s, char c)
{
        int len = strlen(s);
        s[len] = c;
        s[len+1] = '\0';
}

int dirExists(const char *path)
{
    struct stat info;

    if(stat( path, &info ) != 0)
        return 0;
    else if(info.st_mode & S_IFDIR)
        return 1;
    else
        return 0;
}

void print_usage(){
char command[500];
sprintf(command,"echo \"Please refer to the Manual: \"");
system(command);
sprintf(command,"less ./reference_C/MANUAL*.tex");
system(command);
exit(0);
}

void strrev(char *str)
{
        if( str == NULL )
                return;

        char *end_ptr = &str[strlen(str) -1 ];
        char temp;
        while( end_ptr > str )
        {
                temp = *str;
                *str++ = *end_ptr;
                *end_ptr-- = temp;
        }
}


double init();
double digest();

int grid();
int cgrid();
int bgrid();

int main(int argc, char *argv[])
{
      int option,i,j,k,temp,a,chrom,BIN,flag,flag1,flag2;
      int max_iter,corenodes,REwidth,resolution,iter;
      int tot,tot2,tot3,tot5,tot6,totfragments[25],totbaits[25],nchromosomes;
      int save_a,save_b,qc;
      double *genome;

      char *file,*file2,*file3,*file4,*filenorma,*filelogs,*fileref,*fileref2,*filerefb,*fileref3,*fileref4;
      char *fileoutput1,*fileoutput2,*fileoutput3,*fileoutput4;

      double PVAL,FDR,map_thr,weight_threshold,map_threshold;
      double tolerance,conv_check;

      char QNAME1[400],QNAME2[400],CIGAR1[200],CIGAR2[200],SEQ1[400],SEQ2[400];
      char command[500],heatmap[500],graphics[500],end[500];
      char zom[5],chromosome[5],chromosome2[5],iteration[5],sbait[5],wind[5];
      char stamp[200],stamp2[200],schr[10],suffix[200],STRAND1[12],STRAND2[12];
      char *subbuff;
      char lyne[121];
      char *subbuff2;
      char *item,*item1,*item2;
      int FRAG1,FRAG2,POS1,pos2b,POS2,QUAL1,QUAL2,d3,d6,d7,CHR1,CHR2,TRIM;
      int position,position2,position3,enrichment,region,SCALE;
      int area1,area2,window;

      long int tot4; 
      double d1,d2,d4,d5,d1p,d2p;
      double mappability,fdist,gdist,reg_size,reg_size2,mid,mid2,save;
      char enzyme_seq,baitname[100];

      typedef struct {

      int num_chromosome; char name_chromosome[20]; } chromap;
      chromap record[30]; 


int
number_for_key(char *key)
{
    int i = 0;
    char *name = record[i].name_chromosome;
    while (name) {
        if (strcmp(name, key) == 0)
            return record[i].num_chromosome;
        name = record[++i].name_chromosome;
    }
    return 0;
}

      folder="";
      indir="";
      outdir="";
      outdir2="";
      outdir3="";
      folderfrag="";
      folderwrite="";
      folderrefC="";
      folderrefmap="";
      folderN="";
      ALLELE="no";
      PWEIGHT="no";
      lane="DATA";
      BAITFLAG="no";
      resolution=RESOLUTION;
      map_threshold=0.20;
      weight_threshold=0.25;
      outprefix="output";
      tolerance=(double)(0.00001);
      PVAL=(double)(0.00001);
      FDR=(double)(0.01);
      max_iter=100;
      corenodes=12;
      REwidth=3000;
      stats="none";
      blunt="no";
      CUT=1;
      ENZYME="AAGCTT";
      TRIM=100;
      SCALE=1;
      GRID=0;
      PATH="PATH";
      PATH2="PATH2";
      PATHBOWTIE="PATHBOWTIE";
      window=500000;
      chrom_start=1;
      chrom_end=25;
      TRANS="no";
      REFERENCE="hg19";

      if(argc ==1){
	fprintf(stderr, "This program needs arguments .... \n\n");
	print_usage(); exit(0);}

      while((option = getopt(argc, argv, "E:C:l:f:F:D:hVv:S:s:e:r:c:m:w:t:M:J:W:G:B:R:L:P:b:N:Z:A:T:Q:H:U:p:I:O:X:o:")) != -1)

	switch (option)
	  {
        case 'E' :           ENZYME=optarg;                 break;
        case 'C' :           CUT=atoi(optarg);              break;
        case 'l' :           lane=optarg;                   break;
        case 'f' :           folder=optarg;                 break;       
        case 'F' :           folderfrag=optarg;             break;       
        case 'D' :           folderwrite=optarg;            break;
        case 'h' :           print_usage();                 break;
        case 'V' :	     printf("%s\t%s\n",PACKAGE, VERSION);              exit(0);             break;
        case 'v' :           TRANS=optarg;                  break;
	case 'S' :           stats=optarg;                  break;
        case 's' :           chrom_start=atoi(optarg);      break;
        case 'e' :           chrom_end=atoi(optarg);        break;
        case 'r' :           resolution=atoi(optarg);       break;
        case 'c' :           weight_threshold=atof(optarg); break;
        case 'm' :           map_threshold=atof(optarg);    break;
        case 'w' : 	     window=atoi(optarg);           break;
        case 't' :           tolerance=atof(optarg);        break;
        case 'M' :           max_iter=atoi(optarg);         break;
        case 'J' :           corenodes=atoi(optarg);        break;
        case 'W' :           REwidth=atoi(optarg);          break;
        case 'G' :           GRID=atoi(optarg);             break;
        case 'B' :           BAITFLAG=optarg;               break;
        case 'R' :           PATH=optarg;                   break;
        case 'L' :           PATH2=optarg;                  break;
        case 'P' :           PATHBOWTIE=optarg;             break;
        case 'b' :           blunt=(optarg);              break; 
	case 'N' :           folderN=optarg;                    break;
        case 'Z' :           SCALE=atoi(optarg);            break;
        case 'A' :           ALLELE=optarg;                 break;
        case 'T' :           TRIM=atoi(optarg);             break;
	case 'Q' :           folderrefmap=optarg;	    break;
	case 'H' :	     REFERENCE=optarg;		    break;
	case 'U' :	     FDR=atof(optarg);		    break;
        case 'p' :	     PVAL=atof(optarg);		    break;
        case 'I' :           indir=optarg;                  break;
        case 'O' :           outdir=optarg;                 break;
        case 'X' :           PWEIGHT=optarg;                break;
        case 'o' :	     outprefix=optarg;		    break;
	  }
sprintf(zom,"%d",SCALE);
sprintf(wind,"%d",window);
sprintf(stamp,concat(concat(concat("/NewAnalysis_Z",zom),"_W"),wind)); 
    

 init(lane,folder,folderN,PATHBOWTIE,chrom_start,chrom_end,stamp,stamp2,stats);

WIN_SIZE=strlen(ENZYME);
if(CUT*2>strlen(ENZYME)){CUT=strlen(ENZYME)-CUT;
printf("\n You provided a -C parameter that needs to be complemented: %d should read %d .\n",strlen(ENZYME)-CUT,CUT);}

char str1[WIN_SIZE-CUT],str2[WIN_SIZE-CUT];

strncpy(str1,ENZYME,WIN_SIZE-CUT);
strncpy(str2,ENZYME,WIN_SIZE-CUT);

str1[WIN_SIZE-CUT]='\0';
str2[WIN_SIZE-CUT]='\0';

printf("%d\t%s\n",CUT,str2);
strrev(str1);
for (i=0;i<strlen(str1);i++){
if(str1[i]=='A'){str1[i]='T';}
else if(str1[i]=='T'){str1[i]='A';}
else if(str1[i]=='C'){str1[i]='G';}
else if(str1[i]=='G'){str1[i]='C';}
}

if(strcasecmp(blunt,"no")!=0){sprintf(JUNCTION,"%s",concat(concat(str2,str1),""));}else{sprintf(JUNCTION,"%s",str2);}

mkdir(folderfrag,0777);
folderref=folderfrag;

BAITS="GRIDPOINTS_CHR";
if(strcasecmp(BAITFLAG,"no")!=0){BAITS="FRAGMENTS_CHR";}

printf(" %s genome: \n ",REFERENCE);
fileoutput1=concat(concat(PATHBOWTIE,REFERENCE),".chrom.sizes");
printf("%s\n",fileoutput1);
if(fopen(fileoutput1, "rt" )==NULL){ 
sprintf(command,"wget -q -P %s \'http://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips/%s.chrom.sizes\'",PATHBOWTIE,REFERENCE,REFERENCE); 
system(command); 
}

/*sprintf(command,"cat %s/%s.chrom.sizes|grep -v \"_\"|wc -l > %s/chromosomes",PATHBOWTIE,REFERENCE,folder);
printf("%s",command);
system(command);*/


FILE *fchr = fopen (concat(concat(PATHBOWTIE,REFERENCE),".chrom.sizes"), "r");
chrom_start=1;
while (!feof(fchr)){fgets(lyne,120,fchr);
if(strstr(lyne, "_")== NULL){
allchrom++;
item = strtok(lyne,"\t");
strcpy(record[allchrom].name_chromosome,item);
record[allchrom].num_chromosome=allchrom;
printf("**%s**%s\t%d\t%d\n",item,record[allchrom].name_chromosome,record[allchrom].num_chromosome,number_for_key(item));
}
}

chrom_end=allchrom;

        printf("       #################### %s %s ######################\n",PACKAGE,VERSION);
        printf("##                       Hi-Res-Hi-C                            ##\n");
        printf("##                                                              ##\n");
        printf("##          Algorithm by Gabriele Migliorini,                   ##\n");
        printf("##                 Code version %s                              ##\n",VERSION);
        printf("##                       GPLv3                                  ##\n");
        printf(" ################################################################ \n");
        printf("                                                                  \n");

        printf("ENZYME          : %s\n",ENZYME);                printf("CUT_POINT         : %d\n",CUT);
        printf("LANE            : %s\n",lane);                  printf("Folder            : %s\n",folder);
        printf("FolderRef       : %s\n",folderref);             printf("BAIT FILE         : %s\n",BAITFLAG);
        printf("ChromosomeStart : %d\n",chrom_start);           printf("ChromosomeEnd     : %d\n",chrom_end);
        printf("Map threshold   : %lf\n",map_threshold);         printf("blunt           : %s\n",blunt);
        printf("Weight threshold: %lf\n",weight_threshold);     printf("Tolerance         : %lf\n",tolerance);
        printf("Max number of it: %d\n",max_iter);              printf("Minimum Resolution: %d\n",resolution);
        printf("JUNCTION        : %s\n",JUNCTION);              printf("Output Prefix     : %s\n",outprefix);
        printf("Zooming factor  : %d\n",SCALE);                 printf("Window          : %d\n",window);
        printf("Statistics      : %s\n",stats);                 printf("CONFIG:           : %s\n",folderN);
        printf("READ LENGTH     : %d\n",TRIM);                  printf("ALLELE:           : %s\n",ALLELE);
        printf("NUMBER OF CORES : %d\n",corenodes);             printf("NGS OUTPUT        : %s\n",folderN);
        printf("INTER-CHROMOSOME ANALYSIS     : %s\n",TRANS);
        printf("PATHR	: %s\n",PATH); printf("PATHL   : %s\n",PATH2);
        printf("Input Sam : %s\n",indir);             printf("Output : %s\n",outdir);

printf("\n %d Chromosomes will be considered \n",chrom_end);
 
if(strlen(folderN)==0){folderN=folder; }else{ sprintf(zom,"%d",SCALE);

      sprintf(stamp,concat(concat(concat("/NewAnalysis_Z",zom),"_W"),wind));
      mkdir(concat(concat(folder,lane),stamp),0777); mkdir(concat(concat(concat(concat(folder,lane),stamp),"/weights"),stamp2),0777);
      mkdir(concat(concat(concat(concat(concat(folder,lane),stamp),"/"),stats),stamp2),0777);
      mkdir(concat(concat(concat(concat(concat(concat(folder,lane),stamp),"/"),stats),stamp2),"/fdr"),0777);
      mkdir(concat(concat(concat(folder,lane),stamp),"/graphics"),0777);
}

sprintf(command,"rm -r -f %s",folderref);
system(command);

      folderref=concat(concat(concat(concat(concat(folder,lane),stamp),"/reference_HiC"),stamp2),"/");
      mkdir(concat(concat(concat(folder,lane),stamp),"/logs"),0777);
      mkdir(folderref,0777);
      filelogs=concat(concat(concat(concat(concat(folder,lane),stamp),"/logs/log"),stamp2),".log");
      printf("%s\n%s\n",concat(concat(concat(concat(concat(folder,lane),stamp),"/logs/log"),stamp2),".log"),filelogs); 

if(folderN==folder){
   sprintf(command,"cat %s/new.c|tail -n 400|head -n 398 > %s%s.sh",folderrefC,folderN,lane);
   system(command);
   sprintf(command,"sed -i 's/TRIMLENGTH/%d/g' %s%s.sh",TRIM,folderN,lane);
   system(command);
   sprintf(command,"sed -i 's/CHROMEND/%d/g' %s%s.sh",allchrom,folderN,lane);
   system(command);
if(strcasecmp(blunt,"no")!=0){
   sprintf(command,"sed -i 's/rm /rm /g' %s%s.sh",folderN,lane);}
   system(command);
   sprintf(command,"sed -i 's|FLOWCELL|%s|g' %s%s.sh",lane,folderN,lane);
   system(command);
   sprintf(command,"sed -i 's|NODES|%d|g' %s%s.sh",corenodes,folderN,lane);
   system(command);
   sprintf(command,"sed -i 's|ALLELE|%s|g' %s%s.sh",ALLELE,folderN,lane);
   system(command);
   sprintf(command,"sed -i 's|REFERENCEFOLDER|%s|g' %s%s.sh",folderref,folderN,lane);
   system(command);
   sprintf(command,"sed -i 's|REFERENCEFOLDEC|%s|g' %s%s.sh",folderrefC,folderN,lane);
   system(command);
   sprintf(command,"sed -i 's|HOMEFOLDER|%s|g' %s%s.sh",folderwrite,folderN,lane);
   system(command);
   sprintf(command,"sed -i 's|DITAG|%s|g' %s%s.sh",JUNCTION,folderN,lane);
   system(command);
   sprintf(command,"sed -i 's|NEW|%s|g' %s%s.sh",lane,folderN,lane);
   system(command);
   sprintf(command,"sed -i 's|REFBOWTIE|%s|g' %s%s.sh",PATHBOWTIE,folderN,lane);
   system(command);
   sprintf(command,"sed -i 's|STARTPATHL|%s|g' %s%s.sh",PATH,folderN,lane);
   system(command);
   sprintf(command,"sed -i 's|STARTPATHR|%s|g' %s%s.sh",PATH2,folderN,lane);
   system(command);
   sprintf(command,"sed -i 's|CUTTER|%s|g' %s%s.sh",ENZYME,folderN,lane);
   system(command);

   if(strcasecmp(BAITFLAG,"no")!=0){BAITS="BAITS_CHR";

   sprintf(command,"sed -i 's|BAITS|%s|g' %s%s.sh",BAITFLAG,folderN,lane);}else{
   sprintf(command,"sed -i 's|BAITS|%s|g' %s%s.sh",ENZYME,folderN,lane);  }
   system(command);
   sprintf(command,"sed -i 's|REFERENCE|%s|g' %s%s.sh",REFERENCE,folderN,lane);
   system(command);

   sprintf(command,"chmod a+x %s%s.sh",folderN,lane);

   system(command);

sprintf(command,"sed -i 's/TITLE/%s/g' %s%s.sh",lane,folderN,lane);
system(command);

sprintf(command,"sed -i '1i #!/bin/bash ' %s%s.sh",folderN,lane);
system(command);


printf("\n\nrunning the NGS shell script: %s%s.sh. \nCheck the NGS folders for completion.\n\n",folderN,lane);
sprintf(command,"%s%s.sh\n\n",folderN,lane);
printf("%s",command);
system(command);
/*sprintf(command,"mv %s%s.sh %s\n\n",folderN,lane,folderN);
printf("%s",command);
system(command);*/
   }else{    printf("You are on stage II");

      fileref=concat(folderref,"/GRIDPOINTS_CHR");
      filerefb=concat(folderref,"/ALLFRAGMENTS_CHR");

      fileref2=concat(folderref,"/grid_chromosomes");
        printf("                                                                  \n");	      
        printf("       #################### %s %s ######################\n",PACKAGE,VERSION);
        printf("##                       Hi-Res-Hi-C                            ##\n");
        printf("##                                                              ##\n");
        printf("##          Algorithm by Gabriele Migliorini,                   ##\n");
	printf("##                 Code version %s                              ##\n",VERSION);
	printf("##                       GPLv3                                  ##\n");
	printf(" ################################################################ \n");
	printf("                                                                  \n");

        printf("ENZYME          : %s\n",ENZYME);                printf("CUT_POINT         : %d\n",CUT);
        printf("LANE            : %s\n",lane);   	        printf("Folder            : %s\n",folder);
	printf("FolderRef       : %s\n",folderref);             printf("BAIT FILE         : %s\n",BAITFLAG);
        printf("ChromosomeStart : %d\n",chrom_start);           printf("ChromosomeEnd     : %d\n",chrom_end);
	printf("Map threshold   : %lf\n",map_threshold);         printf("Blunt           : %s\n",blunt);                  
        printf("Weight threshold: %lf\n",weight_threshold);	printf("Tolerance         : %lf\n",tolerance);
	printf("Max number of it: %d\n",max_iter);  	        printf("Minimum Resolution: %d\n",resolution);
	printf("JUNCTION        : %s\n",JUNCTION);              printf("Output Prefix     : %s\n",outprefix);
        printf("Zooming factor  : %d\n",SCALE);                 printf("Window            : %d\n",window);
        printf("Statistics      : %s\n",stats);                 printf("CONFIG:           : %s\n",folderN);
        printf("READ LENGTH     : %d\n",TRIM);                  printf("ALLELE:           : %s\n",ALLELE);
        printf("NUMBER OF CORES : %d\n",corenodes);             printf("NGS OUTPUT        : %s\n",folderN);
        printf("INTER-CHROMOSOME ANALYSIS     : %s\n",TRANS);

        printf("Input Sam : %s\n",indir);             printf("Output : %s\n",outdir3);
        
        if(GRID!=0){printf("\nGrid type         : Uniform bins of size %d\n",SCALE*resolution);}
	if(GRID==0&&SCALE!=1){printf("\nGrid type      : Restriction metafragments of size: %d times %s average cut frequency\n",SCALE,ENZYME);}
	if(GRID==0&&SCALE==1){printf("\nGrid type      : Single fragment resolution: Restriction fragment used is: %s\n",ENZYME);}

        printf("\n\n\n\n\n .. UCSC Reference Genome: %s : Enumerating baited/non-baited fragments.\n\n\n",REFERENCE);    

enzyme = malloc (sizeof *enzyme * (chrom_end-chrom_start));
enriched = malloc (sizeof *enriched * (chrom_end-chrom_start));
gridpoints = malloc (sizeof *gridpoints * (chrom_end-chrom_start));
baitpoints = malloc (sizeof *baitpoints * (chrom_end-chrom_start));

if (enzyme == NULL||enriched == NULL||gridpoints == NULL)
{
puts("\nFailure to allocate room for pointers");
exit(0);
 }
nchromosomes=0;
tot=0; tot2=0; tot3=0; tot6=0;
fileoutput1=concat(concat(outdir3,"_pairs_ontarget"),".hic");
fileoutput2=concat(concat(outdir3,"_pairs_bfide"),".hic");
fileoutput3=concat(concat(outdir3,"_pairs_offgrid"),".hic");
fileoutput4=concat(concat(outdir3,"_pairs_pruned"),".hic");

FILE *peakyout = fopen(concat(concat(outdir,lane),".baited_frags.tsv"), "w+");
FILE *peakyout2 = fopen(concat(concat(outdir,lane),".tsv"), "w+");

if(strcasecmp(BAITFLAG,"no")!=0){
fprintf(peakyout2,"chrchrom\tchromStart\tchromEnd\tID\n");}

totfragments[0]=0;
totbaits[0]=0;

for(chrom=chrom_start;chrom<=chrom_end;chrom++){
sprintf(chromosome,"%d",chrom);
digest(chrom,allchrom,record[chrom].name_chromosome,rf,enzyme,gp,gridpoints,ENZYME,CUT,SCALE,GRID,PATHBOWTIE,REFERENCE,resolution); 
/*if(gp[chrom]==0){gp[chrom]+=1;}*/
for(i=1;i<=gp[chrom];i++){  enriched[chrom][i]=0; }
for(i=1;i<=rf[chrom];i++){  baitpoints[chrom][i]=0; }
tot5=0;
if(strcasecmp(BAITFLAG,"no")!=0){

fchic=fopen(BAITFLAG,"r+");

while (!feof(fchic)){fscanf(fchic,"%s %d %d %s",&chromosome2,&position,&position2,&baitname);
if(strcasecmp(chromosome2,record[chrom].name_chromosome)==0){
for(i=cgrid(chrom,(int)(position+1),gp,gridpoints);i<=cgrid(chrom,(int)(position2-2*strlen(ENZYME)),gp,gridpoints);i++){
if(enriched[chrom][i]==0){
tot5++;tot6++;
sprintf(stringarray[tot6],"%s",&baitname);
bp[chrom]++; baitpoints[chrom][tot5]=gridpoints[chrom][i]+RESOLUTION*GRID;
printf("Baited_Bin_number:%d, Targeted Gene: %s, Position:Chr%d:%d-%d.\n",tot6,stringarray[tot6],chrom,baitpoints[chrom][tot5],baitpoints[chrom][tot5]+resolution*SCALE);
enriched[chrom][i]=tot6;}
}
}
}
printf("\nBaits on this chromosome:%d\tTotal number of Baits: %d\n",tot5,tot6);
fclose(fchic);

sprintf(chromosome,"%d",chrom);

fout2b=fopen(concat(filerefb,chromosome),"w+");

for(i=1;i<=gp[chrom];i++){enriched[chrom][i]=0;}
for(i=1;i<=rf[chrom];i++){
flag=0;
for(j=1;j<=bp[chrom];j++){

if((enzyme[chrom][i]==baitpoints[chrom][j])){flag=1;}
}
fprintf(fout2b,"%s\t%d\t%d\t%d\n",record[chrom].name_chromosome,enzyme[chrom][i],enzyme[chrom][i+1],flag);
}
rewind(fout2b);
}else{
for(i=1;i<=gp[chrom];i++){
baitpoints[chrom][i]=gridpoints[chrom][i]; 
}
tot6+=(gp[chrom]-1);
bp[chrom]=gp[chrom]; 

fout2b=fopen(concat(filerefb,chromosome),"w+");
for(i=1;i<rf[chrom]-1;i++){
fprintf(fout2b,"%s\t%d\t%d\t%d\n",record[chrom].name_chromosome,enzyme[chrom][i],enzyme[chrom][i+1],1);}
rewind(fout2b);
}
nchromosomes++; 
totbaits[nchromosomes]=tot6;
tot+=(gp[chrom]-1);  
totfragments[nchromosomes]=tot;
if((rf[chrom]-1)>tot3){tot3=rf[chrom]-1; BIN=rint((enzyme[chrom][rf[chrom]-1])/(rf[chrom]-1))+1; }

fout = fopen(filelogs,"a+");
fout3 = fopen(fileref2,"a+");

printf("Chromosomes analyzed: %d\tFragments on Chromosome %d: %d\t Total number of Grid Points:%d\t Number of Baited Bins :%d\n",nchromosomes,chrom,rf[chrom],totfragments[nchromosomes],totbaits[nchromosomes]);

while(!feof(fout2b)){fscanf(fout2b,"%s %d %d %d %d",&chromosome2,&position,&position2,&enrichment);
  if(enrichment==1){
   for(i=cgrid(chrom,(int)(position+1),gp,gridpoints);i<=cgrid(chrom,(int)(position2-1),gp,gridpoints);i++){
    if(enriched[chrom][i]==0){ 
    tot2++; enriched[chrom][i]= tot2;
   }
  }
 }
}
fout2 = fopen(concat(fileref,chromosome),"w+");
for(i=1;i<=gp[chrom]-1;i++){
fprintf(fout2,"%s\t%d\t%d\t%d\n",record[chrom].name_chromosome,gridpoints[chrom][i],gridpoints[chrom][i+1],enriched[chrom][i]);
if(enriched[chrom][i]!=0){
fprintf(peakyout,"%s\t%d\t%d\t%s\n",record[chrom].name_chromosome,gridpoints[chrom][i],gridpoints[chrom][i+1],stringarray[enriched[chrom][i]]);
}

fprintf(peakyout2,"%s\t%d\t%d\t%d\n",record[chrom].name_chromosome,gridpoints[chrom][i],gridpoints[chrom][i+1],totfragments[chrom-1]+i);
}
fclose(fout2);
fclose(fout2b);
fprintf(fout3,"%s\t%d\t%d\n",record[chrom].name_chromosome,1,gridpoints[chrom][gp[chrom]-1]);
fclose(fout3);
}
fclose(peakyout);
fclose(peakyout2);
printf("closing files\n");
if(tot2!=tot6){printf("check baited input fragments\n");exit;}

fout4 = fopen(fileoutput2,"w+");
fout4b =fopen(fileoutput1,"w+");
fout5 = fopen(fileoutput3,"w+");
fout6 = fopen(fileoutput4,"w+");

totfragments[0]=0;
map=malloc(tot * sizeof(double));
map_count=malloc(tot * sizeof(int));
hic=malloc((tot2+1) * sizeof(int *));    
chromgrid=malloc(tot * sizeof(int));
chrombait=malloc(tot2 * sizeof(int));
 
for(i=1;i<tot; i++){  *(map+i)=0.0; *(map_count+i)=0; *(chromgrid+i)=0;}

for(i=1;i<tot2+1; i++){ *(chrombait+i)=0;
hic[i] = malloc(sizeof *hic[i] * tot); }

for(i=1;i<tot2; i++){
for(j=1;j<tot;j++){
  *(hic[i]+j)=0; }}

printf("Reading baits position on chromosome %d\t Enriched fragments: %d\n",chrom-1,(int)(tot2));
 
tot4=(long)tot*(long)tot2;

printf("\nYou allocated a %d x %d matrix with %li entries. \n\n\n",tot,tot2,tot4);

printf("\nThe average Bin Size on the chromosomes under study is: %d\n",BIN*SCALE);

fprintf(fout,"\nYou allocated a %d x %d matrix with %li entries. \n\n\n",tot,tot2,tot4);

printf("\n You generated input files from chrom %d\t to chrom %d\n",chrom_start,chrom_end);

for(chrom=chrom_start;chrom<=chrom_end;chrom++){

for(i=1;i<=gp[chrom];i++){ chromgrid[i+totfragments[chrom-1]]=chrom; } sprintf(chromosome,"%d",chrom);
for(i=1;i<=bp[chrom];i++){ chrombait[i+totbaits[chrom-1]]=chrom;    }
printf("%s\n",concat(concat(folderrefmap,"mappability_"),record[chrom].name_chromosome));
fmap = fopen(concat(concat(folderrefmap,"mappability_"),record[chrom].name_chromosome), "r");
if(fmap==NULL){printf("mappability files is empty or missing\nPlease refer to the manual\n\n\n\n\n"); exit(0);}

while (!feof(fmap)){fscanf(fmap,"%s %d %lf",&chromosome2,&position,&mappability);
   if(position<gp[chrom]*SCALE*resolution){temp=totfragments[chrom-chrom_start]+cgrid(chrom,position,gp,gridpoints)+1;
   map[temp]+=mappability; map_count[temp]++;}
} 
   fclose(fmap);
}
   
for(i=1;i<=tot; i++){
  if(map_count[i]!=0){map[i]=map[i]/map_count[i];}else{map[i]=0.0;}
}

for(chrom=chrom_start;chrom<=chrom_end;chrom++){
  
sprintf(schr,"%d",chrom);    

if(strcasecmp(BAITFLAG,"no")==0){sprintf(suffix,".dedup.sam");}
if(strcasecmp(BAITFLAG,"no")!=0){sprintf(suffix,".dedup.sam");}

FILE *fin = fopen (concat(concat(concat(indir,".pairs."),record[chrom].name_chromosome),".dedup.sam") , "r" );
if(fin==NULL){printf("ditag file is empty or missing\n"); exit(0);}
fprintf(fout,"\nReading Files:%s\n",concat(concat(concat(concat(indir,".pairs."),record[chrom].name_chromosome),suffix),""));
printf("Files Used:%s\n",concat(concat(concat(concat(indir,".pairs."),record[chrom].name_chromosome),suffix),""));
while (!feof(fin))
{

fscanf(fin,"%s %s %d %d %s %s %d %d %d %s %s %d %s %s %s %s",&STRAND1, &chromosome, &POS1, &FRAG1, &STRAND2, &chromosome2, &POS2, &FRAG2, &QUAL1,&CIGAR1,&SEQ1,&QUAL2,&CIGAR2,&SEQ2,&QNAME1,&QNAME2);

CHR1=number_for_key(chromosome);
CHR2=number_for_key(chromosome2);

flag1=0;
flag2=0;

if((CHR1>=chrom_start&&CHR1<=chrom_end)&&(CHR2>=chrom_start&&CHR2<=chrom_end)){

qc=0;
save_a=grid(CHR1,POS1,rf,enzyme);
save_b=grid(CHR2,POS2,rf,enzyme);

 if(strncmp(STRAND1,"0",1)==0&&strncmp(STRAND2,"0",1)==0){flag1=1;flag2=1; if(save_a==save_b){qc=3;} if((save_a==save_b+1&&CHR1==CHR2)||(save_b==save_a+1&&CHR1==CHR2)){qc=6;}}
 if(strncmp(STRAND1,"16",2)==0&&strncmp(STRAND2,"16",2)==0){flag1=0;flag2=0; if(save_a==save_b){qc=3;} if((save_a==save_b+1&&CHR1==CHR2)||(save_b==save_a+1&&CHR1==CHR2)){qc=6;}}
 if(strncmp(STRAND1,"16",2)==0&&strncmp(STRAND2,"0",1)==0){flag1=0;flag2=1; if(save_a==save_b&&CHR1==CHR2&&POS1<POS2){qc=2;}if(save_a==save_b&&CHR1==CHR2&&POS1>=POS2){qc=1;};if((save_a==save_b+1&&CHR1==CHR2)||(save_b==save_a+1&&CHR1==CHR2)){if(POS1<POS2){qc=5;};if(POS1>=POS2){qc=4;}}}                                 
 if(strncmp(STRAND1,"0",1)==0&&strncmp(STRAND2,"16",2)==0){flag1=1;flag2=0; if(save_a==save_b&&CHR1==CHR2&&POS1<POS2){qc=1;}if(save_a==save_b&&CHR1==CHR2&&POS1>=POS2){qc=2;};if((save_a==save_b+1&&CHR1==CHR2)||(save_b==save_a+1&&CHR1==CHR2)){if(POS1<POS2){qc=4;};if(POS1>=POS2){qc=5;}}}

d1=0;d2=0;
 if(strncmp(STRAND1,"0",1)==0){d1=fabs(POS1-enzyme[CHR1][grid(CHR1,POS1,rf,enzyme)+flag1]);
                               d1p=fabs(POS1-baitpoints[CHR1][bgrid(CHR1,POS1,bp,baitpoints)+flag1]);}
 if(strncmp(STRAND1,"16",2)==0){d1=fabs(POS1-enzyme[CHR1][grid(CHR1,POS1-strlen(SEQ1),rf,enzyme)+flag1])-(strlen(ENZYME));
                               d1p=fabs(POS1-baitpoints[CHR1][bgrid(CHR1,POS1-strlen(SEQ1),bp,baitpoints)+flag1]-strlen(ENZYME));}
 if(strncmp(STRAND2,"0",1)==0){d2=fabs(POS2-enzyme[CHR2][grid(CHR2,POS2,rf,enzyme)+flag2]);
                               d2p=fabs(POS2-baitpoints[CHR2][bgrid(CHR2,POS2,bp,baitpoints)+flag2]);}
 if(strncmp(STRAND2,"16",2)==0){d2=fabs(POS2-enzyme[CHR2][grid(CHR2,POS2-strlen(SEQ2),rf,enzyme)+flag2])-(strlen(ENZYME));
                               d2p=fabs(POS2-baitpoints[CHR2][bgrid(CHR2,POS2-strlen(SEQ2),bp,baitpoints)+flag2]-strlen(ENZYME));}
 if(qc==0){
     if(enriched[CHR1][cgrid(CHR1,POS1,gp,gridpoints)]!=0){

hic[enriched[CHR1][cgrid(CHR1,POS1,gp,gridpoints)]][totfragments[CHR2-chrom_start]+cgrid(CHR2,POS2,gp,gridpoints)]+=1.0; 

fprintf(fout4,"%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%f\t%f\t%6.5f\t%6.5f\t%d\t%d\n",STRAND1,CHR1,POS1,save_a,STRAND2,CHR2,POS2,save_b,QUAL1,CIGAR1,SEQ1,QUAL2,CIGAR2,SEQ2,QNAME1,QNAME2,d1,d2,d1p,d2p,enriched[CHR1][cgrid(CHR1,POS1,gp,gridpoints)],enriched[CHR2][cgrid(CHR2,POS2,gp,gridpoints)]);
 
if(((d1p<=d1)&&(d2p>d2))||((d1p>d1)&&(d2p<=d2))){
fprintf(fout4b,"%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%f\t%f\t%6.5f\t%6.5f\n",STRAND1,CHR1,POS1,save_a,STRAND2,CHR2,POS2,save_b,QUAL1,CIGAR1,SEQ1,QUAL2,CIGAR2,SEQ2,QNAME1,QNAME2,d1,d2,d1p,d2p);
}
}else{
fprintf(fout5,"%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%f\t%f\t%6.5f\t%6.5f\t%d\t%d\n",STRAND1,CHR1,POS1,save_a,STRAND2,CHR2,POS2,save_b,QUAL1,CIGAR1,SEQ1,QUAL2,CIGAR2,SEQ2,QNAME1,QNAME2,d1,d2,d1p,d2p,enriched[CHR1][cgrid(CHR1,POS1,gp,gridpoints)],enriched[CHR2][cgrid(CHR2,POS2,gp,gridpoints)]);
     }}else{

fprintf(fout6,"%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%f\t%f\t%6.5f\t%6.5f\t%d\n",STRAND1,CHR1,POS1,save_a,STRAND2,CHR2,POS2,save_b,QUAL1,CIGAR1,SEQ1,QUAL2,CIGAR2,SEQ2,QNAME1,QNAME2,d1,d2,d1p,d2p,qc);
  }
 }
}
fclose(fin); 
}
         
fclose(fout4);
fclose(fout4b);
fclose(fout5);
fclose(fout6);

FILE *tads= fopen(concat(concat(concat(concat(outdir,lane),"_Z"),zom),".tad"), "w+");
FILE *peakY = fopen(concat(concat(concat(concat(outdir,lane),"_Z"),zom),".pkY"), "w+");
FILE *peakY2 = fopen(concat(concat(concat(concat(outdir,lane),"_Z"),zom),".cool"), "w+");
FILE *gff;
 FILE *peakC; 
FILE *wig= fopen(concat(concat(concat(concat(outdir,lane),"_Z"),zom),".wig"), "w+");

pchic=0;
if(strcasecmp(BAITFLAG,"no")!=0){
pchic=1;}

fprintf(peakY,"baitID\tpreyID\tN\tdist\tdens\tdensdist\n");

for(i=1;i<=tot2;i++){
chrom2=chrombait[i];
sprintf(sbait,"%d",i);

if(strcasecmp(stringarray[i],"")!=0){sprintf(sbait,"%s",stringarray[i]);}

FILE *peakC = fopen(concat(concat(concat(concat(concat(concat(outdir2,lane),"_BAIT_"),sbait),"_Z"),zom),".pkC"), "w+");
FILE *gff = fopen(concat(concat(concat(concat(concat(concat(outdir2,lane),"_BAIT_"),sbait),"_Z"),zom),".gff"), "w+");

if(pchic==1){
}else{
fprintf(tads,"\n%s\t%d\t%d\t",record[chrombait[i]].name_chromosome,gridpoints[chrombait[i]][i-totbaits[chrombait[i]-1]],gridpoints[chrombait[i]][i-totbaits[chrombait[i]-1]+1]);}

mid=0;
reg_size=0;
reg_size2=0;
for(j=1;j<tot;j++){
chrom=chromgrid[j];
if(chrombait[i]==chromgrid[j]){
fdist=fabs(gridpoints[chrom][cgrid(chrom,baitpoints[chrom][i-totbaits[chrom-1]],gp,gridpoints)]-gridpoints[chrom][cgrid(chrom,gridpoints[chrom][j-(totfragments[chrom-1])],gp,gridpoints)]);

if(fdist<(double)(window)&&totfragments[chrombait[i]-1]+cgrid(chrombait[i],baitpoints[chrombait[i]][i-totbaits[chrombait[i]-1]]+RESOLUTION,gp,gridpoints)>j){reg_size+=hic[i][j];}
if(fdist<(double)(window)&&totfragments[chrombait[i]-1]+cgrid(chrombait[i],baitpoints[chrombait[i]][i-totbaits[chrombait[i]-1]]+RESOLUTION,gp,gridpoints)<j){reg_size2+=hic[i][j];}

}else{fdist=0.0;}
}

gdist=0;
for(j=1;j<tot;j++){
chrom=chromgrid[j];

if(chrombait[i]==chromgrid[j]){
fdist=fabs(gridpoints[chrom][cgrid(chrom,baitpoints[chrom][i-totbaits[chrom-1]],gp,gridpoints)]-gridpoints[chrom][cgrid(chrom,gridpoints[chrom][j-(totfragments[chrom-1])],gp,gridpoints)]);
save=gdist;

if(fdist<(double)(window)&&totfragments[chrombait[i]-1]+cgrid(chrombait[i],baitpoints[chrombait[i]][i-totbaits[chrombait[i]-1]]+RESOLUTION,gp,gridpoints)>j){gdist+=hic[i][j];}
if(fdist<(double)(window)&&totfragments[chrombait[i]-1]+cgrid(chrombait[i],baitpoints[chrombait[i]][i-totbaits[chrombait[i]-1]]+RESOLUTION,gp,gridpoints)==j){gdist=0;}
if(fdist<(double)(window)&&totfragments[chrombait[i]-1]+cgrid(chrombait[i],baitpoints[chrombait[i]][i-totbaits[chrombait[i]-1]]+RESOLUTION,gp,gridpoints)<j){gdist+=hic[i][j];}
if(pchic==1){

fprintf(peakC,"%d\t%d\n",gridpoints[chrom][j-(totfragments[chrom-1])],hic[i][j]);
if(fdist<=(double)(window)){fprintf(wig,"%s\t%d\t%d\t%d\n",record[chrom].name_chromosome,gridpoints[chrom][j-(totfragments[chrom-1])],gridpoints[chrom][j-(totfragments[chrom-1])+1],hic[i][j]);
}}else{
fprintf(tads,"%d\t",hic[i][j]);
}
}else{
fdist=100000000;}
if(hic[i][j]>0){
if(j>=i){fprintf(peakY2,"%d\t%d\t%d\n",i,j,hic[i][j]);}
if(chrombait[i]!=chromgrid[j]){
fprintf(peakY,"%d\t%d\t%d\t%f\t%f\t%f\n",totfragments[chrombait[i]-1]+cgrid(chrombait[i],baitpoints[chrombait[i]][i-totbaits[chrombait[i]-1]]+(1-pchic)*RESOLUTION,gp,gridpoints),j,hic[i][j],fdist,1,1);
}else{

if(fdist<(double)(window)){if((totfragments[chrombait[i]-1]+cgrid(chrombait[i],baitpoints[chrombait[i]][i-totbaits[chrombait[i]-1]]+RESOLUTION,gp,gridpoints)>j)){
fprintf(peakY,"%d\t%d\t%d\t%f\t%f\t%f\n",totfragments[chrombait[i]-1]+cgrid(chrombait[i],baitpoints[chrombait[i]][i-totbaits[chrombait[i]-1]]+(1-pchic)*RESOLUTION,gp,gridpoints),j,hic[i][j],fdist,((reg_size-gdist)/reg_size),((reg_size-gdist)/reg_size)*window);
}else if((totfragments[chrombait[i]-1]+cgrid(chrombait[i],baitpoints[chrombait[i]][i-totbaits[chrombait[i]-1]]+RESOLUTION,gp,gridpoints)<j)){
fprintf(peakY,"%d\t%d\t%d\t%f\t%f\t%f\n",totfragments[chrombait[i]-1]+cgrid(chrombait[i],baitpoints[chrombait[i]][i-totbaits[chrombait[i]-1]]+(1-pchic)*RESOLUTION,gp,gridpoints),j,hic[i][j],fdist,(gdist/reg_size2),(gdist/reg_size2)*window);
}
}
}
}
}

fclose(peakC);
fclose(gff);
}

fclose(peakY);
fclose(peakY2);
if(pchic==1){fclose(wig);}else{
fclose(tads);}

sprintf(command,"touch %s/%s/%s_domain_caller.sh",folder,lane,lane);
system(command);
sprintf(heatmap,outdir);
sprintf(graphics,concat(outdir,"/../graphics"));

sprintf(command,"echo '#!/bin/bash ' >> %s/%s/%s_domain_caller.sh",folder,lane,lane);
system(command);
sprintf(command,"echo 'HM=%s' >> %s/%s/%s_domain_caller.sh",heatmap,folder,lane,lane);
system(command);
sprintf(command,"echo 'GRAPHICS=%s' >> %s/%s/%s_domain_caller.sh",graphics,folder,lane,lane);
system(command);
if(strcasecmp(BAITFLAG,"no")==0){sprintf(command,"echo 'RES=%d' >> %s/%s/%s_domain_caller.sh",SCALE*resolution,folder,lane,lane);}else{
sprintf(command,"echo 'RES=%d' >> %s/%s/%s_domain_caller.sh",1,folder,lane,lane);}
system(command);

for(chrom=chrom_start;chrom<=chrom_end;chrom++){

sprintf(command,"echo 'CHROM=%d' >> %s/%s/%s_domain_caller.sh",chrom,folder,lane,lane);
system(command);
if(strcasecmp(BAITFLAG,"no")==0){
sprintf(command,"cat %sgraphics.sh >> %s/%s/%s_domain_caller.sh",folderrefC,folder,lane,lane);
 system(command);
}

sprintf(command,"chmod a+x %s/%s/%s_domain_caller.sh",folder,lane,lane);
sprintf(command,"%s/%s/%s_domain_caller.sh",folder,lane,lane);
}

if(strcasecmp(BAITFLAG,"no")!=0){
 
sprintf(command,"sort -k1,1 -k2,2n %s >  %s",concat(concat(concat(concat(outdir,lane),"_Z"),zom),".wig"),concat(concat(concat(concat(outdir,lane),"_Z"),zom),".wig1"));
printf("%s\n",command);
system(command);

sprintf(command,"LC_COLLATE=C;cat %s |gawk '{a=$1\"_\"$2\"_\"$3;c[a]+=$4;r[NR]=a}END{for(i=1;i<=NR;i++){print r[i] \"_\" c[r[i]]}}'|uniq|sed 's/_/\t/g'|sort -k1,1 -k2,2n >  %s",concat(concat(concat(concat(outdir,lane),"_Z"),zom),".wig1"),concat(concat(concat(concat(outdir,lane),"_Z"),zom),".wig2"));
printf("%s\n",command);
system(command);
sprintf(command,"mkdir %s ",concat(outdir,"/bins/"));
system(command);
printf("%s\n",command);
sprintf(command,"sed -i 's/chr//' %s",concat(concat(outdir,lane),".tsv"));
system(command);
printf("%s\n",command);

sprintf(command,"R --vanilla < %sPEAKY.R --args %s %s %s %s %d ",folderrefC,outdir,concat(concat(concat(concat(outdir,lane),"_Z"),zom),".pkY"),concat(outdir,"/bins/"),concat(concat(outdir,lane),".tsv"),window);
printf("%s\n",command);
system(command);

sprintf(command,"sort -k4,4 -k5,5n %s >  %s",concat(concat(concat(concat(outdir,lane),"_Z"),zom),".pkY_fdr.wig"),concat(concat(concat(concat(outdir,lane),"_Z"),zom),".pkY_fdr.wig1"));
printf("%s\n",command);
system(command);

sprintf(command," LC_COLLATE=C;cat %s |gawk '{a=$4\"_\"$5\"_\"$6;c[a]+=$NF;r[NR]=a;d[NR]=$3}END{for(i=1;i<=NR;i++){if(sqrt(d[i]*d[i])<2000000){print \"chr\"r[i] \"_\" c[r[i]]}}}'|uniq|sed 's/_/\t/g' |sort -k1,1 -k2,2n >  %s",concat(concat(concat(concat(outdir,lane),"_Z"),zom),".pkY_fdr.wig1"),concat(concat(concat(concat(outdir,lane),"_Z"),zom),".pkY_fdr.wig2"));

/*sprintf(command,"cat %s |gawk '{a=$4\"_\"$5\"_\"$6;c[a]+=$NF;r[NR]=a;d[NR]=$3}END{for(i=1;i<=NR;i++){print \"chr\"r[i] \"_\" c[r[i]]}}'|uniq|sed 's/_/\t/g' >  %s",concat(concat(concat(concat(outdir,lane),"_Z"),zom),".pkY_fdr.wig1"),concat(concat(concat(concat(outdir,lane),"_Z"),zom),".pkY_fdr.wig2"));*/
printf("%s\n",command);
system(command);

sprintf(command,"%sbedGraphToBigWig %s %s%s.chrom.sizes %s",folderrefC,concat(concat(concat(concat(outdir,lane),"_Z"),zom),".pkY_fdr.wig2"),PATHBOWTIE,REFERENCE,concat(concat(outdir,lane),"_fdr.bw"));
printf("%s\n",command);
system(command);
sprintf(command,"%sbedGraphToBigWig %s %s%s.chrom.sizes %s",folderrefC,concat(concat(concat(concat(outdir,lane),"_Z"),zom),".wig2"),PATHBOWTIE,REFERENCE,concat(concat(concat(concat(outdir,lane),"_Z"),zom),".bw"));
printf("%s\n",command);
system(command);
}
sprintf(command,"sed 's|TITLE|%s|g' %sNGS2.sh > %s%s_2.sh",lane,folderrefC,folder,lane);
system(command);

sprintf(command,"sed 's|RES|%s|g' %sNGS2.sh > %s%s_2.sh",zom,folderrefC,folder,lane);
system(command);


 sprintf(command,"sed -i 's|HOMEFOLDER|%s|g' %s%s_2.sh",folderwrite,folder,lane);
   system(command);

sprintf(command,"sed -i 's/TRIMLENGTH/%d/g' %s%s_2.sh",TRIM,folder,lane);
system(command);

sprintf(command,"sed -i 's|BAITPATH|%s/%s/|g' %s%s_2.sh",stamp,stats,folder,lane);
system(command);

sprintf(command,"sed -i 's|FLOWCELL|%s|g' %s%s_2.sh",lane,folder,lane);
system(command);

sprintf(command,"sed -i 's|REFBOWTIE|%s|g' %s%s_2.sh",PATHBOWTIE,folder,lane);
printf("%s\n",command);
system(command);

sprintf(command,"sed -i 's|REFERENCEFOLDEC|%s|g' %s%s_2.sh",folderrefC,folder,lane);
printf("%s\n",command);
system(command);

sprintf(command,"sed -i 's|REFERENCE|%s|g' %s%s_2.sh",REFERENCE,folder,lane);
printf("%s\n",command);
system(command);

sprintf(command,"sed -i '1i #!/bin/bash ' %s%s_2.sh",folder,lane);
system(command);

sprintf(command,"chmod a+x %s%s_2.sh",folder,lane);
system(command);
sprintf(command,"%s%s_2.sh\n\n",folder,lane);
printf("%s",command);
system(command);


printf("\n\nICeCAP - Iterative normalization of C-Hi-C biases \n  \n");
   }
}

double digest(chrom,chrom_end,tmpname,rf,enzyme,gp,gridpoints,ENZYME,CUT,SCALE,GRID,folderref,REFERENCE,resolution)
     int **enzyme; int rf[CHROMOSOMES]; int **gridpoints; int gp[CHROMOSOMES]; int chrom; int chrom_end; int CUT; int resolution; char *tmpname; char *ENZYME; int SCALE; int GRID; char *REFERENCE; char *folderref;
    {

    FILE *reference_genome;
    char c,str[20],*chr;      
    int win_size=0;
    int pos=0;
    int diff=0;
    
    int up,down,middle,middle2;
    int saved;
    int quality;
    int max,cuts,cuts2;
    int start,end;
    int res,b;

    long unsigned int sizeofchrom[25];
    unsigned int dist;    
    double tmp3,grid,slope,smooth,average;
    char window1[WIN_SIZE];
    char schrom[25],command[500],filename[500];

sprintf(schrom,"%d",chrom);        
max=0;
sprintf(filename,"%s/%s.fa",folderref,tmpname);
printf("Reading: %s\n",filename);

if(fopen(filename, "rt" )==NULL){ 

sprintf(command,"wget  -P %s \'http://hgdownload.cse.ucsc.edu/goldenpath/'%s'/chromosomes/%s.fa.gz\'",folderref,REFERENCE,tmpname);

printf("%s\n",command);
system(command);
sprintf(command,"gunzip -f %s/%s.fa.gz",folderref,tmpname);

printf("%s\n",command);
system(command); 
}

      sprintf(command,"%s/%s.fa",folderref,tmpname);

      reference_genome=fopen(command, "rt" );
      printf("\n*Digesting reference chromosome %s \n\n",tmpname);    

      if(strcasecmp(BAITFLAG,"no")!=0){printf("Bait \t\t Gene\n");}

      gp[chrom]=0; rf[chrom]=0; 
      while((c=fgetc(reference_genome))!=EOF)
      {  

      if(c=='>')
             {
            while((c=fgetc(reference_genome))!=EOF && c!='\n' )
                {
                }

            continue;

            }
            if(isspace(c)) continue;

            window1[win_size++]=c;
            window1[WIN_SIZE]='\0'; 
            if(win_size==WIN_SIZE)
            {
	    if((strcasecmp(window1,ENZYME) ==0&&GRID==0)||(pos%(resolution) ==0&&GRID==1)){cuts++;}	               
            memmove(window1,&window1[WIN_SHIFT1],WIN_SIZE-WIN_SHIFT1);     pos+=WIN_SHIFT1;   	      	
            win_size=WIN_SIZE-WIN_SHIFT1;
            }
	}

    pos=0;
    cuts2 = (int)(cuts/SCALE);

    win_size=0; window1[0]='\0';
    enzyme[chrom] = malloc(sizeof *enzyme[chrom] * (int) cuts);
    gridpoints[chrom] = malloc(sizeof *gridpoints[chrom] * (int) cuts2);
    enriched[chrom] = malloc(sizeof *enriched[chrom] * (int) cuts2);
    baitpoints[chrom] = malloc(sizeof *baitpoints[chrom] * (int) cuts);

    cuts=0; cuts2=0;
    reference_genome=fopen(command, "rt" );


    while((c=fgetc(reference_genome))!=EOF)
        {
        if(c=='>')
             {
            while((c=fgetc(reference_genome))!=EOF && c!='\n' )
                {
                }

            continue;

            }
            if(isspace(c)) continue;

            window1[win_size++]=c;
            window1[WIN_SIZE]='\0'; 

            if(win_size==WIN_SIZE)
            {
	      if((strcasecmp(window1,ENZYME) ==0&&GRID==0)||(pos%(resolution) ==0&&GRID==1)){
		cuts++;if(cuts%SCALE==0){gp[chrom]++;   *(gridpoints[chrom]+gp[chrom])=pos+CUT+RESOLUTION*GRID; }
                         rf[chrom]++;   *(enzyme[chrom]+rf[chrom])=pos+CUT;                 
	    }
            memmove(window1,&window1[WIN_SHIFT1],WIN_SIZE-WIN_SHIFT1);      pos+=WIN_SHIFT1;
            win_size=WIN_SIZE-WIN_SHIFT1;
	    }  
	}
    return (0);	
    }

double init(lane,folder,folderN,PATHBOWTIE,chrom_start,chrom_end,stamp,stamp2,stats)
     char *lane;char *folder; char *folderN; char *PATHBOWTIE; int *chrom_start; int *chrom_end; char *stamp; char *stamp2; char *stats;
{
  char end[1];

sprintf(stamp2,"_%s_%d_%d",outprefix,chrom_start,chrom_end);

      if(strlen(folder)==0||dirExists(folder)==0){printf("\n\n Please provide mandatory -f field: \n\n Refer to the Manual.Thank you. \n\n "); exit(0);}
      sprintf(end,"%c\0",folder[strlen(folder)-1]);
      if(strcasecmp(end,"/")!=0){printf("Adding backslash at the end of the folder path..\n"); folder=concat(folder,"/");}

      if((strlen(PATHBOWTIE)==0||dirExists(PATHBOWTIE)==0)&&strlen(folderN)==0){printf("\n\n Please provide mandatory -P field: \n\n Refer to the Manual.Thank you. \n\n "); exit(0);}
      sprintf(end,"%c\0",PATHBOWTIE[strlen(PATHBOWTIE)-1]);
      if(strcasecmp(end,"/")!=0){printf("Adding backslash at the end of the bowtie path..\n");PATHBOWTIE=concat(PATHBOWTIE,"/");}

      if(strlen(folderrefmap)==0){folderrefmap=concat(folder,"reference_map/");}
      sprintf(end,"%c\0",folderrefmap[strlen(folderrefmap)-1]);
      if(strcasecmp(end,"/")!=0){printf("Adding backslash at the end of the mapability data folder path..\n");folderrefmap=concat(folderrefmap,"/");}
     
      if(strlen(folderwrite)==0){folderwrite=folder;}else{mkdir(folderwrite,0777);}
      if(dirExists(folderwrite)==0){printf("\n\n Please provide a pre-existing -D field: \n\n Refer to the Manual.Thank you. \n\n "); exit(0);} 

      folderrefC=concat(folder,"reference_C/");

      if(strlen(indir)==0){indir=concat(concat(concat(concat(folderwrite,lane),"/data/"),lane),"");}else{
 
      sprintf(end,"%c\0",indir[strlen(indir)-1]);
      if(strcasecmp(end,"/")!=0){printf("Adding backslash at the end of the input data folder path..\n");indir=concat(indir,"/");}
                              folder=indir;
                              folderwrite=indir;
                              indir=concat(concat(concat(concat(folderwrite,lane),"/data/"),lane),"");}

      sprintf(end,"%c\0",folderwrite[strlen(folderwrite)-1]);
      if(strcasecmp(end,"/")!=0){printf("Adding backslash at the end of the NGS output folder path..");folderwrite=concat(folderwrite,"/");}

      if(strlen(folderfrag)==0){folderfrag=concat(concat(concat(concat(concat(folderwrite,"/"),lane),"/reference_HiC"),stamp2),"/");}
      if(strlen(outdir)==0){mkdir(concat(folderwrite,lane),0777);
                            mkdir(concat(concat(folderwrite,lane),stamp),0777);
                            mkdir(concat(concat(concat(concat(concat(folderwrite,lane),stamp),"/"),stats),"/"),0777);
                            outdir3=concat(concat(concat(concat(folderwrite,lane),"/data/"),lane),"");

                            mkdir(concat(concat(concat(concat(concat(folderwrite,lane),stamp),"/"),stats),"/PKY/"),0777);
                            outdir=concat(concat(concat(concat(concat(folderwrite,lane),stamp),"/"),stats),"/PKY/");
                            mkdir(concat(concat(concat(concat(concat(folderwrite,lane),stamp),"/"),stats),"/PKC/"),0777);
                            outdir2=concat(concat(concat(concat(concat(folderwrite,lane),stamp),"/"),stats),"/PKC/");

}
return(0);
}

int grid(chr,pos,rf,enzyme)
     int **enzyme; int rf[CHROMOSOMES]; int chr; int pos;
{
  int flag;
  int up,down,middle,save;

  up=rf[chr];
  down=1;
  middle=(int)((up+down)/2);
  flag=0;

  while(flag==0){
    if(pos>enzyme[chr][middle]){down=middle;};
    if(pos<enzyme[chr][middle]){up=middle;};
    if(pos==enzyme[chr][middle]){up=middle;flag=1;save=middle;}
    middle=(int)((up+down)/2);
    if(up==(down+1)){flag=2;save=down;}
  }

  if(up==(down+2)){flag=1;save=down;}
    return(save);
}


int cgrid(chr,pos,gp,gridpoints)
     int **gridpoints; int gp[CHROMOSOMES]; int chr; int pos;
{
  int flag;
  int up,down,middle,save;

  up=gp[chr];
  down=1;
  middle=(int)((up+down)/2);
  flag=0;

  while(flag==0){
    if(pos>gridpoints[chr][middle]){down=middle;};
    if(pos<gridpoints[chr][middle]){up=middle;};
    if(pos==gridpoints[chr][middle]){up=middle;flag=1;save=middle;}
    middle=(int)((up+down)/2);
    if(up==(down+1)){flag=2;save=down;} 
    if(up==down){flag=2;save=down;}
  }

  if(up==(down+2)){flag=1;save=down;}
    return(save);
}

 int bgrid(chr,pos,bp,baitpoints)
      int **baitpoints; int bp[CHROMOSOMES]; int chr; int pos;
 {
   int flag;
   int up,down,middle,save;
   up=bp[chr];
   down=1;
   middle=(int)((up+down)/2);
   flag=0;
 
   while(flag==0){
     if(pos>baitpoints[chr][middle]){down=middle;};
     if(pos<baitpoints[chr][middle]){up=middle;};
     if(pos==baitpoints[chr][middle]){up=middle;flag=1;save=middle;}
     middle=(int)((up+down)/2);
if(up==(down+1)){flag=2;save=down;}
if(up==down){flag=2;save=down;}
}
     if(up==(down+2)){flag=1;save=down;}
     return(save);
 }


/* #TITLE
## setting the run parameters ##
################################

########################
# zcat
if [ -n "$ZCAT_PATH" ]; then
    zcat=$ZCAT_PATH"/"
else
    zcat=/bin/zcat
fi

# gawk
if [ -n "$GAWK_PATH" ]; then
    gawk=$GAWK_PATH"/"
else
    gawk=/usr/bin/gawk
fi

# split
if [ -n "$SPLIT_PATH" ]; then
    split=$SPLIT_PATH"/"
else
    split=/usr/bin/split
fi

## input and output folders and file prefixes
#############################################

# reference folder for indexes, fragment lists etc
if [ -n "$HOME_PATH" ]; then
    ref_path=$HOME_PATH"/"
else
    ref_path=REFBOWTIE
fi

if [ -n "$HOME_PATH" ]; then
    ref_path=$HOME_PATH"/"
else
ref_path2=REFERENCEFOLDER
fi

# reference source folder for source etc
if [ -n "$HOME_PATH" ]; then
    ref_pathC=$HOME_PATH"/"
else
    ref_pathC=REFERENCEFOLDEC
fi

# home folder for tools
if [ -n "$HOME_PATH" ]; then
    home_path=$HOME_PATH"/"
else
    home_path=HOMEFOLDER
fi

ALLSPEC="ALLELE"
LANE="FLOWCELL"
STRING="DITAG"
SAVE="STARTPATHL"
SAVE2="STARTPATHR"
OLIGOS="BAITS"
ENZYME="CUTTER"
REF="REFERENCE"
CORE=NODES
TRIM=TRIMLENGTH
CHREND=CHROMEND

exec_cmd()
{
    echo $*
    if [ -z "$DRY_RUN" ]; then
        eval "$@" ##|| die 'Error'
    fi
}

# pruning
if [ -n "$PRUNE" ]; then
    prune=$PRUNE
else
    prune=$ref_pathC"prune.pl"
fi

    flash=$ref_pathC"gflash"
    trimmer=$ref_pathC"trim_galore --trim1 --phred33 --paired -q 20 -a AGATCGGAAGAGC -a2 AGATCGGAAGAGC --length 10 --stringency 1" 
#############################################

if [ -n "$GENOME_REF_FASTA" ]; then
    genome_ref_fasta=$GENOME_REF_FASTA
else
    genome_ref_fasta=$ref_path$REF
fi

file=$LANE
trim=$TRIM

mkdir -p $home_path$LANE
mkdir -p $home_path$LANE/BOWTIE/
mkdir -p $home_path$LANE/BOWTIE/$TRIM
mkdir -p $home_path$LANE/data/

if [ -n "$SHIFTREV_PATH" ]; then
shiftrev=$SHIFTREV_PATH"/"
else
shiftrev=$ref_pathC"/"shiftrev
fi

# scratch folder for temporary files
if [ -n "$WORK_PATH" ]; then
    work_path=$WORK_PATH"/"
else
    work_path=$home_path$LANE/BOWTIE/$TRIM/
fi

# output folder

if [ -n "$OUT_PATH" ]; then
   out_path=$OUT_PATH"/"
else
   out_path=$home_path$LANE/data/
#  echo "No output folder specified; quitting."
#    exit
fi

## other run parameters & files
###############################

# input folder
if [ -n "$IN_PATH" ]; then
    in_path=$IN_PATH"/"
else

    in_path=$SAVE
#    echo "No input folder specified; quitting."
#    exit
fi

# input folder
if [ -n "$IN_PATH" ]; then
    in_path2=$IN_PATH"/"
else

    in_path2=$SAVE2
#    echo "No input folder specified; quitting."
#    exit
fi

# input cell_line/replicate name prefix
if [ -n "$NAME" ]; then
    name=$NAME"/"
else
    name=$LANE
#    echo "No input file prefix specified; quitting."
#    exit
fi

# print all input parameters:

echo -e "Input parameters:"
echo -e "\tHOME_PATH = < $home_path >"
echo -e "\tWORK_PATH = < $work_path >"
echo -e "\tOUT_PATH = < $out_path >"
echo -e "\tIN_PATH = < $in_path >"
echo -e "\tIN_PATH = < $in_path2 >"
echo -e "\tGENOME_REF_FASTA = < $genome_ref_fasta >"
echo -e ""

echo -e "Input parameters:" >> $work_path"HiCinC.txt"
echo -e "\tHOME_PATH = < $home_path >"  >> $work_path"HiCinC.txt"
echo -e "\tSCRATCH = < $work_path >"  >> $work_path"HiCinC.txt"
echo -e "\tOUT_PATH = < $out_path >"  >> $work_path"HiCinC.txt"
echo -e "\tIN_PATH = < $in_path >"  >> $work_path"HiCinC.txt"
echo -e "\tIN_PATH = < $in_path2 >"  >> $work_path"HiCinC.txt"
echo -e "\tGENOME_REF_FASTA = < $genome_ref_fasta>"  >> $work_path"HiCinC.txt"
echo -e ""  >> $work_path"HiCinC.txt"

#####################################################

######################
## START PROCESSING ##
######################
echo -e "Starting the data processing."

var=0

for path in $in_path $in_path2


do

var=$((var+1))
#####################################################
## STEP 0: .gz to .fastq conversion (if required) ###
#####################################################

################################
## STEP 1: alignment          ##
################################
FILES=$(find $path -name "*"$var".fastq*"|sort)

if [ "$FILES" == "" ]; then 
exit
fi

TOTAL_LINES=""
for FILE in $FILES
do

if [ ${FILE: -3} == ".gz" ]; then
           echo -e "\t....$FILE "
cat $FILE >> $work_path"/"$LANE"_"$var".fastq.gz"
TOTAL_LINES=$((TOTAL_LINES+$(zcat $FILE|wc -l)))
else
           echo -e "\t....$FILE >"
cat $FILE >> $work_path"/"$LANE"_"$var".fastq"
TOTAL_LINES=$((TOTAL_LINES+$(wc -l < $FILE)))
fi

done

echo $TOTAL_LINES

NUM=$(expr '(' $TOTAL_LINES / $CORE ')' )
SPLIT_NUM=$(expr '(' $NUM % 4 ')')
NUM2=$(expr '(' $NUM - $SPLIT_NUM ')')
NUM2=$((NUM2+4))

echo $NUM2
echo "$split -d -l $NUM2 $LANE $work_path $LANE $var _split_"

if [ ${FILE: -3} == ".gz" ]; then
zcat $work_path"/"$LANE"_"$var".fastq.gz"|$split -d -l $NUM2 - $work_path"/"$LANE"_"$var"_split_"
else
$split -d -l $NUM2 $work_path"/"$LANE"_"$var".fastq" $work_path"/"$LANE"_"$var"_split_"
fi

done

#cat file.fastq | paste - - - - | sort -T $work_path -k1,1 -t " " | tr "\t" "\n" > file_sort -T $work_pathed.fastq
parallel 'var=`echo {}`;var2=`echo {}|sed 's/_1_split/_2_split/'`; var3=`echo {}|sed 's/_1_split//'`; '" $trimmer "' $var $var2 -o '"$work_path"' ' ::: $work_path*_1_split_*

parallel 'str=`echo '"$STRING"'`; var=`echo {}`;var2=`echo {}|sed 's/_1_split/_2_split/'|sed 's/_val_1/_val_2/'`; var3=`echo {}|sed 's/_1_split//'|sed 's/_val_1.fq//'|sed 's/-/_/g'`; '" $flash "' $var $var2 $str -o $var3' ::: $work_path*_1_split_*_val_1.fq

rm $work_path*split*fq

parallel ' file3=`echo {}`;file2=`echo {}|sed 's/\.fastq//g'`; bowtie2 -p 1 --very-sensitive -x '"$genome_ref_fasta"' -U $file3 -S $file2.bowtiesam >> '"$work_path"'bowtie_out.txt' :::  $work_path*.gflashed.fastq

#rm $work_path*flashed.fastq 
rm $work_path*fq

parallel ' file=`echo {}`; echo "END OF FILE" >> $file ' :::  $work_path*.gflashed.bowtiesam

parallel ' file=`echo {}`; file2=`echo {}|sed 's/\.bowtiesam//g'`; '"$shiftrev"' $file $file2 '"$ALLSPEC"'' :::  $work_path*.gflashed.bowtiesam

fmapped=`cat $work_path*bowtiesam|grep "PE12:"|sed 's/:PE/\t/'|gawk '{if($3==0||$3==16){print $1}}'|sort -T $work_path|uniq|wc -l`
ftotal=`cat $work_path*bowtiesam|grep "PE12:"|sed 's/:PE/\t/'|gawk '{print $1}'|sort -T $work_path|uniq|wc -l`
nf1mapped=`cat $work_path*bowtiesam|grep "PE1:"|sed 's/:PE/\t/'|gawk '{if($3==0||$3==16){print $1}}'|sort -T $work_path|uniq|wc -l`
nf1total=`cat $work_path*bowtiesam|grep "PE1:"|sed 's/:PE/\t/'|gawk '{print $1}'|sort -T $work_path|uniq|wc -l`
nf2mapped=`cat $work_path*bowtiesam|grep "PE2:"|sed 's/:PE/\t/'|gawk '{if($3==0||$3==16){print $1}}'|sort -T $work_path|uniq|wc -l`
nf2total=`cat $work_path*bowtiesam|grep "PE2:"|sed 's/:PE/\t/'|gawk '{print $1}'|sort -T $work_path|uniq|wc -l`

tot=`cat $work_path*_1_split*trimming_rep*|grep "Total reads processed"|sed 's/\,//g'|gawk '{a=a+$4;print a}'|tail -n 1`
#trim1=`cat $work_path*_1_*trimming_rep*|grep "Reads with adapters"|sed 's/\,//g'|gawk '{a=a+$4;print a}'|tail -n 1`
mapped=$(( fmapped+nf1mapped))
total=$(( ftotal+nf1total))
trim2=$(( total-mapped ))

echo -e "\n Raw reads: $NAME \t" >> $out_path$LANE.COUNTS.txt
echo $total >> $out_path$LANE.COUNTS.txt
echo -e "\n Mapped reads: $NAME \t" >> $out_path$LANE.COUNTS.txt
echo $mapped >> $out_path$LANE.COUNTS.txt

echo -e "total\t$tot"  >> $work_path/$LANE.read1.mapstat
echo -e "mapped\t$total"  >> $work_path/$LANE.read1.mapstat
echo -e "local\t$mapped"  >> $work_path/$LANE.read1.mapstat
echo -e "global\t$trim2"  >> $work_path/$LANE.read1.mapstat

#trim1b=`cat $work_path*_2_*trimming_rep*|grep "Reads with adapters"|sed 's/\,//g'|gawk '{a=a+$4;print a}'|tail -n 1`

mapped=$(( fmapped+nf2mapped))
total=$(( ftotal+nf2total))
trim2b=$(( total-mapped ))

echo -e "total\t$tot"  >> $work_path/$LANE.read2.mapstat
echo -e "mapped\t$total"  >> $work_path/$LANE.read2.mapstat
echo -e "local\t$mapped"  >> $work_path/$LANE.read2.mapstat
echo -e "global\t$trim2b"  >> $work_path/$LANE.read2.mapstat

ndigested=`cat $work_path*bowtiesam|grep "PE12:"|sed 's/:PE/\t/'|gawk '{if($3==0||$3==16){print $1}}'|sort -T $work_path|uniq -c|gawk '{if($1==1){print}}'|wc -l`
fmapped=`cat $work_path*bowtiesam|grep "PE12:"|sed 's/:PE/\t/'|gawk '{if($3==0||$3==16){print $1}}'|sort -T $work_path|uniq|wc -l`
ftotal=`cat $work_path*bowtiesam|grep "PE12:"|sed 's/:PE/\t/'|gawk '{print $1}'|sort -T $work_path|uniq|wc -l`
digested=$((fmapped-ndigested))

echo -e "ftotal\t$ftotal"  >> $work_path/$LANE.cflashed.mapstat
echo -e "fmapped\t$fmapped"  >> $work_path/$LANE.cflashed.mapstat
echo -e "local\t$digested"  >> $work_path/$LANE.cflashed.mapstat
echo -e "global\t$ndigested"  >> $work_path/$LANE.cflashed.mapstat

ndigested1=`cat $work_path*bowtiesam|grep "PE1:"|sed 's/:PE/\t/'|gawk '{if($3==0||$3==16){print $1}}'|sort -T $work_path|uniq -c|gawk '{if($1==1){print}}'|wc -l`
ndigested2=`cat $work_path*bowtiesam|grep "PE2:"|sed 's/:PE/\t/'|gawk '{if($3==0||$3==16){print $1}}'|sort -T $work_path|uniq -c|gawk '{if($1==1){print}}'|wc -l`

nf1mapped=`cat $work_path*bowtiesam|grep "PE1:"|sed 's/:PE/\t/'|gawk '{if($3==0||$3==16){print $1}}'|sort -T $work_path|uniq|wc -l`
nf1total=`cat $work_path*bowtiesam|grep "PE1:"|sed 's/:PE/\t/'|gawk '{print $1}'|sort -T $work_path|uniq|wc -l`
nf2mapped=`cat $work_path*bowtiesam|grep "PE2:"|sed 's/:PE/\t/'|gawk '{if($3==0||$3==16){print $1}}'|sort -T $work_path|uniq|wc -l`
nf2total=`cat $work_path*bowtiesam|grep "PE2:"|sed 's/:PE/\t/'|gawk '{print $1}'|sort -T $work_path|uniq|wc -l`
digested=$((ndigested1+ndigested2))
nftotal=$((nf1total+nf2total))
nfmapped=$((nf1mapped+nf2mapped))
ndigested=$((nfmapped-digested))

echo -e "ftotal\t$nftotal"  >> $work_path/$LANE.nflashed.mapstat
echo -e "fmapped\t$nfmapped"  >> $work_path/$LANE.nflashed.mapstat
echo -e "local\t$digested"  >> $work_path/$LANE.nflashed.mapstat
echo -e "global\t$ndigested"  >> $work_path/$LANE.nflashed.mapstat

cmd="R CMD BATCH --no-save --no-restore \"--args picDir='$out_path' bwtDir='$work_path/' sampleName='$LANE' rmMulti='1' rmSingle='1'\" $ref_pathC/plot_mapping_portion.R $out_path/plot_mapping_portion.Rout"

exec_cmd " $cmd"

cmd="R CMD BATCH --no-save --no-restore \"--args picDir='$out_path' bwtDir='$work_path/' sampleName='$LANE' rmMulti='1' rmSingle='1'\" $ref_pathC/plot_mapping_portion2.R $out_path/plot_mapping_portion2.Rout"

exec_cmd " $cmd"

#######################################################
## STEP 2:  (pair-matching and de-duplication)       ##
#######################################################

file=$LANE

for filename in $work_path"/"*.gflashed.sam; do
cat $filename|grep -v '^[[:space:]]*@' >> $work_path$file.matefixed.sam
done

################################################################################################################################################

echo -e "\n Sorted&Matefixed Reads \t" >> $out_path$file.COUNTS.txt
wc -l $work_path$file.matefixed.sam  >> $out_path$file.COUNTS.txt

##############################################################################################################################################

cat $work_path$file.matefixed.sam | gawk '{if($9>29&&$12>29){print $0}}'>$work_path$file.matefixed.MQ30.sam

echo "QC reads (MQ>=30)" >> $out_path$file.COUNTS.txt

wc -l $work_path$file.matefixed.MQ30.sam >> $out_path$file.COUNTS.txt

########################################################################################################
list=`cat $work_path$file.matefixed.MQ30.sam |gawk '{print $2}'|sort -V|uniq`

for i in $list;
do

cat $work_path$file.matefixed.MQ30.sam |gawk '{if($2=="'$i'"){print $0  }}'  >  $out_path$file.pairs.$i.sam

counter=$(wc -l < $out_path$file.pairs.$i.sam)

cat $out_path$file.pairs.$i.sam $out_path$file.pairs.$i.sam |gawk '{if(NR<='"$counter"'){split($16,a,":PE"); b=a[1]; cb[b]=cb[b]+1; cont[b]=cont[b]"_"_$2"_"$3"_"$6"_"$7;}else{split($16,a,":PE"); b=a[1]; flag[cont[b]]++; if(flag[cont[b]]<=cb[b]){print $0 }}}'|gawk '{for(i=1;i<=NF;i++){printf("%s\t",$i);}printf("\n");print $5 "\t" $6 "\t" $7 "\t" $8 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $12 "\t" $13 "\t" $14 "\t" $9 "\t" $10 "\t" $11 "\t" $16 "\t" $15;}' > $out_path$file.pairs.$i.dedup.sam

done

echo "Deduplicated Valid Paired reads " >> $out_path$file.COUNTS.txt
wc -l $out_path$file.pairs.chr_*.dedup.sam >> $out_path$file.COUNTS.txt

##########################################################################################################

pairtotal1=`cat $work_path$file.matefixed.sam |wc -l`

Unmapped_pairs=`cat $work_path*unmapped.sam|grep -v chr|wc -l`

Pairs_with_singleton=`cat $work_path*unmapped.sam|grep chr|wc -l`

Total_pairs_processed=$((pairtotal1 +Pairs_with_singleton+ Unmapped_pairs))

Low_qual_pairs=`cat $work_path$file.matefixed.sam | gawk '{if($9<30||$12<30){print $0}}'|wc -l`

dedupl=`cat $out_path$file.pairs.chr_*.dedup.sam|wc -l`

Multiple_pairs_alignments=$((pairtotal1-Low_qual_pairs-dedupl/2))

Not_Reported_pairs=$((Multiple_pairs_alignments+Low_qual_pairs+Pairs_with_singleton))

Reported_pairs=$((pairtotal1-Multiple_pairs_alignments-Low_qual_pairs))

frac1=`gawk 'BEGIN{print '"$Unmapped_pairs"'/'"$Total_pairs_processed"*100'}'`
frac2=`gawk 'BEGIN{print '"$Low_qual_pairs"'/'"$Total_pairs_processed"*100'}'`
frac3=`gawk 'BEGIN{print '"$Reported_pairs"'/'"$Total_pairs_processed"*100'}'`
frac4=`gawk 'BEGIN{print '"$Multiple_pairs_alignments"'/'"$Total_pairs_processed"*100'}'`
frac5=`gawk 'BEGIN{print '"$Pairs_with_singleton"'/'"$Total_pairs_processed"*100'}'`

echo -e "Total_pairs_processed\t$Total_pairs_processed\t100" >> $work_path/$LANE.pairs.pairstat
echo -e "Unmapped_pairs\t$Unmapped_pairs\t$frac1"  >> $work_path/$LANE.pairs.pairstat
echo -e "Low_qual_pairs\t$Low_qual_pairs\t$frac2"  >> $work_path/$LANE.pairs.pairstat
echo -e "Unique_paired_alignments\t$Reported_pairs\t$frac3" >> $work_path/$LANE.pairs.pairstat
echo -e "Duplicate_pairs\t$Multiple_pairs_alignments\t$frac4"  >> $work_path/$LANE.pairs.pairstat
echo -e "Pairs_with_singleton\t$Pairs_with_singleton\t$frac5"  >> $work_path/$LANE.pairs.pairstat
echo -e "Reported_pairs\t$Reported_pairs\t$frac3"  >> $work_path/$LANE.pairs.pairstat

cmd="R CMD BATCH --no-save --no-restore \"--args picDir='$out_path' bwtDir='$work_path/' sampleName='$LANE' rmMulti='1' rmSingle='1'\" $ref_pathC/plot_pairing_portion.R $out_path/plot_pairing_portion.Rout"
exec_cmd " $cmd"

*/
