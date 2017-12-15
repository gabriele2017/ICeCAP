/*    ###################################################################
     ##                                                                   ##
     ##       #####################  ICeCap ##########################    ##
     ##       ##                   Hi-Resolution-Hi-C @ICR          ##    ##
     ##       ##                                                    ##    ##
     ##       ##          Algorithm by Gabriele Migliorini,         ##    ##
     ##       ########################################################    ##
     ##                                                                   ##
     ################################################################## */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h> 

#define PACKAGE "ICeCap"
#define VERSION "0.01"

#define CHROMOSOMES 25
#define WIN_SHIFT1 1
#define RESOLUTION 1000

#define maxim(a,b)(((a)>(b)) ? (a) : (b))
#define minim(a,b)(((a)<(b)) ? (a) : (b))

double expit(double x){ return(1/(1+exp(-x)));}

int rf[CHROMOSOMES],reg[CHROMOSOMES];
int **enzyme,**enriched;
int *store,*store2,*store2b,*store3,*num_cells;
int *map_count;
int *avg_norma,*avg_norma2;
int WIN_SIZE;
int chrom,chrom_start,chrom_end,GRID;
double *weight,*weight2;
double *avg,*avg2,*dilute,*logistic,*map,*norma;
double **hic;

int CUT;
char *ENZYME; 
char JUNCTION[20];
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

double digest();
int grid();

int main(int argc, char *argv[])
    {
      int option,i,j,ii,jj,k,temp,a,chrom;
      int BIN,fflag;
      int max_iter,corenodes,REwidth,resolution,window,iter;
      double *genome;
      int tot,tot2,tot3,totfragments[25],upstore[25],downstore[25],nchromosomes;
      double PVAL,FDR,map_thr,weight_threshold,map_threshold,trimming;
      double tolerance,conv_check;
      char *file,*file2,*file3,*file4,*filenorma,*filelogs,*fileref,*fileref2,*fileref3,*fileref4;
      char *lane,*folderN,*folder,*tmp,*folderref,*folderrefmap,*folderfrag,*folderrefC,*folderwrite,*BAITS,*BAITFLAG,*PATH,*PATH2;
      char *ALLELE,*PWEIGHT,*PATHBOWTIE,*verbose,*stats,*outprefix,*path;
      char *outdir,*outdir2,*indir,*REFERENCE,*folder1,*TRANS;
      char *fileoutput1,*fileoutput2,*fileoutput3,*fileoutput4;

      char QNAME[100],CIGAR[100],SEQ[100],QSEQ[100];
      char command[500],heatmap[500],graphics[500],end[500];
      char chromosome[5],chromosome2[5],iteration[5];
      char stamp[100],stamp2[100],schr[10],suffix[100];
      int FLAG,POS1,POS2,QUAL,d3,BOH,CHR1,CHR2,TRIM;
      int position,position2,position3,region,SCALE;
      int area1,area2;
      long int tot4; 
      double d1,d2,d4,d5;
      double mappability;
      char enzyme_seq;
      FILE *fmap;
      FILE *fout3;
      FILE *fout3b,*fout3c,*fout3d,*fout3e;

      folder="";
      indir="";
      outdir="";
      outdir2="";
      folderfrag="";
      folderwrite="";
      folderrefmap="";
      folderN="yes";
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
      REwidth=800;
      stats="none";
      verbose="no";
      CUT=1;
      ENZYME="AAGCTT";
      TRIM=100;
      SCALE=1;
      GRID=0;
      PATH="PATH";
      PATH2="PATH2";
      PATHBOWTIE="PATHBOWTIE";
      trimming=(double)(0.01);
      chrom_start=1;
      chrom_end=25;
      TRANS="no";
      REFERENCE="hg19";

      if(argc ==1){
	fprintf(stderr, "This program needs arguments .... \n\n");
	print_usage(); exit(0);}

      while((option = getopt(argc, argv, "E:C:l:f:F:D:hVv:S:s:e:r:w:m:c:t:M:J:W:G:B:R:L:P:b:N:Z:A:T:Q:H:U:p:I:O:X:o:")) != -1)

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
        case 'w' :           weight_threshold=atof(optarg); break;
        case 'm' :           map_threshold=atof(optarg);    break;
        case 'c' : 	     trimming=atof(optarg);         break;
        case 't' :           tolerance=atof(optarg);        break;
        case 'M' :           max_iter=atoi(optarg);         break;
        case 'J' :           corenodes=atoi(optarg);        break;
        case 'W' :           REwidth=atoi(optarg);          break;
        case 'G' :           GRID=atoi(optarg);             break;
        case 'B' :           BAITFLAG=optarg;               break;
        case 'R' :           PATH=optarg;                   break;
        case 'L' :           PATH2=optarg;                   break;
        case 'P' :           PATHBOWTIE=optarg;             break;
        case 'b' :           verbose=(optarg);              break; 
	case 'N' :           folderN="";                    break;
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

sprintf(stamp,"");
/*_HiCinC_%s_Grid_%d_Zoom_%d_Chrom_%d_%d",ENZYME,GRID,SCALE,chrom_start,chrom_end);*/
sprintf(stamp2,"_%s_%d_%d",outprefix,chrom_start,chrom_end);
   
      if(strlen(folder)==0||dirExists(folder)==0){printf("\n\n Please provide mandatory -f field: \n\n Refer to the Manual.Thank you. \n\n "); exit(0);}
      sprintf(end,"%c\0",folder[strlen(folder)-1]);
      if(strcasecmp(end,"/")!=0){printf("Adding backslash at the end of the folder path..\n");folder=concat(folder,"/");}

      if(strlen(PATHBOWTIE)==0||dirExists(PATHBOWTIE)==0){printf("\n\n Please provide mandatory -P field: \n\n Refer to the Manual.Thank you. \n\n "); exit(0);}
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
                            if(strlen(folderN)==0){
                            outdir2=concat(concat(folderwrite,"/"),lane);
                            mkdir(concat(concat(concat(concat(folderwrite,lane),stamp),"/HiC_contacts_chr"),stamp2),0777);
                            outdir=concat(concat(concat(concat(concat(concat(folderwrite,lane),"/HiC_contacts_chr"),stamp2),"/"),lane),"");
                            }}else{
                            mkdir(outdir,0777);
                            if(strlen(folderN)==0){mkdir(concat(concat(outdir,"/HiC_contacts_chr"),stamp2),0777);
                            outdir2=concat(concat(outdir,"/"),lane); 
                            outdir=concat(concat(concat(concat(concat(outdir,"/HiC_contacts_chr"),stamp2),"/"),lane),""); 
                            }}

WIN_SIZE=strlen(ENZYME);
if(CUT*2>strlen(ENZYME)){CUT=strlen(ENZYME)-CUT;
printf("\n You provided a -C parameter that needs to be complemented: %d should read %d .\n",strlen(ENZYME)-CUT,CUT);}

char str1[WIN_SIZE-CUT],str2[WIN_SIZE-CUT];

strncpy(str1,ENZYME,WIN_SIZE-CUT);
strncpy(str2,ENZYME,WIN_SIZE-CUT);
str1[WIN_SIZE-CUT]='\0';

strrev(str1);
for (i=0;i<strlen(str1);i++){
if(str1[i]=='A'){str1[i]='T';}
else if(str1[i]=='T'){str1[i]='A';}
else if(str1[i]=='C'){str1[i]='G';}
else if(str1[i]=='G'){str1[i]='C';}
}

sprintf(JUNCTION,"%s",concat(concat(str2,str1),""));

mkdir(folderfrag,0777);
folderref=folderfrag;

BAITS="ALLFRAGMENTS_CHR";
if(strcasecmp(BAITFLAG,"no")!=0){BAITS="BAITS_CHR";}

if(strlen(folderN)!=0){
printf(" %s genome: \n ",REFERENCE);
fileoutput4=concat(concat(PATHBOWTIE,REFERENCE),".chrom.sizes");

    if(fopen(fileoutput4, "rt" )==NULL){

sprintf(command,"wget -q -P %s \'http://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips/%s.chrom.sizes\'",PATHBOWTIE,REFERENCE,REFERENCE);

system(command);

}

sprintf(command,"cat %s/%s.chrom.sizes|grep -v \"_\"|wc -l > %s/chromosomes",PATHBOWTIE,REFERENCE,PATHBOWTIE);

system(command);

FILE *fchr = fopen (concat(PATHBOWTIE,"chromosomes"), "r");

chrom_start=1;

while (!feof(fchr)) {fscanf(fchr,"%d",&chrom_end);}
sleep(1);

printf("\n %d Chromosomes will be considered \n",chrom_end);
folderN=folder;
sprintf(command,"rm -r -f %s",folderref);
system(command);

mkdir(folderref,0777);
}else{
sprintf(command,"cat %s/%s.chrom.sizes|grep -v \"_\"|grep -v \"M\"|wc -l > %s/chromosomes",PATHBOWTIE,REFERENCE,PATHBOWTIE);

system(command);

sprintf(command,"cat %s/%s.chrom.sizes|grep -v \"_\"|wc -l > %s/chromosomes",PATHBOWTIE,REFERENCE,PATHBOWTIE);

system(command);

FILE *fchr = fopen (concat(PATHBOWTIE,"chromosomes"), "r");

chrom_start=1;

while (!feof(fchr)) {fscanf(fchr,"%d",&chrom_end);}

folderref=concat(concat(concat(concat(concat(folder,lane),stamp),"/reference_HiC"),stamp2),"/");

      mkdir(concat(concat(folder,lane),stamp),0777);
      mkdir(concat(concat(concat(folder,lane),stamp),"/logs"),0777);
      mkdir(folderref,0777);
      filelogs=concat(concat(concat(concat(concat(folder,lane),stamp),"/logs/log_"),stamp2),".log");
      printf("%s\n%s\n",concat(concat(concat(concat(concat(folder,lane),stamp),"/logs/log_"),stamp2),".log"),filelogs);
}

if(strcasecmp(BAITFLAG,"no")==0){}else{ 
for(chrom=chrom_start;chrom<=chrom_end;chrom++){
sprintf(command,"gawk '{ if($1==\"chr%d\") {print $1 \"\t\" $2 \"\t\" $3 \"\t\" $4} }' %s > %sBAITS_CHR%d",chrom,BAITFLAG,folderref,chrom);
/*printf("\n Preparing folders \n\n %d %s\n",chrom,command);*/
system(command);
}}

      fileref=concat(folderref,"/ALLFRAGMENTS_CHR");
      fileref2=concat(folderref,"/FRAGMENTS_CHR");
      fileref3=concat(folderref,"/chr");
      fileref4=concat(concat(concat(folderref,"/"),ENZYME),"");



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
	printf("Map threshold   : %lf\n",map_threshold);         printf("Verbose           : %s\n",verbose);                  
        printf("Weight threshold: %lf\n",weight_threshold);	printf("Tolerance         : %lf\n",tolerance);
	printf("Max number of it: %d\n",max_iter);  	        printf("Minimum Resolution: %d\n",resolution);
	printf("JUNCTION        : %s\n",JUNCTION);              printf("Output Prefix     : %s\n",outprefix);
        printf("Zooming factor  : %d\n",SCALE);                 printf("Trimming          : %lf\n",trimming);
        printf("Statistics      : %s\n",stats);                 printf("CONFIG:           : %s\n",folderN);
        printf("READ LENGTH     : %d\n",TRIM);                  printf("ALLELE:           : %s\n",ALLELE);
        printf("NUMBER OF CORES : %d\n",corenodes);             printf("NGS OUTPUT        : %s\n",folderN);
        printf("INTER-CHROMOSOME ANALYSIS     : %s\n",TRANS);

        printf("Input Sam : %s\n",indir);             printf("Output : %s\n",outdir);
        
        if(GRID!=0){printf("\nGrid type         : Uniform bins of size %d\n",SCALE*resolution);}
	if(GRID==0&&SCALE!=1){printf("\nGrid type      : Restriction metafragments of size: %d times %s average cut frequency\n",SCALE,ENZYME);}
	if(GRID==0&&SCALE==1){printf("\nGrid type      : Single fragment resolution: Restriction fragment used is: %s\n",ENZYME);}
        printf("\n\n\n\n\n .. Reference Genome from UCSC. Enumerating baited/non-baited fragments.Please be patient..\n\n\n");


	/*      CODE   */
enzyme = malloc (sizeof *enzyme * (chrom_end-chrom_start));
enriched = malloc (sizeof *enriched * (chrom_end-chrom_start));

if (enzyme == NULL||enriched ==NULL)
{
puts("\nFailure to allocate room for pointers");
exit(0);
}

tot=0; nchromosomes=0; tot2=0; tot3=0;
totfragments[0]=0;

if(strlen(folderN)!=0){fout3e = fopen(concat(fileref4,""),"w+");}

for(chrom=chrom_start;chrom<=chrom_end;chrom++){

digest(chrom,chrom_end,rf,enzyme,ENZYME,CUT,SCALE,GRID,PATHBOWTIE,REFERENCE,resolution); nchromosomes++; tot+=(rf[chrom]-1);  totfragments[nchromosomes]=tot;
if((rf[chrom]-1)>tot3){tot3=rf[chrom]-1; BIN=rint((enzyme[chrom][rf[chrom]-1])/(rf[chrom]-1))+1; }
printf("Chromosomes analyzed: %d\tFragments on Chromosome %d: %d\t Total number of Fragments:%d\n",nchromosomes,chrom,rf[chrom],totfragments[nchromosomes]);
sprintf(chromosome,"%d",chrom);

fout3b = fopen(concat(fileref,chromosome),"w+");
fout3c = fopen(concat(fileref2,chromosome),"w+");
fout3d = fopen(concat(fileref3,chromosome),"w+");

for(i=1;i<=rf[chrom]-1;i++){
fprintf(fout3b,"chr%d\t%d\t%d\n",chrom,enzyme[chrom][i],enzyme[chrom][i+1]);
if(strlen(folderN)!=0){fprintf(fout3e,"chr%d\t %d\t%d\t%d\t0\t+\n",chrom,enzyme[chrom][i],enzyme[chrom][i]+WIN_SIZE,enzyme[chrom][i]+CUT);}
if(strcasecmp(BAITFLAG,"no")==0){
tmp=concat(fileref2,chromosome);
fprintf(fout3c,"chr%d\t%d\t%d\t%d\t%d\n",chrom,enzyme[chrom][i],enzyme[chrom][i+1],(enzyme[chrom][i]+enzyme[chrom][i+1])/2,region);}
}

    fprintf(fout3d,"chr%d\t%d\t%d\n",chrom,1,enzyme[chrom][rf[chrom]-1]);
    fclose(fout3b);
    fclose(fout3d);

for(i=1;i<=rf[chrom];i++){  enriched[chrom][i]=0; }
FILE *fchic = fopen(concat(concat(folderref,BAITS),chromosome), "r"); 
if(fchic==NULL){printf("BAIT file is empty or missing\n"); exit(0);}
    while (!feof(fchic)){fscanf(fchic,"%s %d %d %d %d",&chromosome,&position,&position2,&region); 
    for(i=grid(chrom,(int)(position+1),rf,enzyme);i<=grid(chrom,(int)(position2-1),rf,enzyme);i++){
    if(enriched[chrom][i]==0){ tot2++; enriched[chrom][i]= tot2; 
    fprintf(fout3c,"chr%d\t%d\t%d\t%d\t%d\n",chrom,enzyme[chrom][i],enzyme[chrom][i+1],(enzyme[chrom][i]+enzyme[chrom][i+1])/2,region);}
     }
    }

fclose(fchic);
fclose(fout3c);

printf("Reading baits position on chromosome %d\t Enriched fragments: %d\n",chrom,(int)(tot2));
 
 }

if(strlen(folderN)!=0){fclose(fout3e);}

      if(strlen(folderN)==0){

  fout3 = fopen(filelogs,"a+");

      mkdir(concat(concat(concat(concat(folder,lane),stamp),"/weights"),stamp2),0777);
      mkdir(concat(concat(concat(concat(folder,lane),stamp),"/stats"),stamp2),0777);
      mkdir(concat(concat(concat(concat(concat(folder,lane),stamp),"/stats"),stamp2),"/fdr"),0777);
      mkdir(concat(concat(concat(folder,lane),stamp),"/graphics"),0777);

      fprintf(fout3,"\n\n\n\n\n\n\n\n\n");	      
      fprintf(fout3,"       #################### %s %s ######################\n",PACKAGE,VERSION);
      fprintf(fout3,"##                       Hi-Res-Hi-C                            ##\n");
      fprintf(fout3,"##                                                              ##\n");
      fprintf(fout3,"##          Algorithm by Gabriele Migliorini,                   ##\n");
      fprintf(fout3,"##                 Code version %s                              ##\n",VERSION);
      fprintf(fout3,"##                              2015                            ##\n");
      fprintf(fout3," ################################################################ \n");
      fprintf(fout3,"                                                                  \n");

      fprintf(fout3,"ENZYME          : %s\n",ENZYME);                fprintf(fout3,"CUT_POINT         : %d\n",CUT);
      fprintf(fout3,"LANE            : %s\n",lane);   	             fprintf(fout3,"Folder            : %s\n",folder);
      fprintf(fout3,"FolderRef       : %s\n",folderref);             fprintf(fout3,"BAIT FILE         : %s\n",BAITFLAG);
      fprintf(fout3,"ChromosomeStart : %d\n",chrom_start);           fprintf(fout3,"ChromosomeEnd     : %d\n",chrom_end);
      fprintf(fout3,"Map threshold   : %lf\n",map_threshold);         fprintf(fout3,"Verbose           : %s\n",verbose);                  
      fprintf(fout3,"Weight threshold: %lf\n",weight_threshold);	     fprintf(fout3,"Tolerance         : %lf\n",tolerance);
      fprintf(fout3,"Max number of it: %d\n",max_iter);  	     fprintf(fout3,"Minimum Resolution: %d\n",resolution);
      fprintf(fout3,"JUNCTION        : %s\n",JUNCTION);              fprintf(fout3,"Output Prefix     : %s\n",outprefix);
      fprintf(fout3,"Zooming factor  : %d\n",SCALE);                 fprintf(fout3,"Trimming          : %lf\n",trimming);
      fprintf(fout3,"Statistics      : %s\n",stats);                 fprintf(fout3,"NGS step:         : %s\n",folderN);
      fprintf(fout3,"READ LENGTH     : %d\n",TRIM);                  fprintf(fout3,"ALLELE:           : %s\n",ALLELE);
      fprintf(fout3,"NUMBER OF CORES : %d\n",corenodes);             fprintf(fout3,"NGS OUTPUT        : %s\n",folderN);
      fprintf(fout3,"INTER-CHROMOSOME ANALYSIS     : %s\n",TRANS);

      fprintf(fout3,"Input Sam : %s\n",indir);             fprintf(fout3,"Output : %s\n",outdir);


      if(GRID!=0){fprintf(fout3,"Grid type         : Uniform bins of size %d\n",SCALE*resolution);}
      if(GRID==0&&SCALE==1){fprintf(fout3,"Grid type   : Restriction fragments %s\n",ENZYME);}
      if(GRID==0&&SCALE!=1){fprintf(fout3,"Grid type   : Restriction metafragments of size: %d times %s average cut frequency\n",SCALE,ENZYME);}


        mkdir(concat(concat(concat(concat(folder,lane),stamp),"/"),concat(concat(lane,stamp2),"_cis")),0777);
        mkdir(concat(concat(concat(concat(folder,lane),stamp),"/"),concat(concat(lane,stamp2),"_self")),0777);
        mkdir(concat(concat(concat(concat(folder,lane),stamp),"/"),concat(concat(lane,stamp2),"_trans")),0777);

	file=concat(concat(concat(concat(concat(folder,lane),stamp),"/weights"),stamp2),"/BAREWEIGHTS");
	file2=concat(concat(concat(concat(concat(folder,lane),stamp),"/weights"),stamp2),"/WEIGHTS");
	file3=concat(concat(concat(concat(concat(folder,lane),stamp),"/weights"),stamp2),"/WEIGHTS_");
        file4=concat(concat(concat(concat(concat(folder,lane),stamp),"/stats"),stamp2),"/NORMALIZATION_");
        filenorma=concat(concat(concat(concat(concat(folder,lane),stamp),"/stats"),stamp2),"/NORMALIZATION");
        
        fileoutput1=concat(concat(concat(concat(folder,lane),stamp),"/"),concat(concat(lane,stamp2),"_"));
        fileoutput2=concat(concat(concat(concat(folder,lane),stamp),"/stats"),stamp2);

  if(GRID!=0){
          fprintf(fout3,"Grid type   : Uniform bins of size %d\n",SCALE*resolution);}
        if(GRID==0&&SCALE==1){
          fprintf(fout3,"Grid type   : Single fragment resolution: Restriction fragment used is: %s\n",ENZYME);}

        if(GRID==0&&SCALE!=1){
          fprintf(fout3,"Grid type   : Restriction metafragments of size: %d times %s average cut frequency\n",SCALE,ENZYME);}

        if(GRID==0){printf("\n\nYou selected a fragment based Hi-C map with a meta-fragment Zooming factor of x:%d\n",SCALE);
          fprintf(fout3,"\n\nYou selected a fragment based Hi-C map with a meta-fragment Zooming factor of x:%d\n",SCALE);}
        if(GRID==1){printf("\n\nYou selected a  uniform-grid  Hi-C map with a multi-kbp-bin Zooming factor of x:%d\n",SCALE);
          fprintf(fout3,"\n\nYou selected a  uniform-grid  Hi-C map with a multi-kbp-bin Zooming factor of x:%d\n",SCALE);}
 

      }



tot4=(long)tot*(long)tot2;

printf("\nYou allocated a %d x %d matrix with %li entries. \n\n\n",tot,tot2,tot4);

printf("\nThe average Bin Size on the chromosomes under study is: %d\n",BIN);

if(strlen(folderN)==0){fprintf(fout3,"\nYou allocated a %d x %d matrix with %li entries. \n\n\n",tot,tot2,tot4);}

 if(strlen(folderN)!=0){ 
   printf("\n You generated input files from chrom %d\t to chrom %d\n",chrom_start,chrom_end);

   sprintf(command,"sed 's|TITLE|%s|g' %sNGS.sh > %s%s.sh",lane,folderrefC,folderN,lane);
   system(command);
   sprintf(command,"sed -i 's/TRIMLENGTH/%d/g' %s%s.sh",TRIM,folderN,lane);
   system(command);
   if(strcasecmp(verbose,"no")!=0){
   sprintf(command,"sed -i 's/rm /#rm /g' %s%s.sh",folderN,lane);}
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
   sprintf(command,"sed -i 's|BAITS|%s|g' %s%s.sh",BAITFLAG,folderN,lane);
   system(command);
   sprintf(command,"sed -i 's|REFERENCE|%s|g' %s%s.sh",REFERENCE,folderN,lane);
   system(command);

   sprintf(command,"chmod a+x %s%s.sh",folderN,lane);
   sleep(30);
   system(command);

sprintf(command,"sed -i 's/TITLE/%s/g' %s%s.sh",lane,folderN,lane);
system(command);

sprintf(command,"sed -i '1i #!/bin/bash ' %s%s.sh",folderN,lane);
system(command);


printf("\n\nrunning the NGS shell script: %s%s.sh. \nCheck the NGS folders for completion.\n\n",folderN,lane);
sprintf(command,"%s%s.sh\n\n",folderN,lane);
printf("%s",command);
system(command);
exit(0); 

}

store=malloc(tot2 * sizeof(int));
store2=malloc(tot2 * sizeof(int));
store2b=malloc(tot2 * sizeof(int));
store3=malloc(tot * sizeof(int));
map_count=malloc(tot * sizeof(int));
num_cells=malloc(tot * sizeof(int));
avg_norma=malloc(tot3 * sizeof(int));
avg_norma2=malloc(tot3 * sizeof(int));
norma=malloc(tot3 * sizeof(double));
avg=malloc(tot3 * sizeof(double));
avg2=malloc(tot3 * sizeof(double));
dilute=malloc(tot3 * sizeof(double));
weight=malloc(tot * sizeof(double));
weight2=malloc(tot * sizeof(double));
logistic=malloc(tot3 * sizeof(double));
map=malloc(tot * sizeof(double));
hic=malloc((tot2+1) * sizeof(double *));

if (store==NULL||store2==NULL||store2b==NULL||store3==NULL||map_count==NULL||num_cells==NULL||avg_norma==NULL||avg_norma2==NULL||norma==NULL||avg==NULL||avg2==NULL||dilute==NULL||weight==NULL||weight2==NULL||map==NULL||hic==NULL)
{
puts("\nFailure to allocate room for pointers");
exit(0); 
}

for(i=1;i<tot2; i++){
  *(store+i)=0;
  *(store2+i)=0;
  *(store2b+i)=0;}

for(i=1;i<tot2+1; i++){
hic[i] = malloc(sizeof *hic[i] * tot); }

for(i=1;i<tot2; i++){
for(j=1;j<=tot;j++){
  *(hic[i]+j)=(double)(0.0); }}

for(i=1;i<tot; i++){
  *(weight+i)=(double)(0.0);
  *(weight2+i)=(double)(0.0);
  *(num_cells+i)=0;
  *(map+i)=0.0;
  *(map_count+i)=0;
  *(store3+i)=0;}

for(i=1;i<tot3; i++){
  *(avg+i)=0.0;
  *(avg2+i)=0.0;
  *(dilute+i)=0.0;
  *(avg_norma+i)=0;
  *(avg_norma2+i)=0;
  *(logistic+i)=0.0;
}

tot2=0; 
for(chrom=chrom_start;chrom<=chrom_end;chrom++){  

for(i=1;i<=rf[chrom];i++){  enriched[chrom][i]=0; }

sprintf(chromosome,"%d",chrom);
fmap = fopen(concat(concat(folderrefmap,"mappability_chr"),chromosome), "r");  
if(fmap==NULL){printf("mappability files is empty or missing\nPlease refer to the manual\n\n\n\n\n"); exit(0);}

while (!feof(fmap)){    fscanf(fmap,"%s %d %lf ",&chromosome2,&position,&mappability);
             temp=totfragments[chrom-chrom_start]+grid(chrom,position,rf,enzyme);
             map[temp]+=mappability; map_count[temp]++;	  
}
fclose(fmap);

downstore[chrom]=tot2+1;
FILE *fchic = fopen (concat(concat(folderref,BAITS),chromosome), "r");

while (!feof(fchic)) {fscanf(fchic,"%s %d %d %d",&chromosome,&position,&position2,&region); 
    /*printf("%s\t%d\t%d\t%d\n",chromosome,position,position2,region);*/

for(i=grid(chrom,(int)(position+1),rf,enzyme);i<=grid(chrom,(int)(position2-1),rf,enzyme);i++){
      if(    enriched[chrom][i]==0){ tot2++; enriched[chrom][i]= tot2;   store[tot2]=enzyme[chrom][i];   store2[tot2]=chrom; store2b[tot2]=region;
      }
    }
  }

reg[chrom]=region;
upstore[chrom]=tot2;
fclose(fchic);

printf("Capture Regions on chromosome %d : %d\n",chrom,reg[chrom]-reg[chrom-1]);
printf("Reading mapability track on chromosome %i\n",chrom);
 
}

for(i=1;i<=tot; i++){
  if(map_count[i]!=0){map[i]=map[i]/map_count[i];}else{map[i]=0.0;}
}

for(chrom=chrom_start;chrom<=chrom_end;chrom++){
  
sprintf(schr,"%d",chrom);    

if(strcasecmp(BAITFLAG,"no")!=0){sprintf(suffix,".sam.bfide_ontarget");}
if(strcasecmp(BAITFLAG,"no")==0){sprintf(suffix,".sam.bfide");}

FILE *fin = fopen (concat(concat(concat(concat(indir,".pairs.chr_"),schr),suffix),"") , "r" );
if(fin==NULL){printf("ditag file is empty or missing\n"); exit(0);}
fprintf(fout3,"\nReading Files:%s\n",concat(concat(concat(concat(indir,".pairs.chr_"),schr),suffix),""));
printf("Files Used:%s\n",concat(concat(concat(concat(indir,".pairs.chr_"),schr),suffix),""));
while (!feof(fin))
{

fscanf(fin,"%s %d %d %d %d %s %d %d %d %s %s %lf %lf %d %lf %lf", &QNAME, &FLAG, &CHR1, &POS1, &QUAL, &CIGAR, &CHR2, &POS2, &BOH,&SEQ,&QSEQ,&d1,&d2,&d3,&d4,&d5);
 
if((CHR1>=chrom_start&&CHR1<=chrom_end)&&(CHR2>=chrom_start&&CHR2<=chrom_end)&&(d1+d2)<REwidth){

  /* printf("%10.20s\t%d\t%5d\t%10d\t%d\t%10.10s\t%5d\t%10d\t%10d\t%100.100s\t%100.100s\n\n",QNAME,FLAG,CHR1,POS1,QUAL,CIGAR,CHR2,POS2,BOH,SEQ,QSEQ);*/

     if(enriched[CHR1][grid(CHR1,POS1,rf,enzyme)]!=0){
     hic[enriched[CHR1][grid(CHR1,POS1,rf,enzyme)]][totfragments[CHR2-chrom_start]+grid(CHR2,POS2,rf,enzyme)]+=1.0; }
 }
 }

fclose(fin);
 }
printf("Data Loaded. Analysis Begins now \n");
fprintf(fout3,"\n\nICE - Iterative normalization of HiC biases \n   iteration number \t  Convergence rate \t User's Tolerance \n\n");
printf("\n\nICeCAP - Iterative normalization of C-Hi-C biases \n  \n");
printf("Iteration number \t Convergence Rate \t Tolerance\n\n");

weights(chrom_start,chrom_end,rf,enzyme,enriched,map_threshold,weight_threshold,totfragments,hic,map,weight2,tot,tot2,file);
weights(chrom_start,chrom_end,rf,enzyme,enriched,map_threshold,weight_threshold,totfragments,hic,map,weight2,tot,tot2,file2);

conv_check=(double)(1.0); iter=0;
filewrite(BIN,resolution,chrom_start,chrom_end,iter,rf,enzyme,enriched,totfragments,hic,tot,tot2,tot3,downstore,upstore,avg,avg2,dilute,logistic,avg_norma,avg_norma2,trimming,stats,TRANS,PVAL,FDR,fileoutput1,fileoutput2,filenorma,folderref,folderrefC);

while(conv_check>tolerance&&iter<max_iter){
iter++; 

for(i=1;i<=tot2; i++){
  if(weight2[totfragments[store2[i]-chrom_start]+grid(store2[i],store[i]+1,rf,enzyme)]!=0&&map[totfragments[store2[i]-chrom_start]+grid(store2[i],store[i]+1,rf,enzyme)]>map_threshold){
for(j=1;j<=tot;j++){
if(weight2[j]!=0&&map[j]>map_threshold){
  *(hic[i]+j)=  *(hic[i]+j)  / (weight2[j]*weight2[totfragments[store2[i]-chrom_start]+grid(store2[i],store[i]+1,rf,enzyme)]);}
 else{  *(hic[i]+j)=(double)(0.0);}
}
}else{ for(k=1;k<=tot;k++){ *(hic[i]+k)=(double)(0.0);}}
} 
 
sprintf(iteration,"%d",iter);
file3=""; file3=concat(file2,iteration); 
file4=""; file4=concat(filenorma,iteration);

weights(chrom_start,chrom_end,rf,enzyme,enriched,map_threshold,weight_threshold,totfragments,hic,map,weight2,tot,tot2,file3);

conv_check=(double)(0.0); for(j=2;j<=tot;j++){if(weight2[j]!=0){conv_check+=fabs((double)(1.0)-weight2[j]);}}

fprintf(fout3,"%d\t\t\t%lf\t\t\t%lf\n",iter,conv_check,tolerance); 

printf("%d\t%lf\t%lf\n",iter,conv_check,tolerance); 
}

filewrite(BIN,resolution,chrom_start,chrom_end,iter,rf,enzyme,enriched,totfragments,hic,tot,tot2,tot3,downstore,upstore,avg,avg2,dilute,logistic,avg_norma,avg_norma2,trimming,stats,TRANS,PVAL,FDR,fileoutput1,fileoutput2,file4,folderref,folderrefC);

if(strcasecmp(stats,"none")!=0){
printf("\n\n Statistical Analysis Begin here \n");

sprintf(iteration,"%d",iter);

for(chrom=chrom_start;chrom<=chrom_end;chrom++){            
sprintf(chromosome,"%d",chrom);
printf(".. Processing Chromosome %d .. Please be patient.. \n",chrom);

sprintf(command,"R --vanilla < %sHIC.R --args %s %s %s %s %s %lf &>>TMP \n",folderrefC,concat(concat(concat(concat(fileoutput2,"/self.sample_it_"),iteration),"_chr_"),chromosome),concat(concat(concat(concat(fileoutput2,"/self.position_it_"),iteration),"_chr_"),chromosome),concat(concat(concat(concat(fileoutput2,"/self_output_it_"),iteration),"_chr_"),chromosome),stats,trimming);

printf("\n\nStatistical Analyisis on Baited-Baited contacts.. \n %s\n .. on chromosome %d\n\n",command,chrom);
system(command);

if(strcasecmp(BAITFLAG,"no")==0){}else{
    sprintf(command,"R --vanilla < %sHIC.R --args %s %s %s %s %s %lf &>>TMP \n",folderrefC,concat(concat(concat(concat(fileoutput2,"/cis.sample_it_"),iteration),"_chr_"),chromosome),concat(concat(concat(concat(fileoutput2,"/cis.position_it_"),iteration),"_chr_"),chromosome),concat(concat(concat(concat(fileoutput2,"/cis_output_it_"),iteration),"_chr_"),chromosome),stats,trimming);
	  system(command);
printf("\n\nStatistical Analyisis on Baited-NonBaited contacts.. \n %s\n .. on chromosome %d\n\n",command,chrom);
if(strcasecmp(PWEIGHT,"no")==0){
sprintf(command,"paste %s  %s |gawk '{print $1\"_\"$2\"_\"$3 \"\t\" $4}' > %s \n",concat(concat(concat(concat(concat(fileoutput2,"/cis_output_it_"),iteration),"_chr_"),chromosome),".pos"),concat(concat(concat(concat(concat(concat(fileoutput2,"/cis_output_it_"),iteration),"_chr_"),chromosome),"."),stats),concat(concat(concat(concat(concat(fileoutput2,"/cis.output_it_"),iteration),"_chr_"),chromosome),".tab"));
}else{
sprintf(command,"paste %s  %s |gawk '{if($3!=0&&$4/$3<1){print $1\"_\"$2\"_\"$3 \"\t\" $4/$3};if($3!=0&&$4/$3>=1){print $1\"_\"$2\"_\"$3 \"\t 1 \" }}' > %s \n",concat(concat(concat(concat(concat(fileoutput2,"/cis_output_it_"),iteration),"_chr_"),chromosome),".pos"),concat(concat(concat(concat(concat(concat(fileoutput2,"/cis_output_it_"),iteration),"_chr_"),chromosome),"."),stats),concat(concat(concat(concat(concat(fileoutput2,"/cis.output_it_"),iteration),"_chr_"),chromosome),".tab"));
}
    system(command);

        sprintf(command,"R --vanilla < %sHIC_FDR.R --args %s model none %lf %lf %s &>>TMP \n",folderrefC,concat(concat(concat(concat(concat(fileoutput2,"/cis.output_it_"),iteration),"_chr_"),chromosome),".tab"),PVAL,FDR,concat(concat(concat(concat(concat(fileoutput2,"/fdr/fdr_cis_output_it_"),iteration),"_chr_"),chromosome),".tab"));
printf("%s\n",command);
system(command);
printf("\n\nFalse Discovery Rate Analysis on Baited-NonBaited contacts.. \n %s\n .. on chromosome %d\n\n",command,chrom);
}
if(strcasecmp(PWEIGHT,"no")==0){
sprintf(command,"paste %s  %s |gawk '{print $1\"_\"$2\"_\"$3 \"\t\" $4}' > %s \n",concat(concat(concat(concat(concat(fileoutput2,"/self_output_it_"),iteration),"_chr_"),chromosome),".pos"),concat(concat(concat(concat(concat(concat(fileoutput2,"/self_output_it_"),iteration),"_chr_"),chromosome),"."),stats),concat(concat(concat(concat(concat(fileoutput2,"/self.output_it_"),iteration),"_chr_"),chromosome),".tab"));
}else{
sprintf(command,"paste %s  %s |gawk '{if($3!=0&&$4/$3<1){print $1\"_\"$2\"_\"$3 \"\t\" $4/$3};if($3!=0&&$4/$3>=1){print $1\"_\"$2\"_\"$3 \"\t 1 \" }}' > %s \n",concat(concat(concat(concat(concat(fileoutput2,"/self_output_it_"),iteration),"_chr_"),chromosome),".pos"),concat(concat(concat(concat(concat(concat(fileoutput2,"/self_output_it_"),iteration),"_chr_"),chromosome),"."),stats),concat(concat(concat(concat(concat(fileoutput2,"/self.output_it_"),iteration),"_chr_"),chromosome),".tab"));
}
printf("%s\n",command);
system(command);

sprintf(command,"R --vanilla < %sHIC_FDR.R --args %s model none %lf %lf %s &>>TMP \n",folderrefC,concat(concat(concat(concat(concat(fileoutput2,"/self.output_it_"),iteration),"_chr_"),chromosome),".tab"),PVAL,FDR,concat(concat(concat(concat(concat(fileoutput2,"/fdr/fdr_self_output_it_"),iteration),"_chr_"),chromosome),".tab"));

system(command);
printf("\n\nFalse Discovery Rate Analysis on Baited-Baited contacts.. \n %s\n .. on chromosome %d\n\n",command,chrom);

sprintf(command,"gawk  -F'[_]' '{split($1,c,\",\");output=\"%s_\"c[1]\"_\"c[2]\"_\"c[3]; ;split($2,a,\",\");split($3,b,\"\t\"); if(b[3] < %lf && b[4]< %lf ){print $1 \"\t\" a[1]\",\"a[2]\",\"a[3]\"\t\"(-1)*log(b[3]) >> output}}' %s |grep -v \"V1\"|sort -r -nk3,3 ",outdir,PVAL,FDR,concat(fileoutput2,"/fdr/*fdr.tab"));

printf("\n\nGenerating significant interactions. You can plot results with, e.g., Wash-U Browser \n\n%s\n",command);
    system(command);
sprintf(command,"gawk  -F'[_]' '{split($1,c,\",\"); ;split($2,a,\",\");split($3,b,\"\t\"); if(b[3]< %lf && b[4]< %lf ){print $1 \"\t\" a[1]\",\"a[2]\",\"a[3]\"\t\"(-1)*log(b[3]) }}' %s |grep -v \"V1\" >> %s ",PVAL,FDR,concat(fileoutput2,"/fdr/*fdr.tab"),concat(outdir2,".WashUformat"));
   system(command);

printf("Statistical Analysis Ends Here. Goodbye.\n"); 
 }


if(strlen(folderN)==0){

sprintf(command,"touch %s/%s/%s_domain_caller.sh",folder,lane,lane);
system(command);
sprintf(heatmap,concat(concat(concat(concat(folder,lane),stamp),"/"),concat(concat(lane,stamp2),"_self")));
sprintf(graphics,concat(concat(concat(concat(folder,lane),stamp),"/"),"/graphics"));

sprintf(command,"echo '#!/bin/bash ' >> %s/%s/%s_domain_caller.sh",folder,lane,lane);
system(command);
sprintf(command,"echo 'HM=%s' >> %s/%s/%s_domain_caller.sh",heatmap,folder,lane,lane);
system(command);
sprintf(command,"echo 'GRAPHICS=%s' >> %s/%s/%s_domain_caller.sh",graphics,folder,lane,lane);
system(command);
if(strcasecmp(BAITFLAG,"no")==0){sprintf(command,"echo 'RES=%d' >> %s/%s/%s_domain_caller.sh",SCALE*resolution,folder,lane,lane);}else{
sprintf(command,"echo 'RES=%d' >> %s/%s/%s_domain_caller.sh",SCALE*resolution,folder,lane,lane);}
system(command);

for(chrom=chrom_start;chrom<=chrom_end;chrom++){

  sprintf(command,"echo 'CHROM=%d' >> %s/%s/%s_domain_caller.sh",chrom,folder,lane,lane);
  system(command);
  sprintf(command,"echo 'INT=%d' >> %s/%s/%s_domain_caller.sh",reg[chrom]-reg[chrom-1],folder,lane,lane);
  system(command);

if(strcasecmp(BAITFLAG,"no")==0){
sprintf(command,"cat %sgraphics.sh >> %s/%s/%s_domain_caller.sh",folderrefC,folder,lane,lane);
}else{
sprintf(command,"cat %scapturegraphics.sh >> %s/%s/%s_domain_caller.sh",folderrefC,folder,lane,lane);}
system(command);
}
sprintf(command,"chmod a+x %s/%s/%s_domain_caller.sh",folder,lane,lane);
sprintf(command,"%s/%s/%s_domain_caller.sh",folder,lane,lane);
 }
 }
exit(0);
}

filewrite(BIN,resolution,chrom_start,chrom_end,iteration,rf,enzyme,enriched,totfragments,hic,tot,tot2,tot3,downstore,upstore,avg,avg2,dilute,logit,avg_norma,avg_norma2,trimming,stats,TRANS,PVAL,FDR,fileout1,fileout2,fileout3,folderref,folderrefC)
int BIN; int resolution;int chrom_start;int chrom_end; int iteration;  int rf[CHROMOSOMES];int **enzyme; int **enriched;  int *totfragments;  double **hic;int tot; int tot2; int tot3; int *downstore; int *upstore; int *avg_norma;int *avg_norma2; double *avg;double *avg2;double *dilute; double *logit; double trimming; char *stats; char *TRANS; double PVAL; double FDR;char *fileout1;char *fileout2;char *fileout3; char *folderref; char *folderrefC;
{

char itera[100];
sprintf(itera,"%d",iteration);
FILE *fout6 = fopen (concat(concat(concat(fileout2,"/stats_self_it"),itera),".sh"), "w+");
FILE *fout6b = fopen (concat(concat(concat(fileout2,"/stats_cis_it"),itera),".sh"), "w+");

if(fout6==NULL||fout6b==NULL){printf("stat file is empty or missing\n"); exit(0);}
int chrom,i,j,k,kk,lg,bin,allbin,allnonzerobin,coverage;
int position,position_near,position2,down,up,total,avg_size;
char viewpoint[100];
char chromosome[5];
double sum,distance,eta;

/**********************************/
double alpha = 34.1158;
double beta = -2.58689;
double gamma = -17.1348;
double delta = -7.0761;

bin=0;sum=0;

total=0;

/************************************************************/

for (i=0;i<tot3;i++){norma[i]=0; }
for(i=0;i<tot3;i++){dilute[i]=0;avg[i]=0;avg2[i]=0; avg_norma[i]=0; avg_norma2[i]=0;}

for(i=1;i<=tot2;i++){

position=store[i]+1;
chrom=store2[i];
down=1;
up=enzyme[chrom][rf[chrom]-1];
for(j=grid(chrom,down,rf,enzyme);j<=grid(chrom,up,rf,enzyme);j++){
bin=minim(tot3-1,rint(fabs(position-(enzyme[chrom][j]+enzyme[chrom][j+1])/2)/BIN));
if(hic[i][totfragments[chrom-chrom_start]+j]>0){
avg[bin]=avg[bin]+hic[i][totfragments[chrom-chrom_start]+j];
avg_norma[bin]++;
}else{avg_norma2[bin]++;}
}
}
for(k=0;k<tot3-1;k++){
down=maxim(0,rint(0.005*k+1));
up=minim(rint((0.005*k+1)),rint(0.005*(tot3)+1));
for(kk=k-down;kk<=minim(k+up,tot3-1);kk++){
norma[k]+=(avg[kk])/((tot3-kk)*(up+down+1));
if((avg_norma[kk]+avg_norma2[kk])!=0){dilute[k]+=((avg_norma[kk]*(tot3-kk))/((avg_norma[kk]+avg_norma2[kk])))/(up+down+1); }
}
}

total=0;
  for(i=0;i<50000;i++){
  for(j=i+1;j<50000;j++){
  avg_size=(j-i)*BIN;
  sum+=expit(alpha+beta*log(avg_size)); total=total+1; }
}

for(i=0;i<tot3;i++){
   logistic[i]=0;
   distance=i*BIN;
   eta=expit(alpha+beta*log(distance));
   logistic[i]=log((expit(delta)-expit(gamma))*eta+expit(gamma))-log((expit(delta)-expit(gamma))*sum/total+expit(gamma));}

FILE *fout2 = fopen (fileout3, "w");
if(fout2==NULL){printf("normalization file cannot be opened or folder is missing\n"); exit(0);}
coverage=0;
for(i=0;i<=tot3;i++){
if(dilute[i]/dilute[1]>0.01){fprintf(fout2,"%d\t%lf\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\n",i*BIN,norma[i],avg[i],dilute[i]/dilute[1],avg_norma[i],avg_norma2[i],norma[i]/(dilute[i]/dilute[1]),logistic[i],exp(logistic[i]));
coverage++;}
if(dilute[i]!=0){norma[i]=(norma[i])/(dilute[i]/dilute[1]);}
}
fclose(fout2);
/*************************************************/
  for(i=1;i<tot2; i++){

  if(iteration==0){  sprintf(viewpoint,"frag_%d_%d_%d_%d_%d.in",i,store2[i],store2b[i],store[i],enzyme[store2[i]][grid(store2[i],store[i]+1,rf,enzyme)+1]);}
  if(iteration!=0){  sprintf(viewpoint,"frag_%d_%d_%d_%d_%d.fin",i,store2[i],store2b[i],store[i],enzyme[store2[i]][grid(store2[i],store[i]+1,rf,enzyme)+1]);}

 sprintf(chromosome,"%d",store2[i]);
 sprintf(itera,"%d",iteration);

 FILE *fout = fopen (concat(concat(fileout1,"cis/"),viewpoint), "w+");
 FILE *fout2b = fopen (concat(concat(fileout1,"self/"),viewpoint), "w+");
 FILE *fout3b = fopen (concat(concat(fileout1,"trans/"),viewpoint), "w+");
 FILE *fout4 = fopen (concat(concat(concat(concat(concat(fileout2,"/self.position_it_"),itera),"_chr_"),chromosome),""), "a+");
 FILE *fout5 = fopen (concat(concat(concat(concat(concat(fileout2,"/self.sample_it_"),itera),"_chr_"),chromosome),""), "a+");
 FILE *fout4b = fopen (concat(concat(concat(concat(concat(fileout2,"/cis.position_it_"),itera),"_chr_"),chromosome),""), "a+");
 FILE *fout5b = fopen (concat(concat(concat(concat(concat(fileout2,"/cis.sample_it_"),itera),"_chr_"),chromosome),""), "a+");

 if(fout==NULL||fout2b==NULL||fout3b==NULL||fout4==NULL||fout5==NULL||fout4b==NULL||fout5b==NULL){printf("stat file is empty or missing\n"); exit(0);}
 for(chrom=chrom_start;chrom<=chrom_end;chrom++){

 allbin=0; allnonzerobin=0;
 if(store2[i]==chrom){

 for(j=totfragments[chrom-chrom_start]+1;j<totfragments[chrom-chrom_start+1];j++){
         lg=rint(fabs((double)(store[i]-enzyme[chrom][j-totfragments[chrom-chrom_start]]))/BIN);
	 allbin++; if(hic[i][j]!=0){allnonzerobin++;}}
         for(j=totfragments[chrom-chrom_start]+1;j<totfragments[chrom-chrom_start+1];j++){
               lg=rint(fabs((double)(store[i]-enzyme[chrom][j-totfragments[chrom-chrom_start]]))/BIN);
               if(enriched[chrom][j-totfragments[chrom-chrom_start]]==0){
	       fprintf(fout,"%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\n",store2[i],store[i],chrom,enzyme[chrom][j-totfragments[chrom-chrom_start]],hic[i][j],norma[lg],(double)(allnonzerobin)/allbin);
	       }else if(enriched[chrom][j-totfragments[chrom-chrom_start]]!=0){
               fprintf(fout2b,"%d %d %d\t%d\t%lf\t%lf\t%lf\n",store2[i],store[i],chrom,enzyme[chrom][j-totfragments[chrom-chrom_start]],hic[i][j],norma[lg],(double)(allnonzerobin)/allbin);}
	       if(enriched[chrom][j-totfragments[chrom-chrom_start]]==0){
               if(allbin!=0&&(double)(allnonzerobin)/allbin>trimming&&lg>1&&norma[lg]>0&&lg<coverage){
	       fprintf(fout4b,"chr%d,%d,%d\tchr%d,%d,%d\t%lf\n",store2[i],store[i],enzyme[store2[i]][grid(store2[i],store[i]+1,rf,enzyme)+1],chrom,enzyme[chrom][j-totfragments[chrom-chrom_start]],enzyme[chrom][j-totfragments[chrom-chrom_start]+1],exp(logistic[lg]));
	       fprintf(fout5b,"%lf\t%lf\t%lf\t%lf\n",hic[i][j]/norma[lg],logistic[lg],norma[lg],(double)(allnonzerobin)/allbin);
	       }}else if(enriched[chrom][j-totfragments[chrom-chrom_start]]!=0&&lg>1&&norma[lg]>0&&lg<coverage){
               fprintf(fout4,"chr%d,%d,%d\tchr%d,%d,%d\t%lf\n",store2[i],store[i],enzyme[store2[i]][grid(store2[i],store[i]+1,rf,enzyme)+1],chrom,enzyme[chrom][j-totfragments[chrom-chrom_start]],enzyme[chrom][j-totfragments[chrom-chrom_start]+1],exp(logistic[lg]));
	       fprintf(fout5,"%lf\t%lf\t%lf\t%lf\n",hic[i][j]/norma[lg],logistic[lg],norma[lg],(double)(allnonzerobin)/allbin);
             
	   }}}
           else{ 
for(j=totfragments[chrom-chrom_start]+1;j<totfragments[chrom-chrom_start+1];j++){
if(strcasecmp(TRANS,"no")!=0){fprintf(fout3b,"%d\t%d\t%d\t%d\t%lf\n",store2[i],enzyme[store2[i]][grid(store2[i],store[i]+1,rf,enzyme)],chrom,enzyme[chrom][j-totfragments[chrom-chrom_start]],hic[i][j]);
}
}
}
}
 fclose(fout);
 fclose(fout2b);
 fclose(fout3b);
 fclose(fout4);
 fclose(fout4b);
 fclose(fout5);
 fclose(fout5b);
  }

for(chrom=chrom_start;chrom<=chrom_end;chrom++){            
sprintf(chromosome,"%d",chrom);
sprintf(itera,"%d",iteration);
fprintf(fout6,"R --vanilla < %sHIC.R --args %s %s %s %s %lf \n",folderrefC,concat(concat(concat(concat(fileout2,"/self.sample_it_"),itera),"_chr_"),chromosome),concat(concat(concat(concat(fileout2,"/self.position_it_"),itera),"_chr_"),chromosome),concat(concat(concat(concat(fileout2,"/self_output_it_"),itera),"_chr_"),chromosome),stats,trimming);
fprintf(fout6,"paste %s  %s |gawk '{if($3!=0&&$4/$3<1){print $1\"_\"$2\"_\"$3 \"\t\" $4/$3};if($3!=0&&$4/$3>=1){print $1\"_\"$2\"_\"$3 \"\t 1 \" }}'  > %s \n",concat(concat(concat(concat(concat(fileout2,"/self_output_it_"),itera),"_chr_"),chromosome),".pos"),concat(concat(concat(concat(concat(concat(fileout2,"/self_output_it_"),itera),"_chr_"),chromosome),"."),stats),concat(concat(concat(concat(concat(fileout2,"/self_output_it_"),itera),"_chr_"),chromosome),".tab"));
fprintf(fout6,"R --vanilla < %sHIC_FDR.R --args %s model none %lf %lf %s \n",folderrefC,concat(concat(concat(concat(concat(fileout2,"/self_output_it_"),itera),"_chr_"),chromosome),".tab"),FDR,FDR,concat(concat(concat(concat(concat(fileout2,"/fdr/fdr_self_output_it_"),itera),"_chr_"),chromosome),".tab"));

fprintf(fout6b,"R --vanilla < %sHIC.R --args %s %s %s %s %lf \n",folderrefC,concat(concat(concat(concat(fileout2,"/cis.sample_it_"),itera),"_chr_"),chromosome),concat(concat(concat(concat(fileout2,"/cis.position_it_"),itera),"_chr_"),chromosome),concat(concat(concat(concat(fileout2,"/cis_output_it_"),itera),"_chr_"),chromosome),stats,trimming);

fprintf(fout6b,"paste %s  %s |gawk '{if($3!=0&&$4/$3<1){print $1\"_\"$2\"_\"$3 \"\t\" $4/$3};if($3!=0&&$4/$3>=1){print $1\"_\"$2\"_\"$3 \"\t 1 \" }}' > %s \n",concat(concat(concat(concat(concat(fileout2,"/cis_output_it_"),itera),"_chr_"),chromosome),".pos"),concat(concat(concat(concat(concat(concat(fileout2,"/cis_output_it_"),itera),"_chr_"),chromosome),"."),stats),concat(concat(concat(concat(concat(fileout2,"/cis_output_it_"),itera),"_chr_"),chromosome),".tab"));

fprintf(fout6b,"R --vanilla < %sHIC_FDR.R --args %s model none %lf %lf %s \n",folderrefC,concat(concat(concat(concat(concat(fileout2,"/cis_output_it_"),itera),"_chr_"),chromosome),".tab"),FDR,FDR,concat(concat(concat(fileout2,"/fdr/fdr_cis_output_it_"),itera),".tab"));
	 }
 fclose(fout6);
 fclose(fout6b);
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


int grid2(chr,pos,store,store2,downstore,upstore,tot2)
     int *store; int *store2; int *downstore; int *upstore; int chr; int pos; int tot2;
{
  int flag;
  int i,up,down,middle,save;
  up=upstore[chr];
  down=downstore[chr];

  middle=(int)((up+down)/2);
  flag=0;

  while(flag==0){
    if(pos>store[middle]){down=middle;};
    if(pos<store[middle]){up=middle;};
    if(pos==store[middle]){up=middle;flag=1;save=middle;}
    middle=(int)((up+down)/2);
    if(up==(down+1)){flag=2;save=down;} 
  }
    if(up==(down+2)){flag=1;save=down;}
return(save);
}

weights(chrom_start,chrom_end,rf,enzyme,enriched,map_threshold,weight_threshold,totfragments,hic,map,weight2,tot,tot2,file)
int **enzyme; int **enriched;  int *totfragments;double map_threshold;double weight_threshold; double *map; int rf[CHROMOSOMES]; char *file;
int chrom_start;int tot; int tot2; int chrom_end; double *weight2; double **hic;
{

int chrom,i,j;
double area1,area2,tot_count1,tot_count2,tmp,tmp2;
area1=0;area2=0; tot_count1=0; tot_count2=0;

for(i=0;i<=tot; i++){
   weight[i]=(double)(0.0);
   *(weight2+i)=(double)(0.0);
   *(num_cells+i)=0;}

FILE *fout = fopen (file, "w");

         for(chrom=chrom_start;chrom<=chrom_end;chrom++){            
         for(i=totfragments[chrom-chrom_start]+1;i<=totfragments[chrom-chrom_start+1];i++){
         store3[i]=chrom;}}

         for(chrom=chrom_start;chrom<=chrom_end;chrom++){
         for(i=totfragments[chrom-chrom_start]+1;i<=totfragments[chrom-chrom_start+1];i++){
	      if(map[i]>map_threshold){	
	      if(enriched[chrom][i-totfragments[chrom-chrom_start]]==0){
              for(j=1;j<tot2;j++){          	
              if(map[totfragments[store2[j]-chrom_start]+grid(store2[j],store[j]+1,rf,enzyme)]>map_threshold&&store2[j]!=chrom)
		{		  
                      weight[i]=weight[i]+hic[j][i];              
                      num_cells[i]=num_cells[i]+1;           
		}   
	      }
              area2+=num_cells[i]; tot_count2+=weight[i];}

	      if(enriched[chrom][i-totfragments[chrom-chrom_start]]!=0){	
              

           for(j=1;j<tot;j++){
              if(map[j]>map_threshold&&store3[j]!=chrom)
                {
                      weight[i]=weight[i]+hic[enriched[chrom][i-totfragments[chrom-chrom_start]]][j];
                      num_cells[i]=num_cells[i]+1;
                }
              }
              area1+=num_cells[i]; tot_count1+=weight[i]; }
	      }
	 }
	 }

	 /*	 fprintf(fout4,"%lf \t %lf\t %lf\t %lf\t%lf\n",area1,area2,tot_count1,tot_count2,tot_count1+tot_count2);*/
	 /*	 if(area1!=0&&area2!=0&&tot_count1>0&&tot_count2>0){*/

         for(chrom=chrom_start;chrom<=chrom_end;chrom++){            
         for(i=totfragments[chrom-chrom_start]+1;i<=totfragments[chrom-chrom_start+1];i++){

	   if(num_cells[i]>0&&map[i]>map_threshold){	    
	   if(enriched[chrom][i-totfragments[chrom-chrom_start]]!=0&&area1!=0){ weight2[i]=(weight[i]/(double)(num_cells[i]))/(tot_count1/area1);          }
	   if(enriched[chrom][i-totfragments[chrom-chrom_start]]==0&&area2!=0){ weight2[i]=(weight[i]/(double)(num_cells[i]))/(tot_count2/area2);          }	   
	   }

if(area1!=0){tmp=area1*(tot_count1/area1);}
if(area2!=0){tmp2=area2*(tot_count2/area2);}
         
	 fprintf(fout,"%d\t%d\t%d\t%d\t%d\t%lf\t%d\t%lf\t%lf\t%lf\t%d\t%lf\t%lf\t%lf\n",i,enzyme[chrom][i-totfragments[chrom-chrom_start]],enriched[chrom][i-totfragments[chrom-chrom_start]],chrom,i-totfragments[chrom-chrom_start],map[i],enriched[chrom][i-totfragments[chrom-chrom_start]],weight[i],area1,area2,num_cells[i],weight2[i],tmp,tmp2); 

         if(weight2[i]<weight_threshold){map[i]=0.0;}

	 }
	 }

fclose(fout);
}


double digest(chrom,chrom_end,rf,enzyme,ENZYME,CUT,SCALE,GRID,folderref,REFERENCE,resolution)
     int **enzyme; int rf[CHROMOSOMES]; int chrom; int chrom_end; int CUT; int resolution; char *ENZYME; int SCALE; int GRID; char *REFERENCE; char *folderref;
    {
    FILE *reference_genome;
    char c,str[20],*chr;      
    int win_size=0;
    int pos=0;
    int diff=0;
    
    int up,down,middle,middle2;
    int saved;
    int quality;
    int max,cuts;
    int start,end;
    int res,b;

    long unsigned int sizeofchrom[25];
    unsigned int dist,dist2;    
    double fdist,tmp,tmp3,grid,slope,smooth,average;
    char window1[WIN_SIZE];
    char schrom[25],command[500],filename[500];

    sprintf(schrom,"%d",chrom);        

max=0;
sprintf(filename,"%s/chr%d.fa",folderref,chrom);
if(chrom==chrom_end-2){sprintf(filename,"%s/chrX.fa",folderref);}
if(chrom==chrom_end-1){sprintf(filename,"%s/chrY.fa",folderref);}
if(chrom==chrom_end){sprintf(filename,"%s/chrM.fa",folderref);}

    if(fopen(filename, "rt" )==NULL){ 

sprintf(command,"wget -nv -P %s \'http://hgdownload.cse.ucsc.edu/goldenpath/'%s'/chromosomes/chr%d.fa.gz\'",folderref,REFERENCE,chrom);

if(chrom==chrom_end-2){sprintf(command,"wget -nv -P %s \'http://hgdownload.cse.ucsc.edu/goldenpath/'%s'/chromosomes/chrX.fa.gz\'",folderref,REFERENCE);}
if(chrom==chrom_end-1){sprintf(command,"wget -nv -P %s \'http://hgdownload.cse.ucsc.edu/goldenpath/'%s'/chromosomes/chrY.fa.gz\'",folderref,REFERENCE);}
if(chrom==chrom_end){sprintf(command,"wget -nv -P %s \'http://hgdownload.cse.ucsc.edu/goldenpath/'%s'/chromosomes/chrM.fa.gz\'",folderref,REFERENCE);}

   system(command);

   sprintf(command,"gunzip -f %s/chr%d.fa.gz",folderref,chrom);

if(chrom==chrom_end-2){sprintf(command,"gunzip -f %s/chrX.fa.gz",folderref);}
if(chrom==chrom_end-1){sprintf(command,"gunzip -f %s/chrY.fa.gz",folderref);}
if(chrom==chrom_end){sprintf(command,"gunzip -f %s/chrM.fa.gz",folderref);}

printf("%s\n",command);

system(command); 

}

      sprintf(command,"%s/chr%d.fa",folderref,chrom);

      if(chrom==chrom_end-2){ sprintf(command,"%s/chrX.fa",folderref,chrom);}
      if(chrom==chrom_end-1){ sprintf(command,"%s/chrY.fa",folderref,chrom);}
      if(chrom==chrom_end){ sprintf(command,"%s/chrM.fa",folderref,chrom);}

      reference_genome=fopen(command, "rt" );
      printf("\n*Digesting reference chromosome %i \n",chrom);    

      rf[chrom]=0;
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
	    if((strcasecmp(window1,ENZYME) ==0&&GRID==0)||(pos%resolution ==0&&GRID==1)){cuts++;}	               
            memmove(window1,&window1[WIN_SHIFT1],WIN_SIZE-WIN_SHIFT1);     pos+=WIN_SHIFT1;   	      	
            win_size=WIN_SIZE-WIN_SHIFT1;
            }
	}

    pos=0;
    cuts = (int)(cuts/SCALE);

    win_size=0; pos=0; window1[0]='\0';
    enzyme[chrom] = malloc(sizeof *enzyme[chrom] * cuts);
    enriched[chrom] = malloc(sizeof *enriched[chrom] * cuts);

    cuts=0;
 rewind(reference_genome);
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
	      if((strcasecmp(window1,ENZYME) ==0&&GRID==0)||(pos%resolution ==0&&GRID==1)){
		cuts++;if(cuts%SCALE==0){rf[chrom]++;   *(enzyme[chrom]+rf[chrom])=pos+CUT; *(enriched[chrom]+rf[chrom])=0; }
	    }
            memmove(window1,&window1[WIN_SHIFT1],WIN_SIZE-WIN_SHIFT1);      pos+=WIN_SHIFT1;
            win_size=WIN_SIZE-WIN_SHIFT1;
	    }  
	}
    return (0);	
}	

