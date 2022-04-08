#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINE_MAX 0x4000


int main(int argc,char *argv[])
{
  int i,j,tot,count,count2,chrom,quality,shift,tmp1,tmp2,icount,jcount;

  char buffer[1000];
  char qname[100],qname2[100],qnames[100],temp[100];
  char head[100],flag[10],flag2[10],rname[25],rname2[25],position[10],position2[10],mapq[10];
  char mapq2[10],cigar[100],cigar2[100],rnext[10],pnext[10],tlen[10],sequence[500];
  char sequence2[500],qual[500],tags[1000];
 
  FILE *output,*output2;
  char o1[500],o2[500],o3[500];
  char *substr;

  char (*vqname)[100];
  char (*vrname2)[25];
  char (*vposition)[10];
  char (*vflag)[10];
  char (*vmapq)[10];
  char (*vcigar)[100];
  char (*vsequence)[500];

tot=500;

vqname= malloc(tot * sizeof *vqname);
vrname2=malloc(tot * sizeof *vrname2);
vposition=malloc(tot * sizeof *vposition);
vflag=malloc(tot* sizeof *vflag);
vmapq=malloc(tot* sizeof *vmapq);
vcigar=malloc(tot* sizeof *vcigar);
vsequence=malloc(tot* sizeof *vsequence);

sprintf(o1,"%s%s",argv[2],".sam");
sprintf(o3,"%s%s",argv[2],".unmapped.sam");
sprintf(o2,"%s",argv[3]);
if(strcasecmp(o2,"no")==0){shift=0;}else{shift=1;}

output=fopen(o1,"w");
output2=fopen(o3,"w");
count=0;
count2=0;
FILE *fin = fopen ( argv[1], "r" );
while(fscanf(fin,"%[^\n]\n", buffer) != EOF){    
 
sscanf (buffer,"%s\n",&head);
if(strncmp(head,"@",1)==0){fprintf(output,"%s\n",buffer);}else{

sscanf (buffer,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%[^\n]\n",qname,flag,rname,position,mapq,cigar,rnext,pnext,tlen,sequence,qual,tags);

count++;
strcpy(qnames,qname);
substr=strtok(qnames,"P");
strcpy(vqname[count2],temp);
strcpy(vrname2[count2],rname2);
strcpy(vposition[count2],position2);
strcpy(vflag[count2],flag2);
strcpy(vmapq[count2],mapq2);
strcpy(vcigar[count2],cigar2);
strcpy(vsequence[count2],sequence2);

if(strcasecmp(qname2,substr)==0){
count2++;
}else{
for(i=0;i<=count2;i++){
icount=0;
if(strncmp(vflag[i],"16",2)==0){
    char* c=vcigar[i];

        while(*c!=0)
            {
            char* p2;
                if(!isdigit(*c)) {       fprintf(stderr,"bad cigar string: %s\n",*cigar); return 0; }
            int n=strtol(c,&p2,10);
                if(n<1){                 fprintf(stderr,"bad cigar string: %s\n",*cigar); return 0; }
                switch(*p2)
                {
                case 'M':{ while(n>0) { icount++; --n; } break; }
                case 'D':{ while(n>0) { icount++; --n; } break; }
                case 'I':{ while(n>0) { --n; } break; }
                case 'P':{ break;}default: { fprintf(stderr,"unsupported operator %c\n",*p2);  return 0; }} c=p2+1; }}

for(j=i+1;j<=count2;j++){

if((strncmp(vflag[j],"16",2)==0||strncmp(vflag[j],"0",1)==0)&&strlen(vsequence[j])>20&&strlen(vrname2[j])<=5&&(strncmp(vflag[i],"16",2)==0||strncmp(vflag[i],"0",1)==0)&&strlen(vsequence[i])>20&&strlen(vrname2[i])<=5){

jcount=0;
if(strncmp(vflag[j],"16",2)==0){
    char* c=vcigar[j];

        while(*c!=0)
            {
            char* p2;
                if(!isdigit(*c)) {       fprintf(stderr,"bad cigar string: %s\n",*cigar); return 0; }
            int n=strtol(c,&p2,10);
                if(n<1){                fprintf(stderr,"bad cigar string: %s\n",*cigar); return 0; }
                switch(*p2)
                {
                case 'M':{ while(n>0) { jcount++; --n; } break; }
                case 'D':{ while(n>0) { jcount++; --n; } break; }
                case 'I':{ while(n>0) { --n; } break; }
                case 'P':{ break;}default: { fprintf(stderr,"unsupported operator %c\n",*p2);  return 0; }} c=p2+1; }

}

tmp1=atoi(strtok(vrname2[i],"chr"));
tmp2=atoi(strtok(vrname2[j],"chr"));

if(strncmp(strtok(vrname2[i],"chr"),"X",1)==0){ tmp1=23; }
if(strncmp(strtok(vrname2[i],"chr"),"Y",1)==0){ tmp1=24; }

if(strncmp(strtok(vrname2[j],"chr"),"X",1)==0){ tmp2=23; }
if(strncmp(strtok(vrname2[j],"chr"),"Y",1)==0){ tmp2=24; }

if(tmp1<tmp2){
fprintf(output,"%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",vflag[i],vrname2[i],atoi(vposition[i])+icount,i,vflag[j],vrname2[j],atoi(vposition[j])+jcount,j,vmapq[i],vcigar[i],vsequence[i],vmapq[j],vcigar[j],vsequence[j],vqname[i],vqname[j]);
}else if(tmp1>tmp2){
fprintf(output,"%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",vflag[j],vrname2[j],atoi(vposition[j])+jcount,j,vflag[i],vrname2[i],atoi(vposition[i])+icount,i,vmapq[j],vcigar[j],vsequence[j],vmapq[i],vcigar[i],vsequence[i],vqname[j],vqname[i]);
}else if(tmp1==tmp2){
if(atoi(vposition[i])+icount<atoi(vposition[j])+jcount){
fprintf(output,"%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",vflag[i],vrname2[i],atoi(vposition[i])+icount,i,vflag[j],vrname2[j],atoi(vposition[j])+jcount,j,vmapq[i],vcigar[i],vsequence[i],vmapq[j],vcigar[j],vsequence[j],vqname[i],vqname[j]);
 }else{
fprintf(output,"%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",vflag[j],vrname2[j],atoi(vposition[j])+jcount,j,vflag[i],vrname2[i],atoi(vposition[i])+icount,i,vmapq[j],vcigar[j],vsequence[j],vmapq[i],vcigar[i],vsequence[i],vqname[j],vqname[i]);
 }
 }
 }else{
fprintf(output2,"%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",vflag[i],vrname2[i],atoi(vposition[i])+icount,i,vflag[j],vrname2[j],atoi(vposition[j])+jcount,j,vmapq[i],vcigar[i],vsequence[i],vmapq[j],vcigar[j],vsequence[j],vqname[i],vqname[j]);

 }
 }
 }
count2=0;

 }

strcpy(qname2,substr);
strcpy(temp,qname);
strcpy(rname2,rname);
strcpy(position2,position);
strcpy(flag2,flag);
strcpy(mapq2,mapq);
strcpy(cigar2,cigar);
strcpy(sequence2,sequence);

 }
 }
}

typedef void(*split_fn)(const char *, size_t, void *);
 
void split(const char *str, char sep, split_fn fun, void *data)
{
    unsigned int start = 0, stop;
    for (stop = 0; str[stop]; stop++) {
        if (str[stop] == sep) {
            fun(str + start, stop - start, data);
            start = stop + 1;
        }
    }
    fun(str + start, stop - start, data);
}


