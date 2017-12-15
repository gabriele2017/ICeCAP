#include <stdio.h>
#include <string.h>

#define LINE_MAX 0x4000
int main(int argc,char *argv[])
{
  char record[100],interval[5],read[5],*chr,*int_1,*int_2;                 /* array to hold each "record" */ 
  int count,chrom,quality,position2,shift;
  char buffer[1000];
    
  char head[100],qname[100],flag[10],rname[10],rname2[10],position[10],mapq[10],cigar[100],rnext[10],pnext[10],tlen[10],sequence[500],qual[500],tags[1000];

 
  FILE *output;
  char o1[100],o2[100],o3[100];

sprintf(o1,"%s%s",argv[2],".sam");
sprintf(o2,"%s",argv[3]);
if(strcasecmp(o2,"no")==0){shift=0;}else{shift=1;}

output=fopen(o1,"w");

FILE *fin = fopen ( argv[1], "r" );

while(fscanf(fin,"%[^\n]\n", buffer) != EOF){    

sscanf (buffer,"%s\n",&head);
sscanf (buffer,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%[^\n]\n",qname,flag,rname,position,mapq,cigar,rnext,pnext,tlen,sequence,qual,tags);

if(strncmp(head,"@",1)==0){fprintf(output,"%s\n",buffer); }else{
strcpy(rname2,rname);
            memmove(rname2,"   ",3);
            chrom=atoi(rname2);

            if(strcmp(rname2,"   X")==0){ chrom=23;}
            if(strcmp(rname2,"   Y")==0){ chrom=24;}
            if(strcmp(rname2,"   M")==0){ chrom=25;}

        quality=atoi(mapq);
        count=0;
        if((strcmp(flag,"16") ==0 ||strcmp(flag,"0") ==0)&&chrom!=0){
        char* seq=sequence;
        char* c=cigar;

        while(*c!=0)
            {
            char* p2;
                if(!isdigit(*c)) {       fprintf(stderr,"bad cigar string: %s\n",*cigar); return 0; }
            int n=strtol(c,&p2,10);
                if(n<=0){                fprintf(stderr,"bad cigar string: %s\n",*cigar); return 0; }
                switch(*p2)
                {
                case 'M':{ while(n>0) { count++; --n; } break; }
                case 'D':{ while(n>0) { count++; --n; } break; }
                case 'I':{ while(n>0) { --n; } break; }
                case 'P':{ break;}default: { fprintf(stderr,"unsupported operator %c\n",*p2);  return 0; }} c=p2+1; }

        /* locate the position of the read, given the fwd/bkw direction and cigar string info 

        printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",qname,flag,rname,position,mapq,cigar,rnext,pnext,tlen,sequence,qual);
*/
/*printf("%s\n%d\n%s\n%d\n%d\n%s\n%s\n%d\n%d\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",name0,f1,name1,f2,f3,name2,c6,f6,f7,name3,name4,tag1,tag2,tag3,tag4,tag5,tag6,tag7,tag8,tag9);}
   */

if(strcmp(flag,"16") ==0 ){position2=atoi(position)+shift*count; }else{position2=atoi(position);}

fprintf(output,"%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",qname,flag,rname,position2,mapq,cigar,rnext,pnext,tlen,sequence,qual,tags);
/*fprintf(output,"%s\tGM:i:%d\n",buffer,count);*/ }
 }
}
}
