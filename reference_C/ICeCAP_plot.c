/* cc -g plot.c -L/usr/X11R6/lib -lX11 -lm */
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <ctype.h>
#include <unistd.h>
#include <dirent.h> 
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/Xproto.h>
#include <string.h>

#define max(a,b) (a>b?a:b)
#define min(a,b) (a<b?a:b)

FILE *heatmap,*tads,*chic,*chic2,*chic3,*chic4,*chic5,*chipseq,*chipseq2,*chipseq3,*chipseq4,*chromhmm,*fourcseq,*genes,*gwas;
char* concat(char *s1, char *s2)
{
    char *result = malloc(strlen(s1)+strlen(s2)+1);
     strcpy(result, s1);
     strcat(result, s2);
     return result;
}

   int * HSL2RGB(double h, double sl, double l)
      {
            double v;
            double r,g,b;
            static int rgb[10];
            int temp1,temp2,temp3; 
            r = l;   // default to gray
            g = l;
            b = l;
            v = (l <= 0.5) ? (l * (1.0 + sl)) : (l + sl - l * sl);

            if (v > 0)
            {
                  double temp;
                  double m;
                  double c;
                  double x;
		  double mm;
 
                  mm = l + l - 1;
                  c = (1-fabs(mm))*sl;
		  m=l-c/2;

                  if(h>270&&h<=360){
                              r = c+m;
                              g = 0+m;
                              b = c*(h-270)/60+m;
                             } 
                       if(h>=0&&h<=60){
                              r = c+m;
                              g = c*h/60+m;
                              b = 0+m;
                             }
                        if(h>60&&h<=120){

                              r = c*(h-60)/60+m;
                              g = c+m;
                              b = 0+m;
                             } 
                         if(h>120&&h<=180){
                              r = 0+m;
                              g = c+m;
                              b = c*(h-120)/60+m;
                             }
                         if(h>180&&h<=240){
                              r = 0+m;
                              g = c*(h-180)/60+m;
                              b = c+m;
                             }
                         if(h>240&&h<=270){
                              r = c*(h-240)/60+m;
                              g = 0+m;
                              b = c+m;
                             }
}

        temp1=(int)(r*255);
        temp2=(int)(g*255);
        temp3=(int)(b*255);

            rgb[1]=temp1;
            rgb[2]=temp2;
            rgb[3]=temp3;
            return rgb;
} 
void RemoveChars(char *s, char c)
{
    int writer = 0, reader = 0;

    while (s[reader])
    {
        if (s[reader]!=c)
        {
            s[writer++] = s[reader];
        }

        reader++;
    }

    s[writer]=0;
}

int main(int argc, char **argv)
{

  int option,i,ii,iii,jj,l,j,f,scale,scale2,scaley,xstart;
  int box1,box2,box3,box4,box5,box6,box7,height,heightsave,heighttag,finalheight,resolution,count;
  int flag,num,viewpointa,viewpointb,check,check2,scalechic,grid,tot,generaw,dist;
  int blackColor,whiteColor,thescreen,thedisplayfont,position;
  int r,g,b,maximum,greyscale;
  int xticks,tmp,bingene,bingene2,taggene,taggene2,intensity,intensity2,buffer;
  int heightgene[10000][100];
  int heighttaggene[10000][100];
  int total[10000],savestart[1000],saveend[1000];

  int *rgb;
  int strindex(char source[], char searchfor[]);
  float brightness,hm,hmm,fourcreads,maximumf,rd,gn,be,dn,all[10000],allb[10000];
  char s[64],red[10],blue[10],green[10];
  char strand[10],start[10],end[10],chromosome[10],*chrom,chroms[10];
  char eS[1000],eE[1000],gS[1000],gT[1000],target[1000],*path,*targetinfo,track[100],track2[100],track3[100],track4[100],track4b[100],track5[100],track5b[100],track6[100],track7b[100],track7[100],track8[100],track10[100];

  const char delims[] = ",";
    const char * fontname = "-*-helvetica-*-r-*-*-34-*-*-*-*-*-*-*";

   XTextProperty textproperty;
   Display *thedisplay;
   GC thecontext;
   XEvent anevent;
   XColor xcolour,draw;
   XPoint v[3];

Colormap thecolormap;
   Window thewindow;

   thedisplay  = XOpenDisplay(NULL);
   blackColor  = BlackPixel(thedisplay,DefaultScreen(thedisplay));
   whiteColor  = WhitePixel(thedisplay,DefaultScreen(thedisplay));
   thescreen   = DefaultScreen(thedisplay);
   thecolormap = DefaultColormap(thedisplay,thescreen);

 /* Create the window */
   thewindow = XCreateSimpleWindow(thedisplay, 
   DefaultRootWindow(thedisplay),0,0, 
   1000,800,0,whiteColor,whiteColor);
   XSelectInput(thedisplay,thewindow,StructureNotifyMask);
   XMapWindow(thedisplay,thewindow);

   /* Label the window */
   XStringListToTextProperty(&(argv[0]),1,&textproperty);
   XSetWMName(thedisplay,thewindow,&textproperty);

   /* Get the context */
   thecontext = XCreateGC(thedisplay,thewindow,0,NULL);
   XSetBackground(thedisplay,thecontext,whiteColor);
   XSetForeground(thedisplay,thecontext,whiteColor);
   /*XSetFont(thedisplay,thecontext,font->fid);
  */ 
char sizes[10000],chpos[12],chpos2[12];
int chr,pos,chr2,pos2,posb,pos2b,pos3,pos3b,pos4,pos5,pos5b,pos6,pos7,pos8,pos9,pos10,excount,delta;
int marginx,marginy;

 /* Erase the display (In the background colour) */

   for (;;) {
      XNextEvent(thedisplay, &anevent);
      if (anevent.type == MapNotify)
         break; }
/* Default Position being the MYC gene locus */
chrom="8";
position=128740000;
buffer=1000000;
targetinfo="MYC";
path=".";

XClearWindow(thedisplay,thewindow); 
XFlush(thedisplay);

maximumf=0;
box1=250000000;
box2=0;
box3=250000000;
box4=0;
delta=10;
char command[200];

  if(argc ==1){
        fprintf(stderr, "This program needs arguments .... \n -c chrom -p position -b buffer -g target\n E.g. -c 8 -p 128740000 -b 1000000 -g MYC \n");}

   while((option = getopt(argc, argv, "c:p:b:g:")) != -1)

        switch (option)
        {
        case 'c' :           chrom=optarg;                  break;
        case 'p' :           position=atoi(optarg);         break;
        case 'b' :           buffer=atoi(optarg);           break;
        case 'g' :           targetinfo=optarg;             break; }

/* make sure any interaction file resolution is such that res/resolution*grid is an integer*/

marginx=100;
resolution=(int)((double)(buffer)/(double)500);
scale=(int)((10000)/(double)(resolution));

scaley=1.0*scale;

box1=marginx;
box3=position-1.0*buffer;
box4=position+1.0*buffer;
pos4=((box4-box3)/resolution);
box2=2*scale+marginx+2*pos4;

marginy=50*(double)(scaley)/(double)scale;
xticks=(double)(box4-box3)/(double)(5*resolution)+1;
char *hmap=concat(targetinfo,"_hm5kb");
heatmap=fopen(hmap,"rt");
if(heatmap != NULL){
while(fscanf(heatmap,"%s %d %s %d %f",&chromosome,&pos,&chromosome,&pos2,&hm)!=EOF)
{
if((hm)>maximumf){maximumf=(hm);}
}
rewind(heatmap);

/*scale;*/
maximum=2;
box1=marginx;
box2=2*scale+marginx+2*pos4;

while(fscanf(heatmap,"%s %d %s %d %f",&chromosome,&pos,&chromosome,&pos2,&hm)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
pos3=((pos-box3)/resolution);
pos4=((pos2-box3)/resolution);

maximumf=2;
brightness=1;
if(hm>=1&&hm<maximumf){hmm=(log(hm)/log(maximumf));}
if(hm>=maximumf){hmm=1;}

xcolour.red=255*255;
xcolour.green=255*255-255*255*pow(hmm,0.5);
xcolour.blue=xcolour.green;

if(hm<1){
hmm=(log(1+fabs(hm-1)));
if(hmm>=1){hmm=0.99;}
xcolour.blue=255*255;
xcolour.green=255*255-255*255*pow(hmm,0.5);
xcolour.red=xcolour.green;
}
xcolour.flags = DoRed | DoGreen | DoBlue;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);


   v[0].y =-fabs(pos3-pos4)*(double)(scaley)/(double)(scale)+scaley+marginy; v[0].x =(pos3+pos4)+2*scale+marginx;
   v[1].y =-fabs(pos3-pos4)*(double)(scaley)/(double)(scale)+2*scaley+marginy;v[1].x =(pos3+pos4)+scale+marginx;
   v[2].y =-fabs(pos3-pos4)*(double)(scaley)/(double)(scale)+scaley+marginy; v[2].x = (pos3+pos4)+0+marginx;
   v[3].y =-fabs(pos3-pos4)*(double)(scaley)/(double)(scale)+0+marginy; v[3].x = (pos3+pos4)+scale+marginx;

  /* if(pos4<pos3){*/
   xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
   XSetForeground(thedisplay,thecontext,blackColor);

   v[0].y =-fabs(pos3-pos4)*(double)(scaley)/(double)(scale)+scaley+marginy; v[0].x =(pos3+pos4)+2*scale+marginx; 
   v[3].y =-fabs(pos3-pos4)*(double)(scaley)/(double)(scale)+2*scaley+marginy;v[1].x =(pos3+pos4)+scale+marginx; 
   v[2].y =-fabs(pos3-pos4)*(double)(scaley)/(double)(scale)+scaley+marginy; v[2].x = (pos3+pos4)+0+marginx;
   v[1].y =-fabs(pos3-pos4)*(double)(scaley)/(double)(scale)+0+marginy; v[3].x = (pos3+pos4)+scale+marginx; 


   if(pos4==pos3){
   xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
   XFillPolygon(thedisplay,thewindow,thecontext,v,3,Complex,CoordModeOrigin);
   XSetForeground(thedisplay,thecontext,blackColor);
}
}
}
fclose(heatmap);

}

/*TADS*/
		
sprintf(chromosome,"%s",chrom); 
char *tad=concat("../CHIP/TADS_chr",chromosome);

tads=fopen(tad,"rt");
if(tads != NULL){

while(fscanf(tads,"%i\t%i\t%i",&chr,&pos,&pos2)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
  pos3=((pos-box3)/resolution);
  pos4=((pos2-box3)/resolution);
  scale2=(pos4-pos3+scale);


xcolour.red =0*255; xcolour.green =0*255; xcolour.blue =0*255;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);

v[0].y =scaley+marginy;   v[0].x =(2*pos3)+0.5*scale+marginx;
v[1].y =scaley+marginy;   v[1].x = (2*pos3)+2*scale2-0.5*scale+marginx;
v[2].y =1.5*scaley-scale2*(double)(scaley)/(double)(scale)+marginy;        v[2].x = (2*pos3)+scale2+marginx;

for (i=0;i<=2;i++) {
      XDrawLine(thedisplay,thewindow,thecontext,
         v[i].x,v[i].y,v[(i+1)%3].x,v[(i+1)%3].y); }

v[0].y =scaley+marginy;   v[0].x =(2*pos3)-0.5*scale+marginx;
v[1].y =scaley+marginy;   v[1].x = (2*pos3)+2*scale2+0.5*scale+marginx;
v[2].y =0.5*scaley-scale2*(double)(scaley)/(double)(scale)+marginy;        v[2].x = (2*pos3)+scale2+marginx;

for (i=0;i<=2;i++) {
      XDrawLine(thedisplay,thewindow,thecontext,
         v[i].x,v[i].y,v[(i+1)%3].x,v[(i+1)%3].y); }
}
}
XFlush(thedisplay);
fclose(tads);


height=marginy+15*scaley;}else{height=marginy+1*scaley;}
XSetForeground(thedisplay,thecontext,blackColor);

char *chics=concat(targetinfo,"_chic");

chic=fopen(chics,"rt");
if(chic != NULL){
char track[100]="CHiCAGO";

XDrawString(thedisplay,thewindow,thecontext,marginx/4,height-scaley,track,strlen(track));
while(fscanf(chic,"%i\t%i\t%i\t%i\t%i",&chr,&pos,&pos2,&flag,&num)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
printf("%d\t%d\t%d\t%d\t%d\n",chr,pos,pos2,flag,num);
  check=(resolution);
  scalechic=check/resolution;
  pos3=(double)((pos-box3)/check)*scalechic;
  pos4=(double)((pos2-box3)/check)*scalechic;
if(flag==0){height=height+2*scaley; xcolour.red =0; xcolour.green =0; xcolour.blue=50000; viewpointb=(2*pos4)+marginx; viewpointa=(2*pos3)+marginx;}
if(flag==1){xcolour.red =50000; xcolour.green =0;
xcolour.blue=0;
xcolour.flags = DoRed | DoGreen | DoBlue;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);
if(viewpointa<(2*pos3)+marginx){
XDrawLine(thedisplay,thewindow,thecontext,viewpointb,height-1*scaley,(2*pos3)+marginx,height-1*scaley);}else
{XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx,height-1*scaley,viewpointa,height-1*scaley);}}
if(viewpointa<(2*pos3)+marginx){
   v[0].y =(height);       v[0].x=(2*pos4)+marginx;
   v[1].y =(height);       v[1].x=(2*pos3)+marginx;
   v[2].y =(height-2*scaley); v[2].x=(2*pos3)+marginx;
   v[3].y =(height-2*scaley); v[3].x=(2*pos4)+marginx;
   xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
   XSetForeground(thedisplay,thecontext,blackColor);
}
else{
   v[0].y =(height);       v[0].x=(2*pos3)+marginx;
   v[1].y =(height);       v[1].x=(2*pos4)+marginx;
   v[2].y =(height-2*scaley); v[2].x=(2*pos4)+marginx;
   v[3].y =(height-2*scaley); v[3].x=(2*pos3)+marginx;
   xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
   XSetForeground(thedisplay,thecontext,blackColor);
}
}
}
XFlush(thedisplay);
fclose(chic);
height=height+5*scaley;
}

XSetForeground(thedisplay,thecontext,blackColor);

XFlush(thedisplay);

char *chics2=concat(targetinfo,"_chic2");
chic2=fopen(chics2,"rt");
if(chic2 != NULL){
char track[100]="ICeCAP";

XDrawString(thedisplay,thewindow,thecontext,marginx/4,height+5*scaley,track,strlen(track));

while(fscanf(chic2,"%i\t%i\t%i\t%i\t%i",&chr,&pos,&pos2,&flag,&num)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
  check=(resolution);
  scalechic=check/resolution;
  pos3=(double)((pos-box3)/check)*scalechic;
  pos4=(double)((pos2-box3)/check)*scalechic;
if(flag==0){height=height+2*scaley; xcolour.red =0; xcolour.green =0; xcolour.blue=50000; viewpointb=(2*pos4)+marginx; viewpointa=(2*pos3)+marginx;}
if(flag==1){xcolour.red =50000; xcolour.green =0;
xcolour.blue=0;
xcolour.flags = DoRed | DoGreen | DoBlue;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);
if(viewpointa<(2*pos3)+marginx){
XDrawLine(thedisplay,thewindow,thecontext,viewpointb,height-1*scaley,(2*pos3)+marginx,height-1*scaley);}else
{XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx,height-1*scaley,viewpointa,height-1*scaley);}}
if(viewpointa<(2*pos3)+marginx){
   v[0].y =(height);       v[0].x=(2*pos4)+marginx;
   v[1].y =(height);       v[1].x=(2*pos3)+marginx;
   v[2].y =(height-2*scaley); v[2].x=(2*pos3)+marginx;
   v[3].y =(height-2*scaley); v[3].x=(2*pos4)+marginx;
   xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
   XSetForeground(thedisplay,thecontext,blackColor);
}
else{
   v[0].y =(height);       v[0].x=(2*pos3)+marginx;
   v[1].y =(height);       v[1].x=(2*pos4)+marginx;
   v[2].y =(height-2*scaley); v[2].x=(2*pos4)+marginx;
   v[3].y =(height-2*scaley); v[3].x=(2*pos3)+marginx;
   xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
   XSetForeground(thedisplay,thecontext,blackColor);
}
}
}
XFlush(thedisplay);
fclose(chic2);
}
XSetForeground(thedisplay,thecontext,blackColor);


XFlush(thedisplay);


char *chics3=concat(targetinfo,"_chic3");

chic3=fopen(chics3,"rt");
if(chic3 != NULL){
char track[100]="";

XDrawString(thedisplay,thewindow,thecontext,marginx/4,height+4*scaley,track,strlen(track));
while(fscanf(chic3,"%i\t%i\t%i\t%i\t%i",&chr,&pos,&pos2,&flag,&num)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
  check=(resolution);
  scalechic=check/resolution;
  pos3=(double)((pos-box3)/check)*scalechic;
  pos4=(double)((pos2-box3)/check)*scalechic;
if(flag==0){height=height+2*scaley; xcolour.red =0; xcolour.green =0; xcolour.blue=50000; viewpointb=(2*pos4)+marginx; viewpointa=(2*pos3)+marginx;}
if(flag==1){xcolour.red =50000; xcolour.green =0;
xcolour.blue=0;
xcolour.flags = DoRed | DoGreen | DoBlue;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);
if(viewpointa<(2*pos3)+marginx){
XDrawLine(thedisplay,thewindow,thecontext,viewpointb,height-0.5*scaley,(2*pos3)+marginx,height-0.5*scaley);}else
{XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx,height-0.5*scaley,viewpointa,height-0.5*scaley);}}
if(viewpointa<(2*pos3)+marginx){
   v[0].y =(height);       v[0].x=(2*pos4)+marginx;
   v[1].y =(height);       v[1].x=(2*pos3)+marginx;
   v[2].y =(height-2*scaley); v[2].x=(2*pos3)+marginx;
   v[3].y =(height-2*scaley); v[3].x=(2*pos4)+marginx;
   xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
   XSetForeground(thedisplay,thecontext,blackColor);
}
else{
   v[0].y =(height);       v[0].x=(2*pos3)+marginx;
   v[1].y =(height);       v[1].x=(2*pos4)+marginx;
   v[2].y =(height-2*scaley); v[2].x=(2*pos4)+marginx;
   v[3].y =(height-2*scaley); v[3].x=(2*pos3)+marginx;
   xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
   XSetForeground(thedisplay,thecontext,blackColor);
}
}
}
XFlush(thedisplay);
fclose(chic3);
}

XSetForeground(thedisplay,thecontext,blackColor);



XFlush(thedisplay);


char *chics4=concat(targetinfo,"_chic4");

chic4=fopen(chics4,"rt");
if(chic4 != NULL){
char track[100]="";

XDrawString(thedisplay,thewindow,thecontext,marginx/4,height-scaley,track,strlen(track));
while(fscanf(chic4,"%i\t%i\t%i\t%i\t%i",&chr,&pos,&pos2,&flag,&num)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
  check=(resolution);
  scalechic=check/resolution;
  pos3=(double)((pos-box3)/check)*scalechic;
  pos4=(double)((pos2-box3)/check)*scalechic;
if(flag==0){height=height+2*scaley; xcolour.red =0; xcolour.green =0; xcolour.blue=50000; viewpointb=(2*pos4)+marginx; viewpointa=(2*pos3)+marginx;}
if(flag==1){xcolour.red =50000; xcolour.green =0;
xcolour.blue=0;
xcolour.flags = DoRed | DoGreen | DoBlue;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);
if(viewpointa<(2*pos3)+marginx){
XDrawLine(thedisplay,thewindow,thecontext,viewpointb,height-1*scaley,(2*pos3)+marginx,height-1*scaley);}else
{XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx,height-1*scaley,viewpointa,height-1*scaley);}}
if(viewpointa<(2*pos3)+marginx){
   v[0].y =(height);       v[0].x=(2*pos4)+marginx;
   v[1].y =(height);       v[1].x=(2*pos3)+marginx;
   v[2].y =(height-2*scaley); v[2].x=(2*pos3)+marginx;
   v[3].y =(height-2*scaley); v[3].x=(2*pos4)+marginx;
   xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
   XSetForeground(thedisplay,thecontext,blackColor);
}
else{
   v[0].y =(height);       v[0].x=(2*pos3)+marginx;
   v[1].y =(height);       v[1].x=(2*pos4)+marginx;
   v[2].y =(height-2*scaley); v[2].x=(2*pos4)+marginx;
   v[3].y =(height-2*scaley); v[3].x=(2*pos3)+marginx;
   xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
   XSetForeground(thedisplay,thecontext,blackColor);
}
}
}
XFlush(thedisplay);
height=height+5*scaley;
fclose(chic4);
}

XSetForeground(thedisplay,thecontext,blackColor);

XFlush(thedisplay);

char *chics5=concat(targetinfo,"_chic5");

chic5=fopen(chics5,"rt");
if(chic5 != NULL){
char track[100]="ICeCAP HT29";

XDrawString(thedisplay,thewindow,thecontext,marginx/4,height+scaley,track,strlen(track));
while(fscanf(chic5,"%i\t%i\t%i\t%i\t%i",&chr,&pos,&pos2,&flag,&num)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
  check=(resolution);
  scalechic=check/resolution;
  pos3=(double)((pos-box3)/check)*scalechic;
  pos4=(double)((pos2-box3)/check)*scalechic;
if(flag==0){height=height+scaley; xcolour.red =0; xcolour.green =0; xcolour.blue=50000; viewpointb=(2*pos4)+marginx; viewpointa=(2*pos3)+marginx;}
if(flag==1){xcolour.red =50000; xcolour.green =0;
xcolour.blue=0;
xcolour.flags = DoRed | DoGreen | DoBlue;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);
if(viewpointa<(2*pos3)+marginx){
XDrawLine(thedisplay,thewindow,thecontext,viewpointb,height-0.5*scaley,(2*pos3)+marginx,height-0.5*scaley);}else
{XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx,height-0.5*scaley,viewpointa,height-0.5*scaley);}}
if(viewpointa<(2*pos3)+marginx){
   v[0].y =(height);       v[0].x=(2*pos4)+marginx;
   v[1].y =(height);       v[1].x=(2*pos3)+marginx;
   v[2].y =(height-scaley); v[2].x=(2*pos3)+marginx;
   v[3].y =(height-scaley); v[3].x=(2*pos4)+marginx;
   xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
   XSetForeground(thedisplay,thecontext,blackColor);
}
else{
   v[0].y =(height);       v[0].x=(2*pos3)+marginx;
   v[1].y =(height);       v[1].x=(2*pos4)+marginx;
   v[2].y =(height-scaley); v[2].x=(2*pos4)+marginx;
   v[3].y =(height-scaley); v[3].x=(2*pos3)+marginx;
   xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
   XSetForeground(thedisplay,thecontext,blackColor);
}
}
}
height=height+5*scaley;
XFlush(thedisplay);
fclose(chic5);
}


maximumf=0.0;
xcolour.red =0; xcolour.green =0; xcolour.blue = 0;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);

/* 4C */

char *ccccb=concat(targetinfo,"_4C");

fourcseq=fopen(ccccb,"rt");
if(fourcseq != NULL){
char track3[100]="Hi-C";
XSetForeground(thedisplay,thecontext,blackColor);
XDrawString(thedisplay,thewindow,thecontext,marginx/4,height-scaley,track3,strlen(track3));
while(fscanf(fourcseq,"%s %d %s %d %f",&chromosome,&pos2,&chromosome,&pos,&fourcreads)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
if(sqrt((pos-position)*(pos-position))>25000){
  if(fourcreads>maximumf){maximumf=fourcreads;}}
 pos3=(int)((double)(pos-box3)/(double)resolution);
  pos4=(int)((double)(pos+5000-box3)/(double)resolution);
  allb[pos3]=allb[pos3]+fourcreads;
}
}
maximumf=10;
rewind(fourcseq);
while(fscanf(fourcseq,"%s %d %s %d %f",&chromosome,&pos2,&chromosome,&pos,&fourcreads)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
 pos3=(int)((double)(pos-box3)/(double)resolution);
 pos4=(int)((double)(pos+5000-box3)/(double)resolution);
 allb[pos3]=allb[pos3]+fourcreads;

if(sqrt((pos-position)*(pos-position))>25000){
if(allb[pos3]>maximumf){maximumf=allb[pos3];}}
}
}

rewind(fourcseq);

while(fscanf(fourcseq,"%s %d %s %d %f",&chromosome,&pos2,&chromosome,&pos,&fourcreads)!=EOF)
 {if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
  check=(pos2-pos);
  scalechic=(int)((double)check/(double)resolution);
  pos3=(int)((double)(pos-box3)/(double)resolution);
  pos4=(int)((double)(pos+5000-box3)/(double)resolution);
if(allb[pos3]<=maximumf){
if(pos3!=pos4){
   v[0].y =(height);       v[0].x=(2*pos3)+marginx;
   v[1].y =(height);       v[1].x=(2*pos4)+marginx;
   v[2].y =(height-15*(allb[pos3]/maximumf)*scaley); v[2].x=(2*pos4)+marginx;
   v[3].y =(height-15*(allb[pos3]/maximumf)*scaley); v[3].x=(2*pos3)+marginx;

   /*xcolour.flags = DoRed | DoGreen | DoBlue;
 *    XAllocColor(thedisplay,thecolormap,&xcolour);
 *       XSetForeground(thedisplay,thecontext,xcolour.pixel);*/
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
   /*XSetForeground(thedisplay,thecontext,blackColor);*/
}else{
 XSetForeground(thedisplay,thecontext,blackColor);
 XDrawLine(thedisplay,thewindow,thecontext,(2*pos4)+marginx,height,(2*pos4)+marginx,height-10*(double)(fourcreads/maximumf)*scaley);

}
}
}
}
fclose(fourcseq);
XDrawLine(thedisplay,thewindow,thecontext,box1,height,box2,height);
 }

/* 4C ViewPoint*/
pos=position;
pos2=position+5000;
check=(pos2-pos);
scalechic=(int)((double)check/(double)resolution);
pos3=(int)((double)(pos-box3)/(double)resolution);
pos4=(int)((double)(pos2-box3)/(double)resolution);

   v[0].y =(height);       v[0].x=(2*pos3)+marginx;
   v[1].y =(height);       v[1].x=(2*pos4)+marginx;
   v[2].y =(height-6*scaley); v[2].x=(2*pos4)+marginx;
   v[3].y =(height-6*scaley); v[3].x=(2*pos3)+marginx;

xcolour.red =0; xcolour.green =255*257; xcolour.blue = 0;

/* xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
 XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
XFlush(thedisplay);

height=height+15*scaley;
maximumf=0.0;

XSetForeground(thedisplay,thecontext,blackColor);
XDrawLine(thedisplay,thewindow,thecontext,box1,height,box2,height);
*/

/* 4C */

char *ccccd=concat(targetinfo,"_4C2");

fourcseq=fopen(ccccd,"rt");
if(fourcseq != NULL){
char track3[100]="4C";
XSetForeground(thedisplay,thecontext,blackColor);
XDrawString(thedisplay,thewindow,thecontext,marginx/4,height-scaley,track3,strlen(track3));
while(fscanf(fourcseq,"%s %d %s %d %f",&chromosome,&pos2,&chromosome,&pos,&fourcreads)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
if(sqrt((pos-position)*(pos-position))>25000){
  if(fourcreads>maximumf){maximumf=fourcreads;}}
 pos3=(int)((double)(pos-box3)/(double)resolution);
  pos4=(int)((double)(pos+5000-box3)/(double)resolution);
  allb[pos3]=allb[pos3]+fourcreads;
}
}
maximumf=20;
rewind(fourcseq);
while(fscanf(fourcseq,"%s %d %s %d %f",&chromosome,&pos2,&chromosome,&pos,&fourcreads)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
 pos3=(int)((double)(pos-box3)/(double)resolution);
  pos4=(int)((double)(pos+5000-box3)/(double)resolution);
  allb[pos3]=allb[pos3]+fourcreads;

if(sqrt((pos-position)*(pos-position))>25000){
  if(allb[pos3]>maximumf){maximumf=allb[pos3];}}
}
}
rewind(fourcseq);

while(fscanf(fourcseq,"%s %d %s %d %f",&chromosome,&pos2,&chromosome,&pos,&fourcreads)!=EOF)
 {if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
  check=(pos2-pos);
  scalechic=(int)((double)check/(double)resolution);
  pos3=(int)((double)(pos-box3)/(double)resolution);
  pos4=(int)((double)(pos+5000-box3)/(double)resolution);
if(allb[pos3]<=maximumf){
if(pos3!=pos4){
   v[0].y =(height);       v[0].x=(2*pos3)+marginx;
   v[1].y =(height);       v[1].x=(2*pos4)+marginx;
   v[2].y =(height-15*(allb[pos3]/maximumf)*scaley); v[2].x=(2*pos4)+marginx;
   v[3].y =(height-15*(allb[pos3]/maximumf)*scaley); v[3].x=(2*pos3)+marginx;

   /*xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);*/
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
   /*XSetForeground(thedisplay,thecontext,blackColor);*/
}else{
 XSetForeground(thedisplay,thecontext,blackColor);
 XDrawLine(thedisplay,thewindow,thecontext,(2*pos4)+marginx,height,(2*pos4)+marginx,height-10*(double)(fourcreads/maximumf)*scaley);

}
}
}
}
fclose(fourcseq);
XDrawLine(thedisplay,thewindow,thecontext,box1,height,box2,height);
 

/* 4C ViewPoint*/
pos=position-53555;
pos2=position-53555+5000;
check=(pos2-pos);
scalechic=(int)((double)check/(double)resolution);
pos3=(int)((double)(pos-box3)/(double)resolution);
pos4=(int)((double)(pos2-box3)/(double)resolution);

   v[0].y =(height);       v[0].x=(2*pos3)+marginx;
   v[1].y =(height);       v[1].x=(2*pos4)+marginx;
   v[2].y =(height-6*scaley); v[2].x=(2*pos4)+marginx;
   v[3].y =(height-6*scaley); v[3].x=(2*pos3)+marginx;

xcolour.red =0; xcolour.green =255*257; xcolour.blue = 0;

 xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
 XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
XFlush(thedisplay);

}
maximumf=0.0;

XSetForeground(thedisplay,thecontext,blackColor);


XFlush(thedisplay);

maximum=0;
intensity=1;

/* PROMOTORS */

xcolour.red =255*257; xcolour.green =0; xcolour.blue = 0;

 xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
sprintf(chromosome,"%s",chrom);
char *prom=concat("./lane/reference_HiC/FRAGMENTS_CHR",chromosome);
chipseq3=fopen(prom,"rt");
if(chipseq3 != NULL){
char track4[100]="Promoter Fragments";


XDrawString(thedisplay,thewindow,thecontext,marginx/4,height-scaley,track4,strlen(track4));
while(fscanf(chipseq3,"%s\t%i\t%i",&chromosome,&pos,&pos2)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
 pos3=((double)(pos-box3)/(double)resolution);
 pos4=((double)(pos2-box3)/(double)resolution);
total[pos3]+=(intensity);
if(total[pos3]>maximum){maximum=total[pos3];}

}
}
rewind(chipseq3);

intensity=1.;
while(fscanf(chipseq3,"%s\t%i\t%i",&chromosome,&pos,&pos2)!=EOF)
{

if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
scalechic=check/resolution;
pos3=((double)(pos-box3)/(double)resolution);
pos4=((double)(pos2-box3)/(double)resolution);

   v[0].y =(height);       v[0].x=(2*(pos3-1))+marginx;
   v[1].y =(height);       v[1].x=(2*(pos4+1))+marginx;
   v[2].y =(height-scaley); v[2].x=(2*(pos4+1))+marginx;
   v[3].y =(height-scaley); v[3].x=(2*(pos3-1))+marginx;

   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);

}
}
fclose(chipseq3);
height=height+8*scaley;
}
XFlush(thedisplay);

intensity2=2;

/* OMOTORS */

xcolour.red =0; xcolour.green =0; xcolour.blue = 255*257;

 xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);

char *prom3=concat("./lane/reference_HiC/FRAGMENTS_",chromosome);
chipseq3=fopen(prom3,"rt");
if(chipseq3 != NULL){
char track4c[100]="5kb Promoter Bins";

XDrawString(thedisplay,thewindow,thecontext,marginx/4,height-scaley,track4c,strlen(track4c));

while(fscanf(chipseq3,"%s\t%i\t%i",&chromosome,&pos,&pos2)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
 pos3=((double)(pos-box3)/(double)resolution);
 pos4=((double)(pos2-box3)/(double)resolution);
total[pos3]+=(intensity);
if(total[pos3]>maximum){maximum=total[pos3];}

}
 }
rewind(chipseq3);

intensity=1.;
while(fscanf(chipseq3,"%s\t%i\t%i\t%i\t%i",&chromosome,&pos,&pos2,&intensity)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
scalechic=check/resolution;
pos3=((double)(pos-box3)/(double)resolution);
pos4=((double)(pos2-box3)/(double)resolution);

   v[0].y =(height);       v[0].x=(2*(pos3-1))+marginx;
   v[1].y =(height);       v[1].x=(2*(pos4+1))+marginx;
   v[2].y =(height-scaley); v[2].x=(2*(pos4+1))+marginx;
   v[3].y =(height-scaley); v[3].x=(2*(pos3-1))+marginx;

   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);

}
}
fclose(chipseq3);
memset(total, 0, sizeof(total)); 
height=height+5*scaley;
}
XFlush(thedisplay);

intensity2=2;
intensity=1.;
XSetForeground(thedisplay,thecontext,blackColor);
XDrawLine(thedisplay,thewindow,thecontext,box1,height,box2,height);

/* RF */

xcolour.green =255*257; xcolour.blue =0; xcolour.red = 0;

 xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
sprintf(chromosome,"%s",chrom);
char *dnase7=concat("MboI_chr",chromosome);
chipseq4=fopen(dnase7,"rt");
if(chipseq4 != NULL){
char track5d[100]="MboI";
XDrawString(thedisplay,thewindow,thecontext,marginx/4,height-scaley,track5d,strlen(track5d));


while(fscanf(chipseq4,"%s\t%i\t%i\t%i\t%i",&chromosome,&pos,&pos2,&intensity)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
 pos3=((double)(pos-box3)/(double)resolution);
 pos4=((double)(pos2-box3)/(double)resolution);
 total[pos3]+=(intensity);
if(total[pos3]>maximum){maximum=total[pos3];}

}
}

rewind(chipseq4);

intensity=0.;
while(fscanf(chipseq4,"%s\t%i\t%i\t%i\t%i",&chromosome,&pos,&pos2,&intensity)!=EOF)
{if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
check=(1000);
scalechic=check/resolution;
pos3=((double)(pos-box3)/(double)resolution);
pos4=((double)(pos2-box3)/(double)resolution);
intensity=(total[pos3]);
XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx,height,(2*pos3)+marginx,height-2*scaley);
}
}
fclose(chipseq4);
memset(total, 0, sizeof(total));
height=height+5*scaley;
}
XFlush(thedisplay);
intensity=2;


/* DNASE */
height=height+5*scaley;

xcolour.blue =0; xcolour.green =0; xcolour.red = 0;

xcolour.flags = DoRed | DoGreen | DoBlue;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);
sprintf(chromosome,"%s",chrom);
char *dnase5=concat("../CHIP/DNASE_chr",chromosome);

chipseq4=fopen(dnase5,"rt");
if(chipseq4 != NULL){
char track5[100]="DNASE";

XDrawString(thedisplay,thewindow,thecontext,marginx/4,height-scaley,track5,strlen(track5));
XSetForeground(thedisplay,thecontext,blackColor);
XDrawString(thedisplay,thewindow,thecontext,marginx/4,height-scaley,track6,strlen(track6));
while(fscanf(chipseq4,"%s\t%i\t%i\t%i\t%i",&chromosome,&pos,&pos2,&intensity)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
 pos3=((double)(pos-box3)/(double)resolution);
 pos4=((double)(pos2-box3)/(double)resolution);
 total[pos3]+=(intensity);
if(total[pos3]>maximum){maximum=total[pos3];}
}
 }
rewind(chipseq4);

intensity=0.;
check=(1000);

 for(pos3=1;pos3<=1000;pos3++){
   pos4=pos3+1;
intensity=(total[pos3]);
rd=1255*257*((double)(maximum-total[pos3])/(double)(maximum));
gn=1255*257*((double)(maximum-total[pos3])/(double)(maximum));
be=1255*257*((double)(maximum-total[pos3])/(double)(maximum));

xcolour.red =rd; xcolour.green =gn; xcolour.blue = be;

xcolour.flags = DoRed | DoGreen | DoBlue;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);
if(pos3!=pos4){
   v[0].y =(height);       v[0].x=(2*pos3)+marginx;
   v[1].y =(height);       v[1].x=(2*pos4)+marginx;
   v[2].y =(height-2*((double)(intensity2))*scaley); v[2].x=(2*pos4)+marginx;
   v[3].y =(height-2*((double)(intensity2))*scaley); v[3].x=(2*pos3)+marginx;
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
}else{
/*XDrawLine(thedisplay,thewindow,thecontext,(2*pos4)+marginx,height,(2*pos4)+marginx,height-5*((double)(intensity)/(double)(maximum))*scaley);
 * */
 }
 }
 }
 
memset(total, 0, sizeof(total));
height=height+5*scaley;
XFlush(thedisplay);
intensity=2;

/* Chipseq */  
sprintf(chromosome,"%s",chrom);
char *ctcf=concat("../CHIP/CTCF_chr",chromosome);

chipseq=fopen(ctcf,"rt");
if(chipseq != NULL){
char track6[100]="Ctcf";
XSetForeground(thedisplay,thecontext,blackColor);
XDrawString(thedisplay,thewindow,thecontext,marginx/4,height-scaley,track6,strlen(track6));
while(fscanf(chipseq,"%s\t%i\t%i\t%i\t%i",&chromosome,&pos,&pos2,&intensity)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
 pos3=((double)(pos-box3)/(double)resolution);
 total[pos3]+=(intensity);
if(total[pos3]>maximum){maximum=total[pos3];}

 }
 }
rewind(chipseq);

intensity=0.;

check=(1000);
scalechic=check/resolution;
 for(pos3=1;pos3<=1000;pos3++){
   pos4=pos3+1;

intensity=(total[pos3]);

rd=1255*257*((double)(maximum-total[pos3])/(double)(maximum));
gn=1255*257*((double)(maximum-total[pos3])/(double)(maximum));
be=1255*257*((double)(maximum-total[pos3])/(double)(maximum));

xcolour.red =rd; xcolour.green =gn; xcolour.blue = be;

xcolour.flags = DoRed | DoGreen | DoBlue;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);

if(pos3!=pos4){
   v[0].y =(height);       v[0].x=(2*pos3)+marginx;
   v[1].y =(height);       v[1].x=(2*pos4)+marginx;
   v[2].y =(height-2*((double)(intensity2))*scaley); v[2].x=(2*pos4)+marginx;
   v[3].y =(height-2*((double)(intensity2))*scaley); v[3].x=(2*pos3)+marginx;
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
}else{
/*XDrawLine(thedisplay,thewindow,thecontext,(2*pos4)+marginx,height,(2*pos4)+marginx,height-5*((double)(intensity)/(double)(maximum))*scaley);
*/
 }}

fclose(chipseq);
memset(total, 0, sizeof(total));
height=height+5*scaley;
}
XFlush(thedisplay);

intensity=1;

/* Chipseq */
sprintf(chromosome,"%s",chrom);
char *hist1=concat("../CHIP/RELA_chr",chromosome);
chipseq=fopen(hist1,"rt");
if(chipseq != NULL){
char track6[100]="RELA";
XSetForeground(thedisplay,thecontext,blackColor);
XDrawString(thedisplay,thewindow,thecontext,marginx/4,height-scaley,track6,strlen(track6));
while(fscanf(chipseq,"%s\t%i\t%i\t%i\t%i",&chromosome,&pos,&pos2,&intensity)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
 pos3=((double)(pos-box3)/(double)resolution);
 pos4=((double)(pos2-box3)/(double)resolution);
 total[pos3]+=(intensity);
if(total[pos3]>maximum){maximum=total[pos3];}

}
}
rewind(chipseq);

intensity=0.;

check=(1000);
scalechic=check/resolution;

 for(pos3=1;pos3<=1000;pos3++){
 pos4=pos3+1;

intensity=(total[pos3]);

rd=1255*257*((double)(maximum-total[pos3])/(double)(maximum));
gn=1255*257*((double)(maximum-total[pos3])/(double)(maximum));
be=1255*257*((double)(maximum-total[pos3])/(double)(maximum));

xcolour.red =rd; xcolour.green =gn; xcolour.blue = be;

xcolour.flags = DoRed | DoGreen | DoBlue;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);

if(pos3!=pos4){
   v[0].y =(height);       v[0].x=(2*pos3)+marginx;
   v[1].y =(height);       v[1].x=(2*pos4)+marginx;
   v[2].y =(height-2*((double)(intensity2))*scaley); v[2].x=(2*pos4)+marginx;
   v[3].y =(height-2*((double)(intensity2))*scaley); v[3].x=(2*pos3)+marginx;
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
}else{
/*XDrawLine(thedisplay,thewindow,thecontext,(2*pos4)+marginx,height,(2*pos4)+marginx,height-5*((double)(intensity)/(double)(maximum))*scaley);
 * */
}
 }
fclose(chipseq);
height=height+5*scaley;
memset(total, 0, sizeof(total));
}
/* Chipseq */
sprintf(chromosome,"%s",chrom);
char *hist2=concat("../CHIP/H3K4me1_chr",chromosome);

chipseq=fopen(hist2,"rt");
if(chipseq != NULL){
char track6[100]="H3K4me1";
XSetForeground(thedisplay,thecontext,blackColor);
XDrawString(thedisplay,thewindow,thecontext,marginx/4,height-scaley,track6,strlen(track6));
while(fscanf(chipseq,"%s\t%i\t%i\t%i\t%i",&chromosome,&pos,&pos2,&intensity)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
 pos3=((double)(pos-box3)/(double)resolution);
 pos4=((double)(pos2-box3)/(double)resolution);
 total[pos3]+=(intensity);
if(total[pos3]>maximum){maximum=total[pos3];}

}
}
rewind(chipseq);

intensity=0.;

check=(1000);
scalechic=check/resolution;

 for(pos3=1;pos3<=1000;pos3++){
 pos4=pos3+1;

intensity=(total[pos3]);

rd=1255*257*((double)(maximum-total[pos3])/(double)(maximum));
gn=1255*257*((double)(maximum-total[pos3])/(double)(maximum));
be=1255*257*((double)(maximum-total[pos3])/(double)(maximum));

xcolour.red =rd; xcolour.green =gn; xcolour.blue = be;

xcolour.flags = DoRed | DoGreen | DoBlue;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);

if(pos3!=pos4){
   v[0].y =(height);       v[0].x=(2*pos3)+marginx;
   v[1].y =(height);       v[1].x=(2*pos4)+marginx;
   v[2].y =(height-2*((double)(intensity2))*scaley); v[2].x=(2*pos4)+marginx;
   v[3].y =(height-2*((double)(intensity2))*scaley); v[3].x=(2*pos3)+marginx;
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
}else{
/*XDrawLine(thedisplay,thewindow,thecontext,(2*pos4)+marginx,height,(2*pos4)+marginx,height-5*((double)(intensity)/(double)(maximum))*scaley);
 * */
}
 }
fclose(chipseq);
height=height+5*scaley;
memset(total, 0, sizeof(total));
 }
/* Chipseq */
sprintf(chromosome,"%s",chrom);
char *hist3=concat("../CHIP/H3K4me3_chr",chromosome);

chipseq=fopen(hist3,"rt");
if(chipseq != NULL){
char track6[100]="H3K4me3";
XSetForeground(thedisplay,thecontext,blackColor);
XDrawString(thedisplay,thewindow,thecontext,marginx/4,height-scaley,track6,strlen(track6));
while(fscanf(chipseq,"%s\t%i\t%i\t%i\t%i",&chromosome,&pos,&pos2,&intensity)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
 pos3=((double)(pos-box3)/(double)resolution);
 pos4=((double)(pos2-box3)/(double)resolution);
 total[pos3]+=(intensity);
if(total[pos3]>maximum){maximum=total[pos3];}

}
}
rewind(chipseq);

intensity=0.;


check=(1000);
scalechic=check/resolution;

 for(pos3=1;pos3<=1000;pos3++){
 pos4=pos3+1;
intensity=(total[pos3]);

rd=1255*257*((double)(maximum-total[pos3])/(double)(maximum));
gn=1255*257*((double)(maximum-total[pos3])/(double)(maximum));
be=1255*257*((double)(maximum-total[pos3])/(double)(maximum));

xcolour.red =rd; xcolour.green =gn; xcolour.blue = be;

xcolour.flags = DoRed | DoGreen | DoBlue;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);

if(pos3!=pos4){
   v[0].y =(height);       v[0].x=(2*pos3)+marginx;
   v[1].y =(height);       v[1].x=(2*pos4)+marginx;
   v[2].y =(height-2*((double)(intensity2))*scaley); v[2].x=(2*pos4)+marginx;
   v[3].y =(height-2*((double)(intensity2))*scaley); v[3].x=(2*pos3)+marginx;
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
}else{
/*XDrawLine(thedisplay,thewindow,thecontext,(2*pos4)+marginx,height,(2*pos4)+marginx,height-5*((double)(intensity)/(double)(maximum))*scaley);
 * */
}
 }

fclose(chipseq);
height=height+5*scaley;
memset(total, 0, sizeof(total));
}
/* Chipseq */
sprintf(chromosome,"%s",chrom);
char *hist4=concat("../CHIP/H3K36me3_chr",chromosome);

chipseq=fopen(hist4,"rt");
if(chipseq != NULL){
char track6[100]="H3K36me3";
XSetForeground(thedisplay,thecontext,blackColor);
XDrawString(thedisplay,thewindow,thecontext,marginx/4,height-scaley,track6,strlen(track6));
while(fscanf(chipseq,"%s\t%i\t%i\t%i\t%i",&chromosome,&pos,&pos2,&intensity)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
 pos3=((double)(pos-box3)/(double)resolution);
 pos4=((double)(pos2-box3)/(double)resolution);
 total[pos3]+=(intensity);
if(total[pos3]>maximum){maximum=total[pos3];}

}
}
rewind(chipseq);

intensity=0.;

check=(1000);
scalechic=check/resolution;

 for(pos3=1;pos3<=1000;pos3++){
 pos4=pos3+1;


intensity=(total[pos3]);

rd=1255*257*((double)(maximum-total[pos3])/(double)(maximum));
gn=1255*257*((double)(maximum-total[pos3])/(double)(maximum));
be=1255*257*((double)(maximum-total[pos3])/(double)(maximum));

xcolour.red =rd; xcolour.green =gn; xcolour.blue = be;

xcolour.flags = DoRed | DoGreen | DoBlue;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);

if(pos3!=pos4){
   v[0].y =(height);       v[0].x=(2*pos3)+marginx;
   v[1].y =(height);       v[1].x=(2*pos4)+marginx;
   v[2].y =(height-2*((double)(intensity2))*scaley); v[2].x=(2*pos4)+marginx;
   v[3].y =(height-2*((double)(intensity2))*scaley); v[3].x=(2*pos3)+marginx;
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
}else{
/*XDrawLine(thedisplay,thewindow,thecontext,(2*pos4)+marginx,height,(2*pos4)+marginx,height-5*((double)(intensity)/(double)(maximum))*scaley);
 * */
}
}
fclose(chipseq);
height=height+5*scaley;
memset(total, 0, sizeof(total));
}
/*XDrawLine(thedisplay,thewindow,thecontext,box1,height,box2,height);
*/
XDrawLine(thedisplay,thewindow,thecontext,box1,height,box2,height);
height=height+1*scaley;


r=255;
g=255;
b=255;

xcolour.red =r; xcolour.green =g; xcolour.blue = b;
xcolour.flags = DoRed | DoGreen | DoBlue;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);

pos3=(double)(position-box3)/((double)(resolution));
pos4=(double)(position-box3)/((double)(resolution));

/*
XDrawString(thedisplay,thewindow,thecontext,(2*pos3)+marginx-strlen(targetinfo),height-scaley,targetinfo,strlen(targetinfo));
XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx,height,(2*pos3)+marginx,height+20*scaley);
*/
XFlush(thedisplay);
height=height+5*scaley;

/* THE GENE TRACK*/ 
maximum=0;
char track9[100]="Ensembl Genes";
XDrawLine(thedisplay,thewindow,thecontext,box1,height,box2,height);
XDrawString(thedisplay,thewindow,thecontext,marginx/4,height-scaley,track9,strlen(track9));
heightsave=height+4*scaley;

XDrawLine(thedisplay,thewindow,thecontext,box1,height,box2,height);

heightsave=height+8*scaley;

box5=(int)(box3/(resolution));
box6=(int)(box4/(resolution));

for(j=0;j<100;j++){for (i=0;i<2*(box6-box5);i++){ heightgene[i][j]=0;}}
for(j=0;j<100;j++){for (i=0;i<2*(box6-box5);i++){ heighttaggene[i][j]=0;}}

tot=0;
flag=0;
genes=fopen("./genes","r");
char *chromosome2=concat("chr",chrom);

while(fscanf(genes,"%s %s %d %d %d %d %d %s %s %s %s",&chromosome,&strand,&pos,&pos2,&posb,&pos2b,&excount,&eS,&eE,&gS,&gT)!=EOF){
if(strncmp(chromosome,chromosome2,5)==0){
if(posb>=box3&&posb<=box4){

RemoveChars(eE, '\"');
RemoveChars(eS, '\"');

const char *s = eS;
const char *e = eE;

if(strcmp(gT, "protein_coding")==0){
xcolour.red =0; xcolour.green =0; xcolour.blue=50000;
xcolour.flags = DoRed | DoGreen | DoBlue;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);
}else{
xcolour.red =50000; xcolour.green =0; xcolour.blue=0;
xcolour.flags = DoRed | DoGreen | DoBlue;
XAllocColor(thedisplay,thecolormap,&xcolour);
XSetForeground(thedisplay,thecontext,xcolour.pixel);
}

/*if((tot/3==tot%3)){
height=height+10*scaley;}
*/
pos7=(double)(pos-box3)/((double)(resolution));
pos8=(double)(pos2-box3)/((double)(resolution));
pos9=(double)(posb-box3)/((double)(resolution));
pos10=(double)(pos2b-box3)/((double)(resolution));

bingene=(int)(pos/(resolution))-box5;
bingene2=(int)(pos2/(resolution))-box5;

taggene=(int)(pos)/(resolution)-box5-strlen(targetinfo);
taggene2=(int)(pos)/(resolution)-box5+strlen(targetinfo);

/*printf("%d\t%d\t%d\t%d\t%d\n",bingene,bingene2,taggene,taggene2,strlen(targetinfo)); */
for(j=0;j<10;j++){
check2=0;
for(i=bingene-10;i<=bingene2+10;i++){
check2+=heightgene[i][j];
}
if(check2==0){break;}}

height=heightsave+j*12*scaley;
for(i=bingene-10;i<=bingene2+10;i++){heightgene[i][j]=1;}

for(jj=0;jj<100;jj++){
check2=0;
for(i=taggene-10;i<=taggene2+10;i++){
check2+=heighttaggene[i+jj][j];
}
if(check2==0){break;}}

heighttag=jj*scale;

for(i=taggene;i<=taggene2+10;i++){heighttaggene[i][j]=1;}

/*if(bingene==bingene2){heightgene[bingene]+=5*scaley;}
 * if(bingene!=bingene2){*/

for(i=bingene-5;i<=bingene2+5;i++){heightgene[i][j]=1;}

/*printf("%d\t%d\t%d\t%d\t%d\t%s\n",pos7,pos8,bingene,bingene2,height,gS);
*/

if(height<900){
if(height>maximum){maximum=height;}

XDrawString(thedisplay,thewindow,thecontext,(2*pos7)+marginx+heighttag-4*strlen(gS),height-1.5*scaley,gS,strlen(gS));
XDrawLine(thedisplay,thewindow,thecontext,(2*pos7)+marginx,height+scaley,(2*pos8)+marginx,height+scaley);
XDrawLine(thedisplay,thewindow,thecontext,(2*pos7)+marginx,height+1.1*scaley,(2*pos8)+marginx,height+1.1*scaley);
XDrawLine(thedisplay,thewindow,thecontext,(2*pos7)+marginx,height+0.9*scaley,(2*pos8)+marginx,height+0.9*scaley);
/* draw here the strand arrow */

count=0;
do {
flag=0;
        size_t field_len = strcspn(s, delims);
        size_t field_len2 = strcspn(e, delims);

         sprintf(chpos,"%.*s\t", (int)field_len, s);
         sprintf(chpos2,"%.*s\n\n", (int)field_len2, e);

s += field_len;
e += field_len2;
count=count+1;
pos5=atoi(chpos);
pos6=atoi(chpos2);
savestart[count]=pos5;
saveend[count]=pos6;

pos3b=(double)(saveend[count-1]-box3)/((double)(resolution));
pos5b=(double)(savestart[count]-box3)/((double)(resolution));
pos3=(double)(pos5-box3)/((double)(resolution));
pos4=(double)(pos6-box3)/((double)(resolution));

if(pos5b>0&&pos3b>0){
if(pos5b-pos3b>7&&strcmp(strand, "+")==0&&count>1){
XDrawLine(thedisplay,thewindow,thecontext,(pos3b+pos5b)+marginx,height+1*scaley+delta,(pos3b+pos5b+delta)+marginx,height+scaley);
XDrawLine(thedisplay,thewindow,thecontext,(pos3b+pos5b)+marginx,height+1*scaley-delta,(pos3b+pos5b+delta)+marginx,height+scaley);
XDrawLine(thedisplay,thewindow,thecontext,(pos3b+pos5b)+marginx+0.1*scaley,height+1*scaley+delta,(pos3b+pos5b+delta)+marginx+0.1*scaley,height+scaley);
XDrawLine(thedisplay,thewindow,thecontext,(pos3b+pos5b)+marginx+0.1*scaley,height+1*scaley-delta,(pos3b+pos5b+delta)+marginx+0.1*scaley,height+scaley);

}
if(pos3b-pos5b>7&&strcmp(strand, "-")==0&&count>1){
XDrawLine(thedisplay,thewindow,thecontext,(pos3b+pos5b)+marginx,height+scaley+delta,(pos3b+pos5b-delta)+marginx,height+scaley);
XDrawLine(thedisplay,thewindow,thecontext,(pos3b+pos5b)+marginx,height+scaley-delta,(pos3b+pos5b-delta)+marginx,height+scaley);
XDrawLine(thedisplay,thewindow,thecontext,(pos3b+pos5b)+marginx-0.1*scaley,height+1*scaley+delta,(pos3b+pos5b-delta)+marginx-0.1*scaley,height+1*scaley);
XDrawLine(thedisplay,thewindow,thecontext,(pos3b+pos5b)+marginx-0.1*scaley,height+1*scaley-delta,(pos3b+pos5b-delta)+marginx-0.1*scaley,height+1*scaley);
/*printf("%s\t%d\t%d\n",gS,pos3b,pos5b);
*/

}}

if(pos5!=0&&pos6!=0){
   xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);

   if(pos3<=pos9&&pos4>=pos9){
           v[0].y =(height)+0*scaley; v[0].x=(2*pos3)+marginx;
           v[1].y =(height)+0*scaley; v[1].x=(2*pos9)+marginx;
           v[2].y =(height)+2*scaley; v[2].x=(2*pos9)+marginx;
           v[3].y =(height)+2*scaley; v[3].x=(2*pos3)+marginx;

XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
if(pos3==pos4){
XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx,height+0*scaley,(2*pos4)+marginx,height+2*scaley);
XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx+0.1*scaley,height+0*scaley,(2*pos4)+marginx+0.1*scaley,height+2*scaley);
}
tmp=min(pos4,pos10);
           v[0].y =(height)+0.*scaley; v[0].x=(2*pos9)+marginx;
           v[1].y =(height)+0.*scaley; v[1].x=(2*tmp)+marginx;
           v[2].y =(height)+2*scaley; v[2].x=(2*tmp)+marginx;
           v[3].y =(height)+2*scaley; v[3].x=(2*pos9)+marginx;

XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);}else if(pos3<pos9&&pos4<pos9){
if(pos3==pos4){
XDrawLine(thedisplay,thewindow,thecontext,(2*pos4)+marginx,height+0*scaley,(2*pos4)+marginx,height+2*scaley);
XDrawLine(thedisplay,thewindow,thecontext,(2*pos4)+marginx+0.1*scaley,height+0*scaley,(2*pos4)+marginx+0.1*scaley,height+2*scaley);
}
           v[0].y =(height)-1*scaley; v[0].x=(2*pos3)+marginx;
           v[1].y =(height)-1*scaley; v[1].x=(2*pos4)+marginx;
           v[2].y =(height)+3*scaley; v[2].x=(2*pos4)+marginx;
           v[3].y =(height)+3*scaley; v[3].x=(2*pos3)+marginx;

XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);}else if(pos3>pos9&&pos4>pos9&&pos3<pos10&&pos4<pos10){
tmp=min(pos4,pos10);
           v[0].y =(height)-1.*scaley; v[0].x=(2*pos3)+marginx;
           v[1].y =(height)-1.*scaley; v[1].x=(2*tmp)+marginx;
           v[2].y =(height)+3*scaley; v[2].x=(2*tmp)+marginx;
           v[3].y =(height)+3*scaley; v[3].x=(2*pos3)+marginx;
if(pos3==pos4){
XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx,height-1*scaley,(2*pos4)+marginx,height+3.0*scaley);
XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx+0.1*scaley,height-1*scaley,(2*pos4)+marginx+0.1*scaley,height+3.0*scaley);
}

XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);}
XFlush(thedisplay);

if(pos3<=pos10&&pos4>=pos10){
           tmp=max(pos3,pos9);
           v[0].y =(height)+0.*scaley; v[0].x=(2*tmp)+marginx;
           v[1].y =(height)+0.*scaley; v[1].x=(2*pos10)+marginx;
           v[2].y =(height)+2.0*scaley; v[2].x=(2*pos10)+marginx;
           v[3].y =(height)+2.0*scaley; v[3].x=(2*tmp)+marginx;

XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
if(pos3==pos10){
XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx,height+0.*scaley,(2*pos3)+marginx,height+2*scaley);
XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx+0.1*scaley,height+0.*scaley,(2*pos3)+marginx+0.1*scaley,height+2*scaley);
}

           v[0].y =(height)-1*scaley; v[0].x=(2*pos10)+marginx;
           v[1].y =(height)-1*scaley; v[1].x=(2*pos4)+marginx;
           v[2].y =(height)+3*scaley; v[2].x=(2*pos4)+marginx;
           v[3].y =(height)+3*scaley; v[3].x=(2*pos10)+marginx;

XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
if(pos10==pos4){
XDrawLine(thedisplay,thewindow,thecontext,(2*pos4)+marginx,height-1*scaley,(2*pos4)+marginx,height+3*scaley);
XDrawLine(thedisplay,thewindow,thecontext,(2*pos4)+marginx+0.1*scaley,height-1*scaley,(2*pos4)+marginx+0.1*scaley,height+3*scaley);
 }
}else if(pos3>pos10&&pos4>pos10){

           v[0].y =(height)-1*scaley; v[0].x=(2*pos3)+marginx;
           v[1].y =(height)-1*scaley; v[1].x=(2*pos4)+marginx;
           v[2].y =(height)+3*scaley; v[2].x=(2*pos4)+marginx;
           v[3].y =(height)+3*scaley; v[3].x=(2*pos3)+marginx;

XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
if(pos3==pos4){
XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx,height-1*scaley,(2*pos4)+marginx,height+3*scaley);
XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx+0.1*scaley,height-1*scaley,(2*pos4)+marginx+0.1*scaley,height+3*scaley);
}
}else if(pos3<pos10&&pos4<pos10&&pos3>pos9&&pos4>pos9){
tmp=max(pos3,pos9);
           v[0].y =(height)+0.*scaley; v[0].x=(2*tmp)+marginx;
           v[1].y =(height)+0.*scaley; v[1].x=(2*pos4)+marginx;
           v[2].y =(height)+2*scaley; v[2].x=(2*pos4)+marginx;
           v[3].y =(height)+2*scaley; v[3].x=(2*tmp)+marginx;

XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
if(pos3==pos4){
XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx,height+0.*scaley,(2*pos4)+marginx,height+2.0*scaley);
XDrawLine(thedisplay,thewindow,thecontext,(2*pos3)+marginx+0.1*scaley,height+0.*scaley,(2*pos4)+marginx+0.1*scaley,height+2.0*scaley);

}}
XFlush(thedisplay);
 }}while (*s++&&*e++);

tot=tot+1;
   XSetForeground(thedisplay,thecontext,blackColor);
   XDrawString(thedisplay,thewindow,thecontext,(2*pos3)+marginx-5*scale,height+scaley,s,strlen(s));
   XFlush(thedisplay);
}
XSetForeground(thedisplay,thecontext,blackColor);
 XFlush(thedisplay);
}
}
 }
fclose(genes);
/*printf("MAX:%d\n",maximum);
THE AXES */

height=maximum+5*scaley;
XFlush(thedisplay);
XSetForeground(thedisplay,thecontext,blackColor);
XDrawLine(thedisplay,thewindow,thecontext,box1,height+0.1*scaley,box2,height+0.1*scaley); 
XDrawLine(thedisplay,thewindow,thecontext,box1,height,box2,height);
XDrawLine(thedisplay,thewindow,thecontext,box1,height-0.1*scaley,box2,height-0.1*scaley);

box5=(int)(box3/(50*resolution));
box6=(int)(box4/(50*resolution));
box7=(int)(box3/(50*resolution));

/*printf("%i\t%i\t%i\t%i\n",box3,box4,box5,box6);
*/
for (i=0;i<(box6-box5);i=i+1) {
ii=((i+box7+1)*(50*resolution));
iii=2*(ii-box3)/resolution+marginx;
XDrawLine(thedisplay,thewindow,thecontext,iii,height,iii,height+scaley);
XDrawLine(thedisplay,thewindow,thecontext,iii+0.1*scaley,height,iii+0.1*scaley,height+scaley);
XDrawLine(thedisplay,thewindow,thecontext,iii-0.1*scaley,height,iii-0.1*scaley,height+scaley);
sprintf(end, "%2.2fMb",(double)(ii)/1000000);
XDrawString(thedisplay,thewindow,thecontext,iii-strlen(end)*0.25*scaley,height+9.5*scaley,end,strlen(end));
}

XFlush(thedisplay); 

 /* Draw some text in black */
height=height+.5*scaley;
ii=(((int)(box6-box5)/2+box5)*(50*resolution));
iii=2*(ii-box3)/(20*resolution)*scale+marginx;

   sprintf(start,"%s, position (Mb)",chromosome2);
   XSetForeground(thedisplay,thecontext,blackColor);
   XDrawString(thedisplay,thewindow,thecontext,iii+strlen(end)*12.*scaley,height+14.5*scaley,start,strlen(start));
   XFlush(thedisplay);

char *hmap2=concat(targetinfo,"_hm2");
heatmap=fopen(hmap2,"rt");
if(heatmap != NULL){
while(fscanf(heatmap,"%s %d %s %d %f",&chromosome,&pos,&chromosome,&pos2,&hm)!=EOF)
{
f=(int)(hm);
if((f)>maximumf){maximumf=(f);}
}
rewind(heatmap);

/*scale;*/

box1=marginx;
box2=2*scale+marginx+2*pos4;

height=height-55*scaley;

while(fscanf(heatmap,"%s %d %s %d %f",&chromosome,&pos,&chromosome,&pos2,&hm)!=EOF)
{
if(pos>=box3&&pos<=box4&&pos2>=box3&&pos2<=box4){
pos3=((pos-box3)/resolution);
pos4=((pos2-box3)/resolution);

brightness=1;
hmm=1-log(hm+1)/log(maximumf+1);
rgb=HSL2RGB(10,0.50,hm);
if(hm>1){
xcolour.red=255*rgb[1];
xcolour.green=0*rgb[2];
xcolour.blue=0*rgb[3];}
if(hm<=1){xcolour.red=0*rgb[1];
xcolour.green=0*rgb[2];
xcolour.blue=255*rgb[3];}



 /* xcolour.red =60*(int)(1.5*(hmm)); xcolour.green =60*(int)(0.5*(hmm)); xcolour.blue =60*(int)(3.7*(hmm));

  printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\n",pos,pos3,pos4,pos2,xcolour.red,xcolour.green,xcolour.blue);
*/
/* Draw a red filled polygon with a black border */

   v[0].y =fabs(pos3-pos4)*(double)(scaley)/(double)(scale)+height+scaley+marginy; v[0].x =(pos3+pos4)+2*scale+marginx; 
   v[1].y =fabs(pos3-pos4)*(double)(scaley)/(double)(scale)+2*scaley+height+marginy;v[1].x =(pos3+pos4)+scale+marginx; 
   v[2].y =fabs(pos3-pos4)*(double)(scaley)/(double)(scale)+height+scaley+marginy; v[2].x = (pos3+pos4)+0+marginx; 
   v[3].y =fabs(pos3-pos4)*(double)(scaley)/(double)(scale)+height+marginy; v[3].x = (pos3+pos4)+scale+marginx; 

   if(pos4<pos3||pos4<pos3){
   xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
   XFillPolygon(thedisplay,thewindow,thecontext,v,4,Complex,CoordModeOrigin);
   XSetForeground(thedisplay,thecontext,blackColor);
   }
 
   v[0].y =fabs(pos3-pos4)*(double)(scaley)/(double)(scale)+scaley+height+marginy; v[0].x =(pos3+pos4)+2*scale+marginx; 
   v[3].y =fabs(pos3-pos4)*(double)(scaley)/(double)(scale)+2*scaley+height+marginy;v[1].x =(pos3+pos4)+scale+marginx; 
   v[2].y =fabs(pos3-pos4)*(double)(scaley)/(double)(scale)+scaley+height+marginy; v[2].x = (pos3+pos4)+0+marginx;
   v[1].y =fabs(pos3-pos4)*(double)(scaley)/(double)(scale)+height+marginy; v[3].x = (pos3+pos4)+scale+marginx; 


   if(pos4==pos3){
   xcolour.flags = DoRed | DoGreen | DoBlue;
   XAllocColor(thedisplay,thecolormap,&xcolour);
   XSetForeground(thedisplay,thecontext,xcolour.pixel);
   XFillPolygon(thedisplay,thewindow,thecontext,v,3,Complex,CoordModeOrigin);
   XSetForeground(thedisplay,thecontext,blackColor);
}
 }
 }
fclose(heatmap);
height=height+4*scaley;
 
XClearWindow(thedisplay,thewindow);

printf("Chromosome:%s\nPosition:%d\nWindow size: %d\nAccessing the Hi-C database, please wait",chrom,position,buffer);
}
sleep(1000);
}

  /* strindex:  return index of t in s, -1 if none */
   int strindex(char s[], char t[])
   {
       int i, j, k;

       for (i = 0; s[i] != '\0'; i++) {
           for (j=i, k=0; t[k]!='\0' && s[j]==t[k]; j++, k++)
               ;
           if (k > 0 && t[k] == '\0')
               return i;
       }
       return -1;


   }

