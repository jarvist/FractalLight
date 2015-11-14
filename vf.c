/* Jarvist Frost 2004-2006
 * Program to create 'video-fractals'
 */

//#include <file.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>

#define MAG (1.4)
#define X_RES 500
#define Y_RES 500

#define X_OFF 0 //offset of newcenter in pixels
#define Y_OFF 0

#define TWIST 3.14/6 //radians twist between zoom's

#define PIXW 0.65 //0.65 //width of sensor pixel in display pixels 

#define BACKGROUND 140

int curpic[X_RES][Y_RES];
int newpic[X_RES][Y_RES];
void outputpic(char *filename);
void inputpic(char *filename);

main()
{
   int x,y,loops;
   char filename[20];

   //fill display with white noise
   srand(123);
   for (x=0;x<X_RES;x++)
     for (y=0;y<Y_RES;y++)
       curpic[x][y]=rand();
   
   inputpic("begin.pgm");


   
   outputpic("first.pgm");
//   zoom();
//   swap();
//   outputpic("test2.pgm");
   
   for(loops=0;loops<150;loops++)
     {
	printf("%d\n",loops);
	sprintf(filename,"pic%.3d.pgm",loops);
	if (loops%10==0)
  	   outputpic(filename);

	zoom(); swap();
     }
   
   outputpic("last.pgm");
}

swap()
{
  int x,y,max=0;
   double light;
   
   for (x=0;x<X_RES;x++)
     for(y=0;y<Y_RES;y++)
       {  
       curpic[x][y]=newpic[x][y];
       if (curpic[x][y]>max)
          max=curpic[x][y];
       }

   for (x=0;x<X_RES;x++)
     for(y=0;y<Y_RES;y++)
       curpic[x][y]*=INT_MAX/max;   
}


zoom()
{
 int x,y;
   double nx,ny,np,dx,dy,r,theta,phi;
   
   for (x=0;x<X_RES;x++)
     for(y=0;y<Y_RES;y++)
       {
	  dx=(0.5+(double)(x-(X_RES/2)))/MAG;
	  dy=(0.5+(double)(y-(Y_RES/2)))/MAG;
	  
	  r=sqrt(dx*dx+dy*dy);
	  theta=atan2(dy,dx);
	  phi=theta+TWIST;
	  
	  ny=r*sin(phi);
	  nx=r*cos(phi);
	  
	  nx=nx+(double)(X_RES/2+X_OFF);
	  ny=ny+(double)(Y_RES/2+Y_OFF);
	  
//	  printf("x: %d nx: %f y: %d ny: %f\n",x,nx,y,ny);
	  
	  np=0; 
	  np+=vo((int)nx+1,(int)ny+1) *(nx-(double)((int)nx))*(ny-(double)((int)ny)); //top-right pixel
	  np+=vo((int)nx,(int)ny+1) *(PIXW-(nx-(double)((int)nx)))*(ny-(double)((int)ny)); //top-left pixel
	  np+=vo((int)nx+1,(int)ny) *(nx-(double)((int)nx))*(PIXW-(ny-(double)((int)ny))); //bot-right pixel
	  np+=vo((int)nx,(int)ny) *(PIXW-(double)(nx-(double)((int)nx)))*(PIXW-(ny-(double)((int)ny))); //bot-left pixel
	  
	  np*=1.0/(PIXW*PIXW); //compensates for size of pixel otherwise 'losing' light from the feedback
	  
	  newpic[x][y]=(int)np;
//	  printf("np: %f\n",np);
       }
   
   
}

int vo(int x, int y) //value of a particular pixel; with bounds checking
{
   if (x<0||x>X_RES) return(BACKGROUND);
   if (y<0||y>Y_RES) return(BACKGROUND);
   
   return (curpic[x][y]);	
}
	

void outputpic(char *filename)
{
   int x,y;
   FILE * f;
   f=fopen(filename,"w");
   
   fprintf(f,"P5 %d %d 255\n",X_RES,Y_RES); //.pgm filetype - binary form
   
   for (x=0;x<X_RES;x++)
      for (y=0;y<Y_RES;y++)
         fprintf(f,"%c",curpic[x][y]/(INT_MAX/255));
   
   fclose(f);
}

void inputpic(char *filename)
{ 
  int x,y;
   int tmp;
   
   FILE * f;
   f=fopen(filename,"r");

   fscanf(f,"P2 500 500\n"); //.pgm filetype

   for (x=0;x<X_RES;x++)
        for (y=0;y<Y_RES;y++)
       {
             fscanf(f,"%d",&tmp);
//	  printf("%d\n",curpic[x][y]);
// printf("tmp: %d\n",tmp);
	  curpic[x][y]=tmp*(INT_MAX/255);
       }
   

   fclose(f);
}
