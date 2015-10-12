#include <stdio.h>
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_xtc.h"

// Created by Veselin Kolev <vesso.kolev@gmail.com>
// 20151011042819


// * This sample program read XTC trajectory file produced
// by GROMACS and display the content of each frame. It uses
// xdrfile library provided by GROMACS project * //

int main(void) {

   char *filename="/storage/vesso/PB/MD/wall/trajout.xtc";
   char *mode="r";
   int natoms,step,i;
   float time;
   int status;
   float prec;
   XDRFILE *t;
   mybool *bRead;
   matrix box;
   rvec *x;

   t=xdrfile_open(filename,mode);
  
   printf("%d\n",read_xtc_natoms(filename,&natoms));
   printf("%d\n",natoms);

   x=(rvec * )calloc(natoms,sizeof(x[0]));

   while (status==0) {

      status=read_xtc(t,natoms,&step,&time,box,x,&prec);

      printf("step %d, time %f\n",step,time);

      for (i=0;i<natoms;i++) {
         printf("%f %f %f\n",x[i][0],x[i][1],x[i][2]);
      }
      printf("\n");
   }

   status=xdrfile_close(t);

   return 0;

}
