/* Copyright (c) Colorado School of Mines, 2021.*/
/* All rights reserved.                       */

/* SUTOOLCSV: $Revision: 1.00 $ ; $Date: 2021/05/15 00:09:00 $        */

#include <stdio.h>
#include <string.h>
#include "math.h"

#include "su.h"
#include "segy.h"

segy tr;

void readkfile(FILE *fpR, cwp_String *names, cwp_String *forms, double *dfield, int *numcases) ;
void writekfile(FILE *fpW, cwp_String *names, cwp_String *forms, double *dfield, int numcasesout) ;
void getCSV(char *textraw, char *textbeg, int maxtext, char rdel, 
            double *dfield, int *nspot, int numcases,   
            int ncount, int *comerr,int *morerr,int *numerr,int *nblank);
void tparse(char *tbuf, char d, char **fields, int *numfields) ; 
void gridset(double *gvals, int *errwarn) ; 
void gridrawxycdpic(double *gvals,double dx,double dy,int *icdp,int *igi,int *igc) ;
void gridicrawxy(double *gvals,int igi,int igc,double *dx,double *dy) ;
void gridicgridxy(double *gvals,int igi,int igc,double *dx,double *dy) ;
void gridiccdp(double *gvals,int igi,int igc,int *icdp) ;
void gridcdpic(double *gvals,int icdp,int *igi,int *igc) ;
void gridrawxygridxy(double *gvals,double dx,double dy,double *tx,double *ty) ;
void gridgridxyrawxy(double *gvals,double dx,double dy,double *tx,double *ty) ;
void gridcheck(double *gvals, int icheck, int *errwarn) ; 

/*********************** self documentation **********************/
char *sdoc[] = {
"                          ",
" SUBINCSV - Various CDP and other binning options                         ",
"                                                                          ",
" subincsv  rfile=in.csv wfile=out.csv <in.su >out.su                      ",
"                                                                          ",
" Parameter overview:                                                      ",
"                                                                          ",
"       rfile= if specified, read a file containing parameter values.      ", 
"       wfile= if specified, write a file with parameter values.           ", 
"                                                                          ",
"       in.su and out.su do not have to be specified (if you just want     ",
"       to specify command line parameters and get the output wfile).      ",
"                                                                          ",
"       Parameters specified on the command line override parameters from  ",
"       the rfile.                                                         ",
"                                                                          ",
"       bintype=        See below for which other parameters are required  ",
"                       for the option numbers here.                       ",
"               30      Compute trace midpoint coordinates (sx+gx)/2 and   ",
"                       (gx+gy)/2 and update cdp,igi,igc keys (grid cell   ",
"                       number and inline and crossline index numbers).    ",
"              -30      Use input cdp (cell) number and update keys        ",
"                       igi,igc,sx,sy,gx,gy. Where igi,igc are grid index  ",
"                       numbers of the cdp (cell) number and sx,sy is the  ",
"                       cell centre in raw coordinates and gx,gy is cell   ",
"                       centre in grid coordinates (which are shifted and  ",
"                       aligned with grid definition, but not scaled by    ",
"                       cell widths). Note that sx,sy are only approximate ",
"                       cell centres since scalco rounds after sin,cosine  ",
"                       computations, but gx,gy are usually much more      ",
"                       precise since they are multiples of cell widths.   ",
"              -31      Use input cdp (cell) number and update igi,igc.    ",
"              -32      Use input igi,igc and update cdp number as well as ",
"                       sx,sy,gx,gy as described in option -30.            ",
"               20      Compute 2D cdp from point numbers.                 ",
"                                                                          ",
"       offset=         By default, bintype=30 recomputes the offset key,  ",
"                       but other bintypes leave it as-is.                 ",
"             =1        Recompute offset key.                              ",
"             =0        Do not recompute offset key.                       ",
"                                                                          ",
"       check=0         Do not print checking details.                     ",
"             1         For grid bintypes, after grid defintion is set,    ",
"                       run some grid functions on the 4 corner points     ",
"                       and print the results. The intention here is to    ",
"                       exercise many functions in case of issues created  ",
"                       by coding or compiler or optimizer errors/changes. ",
"                       But the output may also be useful for users.       ",
"                       For instance, you can see slight differences in    ",
"                       the coordinates of those 4 corners when produced   ",
"                       by different functions, and when run on different  ",
"                       hardware or with different compilers/optimizers.   ",
"                                                                          ",
" Grid parameters (either on command line or in rfile).                    ",
"                                                                          ",
"    grid_xa=  X coordinate of corner A.                                   ",
"    grid_ya=  Y coordinate of corner A.                                   ",
"    grid_xb=  X coordinate for corner B.                                  ",
"    grid_yb=  Y coordinate for corner B.                                  ",
"    grid_xc=  X coordinate for corner C.                                  ",
"    grid_yc=  Y coordinate for corner C.                                  ",
"    grid_wb=  width of cells in A-->B direction.                          ",
"    grid_wc=  width of cells in A-->C direction.                          ",
"                                                                          ",
" Note that corner A coordinates are used exactly, but corner B is reset   ",
" to an exact multiple distance of the cell width in A-->B direction.      ",
" And corner C is reset to a line at right angle to A-->B direction        ",
" through A and also to an exact multiple distance of A-->C cell width.    ",
"                                                                          ",
"                                                                          ",
" Point parameters (either on command line or in rfile).                   ",
"                                                                          ",
"    point_rpb= receiver point base number                                 ",
"    point_rcb= receiver cdp base number                                   ",
"    point_rpi= receiver point increment                                   ",
"    point_rci= receiver cdp increment                                     ",
"    point_spb= source point base number                                   ",
"    point_scb= source cdp base number                                     ",
"    point_spi= source point increment                                     ",
"    point_sci= source cdp increment                                       ",
"                                                                          ",
" Point parameters (only on command line).                                 ",
"                                                                          ",
"    rkey=      key containing receiver point numbers (default is gaps)    ",
"    skey=      key containing source point numbers (default is garad)     ",
"               Note these defaults match the defaults of sugeomcsv.       ",
"                                                                          ",
" ***********************************************************              ",
"   To output this documentation:  subincsv 2> bindoc.txt                  ",
" ***********************************************************              ",
"                                                                          ",
"                                                                          ",
" Seismic Unix has a 240 byte header which already has defined key names.  ",
" In other seismic processing systems the ability to expand trace headers  ",
" and insert intricate grid transform values is both a blessing and a curse",
" (trust me on that). For SU, intricate grid-related values have no keys   ",
" to be stored. Therefore only 3 keys are updated here (cdp, igi, igc).    ",
" Where igi is set to the cell index in the direction from corner A to     ",
" corner B (a mnemonic for igi is index-grid-inline). And where igc is set ",
" to cell index in the direction from corner A to corner C (a mnemonic for ",
" igc is index-grid-crossline).                                            ",
"                                                                          ",
" The input grid definition command line parameters are processed and      ",
" written to an external file. That file follows conventions established   ",
" by SUTOOLCSV and SUGEOMCSV. I call this the K-file (K for Konstants).    ",
" Reading this K-file should allows other SU programs to perform intricate ",
" grid transforms on-the-fly using sx,sy,gx,gy coordinates as well as      ",
" backwards transforms from cdp,igi,igc to cell centre XYs. The grid in    ",
" the K-file will also allow transforms of XYs values in S and R tables in ",
" spreadsheets (eventually). And, once output, the K-file can be re-input  ",
" to this program instead of using command line parameters.                ",
"                                                                          ",
" The grid will be an exact rectangle with 4 corner points A,B,C, and D.   ",
" But you can only specify XYs for corners A,B,C (D is computed herein).   ",
" Corner A coordinates are used exactly as input. Then the direction from  ", 
" corner A to input corner B is determined exactly. After that, corner B   ", 
" coordinates are adjusted to an exact multiple of the cell width you      ", 
" specify for the A-->B direction. Then your input coordinates of corner C ", 
" are used to compute the distance from corner A to corner C. The right    ", 
" angle to A-->B gives direction for output corner C (along line thru A).  ", 
" Corner C is then adjusted to an exact multiple of the cell width you     ", 
" specify for the A-->C direction. Note that this means input corner C is  ", 
" only used to decide how wide the rectangle is, and which side of A-->B   ",
" the output corner C is on. Corner D is computed from the other corners.  ",
"                                                                          ",
" The first cell is centred on corner A and has igi=1 and igc=1.           ",
" Cell centres then increment by their corresponding widths in the A-->B   ",
" and A-->C directions. igi and igc increment by 1 in the same directions. ",
" cdp starts at 1 in the first cell and increments by 1 in the A-->B       ",
" direction until it reaches B, then moves 1 cell in the A-->C direction   ",
" (near corner A) and continues to increment by 1 in the A-->B direction.  ",
"                                                                          ",
" Cells only contain one-half of their boundaries. This ensures that a     ",
" trace midpoint that is exactly between 2 or 4 cell centres is assigned   ",
" to a specific cell. Note: This is why you cannot use proximity to cell   ",
" centres to assign traces to cells. You need to actually compute the cdp  ",
" and igi,igc numbers the way that it is done herein.                      ",
"                                                                          ",
"                                                                          ",
" Warning and advice:                                                      ",
"  Cell boundaries and other grid computations use double precision values ",
"  and are therefore extremely precise. This very precision causes issues. ",
"  When a trace midpoint is very near a cell boundary, it only takes a     ",
"  slight difference in hardware/compiler/optimizer computations for the   ",
"  boundaries to move a bit, and therefore assign some traces to different ",
"  cells. You should expect that. Similarly, trying to reverse or invert a ",
"  grid by exchanging corners A and B and so on, is also not likely to     ",
"  result in exactly the same distribution of traces in the cells.         ",
"                                                                          ",
" ----------------------------------------------------------------------   ", 
" -----------------------------------------------------------------        ",
"                                                                          ",
NULL};

/* Credits:                                                       */
/* Andre Latour                                                   */ 
/*                                                                */
/* Started from sutoolcsv and sugeomcsv                           */
/*                                                                */
/* Trace header fields accessed: sx,sy,gx,gy.cdp,igi,igc,offset   */
/*                                                                */
/**************** end self doc ************************************/

int main(int argc, char **argv) {

   cwp_String Rname=NULL;  /* text file name for values            */
   FILE *fpR=NULL;         /* file pointer for Rname input file    */
   cwp_String Wname=NULL;  /* text file name for output values     */
   FILE *fpW=NULL;         /* file pointer for Wname output file   */

/* Most of following are about 10 times bigger than will be used.  */    
/* But difficult to dynamically allocate since they often will be  */    
/* read-in from C_SU_ records in the input text file.              */    

   cwp_String names[999];   
   cwp_String forms[999];   
   double dfield[999];

   cwp_String gnams[99];   
   double gvals[99];

   int nproct = 0; // number of processed traces

   /* Initialize */
   initargs(argc, argv);
   requestdoc(1);

   getparstring("rfile", &Rname);
   getparstring("wfile", &Wname);

   if(Rname != NULL && Wname != NULL && strcmp(Rname,Wname) == 0) 
     err("**** Error: wfile= output file must be different than rfile= input file.");

   int intraces = 1;
   if(isatty(STDIN_FILENO)==1) { /* do not have input trace file */
     intraces = 0;
     if (Wname == NULL) err("**** Error: wfile= output text file name must be specified when no input traces.");
     if(isatty(STDOUT_FILENO)!=1) { /* have output trace file */
       err("**** Error: Cannot specify output trace file with no input trace file.");
     }
   }
   else {
     if(isatty(STDOUT_FILENO)==1) { /* do not have output trace file */
       err("**** Error: Must have output trace file when input trace file is specified.");
     }
   }

   int bintype;
   if (!getparint("bintype", &bintype)) bintype = -1;
  
   int ioffset;
   if (!getparint("offset", &ioffset)) ioffset = -1;
  
   int icheck;
   if (!getparint("check", &icheck)) icheck = 0;

/* Cycle over rfile records to get some C_SU_ parameters? */ 

   int numcases = 0;

   if (Rname != NULL) {
     fpR = fopen(Rname, "r");
     readkfile(fpR,names,forms,dfield,&numcases);
   } /* end of  else for if (Rname == NULL) { */

/* ----------------------------------------------------------- */

   int numcasesout = numcases;

   int numgnams = 1;
   gvals[0] = -1.1e308;
   gnams[0] = ealloc1(7,1); 
   strcpy(gnams[0],"bintype");

   if(bintype==-1) {
     gvals[0] = -1.;
     for(int j=0; j<numcases; j++) { 
       if(strcmp(names[j],gnams[0]) == 0) gvals[0] = dfield[j]; /* andre, lower case? */ 
     }
     if(gvals[0] > 0.) bintype = (int) (gvals[0] + 0.1);
     else bintype = (int) (gvals[0] - 0.1);
   }
   else {
     gvals[0] = bintype;
   }

   if(bintype!=30 && bintype!=-30 && bintype!=-31 && bintype!=-32 && bintype!=20) {
     err("**** Error: bintype= option not recognized.");
   }

   if(ioffset==-1) {
     if(bintype==30) ioffset = 1;
     else ioffset = 0;
   }

/* Process and set the grid definition values?                                   */

   if(bintype==30 || bintype==-30 || bintype==-31 || bintype==-32) {

     numgnams = 18;

     for(int i=1; i<numgnams; i++) {
       gvals[i] = -1.1e308;
       gnams[i] = ealloc1(7,1); 
     }

     strcpy(gnams[1],"grid_lf");
     strcpy(gnams[2],"grid_xa");
     strcpy(gnams[3],"grid_ya");
     strcpy(gnams[4],"grid_xb");
     strcpy(gnams[5],"grid_yb");
     strcpy(gnams[6],"grid_xc");
     strcpy(gnams[7],"grid_yc");
     strcpy(gnams[8],"grid_xd");
     strcpy(gnams[9],"grid_yd");
     strcpy(gnams[10],"grid_wb");
     strcpy(gnams[11],"grid_wc");
     strcpy(gnams[12],"grid_nb");
     strcpy(gnams[13],"grid_nc");
     strcpy(gnams[14],"grid_fp");
     strcpy(gnams[15],"grid_lp");
     strcpy(gnams[16],"grid_sb");
     strcpy(gnams[17],"grid_cb");

     for(int i=2; i<12; i++) {    /* read-in corners A,B,C and cell widths */ 
       if(i==8 || i==9) continue; /* do not read-in corner D  */
       if(!getpardouble(gnams[i],gvals+i)) { 
         for(int j=0; j<numcases; j++) { 
           if(strcmp(names[j],gnams[i]) == 0) gvals[i] = dfield[j]; /* andre, lower case? */ 
         }
         if(gvals[i] < -1.e308) {
           gvals[i] = i+100;
           err("**** Error bintype=grid and %s not found.",gnams[i]);
         }
       }
     } /* end of  for(int i=2; i<12; i++) { */

     int errwarn;
     gridset(gvals,&errwarn); 
     if(errwarn == 1) err ("**** Error. The grid_wb cell width must be positive.");
     if(errwarn == 2) err ("**** Error. The grid_wc cell width must be positive.");
     if(errwarn == 3) err ("**** Error. Corner B is within grid_wb cell width of corner A.");
     if(errwarn == -1) warn ("**** Corner C is near A and is reset to A.");

     gridcheck(gvals,icheck,&errwarn); 

   } /* end of  if(bintype==30 || ..... */
   else if(bintype==20) {

     numgnams = 9;

     for(int i=1; i<numgnams; i++) {
       gvals[i] = -1.1e308;
       gnams[i] = ealloc1(9,1); 
     }

     strcpy(gnams[1],"point_rpb"); /* receiver point base number */
     strcpy(gnams[2],"point_rcb"); /* receiver cdp base number   */
     strcpy(gnams[3],"point_rpi"); /* receiver point increment   */
     strcpy(gnams[4],"point_rci"); /* receiver cdp increment     */
     strcpy(gnams[5],"point_spb"); /* source   point base number */
     strcpy(gnams[6],"point_scb"); /* source   cdp base number   */
     strcpy(gnams[7],"point_spi"); /* source   point increment   */
     strcpy(gnams[8],"point_sci"); /* source   cdp increment     */

     for(int i=1; i<9; i++) {     
       if(!getpardouble(gnams[i],gvals+i)) { 
         for(int j=0; j<numcases; j++) { 
           if(strcmp(names[j],gnams[i]) == 0) gvals[i] = dfield[j]; /* andre, lower case? */ 
         }
         if(gvals[i] < -1.e308) {
           if(i<5) {
             gvals[i] = i+100;
             warn("**** Error bintype=point and %s not found.",gnams[i]); /* andre, make it an err */
           }
           else {
             gvals[i] = gvals[i-4];
           }
         }
       }
     } /* end of  for(int i=1; i<9; i++) { */

   } /* end of  if(bintype==20) { */

/* -----------------------------------------------------------    */
/*  If outputting a text file, open it.... */

   if (Wname != NULL) {

     fpW = fopen(Wname, "w");
     if (fpW == NULL) err("**** Error opening the wfile output text file.");

/* Update/add whatever names/values are used by whatever options.                */

     for(int i=0; i<numgnams; i++) { 

       int ifound = 0;
       for(int j=0; j<numcases; j++) { /* reset it, if already in input text file    */ 
         if(strcmp(names[j],gnams[i]) == 0) {
           dfield[j] = gvals[i];  
           ifound = 1;
           break;
         }
       }
       if(ifound == 0) {               /* or add it */         
         dfield[numcasesout] = gvals[i];
         names[numcasesout] = ealloc1(7,1);
         strcpy(names[numcasesout],gnams[i]);
         forms[numcasesout] = ealloc1(5,1);
         strcpy(forms[numcasesout],"%.20g");
         numcasesout++;
       }

     } /* end of  for(int i=0; i<numgnams; i++) { */

     writekfile(fpW,names,forms,dfield,numcasesout);

   } /* end of  if (Wname != NULL) { */

/* andre, should fclose ? */

/* -----------------------------------------------------------    */

   if(intraces==0) return(0);

  if (!gettr(&tr))  err("can't get first trace");

  double dx;
  double dy;
  double tx;
  double ty;
  int icdp;
  int igi;
  int igc;

/* loop over traces   */ 

  do {

    if(ioffset==1) {
      dx = tr.sx;
      dy = tr.sy;
      tx = tr.gx;
      ty = tr.gy;

      if(tr.scalco > 1) { 
        dx *= tr.scalco;
        dy *= tr.scalco;
        tx *= tr.scalco;
        ty *= tr.scalco;
      }
      else if(tr.scalco < 0) { /* note that 1 and 0 retain values as-is */
        dx /= -tr.scalco;
        dy /= -tr.scalco;
        tx /= -tr.scalco;
        ty /= -tr.scalco;
      }
      tr.offset = sqrt((dx-tx)*(dx-tx) + (dy-ty)*(dy-ty));
    }

/* apply bintypes   */ 

    if(bintype==30) {

      dx = 0.5 * (double)(tr.sx + tr.gx);
      dy = 0.5 * (double)(tr.sy + tr.gy);

      if(tr.scalco > 1) { 
        dx *= tr.scalco;
        dy *= tr.scalco;
      }
      else if(tr.scalco < 0) { 
        dx /= -tr.scalco;
        dy /= -tr.scalco;
      }

/* Compute cdp,igi,igc from coordinates.                                */

      gridrawxycdpic(gvals,dx,dy,&icdp,&igi,&igc);

      tr.cdp = icdp;
      tr.igi = igi; 
      tr.igc = igc;

    } /* end of  if(bintype==30) { */
    else if(bintype==-30 || bintype==-31) {
      icdp = tr.cdp;
      gridcdpic(gvals,icdp,&igi,&igc);
      tr.igi = igi;
      tr.igc = igc;
    } 
    else if(bintype==-32) {
      igi = tr.igi;
      igc = tr.igc;
      gridiccdp(gvals,igi,igc,&icdp); 
      tr.cdp = icdp;
    } 

    if(bintype==-30 || bintype==-32) {

      gridicrawxy(gvals,igi,igc,&dx,&dy);
      gridicgridxy(gvals,igi,igc,&tx,&ty);

      if(tr.scalco > 1) { 
        dx /= tr.scalco;
        dy /= tr.scalco;
        tx /= tr.scalco;
        ty /= tr.scalco;
      }
      else if(tr.scalco < 0) { 
        dx *= -tr.scalco;
        dy *= -tr.scalco;
        tx *= -tr.scalco;
        ty *= -tr.scalco;
      }
      tr.sx = dx;
      tr.sy = dy;
      tr.gx = tx;
      tr.gy = ty;
    } 

    puttr(&tr);
    nproct++;

  } while (gettr(&tr));

  warn("Number of traces %d ",nproct);

  return 0;

}

/* -----------------------------------------------------------    */
/* -----------------------------------------------------------    */

void getCSV(char *textraw, char *textbeg, int maxtext, char rdel, 
            double *dfield, int *nspot, int numcases,   
            int ncount, int *comerr,int *morerr,int *numerr,int *nblank) {

  int nbeg = -1;
  int nfield = 0;
  int ineed  = 0;
  int igot;
  double dval;
  for(int n=0; n<maxtext; n++) {                         /* linux \n            windows \r */
    if(textraw[n] == rdel || textraw[n] == '\0' || textraw[n] == '\n' || textraw[n] == '\r') {
      if(nfield == nspot[ineed]) {
        dval = 1.1e308; 
        int nb = -1;
        if(n-nbeg-1 > 0) {
          strncpy(textbeg,textraw+nbeg+1,n-nbeg-1);
          textbeg[n-nbeg-1] = '\0'; /* so sscanf knows where to stop */
          int ib = -1;
          for (int m=0; m<n-nbeg-1; m++) {
            if(textbeg[m] != ' ') {
              nb = m;
              if(ib>-1) {
                *morerr = *morerr + 1;
                if(*morerr<4) {
                  warn("Error at record %d   two-numbers in field (%s)",ncount,textbeg);
                  if(*morerr==3) warn("Have 3 two-numbers in field warnings, no more will be printed.");
                }
                nb=-1;
                break;
              }
            }
            if(textbeg[m] == ' ' && nb>-1) ib = m;
          }  

          if(nb>-1) {
            igot = sscanf(textbeg,"%lf",&dval);  
            if(igot<1) {
              *numerr = *numerr + 1;
              if(*numerr<4) {
                warn("Error at record %d   field-unreadable as a number (%s)",ncount,textbeg);
                if(*numerr==3) warn("Have 3 field-unreadable warnings, no more will be printed.");
              }
            }
          }
        } /* end of  if(n-nbeg-1 > 0) { */
        if(nb<0) {
          *nblank = *nblank + 1;
           dval = 0.; 
        }
        dfield[ineed] = dval;
        ineed++;
        if(ineed>=numcases) break; 
      }
      if(textraw[n] == '\0' || textraw[n] == '\n' || textraw[n] == '\r') {
        if(ineed<numcases) {
          *comerr = *comerr + 1;
          if(*comerr<4) {
            warn("Error at record %d   Not-enough-commas in record to get all values",ncount);
            if(*comerr==3) warn("Have 3 Not-enough-comma warnings, no more will be printed.");
          }
        }
      }
      nbeg = n;
      nfield++;
    }
  } 
  return;
}

/* --------------------------- */
/* expects a string with no blanks and no tabs \t   */
void tparse(char *tbuf, char d, char **fields, int *numfields) { 
  int nbeg = -1;
  *numfields = 0;
  for(int n=0; ; n++) {
    if(tbuf[n] == d || tbuf[n] == '\0') {
      if(n-nbeg-1 > 0) {
        fields[*numfields] = ealloc1(n-nbeg-1,1);
        strncpy(fields[*numfields],tbuf+nbeg+1,n-nbeg-1);
      }
      else {
        fields[*numfields] = ealloc1(4,1);
        fields[*numfields][0] = 'n';
        fields[*numfields][1] = 'u';
        fields[*numfields][2] = 'l';
        fields[*numfields][3] = 'l';
/*      strncpy(fields[*numfields],"null",4);  makes compilor unhappy */
      }
      nbeg = n;
      *numfields = *numfields + 1;
    }
    if(tbuf[n] == '\0') break;
  }
 
}

/* Need an exact rectangle with 4 corner points A,B,C,D                              */
/* Corner A coordinates are used exactly as input.                                   */ 
/* The direction from corner A to corner B is determined by A-->B coordinates.       */ 
/* But then corner B coordinates are adjusted to an exact multiple of cell width WB. */ 
/* Input coordinates of A and C are used to compute the distance from A to C.        */ 
/* Right angle to A-->exactB gives direction for corner C (along line thru A).       */ 
/* Corner C is adjusted to an exact multiple of cell width WC away from A.           */ 
/* Note that this means input corner C is only used to decide how wide the           */
/* rectangle is, and which side of A-->B is output corner C located on.              */
/* Finally, compute corner D just so users can see it.                               */
/*                                                                                   */
/* Note: In the Object Orientated paradigm, the processed grid defintion values      */
/*       would all be hidden (in C++ they would be in private variables and only     */
/*       set-able using private methods). And so on. But SU is supposed to be        */
/*       simple, and a learning tool, so I have left everything exposed.             */
/*       If you change gvals after gridset, woe be to you.                           */

void gridset(double *gvals, int *errwarn) { 

     *errwarn = 0;

     if(gvals[10] <= 0.0) {
       *errwarn = 1;
       return;
     }
     if(gvals[11] <= 0.0) {
       *errwarn = 2;
       return;
     }

/* Reset corner B to be at an exact multiple distance of the B cell width.       */
/* Do not want compilors to optimize these computations, so use explicit (int).  */

     double dab = sqrt((gvals[2] - gvals[4])*(gvals[2] - gvals[4]) 
                     + (gvals[3] - gvals[5])*(gvals[3] - gvals[5])); 

     int nwb = (int) (dab/gvals[10]); 

     double dabwb = nwb * gvals[10];

     if(nwb<1) {
       *errwarn = 3;
       return;
     }

     nwb++; /* reset from number of intervals to number of cells */

     gvals[4] = gvals[2] + dabwb/dab * (gvals[4] - gvals[2]);
     gvals[5] = gvals[3] + dabwb/dab * (gvals[5] - gvals[3]);

/* Compute the input distance from A-->C.                                        */
/* And the exact multiple distance of the C cell width.                          */

     int nwc = 0;
     double dac = sqrt((gvals[2] - gvals[6])*(gvals[2] - gvals[6]) 
                     + (gvals[3] - gvals[7])*(gvals[3] - gvals[7])); 

     nwc = (int) (dac/gvals[11]); 

     if(nwc<1) *errwarn = -1; /*  Corner C is near A and is reset to A */

     nwc++; /* reset from number of intervals to number of cells */

     gvals[12] = nwb; /* set number of cells in A-->B direction */
     gvals[13] = nwc; /* set number of cells in A-->C direction */
     gvals[14] = 1.;  /* set first cdp number (may allow different, eventually) */
     gvals[15] = gvals[14] + nwb*nwc; /* set last cdp number */
     gvals[16] = (gvals[5] - gvals[3]) / dabwb; /* set sine   */
     gvals[17] = (gvals[4] - gvals[2]) / dabwb; /* set cosine */

/* Determine what side of A-->B the input C point is on.                         */

     double det = (gvals[4]-gvals[2]) * (gvals[7]-gvals[3])  /* (xb-xa) * (yc-ya)  */ 
                - (gvals[6]-gvals[2]) * (gvals[5]-gvals[3]); /* (xc-xa) * (yb-ya)  */

     gvals[1] = 1.;
     if(det<0.) gvals[1] = -1.;

/* Reset/set coordinates of corners B,C,D to values computed by gridicrawxy */
/* (sin,cosine were determined by floating point computations and we want   */
/*  the corners to be what gridicrawxy produces in case user chooses to use */
/*  the corner coordinates as range-limits and so on).                      */

     double rx;
     double ry;

     gridicrawxy(gvals,nwb,1,&rx,&ry);   /* get corner B XYs */
     gvals[4] = rx;
     gvals[5] = ry;

     gridicrawxy(gvals,1,nwc,&rx,&ry);   /* get corner C XYs */
     gvals[6] = rx;
     gvals[7] = ry;

     gridicrawxy(gvals,nwb,nwc,&rx,&ry); /* get corner D XYs */
     gvals[8] = rx;
     gvals[9] = ry;

}    

void gridrawxycdpic(double *gvals,double dx,double dy,int *icdp,int *igi,int *igc) {

/* Convert raw (real world) coordinates to cdp and igi,igc indexes.    */
/*                                                                     */
/* Inputs:                                                             */
/*   gvals is grid definition after processing by gridset              */
/*   dx    is x coordinate (raw,real world, after scalco applied)      */
/*   dy    is y coordinate (raw,real world, after scalco applied)      */
/* Outputs:                                                            */
/*   icdp is computed cell number. If input X,Y are outside the grid   */
/*        then icdp is output as -2147483645 but all other output      */
/*        values are still properly computed. So igi,igc will be less  */
/*        than 1 or greater than number of cells in their direction.   */
/*        And rx,ry will be negative or greater than grid extent.      */
/*        Note that rx,ry can be negative anyway when icdp is good     */
/*        since corner A is at the centre of cell 1.                   */
/*   igi  is computed cell index in A-->B direction (first igi is 1).  */
/*   igc  is computed cell index in A-->C direction (first igc is 1).  */

  dx =  dx - gvals[2];
  dy = (dy - gvals[3]) * gvals[1];

/* careful here, we want to rotate back to the axis from raw */

  double rx =  dx*gvals[17] + dy*gvals[16]; 
  double ry = -dx*gvals[16] + dy*gvals[17];

/* Compute the cell index number in b and c directions.     */
/* Add 0.5 because corner A is at the centre of first cell. */
/* Add 1 to start indexes at 1 (more user-friendly).        */

  *igi = (int) floor((rx / gvals[10] + 1.5)); 
  *igc = (int) floor((ry / gvals[11] + 1.5)); 

/* Compute the cell cdp number.                             */

  if(*igi<1 || *igi>gvals[12] || *igc<1 || *igc>gvals[13]) {
    *icdp = -2147483645;
  }
  else { 
    *icdp = gvals[14] + *igi-1 + (*igc-1) * gvals[12] + 0.1;
  }
}




void gridrawxygridxy(double *gvals,double dx,double dy,double *tx,double *ty) {

/* Convert raw (real world) coordinates to grid coordinates.           */
/*                                                                     */
/* Inputs:                                                             */
/*   gvals is grid definition after processing by gridset              */
/*   dx    is x coordinate (raw,real world, after scalco applied)      */
/*   dy    is y coordinate (raw,real world, after scalco applied)      */
/* Outputs:                                                            */
/*   tx    is x coordinate within grid.                                */
/*   ty    is y coordinate within grid.                                */
/*                                                                     */
/* Transforming from raw,real world XYs to grid XYs involves:          */
/*  - subtracting cornerA coordinates                                  */
/*  - rotating to cornerA-->B direction using grid sin,cosine          */
/*  - mirroring (multiplying gridY by -1 if C is on right of A-->B)    */
/*                                                                     */
/* Note: this function does not care if dx,dy are outside grid,        */
/*       it returns tx,ty anyway.                                      */

  dx =  dx - gvals[2];
  dy = (dy - gvals[3]) * gvals[1];

/* careful here, we want to rotate back to the axis from raw */

  *tx =  dx*gvals[17] + dy*gvals[16]; 
  *ty = -dx*gvals[16] + dy*gvals[17];

}


void gridgridxyrawxy(double *gvals,double dx,double dy,double *tx,double *ty) {

/* Convert grid coordinates to raw (real world) coordinates            */
/*                                                                     */
/* Inputs:                                                             */
/*   gvals is grid definition after processing by gridset              */
/*   dx    is x in grid coordiantes                                    */
/*   dy    is y in grid coordinates                                    */
/* Outputs:                                                            */
/*   tx    is x raw (real world) coordinate                            */
/*   ty    is y raw (real world) coordinate                            */
/*                                                                     */
/* Transforming from grid XYs to raw,real world XYs involves:          */
/*  - mirroring (multiplying gridY by -1 if C is on right of A-->B)    */
/*  - rotating from cornerA-->B direction using grid sin,cosine        */
/*  - adding cornerA coordinates                                       */
/*                                                                     */
/* Note: this function does not care if dx,dy are outside grid,        */
/*       it returns tx,ty anyway.                                      */

/* careful here, we want to rotate from axis back to raw      */

  *tx =  dx*gvals[17] - dy*gvals[16];             
  *ty = (dx*gvals[16] + dy*gvals[17]) * gvals[1];

  *tx = *tx + gvals[2];
  *ty = *ty + gvals[3];

}

void gridicgridxy(double *gvals,int igi,int igc,double *dx,double *dy) {

/* Convert grid indexes igi,igc to cell centre in grid XYs.            */
/*                                                                     */
/* Inputs:                                                             */
/*   gvals is grid definition after processing by gridset              */
/*   igi  is computed cell index in A-->B direction (first igi is 1).  */
/*   igc  is computed cell index in A-->C direction (first igc is 1).  */
/* Outputs:                                                            */
/*   dx   is cell centre X (in grid coordinates)                       */
/*   dy   is cell centre Y (in grid coordinates)                       */
/*                                                                     */
/* Note: this function does not care if igi or igc are outside grid,   */
/*       it returns the XYs anyway.                                    */

  *dx = (igi-1)*gvals[10]; 
  *dy = (igc-1)*gvals[11]; 

}

void gridicrawxy(double *gvals,int igi,int igc,double *dx,double *dy) {

/* Convert grid indexes igi,igc to cell centre in raw (real world) XYs.*/
/*                                                                     */
/* Inputs:                                                             */
/*   gvals is grid definition after processing by gridset              */
/*   igi  is computed cell index in A-->B direction (first igi is 1).  */
/*   igc  is computed cell index in A-->C direction (first igc is 1).  */
/* Outputs:                                                            */
/*   dx   is cell centre X (in raw, real world coordinates)            */
/*   dy   is cell centre Y (in raw, real world coordinates)            */
/*                                                                     */
/* Transforming from grid XYs to raw,real world XYs involves:          */
/*  - mirroring (multiplying gridY by -1 if C is on right of A-->B)    */
/*  - rotating from cornerA-->B direction using grid sin,cosine        */
/*  - adding cornerA coordinates                                       */
/*                                                                     */
/* Note: this function does not care if igi or igc are outside grid,   */
/*       it returns the XYs anyway.                                    */

  double rx = (igi-1)*gvals[10]; 
  double ry = (igc-1)*gvals[11]; 
       
/* careful here, we want to rotate away from the axis */

  *dx =  rx*gvals[17] - ry*gvals[16];
  *dy = (rx*gvals[16] + ry*gvals[17]) * gvals[1];

  *dx = *dx + gvals[2];
  *dy = *dy + gvals[3];

}

void gridiccdp(double *gvals,int igi,int igc,int *icdp) {

/* Convert grid indexes igi,igc to cdp (cell) number.                  */
/*                                                                     */
/* Inputs:                                                             */
/*   gvals is grid definition after processing by gridset              */
/*   igi  is computed cell index in A-->B direction (first igi is 1).  */
/*   igc  is computed cell index in A-->C direction (first igc is 1).  */
/* Outputs:                                                            */
/*   icdp is cell cdp number                                           */
/* Note: if igi or igc are outside grid, icdp=-2147483645 is returned  */

  if(igi<1 || igi>gvals[12] || igc<1 || igc>gvals[13]) {
    *icdp = -2147483645;
  }
  else { 
    *icdp = gvals[14] + igi-1 + (igc-1) * gvals[12] + 0.1;
  }

}



void gridcdpic(double *gvals,int icdp,int *igi,int *igc) {

/* Convert cdp (cell) number to grid indexes igi,igc.                  */
/*                                                                     */
/* Inputs:                                                             */
/*   gvals is grid definition after processing by gridset              */
/*   icdp is cell cdp number                                           */
/* Outputs:                                                            */
/*   igi  is computed cell index in A-->B direction (first igi is 1).  */
/*   igc  is computed cell index in A-->C direction (first igc is 1).  */
/* Note: if icdp is outside grid, igi,igc returned as -2147483645      */

  if(icdp<gvals[14] || icdp>gvals[15]) {
    *igi = -2147483645;
    *igc = -2147483645;
    return;
  }

  int ncdp = icdp - gvals[14];
  int nwb  = gvals[12];

  *igi = 1 + ncdp%nwb;
  *igc = 1 + ncdp/nwb;

}


void gridcheck(double *gvals, int icheck, int *errwarn) { 

/* Exercise grid functions using the 4 coorners.                       */
/*                                                                     */
/* Inputs:                                                             */
/*   gvals is grid definition after processing by gridset              */
/*   icheck=0 is no checking.                                          */
/*         =1 get the corner XYs from gvals and run them through       */
/*            various grid functions and print the results.            */
/*                                                                     */
/* The basic idea is to check sin,cosine and left/right coding errors. */
/* But users can also make use of the results if they are confused as  */
/* to what is going on. For extensive testing, start with a simple     */
/* situation where A-->B is roughtly aligned with X-axis with corner B */
/* having a larger X coordinate than A. And put corner C on left side  */
/* near corner A (a larger Y coordinate than A). Then put C on right   */
/* side of A. After that, use the output wfile and exchange A and B    */
/* and C and D (for both the previous runs). All these kinds of tests  */
/* should result in cdp,igi,igc making sense (not negative and so on). */
/* Similarly, the raw-to-grid and grid-to-raw coordinate conversion    */
/* functions should produce understandable results.                    */
/*                                                                     */
/* Outputs:                                                            */
/*   errwarn (always 0 in this version).                               */

  *errwarn = 0;
  if(icheck==0) return;

  int nwb = gvals[12] + 0.1;
  int nwc = gvals[13] + 0.1;
  double rx;
  double ry;
  double tx;
  double ty;
  int jcdp;
  int jigi;
  int jigc;
  int kcdp;
  int kigi;
  int kigc;

  gridicrawxy(gvals,1,1,&rx,&ry);   
  warn("gridicrawxy:     corner A raw  XYs= %.20f %.20f ",rx,ry);
  gridrawxygridxy(gvals,rx,ry,&tx,&ty);  
  warn("gridrawxygridxy: corner A grid XYs= %.20f %.20f ",tx,ty);
  gridicgridxy(gvals,1,1,&tx,&ty);   
  warn("gridicgridxy:    corner A grid XYs= %.20f %.20f ",tx,ty);
  gridgridxyrawxy(gvals,tx,ty,&rx,&ry);   
  warn("gridgridxyrawxy: corner A raw  XYs= %.20f %.20f ",rx,ry);
  gridrawxycdpic(gvals,rx,ry,&jcdp,&jigi,&jigc); 
  gridcdpic(gvals,jcdp,&kigi,&kigc);
  gridiccdp(gvals,jigi,jigc,&kcdp); 
  warn("gridrawxycdpic:          corner A cdp,igi,igc = %d %d %d ",jcdp,jigi,jigc);
  warn("gridcdpic and gridiccdp: corner A cdp,igi,igc = %d %d %d ",kcdp,kigi,kigc);


  gridicrawxy(gvals,nwb,1,&rx,&ry);   /* get corner B XYs */
  warn("gridicrawxy:     corner B raw  XYs= %.20f %.20f ",rx,ry);
  gridrawxygridxy(gvals,rx,ry,&tx,&ty);   /* get grid B XYs */
  warn("gridrawxygridxy: corner B grid XYs= %.20f %.20f ",tx,ty);
  gridicgridxy(gvals,nwb,1,&tx,&ty);   /* get grid B XYs */
  warn("gridicgridxy:    corner B grid XYs= %.20f %.20f ",tx,ty);
  gridgridxyrawxy(gvals,tx,ty,&rx,&ry);   /* get corner B XYs */
  warn("gridgridxyrawxy: corner B raw  XYs= %.20f %.20f ",rx,ry);
  gridrawxycdpic(gvals,rx,ry,&jcdp,&jigi,&jigc); /* get icdp,igi,igc from corner A XYs */
  gridcdpic(gvals,jcdp,&kigi,&kigc); /* get igi,igc from cdp */
  gridiccdp(gvals,jigi,jigc,&kcdp); /* get cpd from igi,igc */
  warn("gridrawxycdpic:          corner B cdp,igi,igc = %d %d %d ",jcdp,jigi,jigc);
  warn("gridcdpic and gridiccdp: corner B cdp,igi,igc = %d %d %d ",kcdp,kigi,kigc);


  gridicrawxy(gvals,1,nwc,&rx,&ry);   /* get corner C XYs */
  warn("gridicrawxy:     corner C raw  XYs= %.20f %.20f ",rx,ry);
  gridrawxygridxy(gvals,rx,ry,&tx,&ty);   /* get grid C XYs */
  warn("gridrawxygridxy: corner C grid XYs= %.20f %.20f ",tx,ty);
  gridicgridxy(gvals,1,nwc,&tx,&ty);   /* get grid C XYs */
  warn("gridicgridxy:    corner C grid XYs= %.20f %.20f ",tx,ty);
  gridgridxyrawxy(gvals,tx,ty,&rx,&ry);   /* get corner C XYs */
  warn("gridgridxyrawxy: corner C raw  XYs= %.20f %.20f ",rx,ry);
  gridrawxycdpic(gvals,rx,ry,&jcdp,&jigi,&jigc); /* get icdp,igi,igc from corner A XYs */
  gridcdpic(gvals,jcdp,&kigi,&kigc); /* get igi,igc from cdp */
  gridiccdp(gvals,jigi,jigc,&kcdp); /* get cpd from igi,igc */
  warn("gridrawxycdpic:          corner C cdp,igi,igc = %d %d %d ",jcdp,jigi,jigc);
  warn("gridcdpic and gridiccdp: corner C cdp,igi,igc = %d %d %d ",kcdp,kigi,kigc);


  gridicrawxy(gvals,nwb,nwc,&rx,&ry); /* get corner D XYs */
  warn("gridicrawxy:     corner D raw  XYs= %.20f %.20f ",rx,ry);
  gridrawxygridxy(gvals,rx,ry,&tx,&ty);   /* get grid D XYs */
  warn("gridrawxygridxy: corner D grid XYs= %.20f %.20f ",tx,ty);
  gridicgridxy(gvals,nwb,nwc,&tx,&ty); /* get grid D XYs */
  warn("gridicgridxy:    corner D grid XYs= %.20f %.20f ",tx,ty);
  gridgridxyrawxy(gvals,tx,ty,&rx,&ry);   /* get corner D XYs */
  warn("gridgridxyrawxy: corner D raw  XYs= %.20f %.20f ",rx,ry);
  gridrawxycdpic(gvals,rx,ry,&jcdp,&jigi,&jigc); /* get icdp,igi,igc from corner A XYs */
  gridcdpic(gvals,jcdp,&kigi,&kigc); /* get igi,igc from cdp */
  gridiccdp(gvals,jigi,jigc,&kcdp); /* get cpd from igi,igc */
  warn("gridrawxycdpic:          corner D cdp,igi,igc = %d %d %d ",jcdp,jigi,jigc);
  warn("gridcdpic and gridiccdp: corner D cdp,igi,igc = %d %d %d ",kcdp,kigi,kigc);

}    





void readkfile(FILE *fpR, cwp_String *names, cwp_String *forms, double *dfield, int *numcasesout) {

  int numcases = 0;

  int maxtext = 10001;
  char textraw[10001]; /* fgets puts a \0 after contents */
  char textbeg[10001]; /* so this size wastes memory but not time */
  char textraw2[10001];
  char textfront[10];  

  int *nspot = NULL;

  char rdel = ',';

  cwp_String Rid  =NULL;  /* rejection id for records             */
  Rid = ealloc1(1,1);
  strcpy(Rid,"K");
        
  int num_names = 0;
  int num_forms = 0;

  int num_c_su_names = 0;
  int num_c_su_forms = 0;
  int num_c_su_setid = 0;
  int read_names = 0;
  int read_forms = 0;

  while (fgets(textraw, maxtext, fpR) != NULL) { /* read a line */

/* Stop this looping? Sometimes, it really will just loop through to last record. */
                   
    if(read_names==-1 && read_forms==-1) break;

/* Remove all blanks and tabs because tparse is not designed to handle them.      */

    int tsize = 0;
    for(int n=0; n<maxtext; n++) { /*   linux \n            windows \r */
      if(textraw[n] == '\0' || textraw[n] == '\n' || textraw[n] == '\r') break;
      if(textraw[n] != ' ' && textraw[n] != '\t') {
        textbeg[tsize] = textraw[n];
        tsize++;
      }
    }

    for(int n=0; n<sizeof(textbeg); n++) textbeg[n] = tolower(textbeg[n]);
    if(strncmp(textbeg,"c_su_setid",10) == 0) num_c_su_setid++;
    if(strncmp(textbeg,"c_su_names",10) == 0) num_c_su_names++;
    if(strncmp(textbeg,"c_su_forms",10) == 0) num_c_su_forms++;

    if(read_names>0) {
      textbeg[tsize] = '\0'; 
      tparse(textbeg, ',', names, &num_names) ; 
      read_names = -1;
    }

    if(strncmp(textbeg,"c_su_names",10) == 0) read_names = 1;

    if(read_forms>0) {
      textbeg[tsize] = '\0';  
      tparse(textbeg, ',', forms, &num_forms) ; 
      read_forms = -1;
    }

    if(strncmp(textbeg,"c_su_forms",10) == 0) read_forms = 1;

  } /* end of while (fgets(textraw,..... */
  
  fseek(fpR, 0L, SEEK_SET); /* reposition file to beginning record */ 

  int ierr = 0;
  if(num_c_su_names>1) {
    ierr = 1;
    warn("**** Error: Text file has more than one C_SU_NAMES parameter record.");
  }
  else if(num_c_su_names==0) {
    ierr = 1;
    warn("**** Error: Text file has no C_SU_NAMES parameter record.");
    warn("****        Add it, or specify names= on command line.");
    }

  if(num_c_su_forms>1) {
    ierr = 1;
    warn("**** Error: Text file has more than one C_SU_FORMS parameter record.");
  }
  else if(num_c_su_forms==0) {
    ierr = 1;
    warn("**** Error: Text file has no C_SU_FORMS parameter record.");
    warn("****        Add it, or specify forms= on command line.");
  }

  if(num_c_su_setid>1) {
    ierr = 1;
    warn("**** Error: Text file has more than one C_SU_SETID record.");
  }
  else if(num_c_su_setid==0) {
    ierr = 1;
    warn("**** Error: Text file has no C_SU_SETID record.");
    warn("****        Add it, or specify setid= on command line.");
  }
  
  if(num_forms != num_names) {
    ierr = 1;
    warn("**** Error: Text file has different number of values on C_SU_NAMES and C_SU_FORMS.");
  }

  if(ierr>0) { /* andre, typo in sutoolcsv? */
    err("**** Error: Text file has wrong, duplicate, or missing C_SU_NAMES or C_SU_FORMS or C_SU_SETID records.");
  }

  nspot = calloc(num_names,sizeof(int));
  if(nspot == NULL) err("**** Unable to allocate memory.");

/* Do not try to get c_su_id from record, it is S,R,X or whatever.      */

  if(strncmp(names[0],"c_su_id",7) == 0) names[0] = "null";

  for (int n=0; n<num_names; n++) {
    if(strncmp(names[n],"null",4) != 0) {
      for (int m=n+1; m<num_names; m++) {
        if(strcmp(names[n],names[m]) == 0) {  
          err("**** Error: Name  %s  exists at least twice in the names list.",names[n]);
        }
      }
    }
  }

/* ----------------------------------------------------- */
/* ----------------------------------------------------- */

  for(int i=0; i<num_names;i++) { 
    if(strncmp(names[i],"null",4) != 0) { /* actually removes c_su_id also */
      names[numcases] = names[i];
      forms[numcases] = forms[i];
      nspot[numcases] = i;
      numcases++; 
    }
  }

  memset(textraw,'\0',10001);
  memset(textraw2,'\0',10001);
  memset(textbeg,'\0',10001);

  int ncount = 0;
  int lenerr = 0;
  int comerr = 0;
  int morerr = 0;
  int numerr = 0;
  int nblank = 0;
  int nextrow = 0;

  while (fgets(textraw, maxtext, fpR) != NULL) { /*read a line*/
    ncount++;
    for(int n=0; n<10; n++) textfront[n] = tolower(textraw[n]);
    if(strncmp(textfront,"c_su",4) == 0 || nextrow==1) {
      nextrow = 0; 
      if(strncmp(textfront,"c_su_names",10) == 0 || 
         strncmp(textfront,"c_su_forms",10) == 0) nextrow = 1;
    }
    else {
      if(strncmp(textraw,Rid,1) == 0) { /* Rid compare is case-sensitive */

        getCSV(textraw, textbeg, maxtext, rdel, 
               dfield, nspot, numcases,
               ncount, &comerr,&morerr,&numerr,&nblank);

        break;

      }  
    }
  }

  if(nblank>0) warn("Total all-blank fields: %d. Assumed zero for all.",nblank);

  if(numerr>0) warn("Total Field-unreadable as a number:        %d (will error-halt SUGEOMCSV)",numerr);
  if(morerr>0) warn("Total Two-numbers in one field:            %d (will error-halt SUGEOMCSV)",morerr);
  if(lenerr>0) warn("Total Record-too-short to get all values:  %d (will error-halt SUGEOMCSV)",lenerr);
  if(comerr>0) warn("Total Not-enough-commas to get all values: %d (will error-halt SUGEOMCSV)",comerr);

  *numcasesout = numcases;

}    


void writekfile(FILE *fpW, cwp_String *names, cwp_String *forms, double *dfield, int numcasesout) {

  char textraw[10001]; /* fgets puts a \0 after contents */
  char textbeg[10001]; /* so this size wastes memory but not time */
  char textraw2[10001];
        
  cwp_String Rid  =NULL;  /* rejection id for records             */
  Rid = ealloc1(1,1);
  strcpy(Rid,"K");
        
  memset(textraw,'\0',10001);
  memset(textraw2,'\0',10001);
  memset(textbeg,'\0',10001);

  int mspot;
  int mleng;

  strcpy(textraw,"C_SU_SETID,");
  strcpy(textraw+11,Rid);
  textraw[12] = '\n';
  textraw[13] = '\0';
  fputs(textraw,fpW);

/* write the forms records                       */ 

  strcpy(textraw,"C_SU_FORMS");
  textraw[10] = ' ';
  textraw[10] = '\n';
  textraw[11] = '\0';
  fputs(textraw,fpW);

  mspot = 8;
  strcpy(textraw,"C_SU_ID,");

  for(int i=0; i<numcasesout; i++) { 
    mleng = strlen(forms[i]);
    strncpy(textraw+mspot,forms[i],mleng);
    mspot += mleng;
    textraw[mspot] = ',';
    mspot++;
  } /* end of  for(int i=0; i<numcases; i++) { */ 
  textraw[mspot-1] = '\n';
  textraw[mspot  ] = '\0';
  fputs(textraw,fpW);

/* write the names records                       */ 

  strcpy(textraw,"C_SU_NAMES");
  textraw[10] = '\n';
  textraw[11] = '\0';
  fputs(textraw,fpW);

  mspot = 8;
  strcpy(textraw,"C_SU_ID,");

  for(int i=0; i<numcasesout; i++) { 
    mleng = strlen(names[i]);
    strncpy(textraw+mspot,names[i],mleng);
    mspot += mleng;
    textraw[mspot] = ',';
    mspot++;
  } /* end of  for(int i=0; i<numcases; i++) { */ 
  textraw[mspot-1] = '\n';
  textraw[mspot  ] = '\0';
  fputs(textraw,fpW);

  strncpy(textraw2,Rid,1);
  int mhere = 1;
  for(int ineed=0; ineed<numcasesout; ineed++) { /* nspot not used, why add comma for no value? */ 
    if(dfield[ineed]<1.e308) {
      sprintf(textbeg,forms[ineed],dfield[ineed]);
      int mfill = strlen(textbeg);
      textraw2[mhere] = ','; /* always comma for output */
      strncpy(textraw2+mhere+1,textbeg,mfill);
      mhere += mfill+1;
    }
    else {
      textraw2[mhere] = ',';
      textraw2[mhere+1] = '*';
      mhere += 2;
    }
  }
  textraw2[mhere]   = '\n';
  textraw2[mhere+1] = '\0';
  fputs(textraw2,fpW);

}    


