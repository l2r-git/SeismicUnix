/* Copyright (c) Colorado School of Mines, 2021.*/
/* All rights reserved.                       */

/* SUTOOLCSV: $Revision: 1.00 $ ; $Date: 2021/05/15 00:09:00 $        */

#include <stdio.h>
#include <string.h>
#include "math.h"

#include "su.h"
#include "segy.h"

segy tr;

void getCSV(char *textraw, char *textbeg, int maxtext, char rdel, 
            double *dfield, int *nspot, int numcases,   
            int ncount, int *comerr,int *morerr,int *numerr,int *nblank);
void tparse(char *tbuf, char d, char **fields, int *numfields) ; 


/*********************** self documentation **********************/
char *sdoc[] = {
"                          ",
" SUBINCSV - Various CDP and other binning options                         ",
"                                                                          ",
" subincsv  rfile=in.csv wfile=out.csv                                     ",
"                                                                          ",
" Parameter overview:                                                      ",
"                                                                          ",
"       rfile= read csv file                                               ", 
"       wfile= write csv file                                              ", 
"                                                                          ",
" ***********************************************************              ",
"   To output this documentation:  subincsv 2> bindoc.txt                  ",
" ***********************************************************              ",
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
/* Trace header fields accessed: sx,sy,gx,gy                      */
/*                                                                */
/**************** end self doc ************************************/

int main(int argc, char **argv) {

   cwp_String Rname=NULL;  /* text file name for values            */
   cwp_String Rid  =NULL;  /* rejection id for records             */
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
   double dabwb = 0.0;
   double ds = 0.0;
   double dc = 0.0;
   int icdpf = 1;
   int nwb = 0;
   int nwc = 0;

   int maxtext = 10001;
   char textraw[10001]; /* fgets puts a \0 after contents */
   char textbeg[10001]; /* so this size wastes memory but not time */
   char textraw2[10001];
   char textfront[10];  
        
   int *nspot = NULL;

   int nproct = 0; // number of processed traces

   /* Initialize */
   initargs(argc, argv);
   requestdoc(1);

   char rdel = ',';

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

   int makecdp = -1;
   char *tcdp;
   getparstring("makecdp", &tcdp);
   if(tcdp != NULL) {
     for (int n=0; n<strlen(tcdp); n++) tcdp[n] = tolower(tcdp[n]);
     if(strlen(tcdp)==4 && strcmp(tcdp,"grid") == 0) makecdp = 3;
     else if(strlen(tcdp)==7 && strcmp(tcdp,"point2d") == 0) makecdp = 2;
     else err("**** Error: makecdp= option not recognized.");
   }
  
/* Cycle over rfile records to get some C_SU_ parameters? */ 

   int numcases = 0;
   int lenid = 0;  /* when id is set, this is length (K,R,X have length 1) */
   int isetid = 1;

   if (Rname == NULL) {
     lenid = 1;
     Rid = ealloc1(1,1);
     strcpy(Rid,"K");
   }
   else {

     fpR = fopen(Rname, "r");

     int num_names = 0;
     int num_forms = 0;

     int num_c_su_names = 0;
     int num_c_su_forms = 0;
     int num_c_su_setid = 0;
     int read_names = 0;
     int read_forms = 0;

     while (fgets(textraw, maxtext, fpR) != NULL) { /* read a line */

/* Stop this looping? Sometimes, it really will just loop through to last record. */
                   
       if(read_names==-1 && read_forms==-1 && lenid>0) break;

/* Remove all blanks and tabs because tparse is not designed to handle them.      */

       int tsize = 0;
       for(int n=0; n<maxtext; n++) { /*   linux \n            windows \r */
         if(textraw[n] == '\0' || textraw[n] == '\n' || textraw[n] == '\r') break;
         if(textraw[n] != ' ' && textraw[n] != '\t') {
           textbeg[tsize] = textraw[n];
           tsize++;
         }
       }

       if(lenid==0) { /* note the id itself is case-sensitive. */
         for(int n=0; n<10; n++) textbeg[n] = tolower(textbeg[n]);
         if(strncmp(textbeg,"c_su_setid",10) == 0) {
           char *ids[99];
           int nids;
           textbeg[tsize] = '\0';
           tparse(textbeg, ',', ids, &nids) ; 
           Rid = ids[1];
           lenid = strlen(Rid);
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

/* Resolve setid options.  */

     if(lenid>3 && Rid[0]=='"' && Rid[lenid-1]=='"') {
       for(int n=0; n<lenid-1; n++) Rid[n] = Rid[n+1];
       lenid -= 2;
     }
     else {
       for(int n=0; n<lenid; n++) Rid[n] = toupper(Rid[n]);
     }

     isetid = 1;
     if(lenid==3) {
       char text3[3]; 
       for(int n=0; n<3; n++) text3[n] = tolower(Rid[n]);
       if(strncmp(text3,"any",3) == 0) lenid = 0;       
     }
     else if(lenid==4) {
       char text4[4]; 
       for(int n=0; n<4; n++) text4[n] = tolower(Rid[n]);
       if(strncmp(text4,"none",4) == 0) {
         lenid = 0;  
         isetid = 0;  
       }
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
         if(lenid<1 || strncmp(textraw,Rid,lenid) == 0) { /* Rid compare is case-sensitive */

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

   } /* end of  else for if (Rname == NULL) { */

/* ----------------------------------------------------------- */

   int numcasesout = numcases;

   if(makecdp==3) {

     for(int i=0; i<20; i++) {
       gvals[i] = -1.1e308;
       gnams[i] = ealloc1(7,1); 
     }

     strcpy(gnams[0],"grid_wb");
     strcpy(gnams[1],"grid_wc");
     strcpy(gnams[2],"grid_xa");
     strcpy(gnams[3],"grid_ya");
     strcpy(gnams[4],"grid_xb");
     strcpy(gnams[5],"grid_yb");
     strcpy(gnams[6],"grid_xc");
     strcpy(gnams[7],"grid_yc");
     strcpy(gnams[8],"grid_xd");
     strcpy(gnams[9],"grid_yd");
     strcpy(gnams[10],"grid_sb");
     strcpy(gnams[11],"grid_cb");
     strcpy(gnams[12],"grid_nb");
     strcpy(gnams[13],"grid_nc");
     strcpy(gnams[14],"grid_fp");
     strcpy(gnams[15],"grid_lp");

     for(int i=0; i<8; i++) { /* do not read-in corner D or any other grid stuff  */ 
       if(!getpardouble(gnams[i],gvals+i)) { 
         for(int j=0; j<numcases; j++) { 
           if(strcmp(names[j],gnams[i]) == 0) gvals[i] = dfield[j]; /* andre, lower case? */ 
         }
         if(gvals[i] < -1.e308) {
           gvals[i] = i+100;
           err("**** Error makedcdp=grid and %s not found.",gnams[i]);
         }
       }
     } /* end of  for(int i=0; i<8; i++) { */

/* Need an exact rectangle with 4 corner points A,B,C,D                              */
/* Corner A coordinates are used exactly as input.                                   */ 
/* The direction from corner A to corner B is determined by A-->B coordinates.       */ 
/* But then corner B coordinates are adjusted to an exact multiple of cell width WB. */ 
/* Input coordinates of A and C are used to compute the distance from A to C.        */ 
/* Right angle to A-->exactB gives direction for corner C (along line thru A).       */ 
/* Corner C is adjusted to an exact multiple of cell width WC away from A.           */ 
/* Compute corner D just so users can see it.                                        */
/* Note that this means input corner C is only used to decide how wide the           */
/* rectangle is, and which side of A-->B is output corner C located on.              */

   if(gvals[0] <= 0.0) err ("**** Error. The grid_wb cell width must be positive.");
   if(gvals[1] <= 0.0) err ("**** Error. The grid_wc cell width must be positive.");

/* Reset corner B to be at an exact multiple distance of the B cell width.       */
/* Do not want compilors to optimize these computations, so use explicit (int).  */

   double dab = sqrt((gvals[2] - gvals[4])*(gvals[2] - gvals[4]) 
                   + (gvals[3] - gvals[5])*(gvals[3] - gvals[5])); 

   nwb = (int) (dab/gvals[0]); 

   dabwb = nwb * gvals[0];

   if(nwb<1) err ("**** Error. Corner B is within grid_wb cell width of corner A.");

   nwb++; /* reset from number of intervals to number of cells */

   gvals[4] = gvals[2] + dabwb/dab * (gvals[4] - gvals[2]);
   gvals[5] = gvals[3] + dabwb/dab * (gvals[5] - gvals[3]);

/* set sine, set cosine, set first cdp number */

   ds = (gvals[5] - gvals[3]) / dabwb; 
   dc = (gvals[4] - gvals[2]) / dabwb; 
   icdpf = 1;                          

/* Compute the input distance from A-->C.                                        */
/* And the exact multiple distance of the C cell width.                          */

   double dac = sqrt((gvals[2] - gvals[6])*(gvals[2] - gvals[6]) 
                   + (gvals[3] - gvals[7])*(gvals[3] - gvals[7])); 

   nwc = (int) (dac/gvals[1]); 

   if(nwc<1) warn ("**** Corner C is so close to corner A that it has been reset to A (forming a 2D grid).");

   double dacwc = nwc * gvals[1];

   nwc++; /* reset from number of intervals to number of cells */

// strcpy(gnams[0],"grid_wb");
// strcpy(gnams[1],"grid_wc");
// strcpy(gnams[2],"grid_xa");
// strcpy(gnams[3],"grid_ya");
// strcpy(gnams[4],"grid_xb");
// strcpy(gnams[5],"grid_yb");
// strcpy(gnams[6],"grid_xc");
// strcpy(gnams[7],"grid_yc");
// strcpy(gnams[8],"grid_xd");
// strcpy(gnams[9],"grid_yd");
// strcpy(gnams[10],"grid_sb");
// strcpy(gnams[11],"grid_cb");
// strcpy(gnams[12],"grid_nb");
// strcpy(gnams[13],"grid_nc");
// strcpy(gnams[14],"grid_fp");
// strcpy(gnams[15],"grid_lp");

/* Determine what side of A-->B the input C point is on.                         */

   double det = (gvals[4]-gvals[2]) * (gvals[7]-gvals[3])  /* (xb-xa) * (yc-ya)  */ 
              - (gvals[6]-gvals[2]) * (gvals[5]-gvals[3]); /* (xc-xa) * (yb-ya)  */

/* Rotate the A to newB vector by +/- 90 degrees and scale it to dacwc length    */
/* then add it to A to get newC. Same for output corner D except add it to newB. */

   if(det>=0.) {
     gvals[6] = gvals[2] + dacwc/dabwb * (gvals[3] - gvals[5]);
     gvals[7] = gvals[3] + dacwc/dabwb * (gvals[4] - gvals[2]);                       
     gvals[8] = gvals[4] + dacwc/dabwb * (gvals[3] - gvals[5]);
     gvals[9] = gvals[5] + dacwc/dabwb * (gvals[4] - gvals[2]);                       
   }
   else {
     gvals[6] = gvals[2] + dacwc/dabwb * (gvals[5] - gvals[3]);
     gvals[7] = gvals[3] + dacwc/dabwb * (gvals[2] - gvals[4]);
     gvals[8] = gvals[4] + dacwc/dabwb * (gvals[5] - gvals[3]);
     gvals[9] = gvals[5] + dacwc/dabwb * (gvals[2] - gvals[4]);
   }

   gvals[10] = ds;
   gvals[11] = dc;
   gvals[12] = nwb;
   gvals[13] = nwc;
   gvals[14] = icdpf;
   gvals[15] = icdpf + nwb*nwc;
   
/* Update/add the grid definition to the output text file.                       */

     for(int i=0; i<16; i++) { 

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

     } /* end of  for(int i=0; i<10; i++) { */

   } /* end of  if(makecdp==3) { */

/* -----------------------------------------------------------    */
/*  If outputting a text file, open it.... */

   if (Wname != NULL) {

     fpW = fopen(Wname, "w");
     if (fpW == NULL) err("**** Error opening the wfile output text file.");

     memset(textraw,'\0',10001);
     memset(textraw2,'\0',10001);
     memset(textbeg,'\0',10001);

     int mspot;
     int mleng;

     if(lenid>0) { 
       strcpy(textraw,"C_SU_SETID,");
       strcpy(textraw+11,Rid);
       textraw[11+lenid] = '\n';
       textraw[12+lenid] = '\0';
     }
     else if(isetid==1) { /* andre, remove ANY,NONE options */
       strcpy(textraw,"C_SU_SETID,ANY");
       textraw[14] = '\n';
       textraw[15] = '\0';
     }
     else {
       strcpy(textraw,"C_SU_SETID,NONE");
       textraw[15] = '\n';
       textraw[16] = '\0';
     }
     fputs(textraw,fpW);

/* write the forms records                       */ 

     strcpy(textraw,"C_SU_FORMS");
     textraw[10] = ' ';
     textraw[10] = '\n';
     textraw[11] = '\0';
     fputs(textraw,fpW);

     mspot = 0;
     if(isetid==1) {
       mspot = 8;
       strcpy(textraw,"C_SU_ID,");
     }

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

     mspot = 0;
     if(isetid==1) {
       mspot = 8;
       strcpy(textraw,"C_SU_ID,");
     }

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

     strncpy(textraw2,Rid,lenid);
     int mhere = lenid;
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
     if(isetid==0) textraw2[0] = ' '; /* rub out the leading comma */
     textraw2[mhere]   = '\n';
     textraw2[mhere+1] = '\0';
     fputs(textraw2,fpW);

   } /* end of  if (Wname != NULL) { */

/* -----------------------------------------------------------    */

   if(intraces==0) return(0);

// strcpy(gnams[0],"grid_wb");
// strcpy(gnams[1],"grid_wc");
// strcpy(gnams[2],"grid_xa");
// strcpy(gnams[3],"grid_ya");
// strcpy(gnams[4],"grid_xb");
// strcpy(gnams[5],"grid_yb");
// strcpy(gnams[6],"grid_xc");
// strcpy(gnams[7],"grid_yc");
// strcpy(gnams[8],"grid_xd");
// strcpy(gnams[9],"grid_yd");
// strcpy(gnams[10],"grid_sb");
// strcpy(gnams[11],"grid_cb");
// strcpy(gnams[12],"grid_nb");
// strcpy(gnams[13],"grid_nc");
// strcpy(gnams[14],"grid_fp");
// strcpy(gnams[15],"grid_lp");

  if (!gettr(&tr))  err("can't get first trace");

/* loop over traces   */ 

  do {

    if(makecdp==3) {

      double sx = tr.sx;
      double sy = tr.sy;
      double gx = tr.gx;
      double gy = tr.gy;

      if(tr.scalco > 1) { /* andre, or the other way around */
        sx /= tr.scalco;
        sy /= tr.scalco;
        gx /= tr.scalco;
        gy /= tr.scalco;
      }
      else if(tr.scalco < 0) { /* note that 1 and 0 retain values as-is */
        sx *= -tr.scalco;
        sy *= -tr.scalco;
        gx *= -tr.scalco;
        gy *= -tr.scalco;
      }

      double dx = (sx+gx) * 0.5 - gvals[2];
      double dy = (sy+gy) * 0.5 - gvals[3];
 
      double rx = dx*dc - dy*ds;
      double ry = dx*ds + dy*dc;


//  warn("andre0005   %f %f %f %f %f %f %f %f %f %f",sx,sy,gx,gy,ds,dc,dx,dy,rx,ry);


      int ib = (int) (rx / gvals[0] + 0.5000000000000001);  
      int ic = (int) (ry / gvals[1] + 0.5000000000000001); 

      tr.cdp = icdpf + ib + ic*nwb;
      tr.igi = 1 + ib; /* decided to start indexes at 1 (more user-friendly)  */
      tr.igc = 1 + ic;

//  warn("andre0006   %d %d %d %d %d ",nwb,nwc,ib,ic,icdp);

    } /* end of  if(makecdp==3) { */

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
