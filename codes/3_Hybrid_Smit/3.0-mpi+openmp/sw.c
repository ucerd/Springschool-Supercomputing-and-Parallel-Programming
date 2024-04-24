/*-------------------------------------------------------------

Copyright (C) 2000 Peter Clote. 
All Rights Reserved.

Permission to use, copy, modify, and distribute this
software and its documentation for NON-COMMERCIAL purposes
and without fee is hereby granted provided that this
copyright notice appears in all copies.


THE AUTHOR AND PUBLISHER MAKE NO REPRESENTATIONS OR
WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE AUTHORS
AND PUBLISHER SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED
BY LICENSEE AS A RESULT OF USING, MODIFYING OR DISTRIBUTING
THIS SOFTWARE OR ITS DERIVATIVES.

-------------------------------------------------------------*/


/*************************************************
	Program: smithWaterman.c
	Peter Clote, 11 Oct 2000

Program for local sequence alignment, using the Smith-Waterman
algorithm and assuming a LINEAR gap penalty.
A traceback is used to determine the alignment, and
is determined by choosing that direction, for
which S(i-1,j-1)+sigma(a_i,b_j), S(i-1,j)+Delta and 
S(i,j-1)+Delta is maximum, i.e.  for which 

                    _
                   |
                   | H(i-1,j-1) + sigma(a_i,b_j)  (diagonal)
H(i,j) =  MAX of   | H(i-1,j)   + delta           (up)
                   | H(i,j-1)   + delta           (left)
                   | 0
                   | _


is a maximum.

*************************************************/

/*begin AMPP
 AMPP: SOME BUGS HAS BEEN SOLVED
       SOME CODE HAS BEEN IMPROVED 
      Code added or re-coded starts with begin AMPP and ends with
      end AMPP
  end AMPP 
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>   	// character handling
#include <stdlib.h>     // def of RAND_MAX 
#include "mpi.h"

/* begin AMPP */
   /* Just a note:                       */
   /* N must be the size of the arrays   */
   /* Here is assume that the two arrays */
   /* have the same size                 */


#define MAX_SEQ 50

#define CHECK_NULL(_check) {\
   if ((_check)==NULL) \
      fprintf(stderr, "Null Pointer allocating memory\n");\
   }
  
/* end AMPP */

#define AA 20           // number of amino acids
#define MAX2(x,y)     ((x)<(y) ? (y) : (x))
#define MAX3(x,y,z)   (MAX2(x,y)<(z) ? (z) : MAX2(x,y))

// function prototypes
void error(char *);		/** error handling */

/* begin AMPP*/
int char2AAmem[256];
int AA2charmem[AA];
void initChar2AATranslation(void);
/* end AMPP */

#define MASTER 0

long usecs (void)
{
  struct timeval t;

  gettimeofday(&t,NULL);
  return t.tv_sec*1000000+t.tv_usec;
}

main(int argc, char *argv[]) {

	// variable declarations
	FILE * in1, *in2, *pam;
	char ch;
	int temp;
	int i,j,k,tempi,tempj,x,y,diag,down,right,DELTA;
	int topskip,bottomskip;
	char *aout,*bout;
	int Aend,Bend,Abegin,Bbegin;
	int max, Max, xMax, yMax;	
		// Max is first found maximum in similarity matrix H
		// max is auxilliary to determine largest of
		// diag,down,right, xMax,yMax are h-coord of Max
	short *a, *b;
	int *hptr;
	int **h, **h_out;
	int sim[AA][AA];		// PAM similarity matrix
	short **xTraceback, **yTraceback;
	short *xTracebackptr, *yTracebackptr;
        int N,B,I;
        int nc;
	int numprocs,myid;
        int l,m;         
	int stride;       
	int num_block_rows, num_block_cols;
	int num_cols, num_rows;
	int num_steps;
	int CHUNK_SIZE = 1;
 
        long t_start,t_end;
	double time;

	MPI_Status status;
	MPI_Request requestr;
	MPI_Request requests;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	

/* begin AMPP */
	/**** Error handling for input file ****/	
	if ( argc != 9) {
	     	fprintf(stderr,"%s error\n",argv[0]);
		exit(1);
	}
        else 
        {
	   /***** Initialization of input file and pair array **/
	    in1   = fopen(argv[1],"r");
    	    in2   = fopen(argv[2],"r");
	    pam   = fopen(argv[3],"r");
	    DELTA = atoi(argv[4]);
            N = atoi(argv[5]);
	    B = atoi(argv[6]);
	    I = atoi(argv[7]);
		 CHUNK_SIZE = atoi(argv[8]);

        } 
/* end AMPP */
	
		/**** Error handling for input****/
   	if (B < 2 && B > N-1) {
            fprintf(stderr,"B is out of range!!!");
        exit(1);
        }
	/**** Error handling for input****/
    	if (I < 1 && I > N-1) {
            fprintf(stderr,"B is out of range!!!");
        exit(1);
        }


//	printf("\n hola");
//	printf("\n hola");
//	printf("\n hola");

	/* begin AMPP */
        CHECK_NULL((aout = (char *) malloc(sizeof(char)*2*N)));
        CHECK_NULL((bout = (char *) malloc(sizeof(char)*2*N)));
        CHECK_NULL((a = (short *) malloc(sizeof(short)*(N+1))));
        CHECK_NULL((b = (short *) malloc(sizeof(short)*(N+1))));
        CHECK_NULL((hptr = (int *) malloc(sizeof(int)*(N+1)*(N+1))));
        CHECK_NULL((h = (int **) malloc(sizeof(int*)*(N+1))));
        /* Mount h[N][N] */
        for(i=0;i<=N;i++) 
           h[i]=hptr+i*(N+1);

        CHECK_NULL((xTracebackptr = (short *) malloc(sizeof(short)*(N+1)*(N+1))));
        CHECK_NULL((xTraceback = (short**) malloc(sizeof(short*)*(N+1))));
        /* Mount xTraceback[N][N] */
        for(i=0;i<=N;i++) 
           xTraceback[i]=xTracebackptr+i*(N+1);

        CHECK_NULL((yTracebackptr = (short *) malloc(sizeof(short)*(N+1)*(N+1))));
        CHECK_NULL((yTraceback = (short**) malloc(sizeof(short*)*(N+1))));
        /* Mount yTraceback[N][N] */
        for(i=0;i<=N;i++) 
           yTraceback[i]=yTracebackptr+i*(N+1);

        initChar2AATranslation();
        /* end  AMPP */


if(myid == MASTER){

	/** read PAM250 similarity matrix **/	
        /* begin AMPP */
        fscanf(pam,"%*s"); 
        /* end  AMPP */
	for (i=0;i<AA;i++)
		for (j=0;j<=i;j++) {
		if (fscanf(pam, "%d ", &temp) == EOF) {
			fprintf(stderr, "PAM file empty\n");
			fclose(pam);
			exit(1);
			}
		sim[i][j]=temp;
		}
	fclose(pam);
	for (i=0;i<AA;i++)
		for (j=i+1;j<AA;j++) 
			sim[i][j]=sim[j][i]; 	// symmetrify



/* begin AMPP */
	/** read first file in array "a" **/	
      	i=0;
        do {
           nc=fscanf(in1,"%c",&ch);
           if (nc>0 && char2AAmem[ch]>=0) 
           {
              a[++i] = char2AAmem[ch]; 
           }
        } while (nc>0 && (i<N));
	a[0]=i;  
        fclose(in1);

	/** read second file in array "b" **/	
	i=0;
        do {
           nc=fscanf(in2,"%c",&ch);
           if (nc>0 && char2AAmem[ch]>=0) 
           {
              b[++i] = char2AAmem[ch]; 
           }
        } while (nc>0 && (i<N));
	b[0]=i;  
        fclose(in2);
/* end AMPP */


}//end of MASTER
       
	/*h_out =  (int **) malloc((N+1)*sizeof(int *));

    	for (l = 0;l < N+1;l++){
    		h_out[l] = (int *) malloc(sizeof(int)*(N+1));
	    	}
*/

	

        MPI_Bcast(a, 1, MPI_SHORT, 0, MPI_COMM_WORLD );
        MPI_Bcast(&a[1], a[0], MPI_SHORT, 0, MPI_COMM_WORLD );

	MPI_Bcast(b, 1, MPI_SHORT, 0, MPI_COMM_WORLD );
        MPI_Bcast(&b[1], b[0], MPI_SHORT, 0, MPI_COMM_WORLD );
	
	MPI_Bcast(sim, AA*AA, MPI_INT, MASTER, MPI_COMM_WORLD);

	
	/* begin AMPP */
        /* You may want to delete yTraceback and xTraceback updates on the following      */
        /* process since we are not interested on it. See comments below. It is up to you.*/
        /* end AMPP */

	

	/*printf("b[0]: %d\n",b[0]);
        for (i=1;i<=b[0];i++)
		printf("\n %d",b[i]);
	
	printf("\n");
	FILE *dat2;

        dat2=fopen("2.dat","w");	
      
        for (i=0;i<AA;i++){
	  for (j=0;j<AA;j++){
				
		  fprintf(dat2,"%d\t",sim[j][i]);
					
	   }
		fprintf(dat2,"\n");
	}

	fclose(dat2);
*/


	MPI_Barrier(MPI_COMM_WORLD);

	/** initialize traceback array **/
	/*Max=xMax=yMax=0;
	for (i=0;i<=a[0];i++)
		for (j=0;j<=b[0];j++) {
			xTraceback[i][j]=-1;
			yTraceback[i][j]=-1;
			}
	*/


	if (myid == MASTER){
    		for (l = 0;l < b[0]+1; l++) {
				h[0][l] = 0;
			}
		}

	for(m = 0;m < a[0]+1;m++){
		h[m][0] = 0;
		}

	
/** compute "h" local similarity array **/
	//for (i=0;i<=a[0];i++) h[i][0]=0;
	//for (j=0;j<=b[0];j++) h[0][j]=0;
	
	
	stride = 0;

	num_block_rows = a[0]%I==0 ? a[0]/I : a[0]/I+1;
	
	num_steps = num_block_rows%numprocs==0 ? num_block_rows/numprocs : num_block_rows/numprocs + 1;

	num_block_cols = b[0]%B==0 ? b[0]/B : b[0]/B+1;
	

/*
if(MASTER==myid){
	printf("num_block_rows:%d\n", num_block_rows);	
	printf("num_steps:%d\n", num_steps);	
	printf("num_block_cols:%d\n", num_block_cols);	
	}
*/	
	t_start=usecs();	

	for(i=0;i<num_steps;i++){
        	if(!((i==num_steps-1) && (myid >= num_block_rows%numprocs) && (num_block_rows%numprocs!=0))){

			for(j=0;j<num_block_cols;j++){
			 	if(j==num_block_cols-1 && b[0]%B!=0)
					num_cols = b[0]%B;
				else
					num_cols = B;

			   if((numprocs*i+myid)!=0){
				MPI_Irecv(&h[(myid)*I+stride][B*j+1], num_cols, MPI_INT, (numprocs*i+myid-1)%numprocs, 0, MPI_COMM_WORLD, &requestr);
							 MPI_Wait(&requestr, &status);
				/*if(i==num_steps-1){
						 printf("\n ovo nece da valja k-1: %d",k-1);
						 printf("\n B*j+1: %d\n",B*j+1);
						 printf("\n num_cols:%d\n",num_cols);
						 MPI_Finalize();
						 return 0;
			 			}
				printf(" %d\t %d\t %d\t myid Recv:%d\t from process: %d\n", (myid)*I+stride,B*j, num_cols, myid, (numprocs*num_steps+myid-1)%numprocs); */
			/*	if(i==2 && myid==1 && j==num_block_cols-1){
					printf("\n num_rows:%d\n",num_rows);
					MPI_Finalize();
					return 0;
				}
			*/
			    }

				if((i==num_steps-1) && (a[0]%I!=0)){
				    if(numprocs*i + myid == num_block_rows-1)
						num_rows = a[0]%I;
					else
				  		num_rows = I;
				}
				else
					num_rows = I;

				int d;
				int n = B;
				int num_waves = ((n / CHUNK_SIZE) * 2);
				int div, mod, dn;
				int xbase, ybase;
				int ix, iy;
				int kk;

				int ii;
				int i_off = stride + myid * I;

				int jj;
				int j_off = j * B;

				for (d = 1; d < num_waves; d++) {
					div = d / ((num_waves) / 2 + 1);
					mod = d % ((num_waves) / 2 + 1);
					dn = d;
					if (div)
						dn = num_waves - d;

					#pragma omp parallel for shared(h) private(k, ii, jj, xbase, ybase, ix, iy, diag, down, right, max)
			      for (k = 0; k < dn; k++) {
						xbase = 1 + i_off + CHUNK_SIZE * (d - k - 1) - (div * ((d - dn) * (CHUNK_SIZE/2)));
						ybase = 1 + j_off + CHUNK_SIZE * (k) + (div * (num_waves / 2 - dn) * CHUNK_SIZE);
						for (ii = 0; ii < CHUNK_SIZE; ii++) {
							for (jj = 0; jj < CHUNK_SIZE; jj++) {
								ix = xbase + ii;
								iy = ybase + jj;
								diag  = h[ix-1][iy-1] + sim[a[ix]][b[iy]];
								down  = h[ix-1][iy] + DELTA;
								right = h[ix][iy-1] + DELTA;
								max   = MAX3(diag, down, right);
								if (max <= 0) {
									h[ix][iy] = 0;
								}
								else if (max == diag) {
									h[ix][iy] = diag;
								}
								else if (max == down) {
									h[ix][iy] = down;
								}
								else {
									h[ix][iy] = right;
								}
								/* printf("NEW: %d %d -- %d %d -- %d -- %d %d %d\n", ix, iy, i, j, h[ix][iy], diag, down, right); */
							}
						}
					}
				}
				/*
				*/

				if((numprocs*i + myid)!= (num_block_rows-1)){
					/* printf("FOO1: %d\n", ix); */
					/* printf("FOO2: %d\n", ix); */
					MPI_Isend(&h[i_off + CHUNK_SIZE][B*j+1], num_cols, MPI_INT, (numprocs*i+myid+1)%numprocs, 0, MPI_COMM_WORLD, &requests);

					/* MPI_Isend(&h[k-1][B*j+1], num_cols, MPI_INT, (numprocs*i+myid+1)%numprocs, 0, MPI_COMM_WORLD, &requests); */

					/*if(i==num_steps-1){
						 printf("\n ovo nece da valja k-1: %d",k-1);
						 printf("\n B*j+1: %d\n",B*j+1);
						 printf("\n num_cols:%d\n",num_cols);
						 MPI_Finalize();
						 return 0;
			 			}
						
			  		
					printf(" %d\t %d\t %d\t myid Send:%d\t to process: %d\n", k-1, B*j, num_cols, myid,(numprocs*num_steps+myid+1)%numprocs); 				if(i==2 && myid==0 && j==num_block_cols-1){			
						MPI_Finalize();
						return 0;
					}*/
				}
				
	    	
				
		      }// for j
		      stride = I*numprocs + stride;
		}//if
	
				
	}//for i	



	
	MPI_Barrier( MPI_COMM_WORLD ); 

	if(myid == MASTER){
                 t_end=usecs();
    		 time = ((double)(t_end-t_start))/1000000;
      		 printf("computation %f\n", time);
	}

	
	if(numprocs!=1){

	if(num_block_rows%numprocs==0 && a[0]%I==0){
	
	 
	for(i=0;i<num_steps;i++)	
		MPI_Gather( &h[I*i*numprocs+myid*I+1][0], I*(b[0]+1), MPI_INT, &h[I*i*numprocs + 1][0], I*(b[0]+1), MPI_INT, MASTER, MPI_COMM_WORLD); 	

	}	
	else{
	
		
		for(i=0;i<num_steps-1;i++)		
			MPI_Gather( &h[I*i*numprocs+myid*I+1][0], I*(b[0]+1), MPI_INT, &h[I*i*numprocs + 1][0], I*(b[0]+1), MPI_INT, MASTER, MPI_COMM_WORLD); 	
	
				
	    if(a[0]%I==0){ 
			if(myid>0 && myid<num_block_rows%numprocs)		    
				MPI_Send(&h[I*(num_steps-1)*numprocs + I*myid + 1][0], I*(b[0]+1), MPI_INT, MASTER, 0, MPI_COMM_WORLD);   
  	    }else{
		
		
		//printf("\n %d vece od 2 %d",num_block_rows%numprocs);	
			

				
		if(num_block_rows%numprocs > 2){

			if(myid>0 && (numprocs*num_steps+myid)%numprocs < numprocs-1)			
	            		MPI_Send(&h[I*(num_steps-1)*numprocs + I*myid + 1][0], I*(b[0]+1), MPI_INT, MASTER, 0, MPI_COMM_WORLD);
			else if(myid!=0 && numprocs*(num_steps-1)+myid == num_block_rows-1)
				MPI_Send(&h[I*(num_steps-1)*numprocs + I*myid + 1][0],  (a[0]%I)*(b[0]+1), MPI_INT, MASTER, 0, MPI_COMM_WORLD);   
		}
		else if(myid!=0 && (numprocs*(num_steps-1)+myid) == (num_block_rows-1)){
			//printf("\n %d vece od 2 %d",num_block_rows%numprocs);	
			MPI_Send(&h[I*(num_steps-1)*numprocs + I*myid + 1][0],  (a[0]%I)*(b[0]+1), MPI_INT, MASTER, 0, MPI_COMM_WORLD);   
	
			
			}	

	   	}// if(a[0]%I==0){ 
	}//if(num_block_rows%numprocs==0 && a[0]%I==0)
	


 	
	if(myid == MASTER){
		
		if(num_block_rows%numprocs > 1){

			if(num_block_rows%numprocs!=0 && a[0]%I==0){
				for(i=1;i<num_block_rows%numprocs;i++)
					MPI_Recv(&h[I*(num_steps-1)*numprocs + i*I + 1][0], I*(b[0]+1), MPI_INT, i, 0, MPI_COMM_WORLD, &status);
			}else if(num_block_rows%numprocs!=0 && a[0]%I!=0){
			
				
				if(num_block_rows%numprocs > 2){

					for(i=1;i<num_block_rows%numprocs-1;i++)
						MPI_Recv(&h[I*(num_steps-1)*numprocs + i*I + 1][0], I*(b[0]+1), MPI_INT, i, 0, MPI_COMM_WORLD, &status);
				
				MPI_Recv(&h[I*(num_steps-1)*numprocs + I*(num_block_rows%numprocs-1) + 1][0], (a[0]%I)*(b[0]+1), MPI_INT, num_block_rows%numprocs-1, 0, MPI_COMM_WORLD, &status);	
				}
				else
					MPI_Recv(&h[I*(num_steps-1)*numprocs + I*(num_block_rows%numprocs-1) + 1][0], (a[0]%I)*(b[0]+1), MPI_INT, num_block_rows%numprocs-1, 0, MPI_COMM_WORLD, &status);

			}else if(num_block_rows%numprocs==0 && a[0]%I!=0){
				if(numprocs > 2){
					for(i=1;i<numprocs-1;i++)
						MPI_Recv(&h[I*(num_steps-1)*numprocs + i*I + 1][0], I*(b[0]+1), MPI_INT, i, 0, MPI_COMM_WORLD, &status);
				}
			
				MPI_Recv(&h[I*(num_steps-1)*numprocs + I*(numprocs-1) + 1][0], (a[0]%I)*(b[0]+1), MPI_INT, numprocs-1, 0, MPI_COMM_WORLD, &status);
			}
		}
	}//if(myid == MASTER)
	}//if(numprocs!=1)

	
	   

	if(myid == MASTER){
                 t_end=usecs();
    		 time = ((double)(t_end-t_start))/1000000;
      		 printf("total %f\n", time);
	}

	
	
        /* begin AMPP */
        /* AMPP parallelization STOPS here. We are not interested on the parallelization     */
        /* of the traceback process, and the generation of the match result.*/
        /* end AMPP */

	MPI_Finalize();
	return 0;

	if(myid==MASTER){

	
	
	FILE *dat;
	
	dat=fopen("h.dat","wt");
	
	for (i=0;i<=a[0];i++){
		
                for (j=0;j<=b[0];j++){
			fprintf(dat,"%d\t ",h[i][j]);

			   	
		} 
		fprintf(dat,"\n");
	}
	
	}

	/* printf("hola"); */
	MPI_Finalize();
	return 0;
	
	


		
	
	// initialize output arrays to be empty -- this is unnecessary
	for (i=0;i<N;i++) aout[i]=' ';
	for (i=0;i<N;i++) bout[i]=' ';
	
	
	// reset to max point to do alignment
	i=xMax; j=yMax;
	x=y=0;
	topskip = bottomskip = 1;
	while (i>0 && j>0 && h[i][j] > 0){
		if (topskip && bottomskip) {
			aout[x++]=AA2charmem[a[i]];
			bout[y++]=AA2charmem[b[j]];
			}
		else if (topskip) {
			aout[x++]='-';
			bout[y++]=AA2charmem[b[j]];
			}
		else if (bottomskip) {
			aout[x++]=AA2charmem[a[i]];
			bout[y++]='-';
			}
		topskip    = (j>yTraceback[i][j]);
		bottomskip = (i>xTraceback[i][j]);
		tempi=i;tempj=j;
		i=xTraceback[tempi][tempj];
		j=yTraceback[tempi][tempj];
	}

	

	// print alignment
	printf("\n");
	printf("\n");
	for (i=x-1;i>=0;i--) printf("%c",aout[i]);	
	printf("\n");
	for (j=y-1;j>=0;j--) printf("%c",bout[j]);	
	printf("\n");
	printf("\n");

	//MPI_Finalize();

	//return 0;

	
	
}

void error(char * s) {
	fprintf(stderr,"%s\n",s);
	exit(1);
}


/* Begin AMPP */
void initChar2AATranslation(void)
{
    int i; 
    for(i=0; i<256; i++) char2AAmem[i]=-1;
    char2AAmem['c']=char2AAmem['C']=0;
    AA2charmem[0]='c';
    char2AAmem['g']=char2AAmem['G']=1;
    AA2charmem[1]='g';
    char2AAmem['p']=char2AAmem['P']=2;
    AA2charmem[2]='p';
    char2AAmem['s']=char2AAmem['S']=3;
    AA2charmem[3]='s';
    char2AAmem['a']=char2AAmem['A']=4;
    AA2charmem[4]='a';
    char2AAmem['t']=char2AAmem['T']=5;
    AA2charmem[5]='t';
    char2AAmem['d']=char2AAmem['D']=6;
    AA2charmem[6]='d';
    char2AAmem['e']=char2AAmem['E']=7;
    AA2charmem[7]='e';
    char2AAmem['n']=char2AAmem['N']=8;
    AA2charmem[8]='n';
    char2AAmem['q']=char2AAmem['Q']=9;
    AA2charmem[9]='q';
    char2AAmem['h']=char2AAmem['H']=10;
    AA2charmem[10]='h';
    char2AAmem['k']=char2AAmem['K']=11;
    AA2charmem[11]='k';
    char2AAmem['r']=char2AAmem['R']=12;
    AA2charmem[12]='r';
    char2AAmem['v']=char2AAmem['V']=13;
    AA2charmem[13]='v';
    char2AAmem['m']=char2AAmem['M']=14;
    AA2charmem[14]='m';
    char2AAmem['i']=char2AAmem['I']=15;
    AA2charmem[15]='i';
    char2AAmem['l']=char2AAmem['L']=16;
    AA2charmem[16]='l';
    char2AAmem['f']=char2AAmem['F']=17;
    AA2charmem[17]='L';
    char2AAmem['y']=char2AAmem['Y']=18;
    AA2charmem[18]='y';
    char2AAmem['w']=char2AAmem['W']=19;
    AA2charmem[19]='w';
}

/* end AMPP*/
