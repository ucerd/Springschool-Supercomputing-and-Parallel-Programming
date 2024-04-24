/*-------------------------------------------------------------

Copyright (C) 2000 Peter Clote.
All Rights Reserved.

Permission to use, copy, modify, and distribute this
software and its documentation for NON-COMMERCIAL purposes
and without fee is hereby granted provided that this
copyright notice appears in all copies.

-------------------------------------------------------------*/

/*

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

*/

/*
AMPP: SOME BUGS HAS BEEN SOLVED
      SOME CODE HAS BEEN IMPROVED
      Code added or re-coded starts with begin AMPP and ends with
      end AMPP
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdlib.h>

#include <omp.h>
#include <sys/times.h>

/* begin AMPP */
/* Just a note:                       */
/* N must be the size of the arrays   */
/* Here we assume that the two arrays */
/* have the same size                 */

#define MAX_SEQ 50

#define CHECK_NULL(_check) {\
   if ((_check)==NULL) \
      fprintf(stderr, "Null Pointer allocating memory\n");\
   }

/* end AMPP */

#define AA 20 // number of amino acids
#define MAX2(x,y)     ((x)<(y) ? (y) : (x))
#define MAX3(x,y,z)   (MAX2(x,y)<(z) ? (z) : MAX2(x,y))

void error(char *);

/* begin AMPP*/
int char2AAmem[256];
int AA2charmem[AA];
void initChar2AATranslation(void);
/* end AMPP */

int main(int argc, char *argv[]) {

	// variable declarations
	FILE * in1, *in2, *pam;
	char ch;
	int temp;
	int i,j,tempi,tempj,x,y,diag,down,right,DELTA;
	int topskip,bottomskip;
	char *aout,*bout;
	int Aend,Bend,Abegin,Bbegin;
	int max, Max, xMax, yMax;
	// Max is first found maximum in similarity matrix H
	// max is auxilliary to determine largest of
	// diag,down,right, xMax,yMax are h-coord of Max
	short *a, *b;
	int *hptr;
	int **h;
	int sim[AA][AA];		// PAM similarity matrix
	short **xTraceback, **yTraceback;
	short *xTracebackptr, *yTracebackptr;
        int N;
        int nc;


/* begin AMPP */
	/**** Error handling for input file ****/
	if (( argc != 5) && (argc!=6)) {
	     	fprintf(stderr,"%s protein1 protein2 PAM gapPenalty [N]\n",argv[0]);
		exit(1);
	}
        else if (argc==5)  /* Maximum size of the proteins, they should have equal size */
        {
	   /***** Initialization of input file and pair array **/
	    in1   = fopen(argv[1],"r");
    	    in2   = fopen(argv[2],"r");
	    pam   = fopen(argv[3],"r");
	    DELTA = atoi(argv[4]);
            N = MAX_SEQ;
        } else
        {
	    in1   = fopen(argv[1],"r");
    	    in2   = fopen(argv[2],"r");
	    pam   = fopen(argv[3],"r");
	    DELTA = atoi(argv[4]);
            N     = atoi(argv[5]);
        }
/* end AMPP */

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
			sim[i][j]=sim[j][i]; // symmetrify


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

	/* begin AMPP
	 * You may want to delete yTraceback and xTraceback updates on the
	 * following process since we are not interested on it. See comments below.
	 * It is up to you.  * end AMPP */

	/** initialize traceback array **/
	Max=xMax=yMax=0;
	for (i=0;i<=a[0];i++)
		for (j=0;j<=b[0];j++) {
			xTraceback[i][j]=-1;
			yTraceback[i][j]=-1;
			}


	/** compute "h" local similarity array **/
	for (i=0;i<=a[0];i++) h[i][0]=0;
	for (j=0;j<=b[0];j++) h[0][j]=0;

	/* int nthreads = omp_get_num_threads(); */
	/* printf("%d\n", nthreads); */

	#define CHUNK_SIZE 250

	int d;
	int n = a[0];
	int num_waves = ((n / CHUNK_SIZE) * 2);
	int div, mod, dn;
	int k;
	int xbase, ybase;
	int ix, iy;

	/* omp_lock_t *lock = (omp_lock_t *) malloc(sizeof(omp_lock_t) * (num_waves + 1)); */
	/* for (i = 0; i < num_waves; i++) */
	/*    omp_init_lock(lock + i); */

	struct timeval start;
	struct timeval end;
	struct timeval elapsed;
	gettimeofday(&start, NULL);

	/* omp_set_lock(lock + 0); */

	for (d = 1; d < num_waves; d++) {
		/* omp_set_lock(lock + d); */
		/* omp_unset_lock(lock + d - 1); */

		div = d / ((num_waves) / 2 + 1);
		mod = d % ((num_waves) / 2 + 1);
		dn = d;
		if (div)
			dn = num_waves - d;

		/* printf("--\n"); */

		#pragma omp parallel for shared(h) private(k, i, j, xbase, ybase, ix, iy, diag, down, right, max)
		for (k = 0; k < dn; k++) {

			/* int tid = omp_get_thread_num(); */
			/* printf("TID: %d\n", tid); */

			xbase = 1 + CHUNK_SIZE * (d - k - 1) - (div * ((d - dn) * (CHUNK_SIZE/2)));
			ybase = 1 + CHUNK_SIZE * (k) + (div * (num_waves / 2 - dn) * CHUNK_SIZE);
			/* if (div) { */
			/*    xbase -= (d - dn); */
			/*    ybase += (num_waves / 2 - dn) * CHUNK_SIZE; */
			/* } */

			/* printf("%d %d - %d, %d, %d\n", xbase, ybase, (d - k - 1), d, dn); */

			for (i = 0; i < CHUNK_SIZE; i++) {
				for (j = 0; j < CHUNK_SIZE; j++) {
					ix = xbase + i;
					iy = ybase + j;

					/* printf("%d %d\n", ix, iy); */

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
				}
			}
		}
	}

	gettimeofday(&end, NULL);
	if (end.tv_usec < start.tv_usec) {
		int nsec = (start.tv_usec - end.tv_usec) / 1000000 + 1;
		start.tv_usec -= 1000000 * nsec;
		start.tv_sec += nsec;
	}
	if (end.tv_usec - start.tv_usec > 1000000) {
		int nsec = (end.tv_usec - start.tv_usec) / 1000000;
		start.tv_usec += 1000000 * nsec;
		start.tv_sec -= nsec;
	}
	elapsed.tv_sec = end.tv_sec - start.tv_sec;
	elapsed.tv_usec = end.tv_usec - start.tv_usec;
	printf("%ld.%06ld\n", elapsed.tv_sec, elapsed.tv_usec);

	/* begin AMPP
	 * AMPP parallelization STOPS here. We are not interested on the
	 * parallelization of the traceback process, and the generation of the
	 * match result.
	 * end AMPP */

	/*
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
	*/

#ifdef VERBOSE
	for (i = 0; i <= a[0]; i++) {
		char first = 1;
		for (j = 0; j <= b[0]; j++) {
			char* padding = " ";
			if (first)
				padding = "";
			printf("%s%d", padding, h[i][j]);
			first = 0;
		}
		printf("\n");
	}
#endif

	return 0;
}

void error(char * s) {
	fprintf(stderr,"%s\n",s);
	exit(1);
}

/* begin AMPP */

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

/* end AMPP */
