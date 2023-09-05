#include <malloc.h>
#include <stdio.h>
#define NR_END 1
#define FREE_ARG char* 

void nrerror(error_text)  /* Numerical Recipe standard error handler */
char error_text[];
{
    void exit();
    
    fprintf(stderr,"Numerical Recipes run-time error...\n") ;
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}

int *ivector(nl, nh)
long nl, nh ;
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

float *vector(nl,nh)  /* allocates a 'float' vector with range [nl..nh]
                                                                      */
int  nl, nh;
{
    float *v;
    
    v=(float *) malloc((unsigned) (nh-nl+1)*sizeof(float)) ;
    if (!v) nrerror("allocation failure in vector()") ;
    return v-nl;
}

double *dvector(nl,nh)  /* allocates a 'float' vector with range [nl..nh]
                                                                      */
int  nl, nh;
{
    double *v;
    
    v=(double *) malloc((unsigned) (nh-nl+1)*sizeof(double)) ;
    if (!v) nrerror("allocation failure in dvector()") ;
    return v-nl;
}


float **matrix(nrl,nrh,ncl,nch)  
/* Allocates a 'float' matrix with range [nrl..nrh][ncl..nch]. */
int  nrl,nrh,ncl,nch ;
{
     int i;
     float **m;
   
     /* allocates pointers to rows */
     m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*)) ;
     if (!m) nrerror("allocation failure 1 in matrix()");
     m -= nrl ;

     /* allocates rows and set pointers to them. */
     for(i=nrl;i<=nrh;i++) {
        m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
        if (!m[i]) nrerror("allocation failure 2 in matrix()");
        m[i] -= ncl;
     }
     /* return pointer to array of pointers to rows. */
     return m ;
}

double **dmatrix(nrl,nrh,ncl,nch)  
/* Allocates a 'double' matrix with range [nrl..nrh][ncl..nch]. */
int  nrl,nrh,ncl,nch ;
{
     int i;
     double **m;
   
     /* allocates pointers to rows */
     m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*)) ;
     if (!m) nrerror("allocation failure 1 in matrix()");
     m -= nrl ;

     /* allocates rows and set pointers to them. */
     for(i=nrl;i<=nrh;i++) {
        m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
        if (!m[i]) nrerror("allocation failure 2 in matrix()");
        m[i] -= ncl;
     }
     /* return pointer to array of pointers to rows. */
     return m ;
}

void free_ivector(v, nl, nh)
int *v ;
long nl, nh ;
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_vector(v,nl,nh)
double  *v;
int nl,nh;
/* Frees a 'float' vector allocated by vector(). */
{
   free((void *) (v+nl));
}

void free_dvector(v,nl,nh)
double  *v;
int nl,nh;
/* Frees a 'float' vector allocated by vector(). */
{
   free((void*) (v+nl)); 
}

void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
int  nrl,nrh,ncl,nch;
/* Frees a matrix allocated with matrix */
{
   int i ;
   
   for(i=nrh;i>=nrl;i--) free((void*) (m[i]+ncl)) ;
   free((void*) (m+nrl)) ;
} 

