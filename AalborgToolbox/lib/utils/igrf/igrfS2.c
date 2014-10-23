#define S_FUNCTION_NAME igrfS2
#include "simstruc.h"

#include <stdio.h>
#define YR ssGetArg(S,0)
#define MN ssGetArg(S,1)

#include "igrf11c.h"

double *Nmaxx;

#define NCOEF 210 
static int  nyrs, Nmax;
double Gh[NCOEF]; /* IGRF coefficients array */

static double absv(double *vec)
{

    register double sum = 0.0;
    register double *d1;

    for (d1 = vec; d1 < (vec + 3); ++d1)
	sum += *d1 * *d1;
    return (sqrt(sum));
}



/* vector product c =  a x b */

static int cross(double *a, double *b, double *c)
{
    /* Function Body */
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return (0);
}



int crossn (double *a, double *b, double *c)
{
    /* Local variables */
    register double *d1;
    register double x;

    /* Function Body */

    cross (a, b, c);
    x = absv (c);

    if (x >= 1e-30)
        for (d1 = c; d1 < c + 3; )
            *d1++ /= x;

    return (0);
}



int com_read (FILE *fo, char *str)
{
    char    c, d;   
    register char   *c1;
    char    ok = 0, done = 1;
    int cnt = 0;

    c1 = str; 

    c = getc(fo);
    
    
    while (c != '\n' && c != EOF && done == 1) {
        if (c == '/') {
            d = getc(fo);
            if (d == '*')
                done = 0;
            else 
             {
                cnt += 2;   
                *c1++ = c;  
                *c1++ = d;  
                c = d;      
            }
        } else
         {
            ++cnt;
            *c1++ = c;
            c = getc(fo);
        }
    }

    *c1 = 0;

    if (done == 0)
        while (c != '\n' && c != EOF)
            c = getc(fo);

    /*** check for all blanks -- if all blanks return -cnt  ***/

    if (cnt > 0) {
        c1 = str;
        while (c1 < str + cnt && !ok) {
            if (*c1++ != ' ')
                ok = 1;
        }
        if (!ok)
            cnt = -cnt;
    }
    
    

    return (cnt);
}



char *skip_in_str (int n, char *str)
{
    register char   *c1;
    register int    i;

    c1 = str;

    while (*c1 == ' ')
        ++c1;

    if (*c1 == 0)
        return (0);

    for (i = 0; i < n; ++i) {

        while (*c1 != ' ' && *c1 != 0)
            ++c1;

        if (*c1 == 0)
            return (0);

        while (*c1 == ' ')
            ++c1;

        if (*c1 == 0)
            return (0);
    }

    return (c1);
}



void setigrf(double year){
	char    gh[4];
	char    str[192];
	char    c, d;   
    char   *c1;
    char    ok = 0, done = 1;
    int cnt = 0;
	FILE *fo;
	char file_name[60] = "./tmp_igrf_coeffs";
	int n = 10, i, m, j;
	static double   f, tmp1, tmp2;
	static double   vdata[24], years[24], f0, h0, w1, w2, gg[196], hh[196]; 
	double *d1, *d2, *d3;
	static int  in;
	
    if (( fo = fopen (file_name,"w")) == NULL) {
        printf("cannot write %s file\n", file_name);
    }

    fwrite(igrf11_1, 1, IGRF_COEFF1_SIZE, fo);
    fwrite(igrf11_2, 1, IGRF_COEFF2_SIZE, fo);
    fwrite(igrf11_3, 1, IGRF_COEFF3_SIZE, fo);
    fwrite(igrf11_4, 1, IGRF_COEFF4_SIZE, fo);


	fclose(fo);
    
	if (( fo = fopen (file_name,"r")) == NULL) {
        printf("cannot read %s file\n", file_name);
    }
    

	while (com_read(fo, str) <= 0)
        ;
    sscanf (str, "%d", &nyrs);

    printf("nyrs is: %i\n", nyrs);
	if (nyrs < 2 || nyrs > 24) {
		
        printf ("Wrong nyrs in igrf.d\n");
	return;
    }

	
    while (com_read(fo, str) <= 0)
        ;
    sscanf (str, "%d", &Nmax);

	printf("Nmax is: %i\n", Nmax);

    if (Nmax < 3 || Nmax > 13) { 
        printf ("Wrong Nmax in igrf.d\n");
	return;
    }

    while (com_read(fo, str) <= 0)
        ;

    c1 = str;
    for (i = 0; i < nyrs; ++i) {
        sscanf (c1, "%lf", &years[i]);
        c1 = skip_in_str (1, c1); 

    }
	
	    in = 0;	
    for (i = 0; i < nyrs - 1; ++i) 
        if (year >= years[i])
            in = i;

    w2 = year - years[in];
    if (w2 < 0.)
        w2 = 0.;
    w1 = 1.;

    if (in + 1 < nyrs - 1) {
        w2 /= years[in+1] - years[in]; 
        w1 = 1.0 - w2;
    }
             
    
    for (d1 = Gh; d1 < Gh + NCOEF; )
        *d1++ = 0.;

    for (n = 0; n < Nmax + 1; ++n) {
        d1 = gg + 14 * n; 
        d2 = hh + 14 * n; 

        for (m = 0; m <= n; ++m) {

            fscanf(fo, "%s %d %d", gh, &i, &j);

            for (d3 = vdata; d3 < vdata + nyrs; ) {
                fscanf(fo, "%lf", d3++);

			}
			
            if (n != i || m != j) {
                printf ("\nWrong gg in igrf.d\n");
                return;
            }
            *d1++ =  w1 * vdata[in] + w2 * vdata[in + 1]; 
            
            fscanf (fo, "%s %d %d", gh, &i, &j);
            for (d3 = vdata; d3 < vdata + nyrs; ) {
                fscanf (fo, "%lf", d3++);
            }
            if (n != i || m != j) {
                printf ("\nWrong gg in igrf.d\n");
                return;
            }
            *d2++ = w1 * vdata[in] + w2 * vdata[in + 1]; /* h parameters */
            
        }
    }


    fclose (fo);



	
	d1 = Gh; 
    *d1++ = 0.0; 

    f0 = -1.;   

    for (n = 1; n <= Nmax; ++n) {
        d2 = gg + 14 * n;
        d3 = hh + 14 * n + 1;
        f0 = f0 * .5 * n; 
        f = f0 * sqrt(2.0) / 2.; 
        *d1++ = *d2++ * f0;
        ++i;
        for (m = 1; m <= n; ++m) {
            tmp1 = (double) (n + m);
            tmp2 = (double) (n - m + 1);
            f = f * tmp1 / tmp2 * sqrt(tmp2 / tmp1);
			*d1++ = *d2++ * f;
            *d1++ = *d3++ * f;
        }
    }

	
/*

    h0 = gg[11] * gg[11] + gg[12] * gg[12] + hh[12] * hh[12];
    dipmom = -sqrt(h0);
    w1 = fabs (gg[11] / dipmom);
    w2 = sqrt(1. - w1 * w1);
    tmp1 = atan(hh[12] / gg[12]);
    Eccdz[0] = w2 * cos(tmp1);
    Eccdz[1] = w2 * sin(tmp1);
    Eccdz[2] = w1;
    Eccdx[0] = 0.0;
    Eccdx[1] = 0.0;
    Eccdx[2] = 1.0;

    crossn(Eccdx, Eccdz, Eccdy);
    crossn(Eccdy, Eccdz, Eccdx);

    lx = -gg[12] * gg[22] + sqrt(3.) * (gg[11] * gg[23] + gg[12] * gg[24]
         + hh[12] * hh[24]);
    ly = -hh[12] * gg[22] + sqrt(3.) * (gg[11] * hh[23] - hh[12] * gg[24]
         + gg[12] * hh[24]);
    lz = gg[11] * 2.0 * gg[22] + sqrt(3.) * (gg[12] * gg[23] + hh[12] * hh[23]);
    tmp2 = (lz * gg[11] + lx * gg[12] + ly * hh[12]) * 0.25 / h0;

    Eccrr[0] = (lx - gg[12] * tmp2) / 3. / h0;
    Eccrr[1] = (ly - hh[12] * tmp2) / 3. / h0;
    Eccrr[2] = (lz - gg[11] * tmp2) / 3. / h0;

	*/


}




void igrfmodel(double *bv, double *geo, double *Nmax_loc)
{
    /* Local variables */
    register double *d1, *d2;
    int imax, nmax, N;
    double  f, h[NCOEF]; /* KV */
    int i, k, m;
    double  s, x, y, z;
    int ihmax, ih, il;
    double  xi[3], rq;
    int ihm, ilm;
    double  srq;

    srq=absv(geo); 
    rq = srq*srq; 
    
    N = (int) *Nmax_loc; 
    if (rq < .8) { 
      printf ("igrf call below surface !!!\n");
    }
    
	
    rq = 1. / rq;
    srq = sqrt(rq);
    if (rq < 0.25)
	nmax = (N - 3) * 4.0 * rq + 3;
    else
	nmax = N;
    
    

    for (d1 = xi, d2 = geo; d1 < xi + 3; )
	*d1++ = *d2++ * rq;
    
    ihmax = nmax * nmax; 
    imax = nmax + nmax - 2; 
    il = ihmax + nmax + nmax;

    d1 = h + ihmax;	
    d2 = Gh + ihmax; 
    for ( ; d1 <= h + il; )
	*d1++ = *d2++;

    for (k = 0; k < 3; k += 2) {
	i = imax;
	ih = ihmax;
	while (i >= k) {
	    il = ih - i - 1;
	    f = 2. / (double) (i - k + 2);
	    x = xi[0] * f;
	    y = xi[1] * f;
	    z = xi[2] * (f + f);
	    i += -2;
	    if (i >= 2) {
		for (m = 3; m <= i + 1; m += 2) {
		    ihm = ih + m;
		    ilm = il + m;
		    h[ilm+1] = Gh[ilm+1]+z*h[ihm+1]+x*(h[ihm+3]-h[ihm-1])-y*(h[ihm+2]+h[ihm-2]);
		    h[ilm] = Gh[ilm]+z*h[ihm]+x*(h[ihm+2]-h[ihm-2])+y*(h[ihm+3]+h[ihm-1]);
		}
		h[il+2] = Gh[il+2]+z*h[ih+2]+x*h[ih+4]-y*(h[ih+3]+h[ih]);
		h[il+1] = Gh[il+1]+z*h[ih+1]+y*h[ih+4]+x*(h[ih+3]-h[ih]);
	    } else if (i == 0) {
		h[il + 2] = Gh[il+2]+z*h[ih+2]+x*h[ih+4]-y*(h[ih+3]+h[ih]);
		h[il+1] = Gh[il+1]+z*h[ih+1]+y*h[ih+4]+x*(h[ih+3]-h[ih]);
	    }
	    h[il] = Gh[il]+z*h[ih]+(x*h[ih+1]+y*h[ih+2])*2.;
	    ih = il;
	}
    }

    s = h[0]*.5+(h[1]*xi[2]+h[2]*xi[0]+h[3]*xi[1])*2.;
    f = (rq+rq)*srq;
    x = f*(h[2]-s*(*(geo+0)))*1e-9;
    y = f*(h[3]-s*(*(geo+1)))*1e-9;
    z = f*(h[1]-s*(*(geo+2)))*1e-9;

    *(bv+0) = x;
    *(bv+1) = y;
    *(bv+2) = z;
    /* *(bv+3) = sqrt(x*x+y*y+z*z); */ 	

}

/****************************************************/


/* ------------------------------------------------------------
   FUNCTION: compute igrf field for cartesian geo
   CALL:     bv = irgfmodel(geo,Nmax)

   input:    geo(3) position vector (ECEF) in earth radii
             (re = 6371.2 km)
             Nmax   order of model 3 < Nmax <= 10

   output:   bv(1:3) magnetic field vector in ECEF coordinate system
             bv(4)   norm of bv(1:3)
             (all in units: nT)

   c-source code from bigrf.c Gh coefficients obtained from setigrf 
   routine in bigrf.c 1996 coefficients based on 1990 igrf model 
   [metoefc@uts.uni-c.dk (Eigil Friis-Christensen)]

   Copyright (c) 1995 Orsted Satellite Project/tb
		 Institute for Electronic Systems
	         Aalborg University, Denmark
--------------------------------------------------------------- */


/* Prototype */

/****************************************************/

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 0);   /* number of continuous states           */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states             */
    ssSetNumInputs(        S, 3);   /* number of inputs                      */
    ssSetNumOutputs(       S, 3);   /* number of outputs                     */
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag               */
    ssSetNumSampleTimes(   S, 0);   /* number of sample times                */
    ssSetNumSFcnParams(    S, 2);   /* number of input arguments             */
    ssSetNumRWork(         S, 0);   /* number of real work vector elements   */
    ssSetNumIWork(         S, 0);   /* number of integer work vector elements*/
    ssSetNumPWork(         S, 0);   /* number of pointer work vector elements*/
}


static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}


static void mdlInitializeConditions(double *x0, SimStruct *S)
{
  double *year;
  year = mxGetPr(YR); 
  Nmaxx = mxGetPr(MN);
  setigrf(*year);
}


static void mdlOutputs(double *y, double *x, double *u, SimStruct *S, int tid)
{
  double coor[3], *order, *res;
  res = y;

  coor[0] = (*(u+0))/(6371.2*1000);
  coor[1] = (*(u+1))/(6371.2*1000);
  coor[2] = (*(u+2))/(6371.2*1000);
  order = Nmaxx;
  igrfmodel(res, coor, order);

}


static void mdlUpdate(double *x, double *u, SimStruct *S, int tid)
{
}


static void mdlDerivatives(double *dx, double *x, double *u, SimStruct *S, int tid)
{
}


#if 0  /* Change to a 1 if VARIABLE_SAMPLE_TIME has been specified */
#define MDL_GET_TIME_OF_NEXT_VAR_HIT
static void mdlGetTimeOfNextVarHit(SimStruct *S)
{
    ssSetTNext(S, <timeOfNextHit>);
}
#endif


static void mdlTerminate(SimStruct *S)
{
}

#ifdef	MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif



