#include "mex.h"

/****************************************************/
#include <math.h>
#include <matrix.h>

#define  DE2RA    0.174532925E-1
#define  E6A      1E-6
#define  PI       3.14159265
#define  PIO2     1.57079633
#define  QO       120.0
#define  SO       78.0
#define  TOTHRD   2./3.
#define  TWOPI    6.2831853
#define  X3PIO2   4.71238898
#define  XJ2      1.082616E-3
#define  XJ3      -.253881E-5
#define  XJ4      -1.65597E-6
#define  XKE      .743669161E-1
#define  XKMPER   6378.135
#define  XMNPDA   1440.0
#define  AE       1.0

/* input from matlab: initFlag,tsince,TLEdata */
/* output to matlab: x,y,z,xdot,ydot,zdot */
/* Input Arguments */
#define	JD_IN	prhs[0]
#define	TLE_IN	prhs[1]
#define NR_IN 2

/* Output Arguments */
#define	POS_OUT	plhs[0]
#define	VEL_OUT	plhs[1]
#define NR_OUT 2


/***** fmod2p *****/
static double fmod2p(double arg) {

  double x;

  x = fmod(arg,TWOPI);
  if(x<0) x=x+TWOPI;  
  return x;
}


/***** actan *****/
static double actan(double s, double c) {

  double res,TEMP;

      res=0;
      if (c == 0) goto lab5;
      if (c > 0) goto lab1;
      res=PI;
      goto lab7;
lab1: if (s == 0) goto lab8;
      if (s > 0) goto lab7;
      res=TWOPI;
      goto lab7;
lab5: if (s == 0) goto lab8;
       if (s > 0) goto lab6;
      res=X3PIO2;
      goto lab8;
lab6: res=PIO2;
      goto lab8;
lab7: TEMP=s/c;
      res=res+atan(TEMP);
lab8: return res; 

}



/***** SGP4 *****/
static void sgp4func(int IFLAG, double TSINCE, double *tle_data, double *res) {

/* variables */
static  double XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,XNDD6O,BSTAR;
static  double X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50;
static  double CK2,CK4,QOMS2T,S;
static  double IEXP, IBEXP,TEMP;

static  double COSIO,THETA2, X3THM1,EOSQ,BETAO2,BETAO,DEL1,A1,AO;
static  double DELO,XNODP,AODP,ISIMP,S4,QOMS24,PERIGE,PINVSQ;
static  double TSI,ETA,ETASQ,EETA,PSISQ,COEF,COEF1,C2,C1,SINIO;
static  double A3OVK2,C3,X1MTH2,C4,C5,THETA4,TEMP1,TEMP2,TEMP3;
static  double XMDOT,X1M5TH,OMGDOT,XHDOT1,XNODOT,OMGCOF,XMCOF,XNODCF;
static  double T2COF,XLCOF,AYCOF,DELMO,SINMO,X7THM1,D2,C1SQ;
static  double D3,D4,T3COF,T4COF,T5COF;
static  double OMGADF,XNODDF,OMEGA,XMP,XMDF,TSQ,XNODE;
static  double TEMPA,TEMPE,TEMPL,DELM,DELOMG,TCUBE,TFOUR,A,E,XL;
static  double BETA,XN,AXN,XLL,AYNL,XLT,AYN,CAPU,i,SINEPW,COSEPW;
static  double TEMP4,TEMP5,TEMP6,EPW,ECOSE,ESINE,ELSQ,PL;
static  double R,RDOT,RFDOT,BETAL,COSU,SINU,U,SIN2U,COS2U,RK,UK;
static  double XNODEK,XINCK,RDOTK,RFDOTK,SINUK,COSUK,SINIK,COSIK;
static  double SINNOK,COSNOK,XMX,XMY,UX,UY,UZ,VX,VY,VZ;


/* extract TLE data */
/* Card 1 */
  EPOCH = *(tle_data+0);
  XNDT2O = *(tle_data+1);
  XNDD6O = *(tle_data+2);
  IEXP = *(tle_data+3); 
  BSTAR = *(tle_data+4);
  IBEXP = *(tle_data+5);
/* Card 2 */
  XINCL = *(tle_data+6);
  XNODEO = *(tle_data+7);
  EO = *(tle_data+8);
  OMEGAO = *(tle_data+9);
  XMO = *(tle_data+10);
  XNO = *(tle_data+11);

 
  CK2=0.5*XJ2*pow(AE,2.0);
  CK4=-0.375*XJ4*pow(AE,4.0);
  QOMS2T=pow(((QO-SO)*AE/XKMPER),4.0);
  S=AE*(1.+SO/XKMPER);


  /* if(XNO <= 0) break; */

  XNDD6O=XNDD6O*pow(10.0,IEXP);


  XNODEO=XNODEO*DE2RA;
  OMEGAO=OMEGAO*DE2RA;
  XMO=XMO*DE2RA;
  XINCL=XINCL*DE2RA;
  TEMP=TWOPI/XMNPDA/XMNPDA;
  XNO=XNO*TEMP*XMNPDA;
  XNDT2O=XNDT2O*TEMP;
  XNDD6O=XNDD6O*TEMP/XMNPDA;

  BSTAR=BSTAR*pow(10.0,IBEXP)/AE;



/*************** from SGP4.FOR ***************/



  if (IFLAG == 0) goto lab100;
    
/*      RECOVER ORIGINAL MEAN MOTION (XNODP) AND SEMIMAJOR AXIS (AODP) */
/*      FROM INPUT ELEMENTS */


    A1=pow((XKE/XNO),TOTHRD);
    COSIO=cos(XINCL);
    THETA2=COSIO*COSIO;
    X3THM1=3.0*THETA2-1.0;
    EOSQ=EO*EO;
    BETAO2=1.0-EOSQ;
    BETAO=sqrt(BETAO2);
    DEL1=1.5*CK2*X3THM1/(A1*A1*BETAO*BETAO2);
    AO=A1*(1.0-DEL1*(0.5*TOTHRD+DEL1*(1.0+134.0/81.0*DEL1)));
    DELO=1.5*CK2*X3THM1/(AO*AO*BETAO*BETAO2);
    XNODP=XNO/(1.0+DELO);
    AODP=AO/(1.0-DELO);


/*      INITIALIZATION */

/*      FOR PERIGEE LESS THAN 220 KILOMETERS, THE ISIMP FLAG IS SET AND
*      THE EQUATIONS ARE TRUNCATED TO LINEAR VARIATION IN SQRT A AND
*      QUADRATIC VARIATION IN MEAN ANOMALY.  ALSO, THE C3 TERM, THE
*      DELTA OMEGA TERM, AND THE DELTA M TERM ARE DROPPED. */

    ISIMP=0;
    if((AODP*(1.0-EO)/AE) < (220.0/XKMPER+AE)) ISIMP=1;

/*      FOR PERIGEE BELOW 156 KM, THE VALUES OF
*      S AND QOMS2T ARE ALTERED */

    S4=S;
    QOMS24=QOMS2T;
    PERIGE=(AODP*(1.0-EO)-AE)*XKMPER;
    if(PERIGE >= 156.0) goto lab10;
    S4=PERIGE-78.0;
    if(PERIGE > 98.0) goto lab9;
    S4=20.0;
lab9: QOMS24=pow(((120.0-S4)*AE/XKMPER),4.0);
    S4=S4/XKMPER+AE;
lab10: PINVSQ=1/(AODP*AODP*BETAO2*BETAO2);
    TSI=1.0/(AODP-S4);
    ETA=AODP*EO*TSI;
    ETASQ=ETA*ETA;
    EETA=EO*ETA;
    PSISQ=fabs(1.0-ETASQ);
    COEF=QOMS24*pow(TSI,4.0);
    COEF1=COEF/pow(PSISQ,3.5);
    C2=COEF1*XNODP*(AODP*(1.0+1.5*ETASQ+EETA*(4.0+ETASQ))
	    +0.75*CK2*TSI/PSISQ*X3THM1*(8.0+3.0*ETASQ*(8.0+ETASQ)));
    C1=BSTAR*C2;
    SINIO=sin(XINCL);
      A3OVK2=-XJ3/CK2*pow(AE,3.0);
      C3=COEF*TSI*A3OVK2*XNODP*AE*SINIO/EO;
      X1MTH2=1.0-THETA2;
      C4=2.0*XNODP*COEF1*AODP*BETAO2*(ETA*
              (2.0+0.5*ETASQ)+EO*(0.5+2.0*ETASQ)-2.0*CK2*TSI/
              (AODP*PSISQ)*(-3.*X3THM1*(1.-2.0*EETA+ETASQ*
              (1.5-0.5*EETA))+0.75*X1MTH2*(2.*ETASQ-EETA*
              (1.+ETASQ))*cos(2.*OMEGAO)));
      C5=2.*COEF1*AODP*BETAO2*(1.+2.75*(ETASQ+EETA)+EETA*ETASQ);
      THETA4=THETA2*THETA2;
      TEMP1=3.*CK2*PINVSQ*XNODP;
      TEMP2=TEMP1*CK2*PINVSQ;
      TEMP3=1.25*CK4*PINVSQ*PINVSQ*XNODP;
      XMDOT=XNODP+0.5*TEMP1*BETAO*X3THM1+0.0625*TEMP2*BETAO*
              (13.-78.*THETA2+137*THETA4);
      X1M5TH=1.-5.*THETA2;
      OMGDOT=-0.5*TEMP1*X1M5TH+0.0625*TEMP2*(7.-114.*THETA2+
              395.*THETA4)+TEMP3*(3.-36.*THETA2+49.*THETA4);
      XHDOT1=-TEMP1*COSIO;
      XNODOT=XHDOT1+(0.5*TEMP2*(4.-19.*THETA2)+2.*TEMP3*(3.-
              7.*THETA2))*COSIO;
      OMGCOF=BSTAR*C3*cos(OMEGAO);
      XMCOF=-TOTHRD*COEF*BSTAR*AE/EETA;
      XNODCF=3.5*BETAO2*XHDOT1*C1;
      T2COF=1.5*C1;
      XLCOF=0.125*A3OVK2*SINIO*(3.+5.*COSIO)/(1.+COSIO);
      AYCOF=0.25*A3OVK2*SINIO;
      DELMO=pow((1.+ETA*cos(XMO)),3.0);
      SINMO=sin(XMO);
      X7THM1=7.*THETA2-1.;
      if(ISIMP == 1) goto lab90;
      C1SQ=C1*C1;
      D2=4.*AODP*TSI*C1SQ;
      TEMP=D2*TSI*C1/3.;
      D3=(17.*AODP+S4)*TEMP;
      D4=0.5*TEMP*AODP*TSI*(221.*AODP+31.*S4)*C1;
      T3COF=D2+2.*C1SQ;
      T4COF=0.25*(3.*D3+C1*(12.*D2+10.*C1SQ));
      T5COF=0.2*(3.*D4+12.*C1*D3+6.*D2*D2+15.*C1SQ*(2.*D2+C1SQ));
      
lab90: IFLAG == 0;

/*      UPDATE FOR SECULAR GRAVITY AND ATMOSPHERIC DRAG */

lab100: XMDF=XMO+XMDOT*TSINCE;
      OMGADF=OMEGAO+OMGDOT*TSINCE;
      XNODDF=XNODEO+XNODOT*TSINCE;
      OMEGA=OMGADF;
      XMP=XMDF;
      TSQ=TSINCE*TSINCE;
      XNODE=XNODDF+XNODCF*TSQ;
      TEMPA=1.-C1*TSINCE;
      TEMPE=BSTAR*C4*TSINCE;
      TEMPL=T2COF*TSQ;
      if(ISIMP != 1) {
	DELOMG=OMGCOF*TSINCE;
	DELM=XMCOF*(pow((1.+ETA*cos(XMDF)),3.0)-DELMO);
	TEMP=DELOMG+DELM;
	XMP=XMDF+TEMP;
	OMEGA=OMGADF-TEMP;
	TCUBE=TSQ*TSINCE;
	TFOUR=TSINCE*TCUBE;
	TEMPA=TEMPA-D2*TSQ-D3*TCUBE-D4*TFOUR;
	TEMPE=TEMPE+BSTAR*C5*(sin(XMP)-SINMO);
	TEMPL=TEMPL+T3COF*TCUBE+TFOUR*(T4COF+TSINCE*T5COF);
      }

      A=AODP*pow(TEMPA,2.0);
      E=EO-TEMPE;
      XL=XMP+OMEGA+XNODE+XNODP*TEMPL;
      BETA=sqrt(1.-E*E);
      XN=XKE/pow(A,1.5);

/*      LONG PERIOD PERIODICS */

      AXN=E*cos(OMEGA);
      TEMP=1./(A*BETA*BETA);
      XLL=TEMP*XLCOF*AXN;
      AYNL=TEMP*AYCOF;
      XLT=XL+XLL;
      AYN=E*sin(OMEGA)+AYNL;


/*      SOLVE KEPLERS EQUATION */

      CAPU=fmod2p(XLT-XNODE);
      TEMP2=CAPU;
      
      for(i=1;i<=10;i++) {
        SINEPW=sin(TEMP2);
	COSEPW=cos(TEMP2);
	TEMP3=AXN*SINEPW;
	TEMP4=AYN*COSEPW;
	TEMP5=AXN*COSEPW;
	TEMP6=AYN*SINEPW;
	EPW=(CAPU-TEMP4+TEMP3-TEMP2)/(1.-TEMP5-TEMP6)+TEMP2;
	if(fabs(EPW-TEMP2) <= E6A) goto lab140;
	TEMP2=EPW;
      }

/*      SHORT PERIOD PRELIMINARY QUANTITIES */

lab140:  ECOSE=TEMP5+TEMP6;
      ESINE=TEMP3-TEMP4;
      ELSQ=AXN*AXN+AYN*AYN;
      TEMP=1.-ELSQ;
      PL=A*TEMP;
      R=A*(1.-ECOSE);
      TEMP1=1./R;
      RDOT=XKE*sqrt(A)*ESINE*TEMP1;
      RFDOT=XKE*sqrt(PL)*TEMP1;
      TEMP2=A*TEMP1;
      BETAL=sqrt(TEMP);
      TEMP3=1./(1+BETAL);
      COSU=TEMP2*(COSEPW-AXN+AYN*ESINE*TEMP3);
      SINU=TEMP2*(SINEPW-AYN-AXN*ESINE*TEMP3);
      U=actan(SINU,COSU);
      SIN2U=2*SINU*COSU;
      COS2U=2*COSU*COSU-1.;
      TEMP=1./PL;
      TEMP1=CK2*TEMP;
      TEMP2=TEMP1*TEMP;

/*      UPDATE FOR SHORT PERIODICS */

      RK=R*(1.-1.5*TEMP2*BETAL*X3THM1)+0.5*TEMP1*X1MTH2*COS2U;
      UK=U-0.25*TEMP2*X7THM1*SIN2U;
      XNODEK=XNODE+1.5*TEMP2*COSIO*SIN2U;
      XINCK=XINCL+1.5*TEMP2*COSIO*SINIO*COS2U;
      RDOTK=RDOT-XN*TEMP1*X1MTH2*SIN2U;
      RFDOTK=RFDOT+XN*TEMP1*(X1MTH2*COS2U+1.5*X3THM1);

/*      ORIENTATION VECTORS */

      SINUK=sin(UK);
      COSUK=cos(UK);
      SINIK=sin(XINCK);
      COSIK=cos(XINCK);
      SINNOK=sin(XNODEK);
      COSNOK=cos(XNODEK);
      XMX=-SINNOK*COSIK;
      XMY=COSNOK*COSIK;
      UX=XMX*SINUK+COSNOK*COSUK;
      UY=XMY*SINUK+SINNOK*COSUK;
      UZ=SINIK*SINUK;
      VX=XMX*COSUK-COSNOK*SINUK;
      VY=XMY*COSUK-SINNOK*SINUK;
      VZ=SINIK*COSUK;

/*      POSITION AND VELOCITY */

      X=RK*UX;
      Y=RK*UY;
      Z=RK*UZ;
      XDOT=RDOTK*UX+RFDOTK*VX;
      YDOT=RDOTK*UY+RFDOTK*VY;
      ZDOT=RDOTK*UZ+RFDOTK*VZ;


/************************************/

      X=X*XKMPER/AE;
      Y=Y*XKMPER/AE;
      Z=Z*XKMPER/AE;
      XDOT=XDOT*XKMPER/AE*XMNPDA/86400;
      YDOT=YDOT*XKMPER/AE*XMNPDA/86400;
      ZDOT=ZDOT*XKMPER/AE*XMNPDA/86400;



/* return parameters */


      *(res+0) = X;
      *(res+1) = Y;
      *(res+2) = Z;
      *(res+3) = XDOT;
      *(res+4) = YDOT;
      *(res+5) = ZDOT;


}

/****************************************************/

double epoc2jd(double epoc)
{
    double temp;
    
    /* Is it the 20th century or the 21th century*/
    if((epoc/1000.0)<50){
        epoc=(epoc/1000.0)+2000.0;
    }
    else{
        epoc=(epoc/1000.0)+1900.0;
    }
    
    temp=floor(epoc);
    epoc=(epoc-temp)*1000.0;
    temp=367.0*temp-floor(7.0*(temp+floor((1.0+9.0)/12.0))/4.0)+
            floor(275.0*(1.0/9.0))+1.0+1721013.5+(0.0/24.0);
    epoc=epoc+temp-1.0;
    return epoc;
}

/* Matlab wrap
 * -----------------------------------------------------------------------*/

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *simjd, tsince, *TLEdata, *pos, *vel, res[6];
  int initFlag;
  
  /* Get the JD */
  simjd=mxGetPr(JD_IN);
  
  /* Get the TLE */
  TLEdata=mxGetPr(TLE_IN);
  
  /* create the output matrixes */
  POS_OUT = mxCreateDoubleMatrix(3,1,mxREAL);
  VEL_OUT = mxCreateDoubleMatrix(3,1,mxREAL);
  
  /* get a pointer to the real data in the output matrix */
  pos = mxGetPr(POS_OUT);
  vel = mxGetPr(VEL_OUT);

  /* Get the epoc from TLE */
  tsince = *(TLEdata+0);
  
  /* Convert epoc to JD */
  tsince = epoc2jd(tsince);
  
  /* Calculate the time since epoc and convert to minutes */
  tsince = (*simjd-tsince)*1440.0;
  
  /* Used in S-function, but must be hardcoded here */
  initFlag = 1;

  /* ts is -1 if tSinceEp is not ready. */
  /* epoch, first element in TLEdata, must be greater than zero */ 
  if((tsince == -1) || (*(TLEdata+0) == 0)) { 
      *(pos+0) = -1;
      *(pos+1) = -1;
      *(pos+2) = -1;
      *(vel+0) = -1;
      *(vel+1) = -1;
      *(vel+2) = -1;
  }
  else {
      /* Run the C-function */
      sgp4func(initFlag,tsince,TLEdata,res);
      *(pos+0) = *(res+0);
      *(pos+1) = *(res+1);
      *(pos+2) = *(res+2);
      *(vel+0) = *(res+3);
      *(vel+1) = *(res+4);
      *(vel+2) = *(res+5);
  }
}
