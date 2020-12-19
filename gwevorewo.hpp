#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include<functional>
#include<thread>
#include<mutex>
#include<fstream>
#include<algorithm>
#include<iomanip>
#include<sstream>
#include<iterator>
#include<chrono>
#include<omp.h>
#include<gsl/gsl_sf_erf.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_sf_expint.h>
#include<gsl/gsl_sf_gamma.h>
#include<limits>
#include<gsl/gsl_sf_lambert.h>

/*//unit conversion
const double sikm=2.99792458E5;
const double pcikm=3.085677581E13; //pc is taken as unity scale
const double sipc=sikm/pcikm;
const double mega=1.0E6;
const double KtoeV=8.621738E-5;
//some factors
const double pi=M_PI;
const double tpi=2.0*pi;
const double pii=1.0/pi;
const double tpii=1.0/tpii;
//cosmological parameters
const double TCMB=2.725*KtoeV;
const double TQCD=1.98E8;
const double TEW=2.E11;
const double Tndec=1.0E6;
const double aQCD=TCMB/TQCD;
const double aTEW=TCMB/TEW;
const double andec=TCMB/Tndec;
const double Omr=8.24E-5;
const double Omm=0.315;
const double Oml=0.685;
const double h=0.7;
const double cli=1.0; //speed of light set as unity
const double H0=h/sikm/mega;
const double GN=1.0;
const double kgiGN=6.674E-11/(mega*sikm*sikm)*1.0E3*pcikm/GN;
const double rhoc=3.*H0*H0/(8.0*pi*GN);
const double Hsor=H0*sqrt(Omr);
*/
const double KtoeV=8.621738E-5;
//some factors
const double pi=M_PI;
const double tpi=2.0*pi;
const double pii=1.0/pi;
const double tpii=1.0/tpii;
//cosmological parameters
const double TCMB=2.725*KtoeV;
const double TQCD=1.98E8;
const double TEW=1.59E11;
const double Tndec=2.7E6;
const double aQCD=TCMB/TQCD;
const double aTEW=TCMB/TEW;
const double andec=TCMB/Tndec;
const double apdec=1.0E200;
const double Omr=8.24E-5;
const double Omm=0.315;
const double Oml=0.685;
const double h=0.677;
const double cli=1.0; //speed of light set as unity
//const double H0=h/sikm/mega;
const double H0t=h*1./(3.089E19); //in u of 1/s
//const double GN=1.0;
//const double kgiGN=6.674E-11/(mega*sikm*sikm)*1.0E3*pcikm/GN;
//const double rhoc=3.*H0*H0/(8.0*pi*GN);
//const double Hsor=H0*sqrt(Omr);

const double Tini=TEW;//3.E11; //in eV
const double ainin=TCMB/Tini;
const double aini=ainin/ainin;

//const double Hin=(H0t*(Omm*pow(aini,-3)+Omr*pow(aini,-4)+Oml));// initial Hubble rate
//const double sidh=;

const double bohi=0.072;//1./12.5;
const double LI=2.*bohi; //fraction of Hubble scale as integral scale

// Settings - boolean

// Lagrangian rather than Eulerian decorrelation rate
const bool lagrangian=false;
// Heaviside step function cutoff to reduce computation times i.e. ignore components in the integration that barely contribute, yet are computational expensive
const bool thetcut=true;
// Cutoff scaler
const double dscalefac=2.;
// Approximative integration in certain boundaries, where an acceptable approximation is given
const bool appint=false;
// Variable Transformation in scaling
const bool simplyi=false;
// Gaussian (false) or Exponantial cutoff (true)
const bool altrate=false;
// True for compressible turbulence
const bool compr=false;
// Pure sweeping effect according to Kraichnan 1965, alternatively (false) approximation due to Kaneda 1993
const bool puresweep=false;
// Magnetosonic wave impact on the decorrelation rate
const bool magson=false;
// additional contribution for testing purposes regarding an changed decorrelation function that also takes sound wave oscillations into account -> keep false
const bool addcon=false;
// build up of the spectrum
const bool forcing=false;
// decorrelation between linear build up and turbulent decay
const bool forcing2=true;

//compressibility spec.
const double initsf=0.;
const double sfsf=1.;
const int sfmodel=0; // 0 0, 1 direct decay, 2 const, 3 const then decay

//const bool altspec=true;

const double Hsor=LI/(2.0*pi);  //initial Hubble rate fixed unit of length as 1;
const double siH=(2.0*pi/LI)*H0t*sqrt(Omm*pow(ainin,-3)+Omr*pow(ainin,-4)+Oml); // s in units of LI*H0=1.;
const double evisi=1./(6.582119E-16);  //hb=1 ev in s^{-1}
const double eViH0=evisi/siH;
//parameters of sim
const std::string direc="datagwrev1";
//const std::string fname="data43eul_0.1_0.13soundwm0okm2s0.9999irightmfm3s.dat";
const std::string fname="data43eul_0.072_0.025helforct9tesg2i2mc.dat";
//const std::string fname="sound0.999G4Fsl.dat";
//const std::string fname="nosoundl3.dat";
//const std::string fname="sound0.9G4Fst8.dat";//"sound0.13-10ttt.dat";
const int dratech=0; //fixes choice of rate 0: gaussian, 1: caprini, 2: delta, 4: coherent
//const double aini=aQCD; //initial scalefactor
//const double amax=3.*aQCD;//andec; //for test purpose
const double co1=1.; //factor for damping rate
const double Omvir=0.013; //initial turb energy
const double Ombir=0.013; //initial mag energy
const bool hrate=true;
const double helbr=1.0; //inital hel vel ratio
const double helvr=1.0; //inital hel mag ratio
const double Lini=2.*pi;//1.0E-2*(aini/aQCD)*(aini/aQCD); //inital integral scale in units of parsec today
//const double kap=0.26667;//1.3824;//0.843625;//0.26667 for 3 // 0.83304 for 2
const double kap=1.3824;//0.843625;//0.26667 for 3 // 0.83304 for 2
const double kap2=0.83304;//0.843625;//0.26667 for 3 // 0.83304 for 2
//const double kap=2.2245;//0.843625;//0.26667 for 3 // 0.83304 for 2
const double enorm=3./(4.*pi*kap); //norm fac of spectr.
const double enorm2=3./(4.*pi*kap2); //norm fac of spectr.
//const double tauLini=Lini/(2.0*sqrt(3.0/2.0*(Omvir+Ombir))); //cascade timescale
//const double tauLini=Lini/(2.0*sqrt(9.0/(8.0*1.3824)*(Omvir+Ombir))); //cascade timescale
const double tauLini=Lini/(sqrt(2.*(Omvir+Ombir))); //cascade timescale
const double taubuild=tauLini;///0.185*0.08; //buildup timescale
//const double taubuild=bohi/Hsor;//tauLini;///0.185*0.08; //buildup timescale
const double abuild=taubuild*Hsor;
const double amaxn=4.*TCMB/Tini;//andec; //for test purpose
//const double amax=amaxn/ainin;//aini*1.85;//+abuild*3.0;//amaxn/ainin;
const double amax=aini*2.;//+abuild*3.0;//amaxn/ainin;
const double gamma1=0.4;//0.4;//0.31; //pow law for integralscale cascade, 0.72 for max hel
const double gamma2=1.2;//1.2;//1.38; //pow law for energy cascade, 0.56 for max hel
const double gamma1h=2./3.;//0.72;
const double gamma2h=2./3.;//0.667;//0.56;
const double gammab=1.0; //buildup pow law
const int Ngw=300;
const double lgw=1.0E-5*(tpi/Lini);
const double ugw=1.0E3*(tpi/Lini);
const double initds=1.0E-6*Lini;
//num parameters
const double arr=1.0E-20;

const double snorm=2.*pi*kap;//4.80018237236; //prefactor for damp timescale
const double norm=snorm*snorm;
//const double Lc=3./4.;//5./12.;//11./6.//3./4. for 3; //norm. of spectra to L=1 initially // 1/2 for 2
const double Lc=5./12.;//11./6.//3./4. for 3; //norm. of spectra to L=1 initially // 1/2 for 2
const double Lc2=1./2.;//11./6.//3./4. for 3; //norm. of spectra to L=1 initially // 1/2 for 2
//const double Lc=1./10.;//11./6.//3./4. for 3; //norm. of spectra to L=1 initially // 1/2 for 2
//const double Kolpot=7./2.;//17./6.;//7./2. for 3 //Kol. law scaling for spectra // 3 for 2
const double Kolpot=17./6.;//7./2. for 3 //Kol. law scaling for spectra // 3 for 2
const double Kolpot2=3.;//7./2. for 3 //Kol. law scaling for spectra // 3 for 2
const int Npsi=3000; //number of points for approximation of spectral integration
const double ls=1.0E-7*(tpi/Lini); //lowest scale for appr.
const double us=1.0E7*(tpi/Lini); // largest scale for appr.
const double infind=-2.0/3.0;//2.0;//4.0/3.0; //pow law ind for appr in inf for inte
const double infindlag=4.0/3.0;//2.0;//4.0/3.0; //pow law ind for appr in inf for inte
const double lowbcoeff=0.0256479; //coeff for pow law at initial for inte
const double initpow=0.;//7.0; //inital pow law for inte appr.     //////// now a factor k^2 is shifted
const double min=1.0E-10;
const double max=1.0E10;
const double fp=1.0E-15; //floating point prec
const double errafa=1.0E-15;

const double g0=3.36;


const double iabserr=1.0E-18;


const double numtol=1.0E-18;

const double tomin=1.0E-50;

const double sltomin=sqrt(-log(tomin));
//const double tomins=tomin*1.0E-5;











// sinpq transformation of order p and q evalueated at t where val is the value of the transformation e.g. x=val and
// dval is the derivative e.g. dx/dt=dval
void sinpqt(int p, int q, double t, double *val, double *dval){
    double sint=sin(pi*t);
    double cost=cos(pi*t);
    if(p==q){
        switch(p){
            case 0:
            *val=t;
            *dval=1;
            return;
            break;
            case 1:
            *val=0.5*(1.-cost);
            *dval=0.5*sint;
          //  std::cout<<"someone home "<<*val<<" "<<*dval<<"\n";
            return;
            break;
            default:
            {
                if(p%2==0){
                    double ipo=t;
                    double inorm=1.;
                    for(int i=1; i<=p/2;i++){
                        double ind=2*(double)i;
                        ipo=-pow(2.,-ind)/(pi*ind)*pow(sint,ind-1)*cost+(ind-1.)/(4.*ind)*ipo;
                        inorm=(ind-1.)/(4.*ind)*inorm;
                    }
                    *val=ipo/inorm;
                    *dval=pow(sin(0.5*pi*t),p)*pow(cos(0.5*pi*t),q)/inorm;
                    return;
                } else {
                    double ipo=(1.-cost)/(2.*pi);
                    double inorm=1./pi;
                    for(int i=1; i<=(p-1)/2;i++){
                        double ind=2*i+1;
                        ipo=-pow(2.,-ind)/(pi*ind)*pow(sint,ind-1)*cost+(ind-1.)/(4.*ind)*ipo;
                        inorm=(ind-1.)/(4.*ind)*inorm;
                    }
                    *val=ipo/inorm;
                    *dval=pow(sin(0.5*pi*t),p)*pow(cos(0.5*pi*t),q)/inorm;
                    return;
                }    
            }
        }
    } else {
       if(q%2==0 && p%2==0){
        double ipo=t;
        double inorm=1;
        if(p>q){
          //  double inorm=1.;
            for(int i=1; i<=q/2;i++){
                double ind=2*(double)i;
                ipo=-pow(2.,-ind)/(pi*ind)*pow(sint,ind-1)*cost+(ind-1.)/(4.*ind)*ipo;
                inorm=(ind-1.)/(4.*ind)*inorm;
            }
            double conq=pow(cos(pi*0.5*t),q+1);
            for(int i=q/2+1;i<=p/2;i++){
                double ind=2*(double)i;
                ipo=-2.0/(pi*(ind+(double)q))*pow(sin(pi*0.5*t),ind-1)*conq+(ind-1.)/(ind+(double)q)*ipo;
                inorm=(ind-1.)/(ind+(double)q)*inorm;
            }
            *val=ipo/inorm;
            *dval=pow(sin(0.5*pi*t),p)*pow(cos(0.5*pi*t),q)/inorm;
            return;
        } else {
           // double inorm=1.;
            for(int i=1; i<=p/2;i++){
                double ind=2*(double)i;
                ipo=-pow(2.,-ind)/(pi*ind)*pow(sint,ind-1)*cost+(ind-1.)/(4.*ind)*ipo;
                inorm=(ind-1.)/(4.*ind)*inorm;
            }
            double sinp=pow(sin(pi*0.5*t),p+1);
            for(int i=p/2+1;i<=q/2;i++){
                double ind=2*(double)i;
                ipo=2.0/(pi*(ind+(double)p))*pow(cos(pi*0.5*t),ind-1)*sinp+(ind-1.)/(ind+(double)p)*ipo;
                inorm=(ind-1.)/(ind+(double)p)*inorm;
            }
            *val=ipo/inorm;
            *dval=pow(sin(0.5*pi*t),p)*pow(cos(0.5*pi*t),q)/inorm;
        }
       } else if(q%2==1 && q%2==1){
        double ipo=(1.-cost)/(2.*pi);
        double inorm=1./pi;
        if(p>q){
         //   double inorm=1.;
            for(int i=1; i<=(q-1)/2;i++){
                double ind=2*(double)i+1.;
                ipo=-pow(2.,-ind)/(pi*ind)*pow(sint,ind-1)*cost+(ind-1.)/(4.*ind)*ipo;
                inorm=(ind-1.)/(4.*ind)*inorm;
            }
            double conq=pow(cos(pi*0.5*t),q+1);
            for(int i=(q-1)/2+1;i<=(p-1)/2;i++){
                double ind=2*(double)i+1.;
                ipo=-2.0/(pi*(ind+(double)q))*pow(sin(pi*0.5*t),ind-1)*conq+(ind-1.)/(ind+(double)q)*ipo;
                inorm=(ind-1.)/(ind+(double)q)*inorm;
            }
            *val=ipo/inorm;
            *dval=pow(sin(0.5*pi*t),p)*pow(cos(0.5*pi*t),q)/inorm;
        } else {
         //   double inorm=1.;
            for(int i=1; i<=(p-1)/2;i++){
                double ind=2*(double)i+1.;
                ipo=-pow(2.,-ind)/(pi*ind)*pow(sint,ind-1)*cost+(ind-1.)/(4.*ind)*ipo;
                inorm=(ind-1.)/(4.*ind)*inorm;
            }
            double sinp=pow(sin(pi*0.5*t),p+1);
            for(int i=(p-1)/2+1;i<=(q-1)/2;i++){
                double ind=2*(double)i+1.;
                ipo=2.0/(pi*(ind+(double)p))*pow(cos(pi*0.5*t),ind-1)*sinp+(ind-1.)/(ind+(double)p)*ipo;
                inorm=(ind-1.)/(ind+(double)p)*inorm;
            }
            *val=ipo/inorm;
            *dval=pow(sin(0.5*pi*t),p)*pow(cos(0.5*pi*t),q)/inorm;
        }
       } else if(q%2==1){
        double ipo=2./pi*sin(0.5*pi*t);
        double inorm=2./pi;
        double conq=pow(cos(pi*0.5*t),2);
        double sinp=pow(sin(pi*0.5*t),p+1);
        for(int i=1;i<=p/2;i++){
            double ind=2*(double)i;
            ipo=-2.0/(pi*(ind+1.))*pow(sin(pi*0.5*t),ind-1)*conq+(ind-1.)/(ind+1.)*ipo;
            inorm=(ind-1.)/(ind+1.0)*inorm;
        }
        for(int i=1;i<=(q-1)/2;i++){
            double ind=2*(double)i+1.;
            ipo=2.0/(pi*(ind+(double)p))*pow(cos(pi*0.5*t),ind-1)*sinp+(ind-1.)/(ind+(double)p)*ipo;
            inorm=(ind-1.)/(ind+(double)p)*inorm;
        }
        *val=ipo/inorm;
        *dval=pow(sin(0.5*pi*t),p)*pow(cos(0.5*pi*t),q)/inorm;
       } else {
        double ipo=2./pi*(1.-cos(0.5*pi*t));
        double inorm=2./pi;
        double conq=pow(cos(pi*0.5*t),2);
        double sinp=pow(sin(pi*0.5*t),p+1);
        for(int i=1;i<=p/2;i++){
            double ind=2*(double)i+1;
            ipo=-2.0/(pi*(ind+1.))*pow(sin(pi*0.5*t),ind-1)*conq+(ind-1.)/(ind+1.)*ipo;
            inorm=(ind-1.)/(ind+1.0)*inorm;
        }
        for(int i=1;i<=(q-1)/2;i++){
            double ind=2*(double)i;
            ipo=2.0/(pi*(ind+(double)p))*pow(cos(pi*0.5*t),ind-1)*sinp+(ind-1.)/(ind+(double)p)*ipo;
            inorm=(ind-1.)/(ind+(double)p)*inorm;
        }
        *val=ipo/inorm;
        *dval=pow(sin(0.5*pi*t),p)*pow(cos(0.5*pi*t),q)/inorm;
       }
    }
}


/*
// calculate Re(erf(z))
double erfz(double a, double b){
    const double reltol=1.0E-10;
    double ta=2.0*a;
    double tab=2*a*b;
    double aq=a*a;
    double ext=exp(-aq);
    double sinf=sin(tab);
    double cosf=cos(tab);
    double sum0=gsl_sf_erf(a)+ext/ta*pii*(1.-cosf);
    double prif=2.0*pii*ext;
    int n=1;
    n2=n*n;
    double f1=prif*(ta*(1.0-cosh((double)n*b)*cosf)+sinh((double)n*b)*sinf)*(exp(-0.25*(double)n2)/((double)n2+4.0*aq));
    while(fabs(f1/sum0)>reltol){
        sum0+=f1;
        n++;
        f1=prif*(ta*(1.0-cosh((double)n*b)*cosf)+sinh((double)n*b)*sinf)*(exp(-0.25*(double)n2)/((double)n2+4.0*aq));
    }
    sum0+=f1;
    return sum0;
}
*/
// solves integrate_0^b exp(-v*z^2)*{cos(k*z),sin(k*z)}
const double maxp=-2.3*200.;



double calCISE(double v, double k, double b, double a){
    if(b<a*numtol) return 0.;
    v*=a*a;
    b=b/a;
    k*=a;
    a=1.;
    double as=sin(k*b);
    double ac=cos(k*b);
    double ai=a-b;
    double Cne=cos(k*a)*(gsl_sf_Ci(k*a)-gsl_sf_Ci(k*ai))+sin(k*a)*(gsl_sf_Si(k*a)-gsl_sf_Si(k*ai));
  //  if(std::isinf(Cne)|| std::isnan(Cne)) std::cout<<"shit in Cne\n";
    double C0=as/k;
    double C1=(1.-ac)/(k*k)+ai*C0;
    double comp=1.;
    std::vector<double>Cmv;
    Cmv.push_back(C0);
    Cmv.push_back(C1);
    int s=2;
    const double reltol=1.0E-4;
    double Csum=Cne;
    int n=1;
    int nfa=1;
    double vmul=1.;
    double compara=log(1./(1.-b));
    bool ftr=false;
    double tole=1.+reltol;
    while(comp>reltol || fabs(Csum)>compara*tole || !ftr){
   //    Csum+=
        if(comp<reltol){
            ftr=true;
        }
       double Cadd=-Cne;//-gsl_sf_choose(tn,);// need to add m=1;
       int tn=2*n;
       double am=pow(a,tn)/(-a);
       for(int m=1;m<=tn;m++){
        Cadd+=-gsl_sf_choose(tn,tn-m)*am;
        if(m<=s){
            Cadd*=Cmv[m-1];
            if(std::isinf(Cadd)){
                std::cout<<Cmv[m-1]<<" "<<m-1<<" "<<Cmv[m-2]<<" "<<Cmv[m-3]<<" "<<k<<"\n";
                exit(1);
            }
        } else {
            double cmul=pow(ai,m-1)/k*as+((double)m-1.)/(k*k)*(pow(a,m-2)-pow(ai,m-2)*ac)-(double)((m-1)*(m-2))/(k*k)*Cmv[m-3];
            Cmv.push_back(cmul);
            if(std::isinf(cmul)) {
                std::cout<<cmul<<" "<<m<<"\n";
                exit(1);
                }
            Cadd*=cmul;
            s++;
        }

        am=am/(-a);
       }
       Cadd*=vmul/((double) nfa);
       comp=fabs(Cadd/Csum);
       if(comp>reltol) ftr=false;
//       if(std::isinf(Cadd)) std::cout<<"ohmy"<<am;
       Csum+=Cadd;
      // if(Csum>compara) ftr=false;
       n++;
       nfa*=n;
       vmul*=-v;
       //if(Cadd>Csum) std::cout<<vmul<<" "<<Csum

    }
    if(fabs(Csum)>compara || std::isinf(Csum) || std::isnan(Csum)) std::cout<<"compara "<<Csum<<" "<<log(1./(1.-b))
        <<" "<<n<<" "<<k<<" "<<a<<" "<<v<<" "<<b<<" "<<vmul<<" "<<Cne<<" "<<C0<<" "<<"\n";
 //   if(true)std::cout<<"shit\n";

    return Csum;
}

// calc. int e^(-v^2*x^2)*cos(k*x)/(a-x) using series of e^(-v^2*x^2)
// trouble with small k

const double signpow=-15;
const double cut=signpow*log(10.);
double cosint(int n, double k, double b){
    double kb=k*b;
    double kb2=kb*kb;
    double lb=log(kb);
    int nmax=std::ceil((cut-lb)/(2.*lb));
    double sum=0.;
    double bpow=pow(kb,n+1);
    double dfac=1.;
    double dsig=1.;
    double dpar=(double)n+1.;
    for(int i=0;i<=nmax;i++){
       sum+=dsig/(dfac*dpar)*bpow;
       dpar+=2.;
       dfac*=((double)i+1.)*((double)i+2.);
       dsig*=-1.;
       bpow*=kb2;
    }
    return sum;
}
double sinint(int n, double k, double b){
    double kb=k*b;
    double kb2=kb*kb;
    double lb=log(kb);
    int nmax=std::ceil((cut-lb)/(2.*lb));
    nmax++;
    double sum=0.;
    double bpow=pow(kb,n+2);
    double dfac=1.;
    double dsig=1.;
    double dpar=(double)n+2.;
    for(int i=0;i<=nmax;i++){
       sum+=dsig/(dfac*dpar)*bpow;
       dpar+=2.;
       dfac*=((double)i+2.)*((double)i+3.);
       dsig*=-1.;
       bpow*=kb2;
    }
    return sum;
}
double calCISEb(double v, double k, double b, double a){
    if(b<a*numtol) return 0.;
    v*=a*a;
    b=b/a;
    k*=a;
    a=1.;
    bool sw=false;
    if(k<1.) sw=true;
  //  double as=sin(k*b);
  //  double ac=cos(k*b);
    double ai=a-b;
    double as=sin(k);
    double ac=cos(k);
    double asi=sin(k*ai);
    double aci=cos(k*ai);
    double Cia=gsl_sf_Ci(k);
    double Sia=gsl_sf_Si(k);
    double Ciai=gsl_sf_Ci(k*ai);
    double Siai=gsl_sf_Si(k*ai);
    double ss0=(aci-ac)/k;
    double cs0=(as-asi)/k;
    double Co0=-ai*Ciai+Cia-cs0;
    double So0=-ai*Siai+Sia-ss0;
    double Cne=ac*(Cia-Ciai)+as*(Sia-Siai);
    double Csum=Cne;
    const double reltol=1.0E-4;
    double comp=1.;
    double compara=log(1./(1.-b));
    double ssnp=ss0;
    double csnp=cs0;
    double Conp=Co0;
    double Sonp=So0;
    double bn=b;
    double vn=v;
    double facu=1.;
    double ssn=bn/k*aci-csnp/k;
    double csn=-bn/k*asi+ssnp/k;
    if(sw){
        ssn=sinint(1,k,b);
        csn=cosint(1,k,b);
    }
    double Con=0.5*(Conp-ssnp/k-bn*ai*Ciai+bn*asi);
    double Son=0.5*(Sonp+csnp/k-bn*ai*Siai-bn*aci);
    ssnp=ssn;
    csnp=csn;
    Conp=Con;
    Sonp=Son;
    //bn*=b;
    //ssn=-bn/k*aci+2.*ssnp/k;
    //csn=bn/k*asi-2.*csnp/k;
    //Con=(2.*(Conp-ssnp/k)-bn*ai*Ciai+bn*asi)/3.;
    //Son=(2.*(Sonp+csnp/k)-bn*ai*Siai-bn*aci)/3.;
    double Cadd=vn*(bn*(ac*Ciai+as*Siai)+2.*(ac*Con+as*Son))/facu;
    comp=fabs(Cadd/Cne);
    Csum+=Cadd;
    double nl=1.;
    double nlo=1.;
    double tole=1.+reltol;
    int imax=ceil(log(reltol)/log(v*b*b));
    while((comp>reltol || fabs(Csum)>compara*tole) && vn*v*bn*b*b>reltol){
        bn*=b;
        vn*=v;
        nl+=1.;
        if(nl>imax) break;
        ssn=bn/k*aci-nl*csnp/k;
        csn=-bn/k*asi+nl*ssnp/k;
        if(sw){
            ssn=sinint(nl,k,b);
            csn=cosint(nl,k,b);
        }
        Con=(nl*(Conp-ssnp/k)-bn*ai*Ciai+bn*asi)/nl;
        Son=(nl*(Sonp+csnp/k)-bn*ai*Siai-bn*aci)/nl;
        ssnp=ssn;
        csnp=csn;
        Conp=Con;
        Sonp=Son;
        nl+=1.;
        double test1=bn/k*aci;
        double test2=nl*csnp/k;
        double test3=-bn/k*asi;
        double test4=nl*ssnp/k;
        ssn=(bn*aci-nl*csnp)/k;
        csn=(-bn*asi+nl*ssnp)/k;
        if(sw){
            ssn=sinint(nl,k,b);
            csn=cosint(nl,k,b);
        }
        Con=(nl*(Conp-ssnp/k)-bn*ai*Ciai+bn*asi)/(nl+1.);
        Son=(nl*(Sonp+csnp/k)-bn*ai*Siai-bn*aci)/(nl+1.);
        ssnp=ssn;
        csnp=csn;
        Conp=Con;
        Sonp=Son;
        nlo+=1.;
        facu*=nlo;
        if((ssn>=pow(b,nl)/((double)nl+1.) || csn>=pow(b,nl)/((double)nl+1.)) && k>1.0E-3){
            std::cout<<"ssn csn shit "<<csn<<" "<<ssn<<" "<<pow(b,nl)/((double)nl+1.)<<" "<<aci<<" "<<asi<<" "<<bn/k<<
            " "<<k<<" "<<nl<<" "<<test1-test2<<" "<<test2<<" "<<test3-test4<<" "<<test4<<" "<<bn/k<<"\n";
            exit(1);
        }
        if(((int)nlo)%2==1){
            Cadd=vn*(bn*(ac*Ciai+as*Siai)+2.*nlo*(ac*Con+as*Son))/facu;
        } else {
            Cadd=-vn*(bn*(ac*Ciai+as*Siai)+2.*nlo*(ac*Con+as*Son))/facu;
        }
        comp=fabs(Cadd/Cne);
        Csum+=Cadd;

    } 
    if(fabs(Csum)>compara*(1.+reltol) || std::isinf(Csum) || std::isnan(Csum)){
        std::cout<<"compara "<<Csum<<" "<<log(1./(1.-b))
        <<" "<<nlo<<" "<<k<<" "<<a<<" "<<v<<" "<<b<<" "<<vn<<" "<<Cne<<" "<<Cadd<<" "<<" "<<Cia<<" "<<Ciai<<" "<<Sia<<
        " "<<Siai<<" "<<bn<<" "<<Son<<" "<<Con<<" "<<ssn<<" "<<csn<<" "<<Co0<<" "<<So0<<" "<<"\n";
        std::cout<<imax<<"\n";
     //   exit(1);
     Csum=compara;
    }

    return Csum;
}



//cal integral for moderate v and small k
double calCISEc(double v, double k, double b, double a){
    if(b<a*numtol) return 0.;
    v*=a*a;
    b=b/a;
    k*=a;
    a=1.;
    double vb=v*b*b;
    double lab=log(vb);
    double svi=1./sqrt(v);
    double Csum=0.;
    const double pref=-0.5;
    double calbval=gsl_sf_gamma_inc_P(0.5,vb);
    double cal0val=gsl_sf_gamma(0.5);
    calbval*=cal0val;
    double calbvalf=gsl_sf_gamma_inc_P(1.,vb);
    double cal0valf=gsl_sf_gamma(1.);
    calbvalf*=calbvalf;
    double alip=pow(v,-0.5);
    double kb=k*b;
   // double kb2=kb*kb;
    double lb=log(kb);
    int nmax=std::ceil((cut-lb)/(2.*lb));
    const double reltol=1.0E-4;
    double comp=1.0;
    double compara=log(1./(1.-b));
    int n=0;
    double cbvn=calbval;
    double c0vn=cal0val;
    double cbvnf=calbvalf;
    double c0vnf=cal0valf;
    double k2=k*k;
    double tole=1.+reltol;
    while((comp>reltol || fabs(Csum)>compara*tole)){
        double ssum=0.; 
        double cbvns;
        double c0vns;
        if(n%2==0){
            cbvns=cbvn;
            c0vns=c0vn;
        } else {
            cbvns=cbvnf;
            c0vns=c0vnf;
        }
        double alips=alip;
        double facu=1.;
        double sign=1.;
        double kpow=1.;
        for(int i=0;i<nmax;i++){
           ssum+=sign/facu*(cbvns)*alips*kpow;
           sign*=-1.;
           facu*=((double)i+1.)*((double)i+2.);
           double fl=0.5*((double)n+1.)+(double)i;
           cbvns=(fl)*cbvns-exp(-vb+(fl)*lab);
           c0vns=(fl)*c0vns;
           alips/=v;
           kpow*=k2;
        }
        comp=fabs(ssum/Csum);
        if(Csum==0.) comp=1.;
        Csum+=ssum;
        n++;
        //adjust cbvn c0vn alip for next step
        double expada=exp(-vb+0.5*((double)n+1.)*lab);
        alip*=svi;
        if(n%2==0){
            cbvn=0.5*((double)n-1.)*cbvn-expada;
            c0vn=0.5*((double)n-1.)*c0vn;
        } else {
            cbvnf=0.5*((double)n-1.)*cbvnf-expada;
            c0vnf=0.5*((double)n-1.)*c0vnf;
        }
    }
    Csum*=pref;
    if(fabs(Csum)>compara*(1.+reltol) || std::isinf(Csum) || std::isnan(Csum)){
        std::cout<<"compara2 "<<Csum<<" "<<log(1./(1.-b))
        <<" "<<k<<" "<<a<<" "<<v<<" "<<b<<" "<<cbvn<<" "<<n<<" "<<nmax<<"\n";
     //   exit(1);
 //    Csum=compara;
    }
    return Csum;
}
//cal integral for moderate v and small k
double calCISEd(double v, double k, double b, double a){
    if(b<a*numtol) return 0.;
    v*=a*a;
    b=b/a;
    k*=a;
    a=1.;
    double vb=v*b*b;
//    double lab=log(vb);
    double sv=sqrt(v);
   // double svi=1./sqrt(v);
    double Csum=0.;
   // const double pref=-0.5;
    double E0=0.5*sqrt(pi/v)*std::erf(b*sv);
    double E1=(1.-exp(-vb))/(2.*v);
// E(n)=(E(n-2)-b^(n-1)/((n-1))*exp(-v*b^2))/v
    double kb=k*b;
    double kb2=kb*kb;
    double lb=log(kb);
    int nmax=std::ceil((cut-lb)/(2.*lb));
    const double reltol=1.0E-4;
    double comp=1.0;
    double compara=log(1./(1.-b));
    int n=0;
    //double Es=E0;
    double Een=E0;
    double Eun=E1;
    double k2=k*k;
    double bnm=1.;
    double b2=b*b;
    double evb=exp(-vb);
    double Cpa;
    double tole=1.+reltol;
    bool saf=false;
    while((comp>reltol/* || fabs(Csum)>compara*tole*/)){
        double ssum=0.;
        double Esr;
        double facu=1.;
        double sign=1.;
        double kpow=1.;
        if(n%2==0){
            Esr=Een;
        } else {
            Esr=Eun;
        }
        double bnmr=bnm*b;
        for(int i=0;i<nmax;i++){
           ssum+=sign/facu*(Esr)*kpow;
           sign*=-1.;
           facu*=((double)i+1.)*((double)i+2.);
           kpow*=k2;
           Esr=(Esr*(2.*((double)i+1.)+(double)n-1.)-bnmr*evb)/(2.*v);
           bnmr*=b2;
        }
       // ssum+=Esr;
        comp=fabs(ssum/Csum);
        if(n==0) comp=1.;
        Cpa=Csum;
     /*   if(n>21){
            std::cout<<"cra "<<Csum<<" "<<ssum<<" "<<comp<<" "<<Een<<" "<<Eun<<" "<<compara<<" "<<nmax<<" "<<saf<<
            " "<<v<<" "<<b<<" "<<k<<"\n";
            if(saf)
            exit(1);
            saf=true;
        }*/
        Csum+=ssum;

        n++;
        if(n>1){
            bnm*=b;
            if(n%2==0){
               Een=(Een*((double)n-1.)-bnm*evb)/(2.*v);
            } else {
               Eun=(Eun*((double)n-1.)-bnm*evb)/(2.*v);
            }
         /*   if(n==2){
                std::cout<<"testc "<<Een<<" "<<v<<" "<<b<<" "<<E1<<" "<<E0<<" "<<compara<<"\n";
               exit(1);
            }*/
        } 
    }
   // Csum*=pref;
    if(/*fabs(Csum)>compara*(1.+reltol) ||*/ std::isinf(Csum) || std::isnan(Csum)){
        std::cout<<"compara3 "<<Csum<<" "<<log(1./(1.-b))
        <<" k "<<k<<" "<<a<<" "<<v<<" "<<b<<" "<<n<<" "<<nmax<<" "<<Een<<" "<<Eun<<" "<<E0<<" "<<E1<<" "
        <<bnm<<" "<<Cpa<<"\n";
      //  exit(1);
 //    Csum=compara;
    }
    return Csum;
}
double calCIr(double v, double k, double b, double a){
    //    return calCISE(v,k,b,a);
    if(v*b*b<0.1){
  //      return 0.;
     //   return calCISEb(v,k,b,a);
 //    std::cout<<"in here\n";
//     return calCIiter(v,k,b,a);
   // return log(1./(1.-(b/a)));
  // } else if(k*a<0.1){
       //or incorparate here
   // return calCISEd(v,k,b,a);
  //  return log(1./(1.-(b/a)));
   }
  // std::cout<<"out here\n";
  // std::cout<<"va "<<v*a*a<<" a "<<a<<" ka "<<k*a<<"\n";

  //  return log(1./(1.-(b/a)));
   // case small k moderate v
   if(b<a*numtol || b==a) return 0.;
    double C0=0.;
    double S0=0.;
    b=b/a;
    v*=a*a;
    k*=a;
 //   b=0.0218789;
 //   k=8.06858;
 //   v=5.4225;
 //   calC0S0(v,k,b,&C0,&S0);
 //   std::cout<<"shitte "<<C0<<" "<<S0<<"\n";
    double b2=b*b;
  //  double v2=v*v;
    double expf=exp(-v*b2);
    double sb=sin(k*b);
    double cb=cos(k*b);
    double lsb=sb*expf;
    double lcb=cb*expf;
    double tv2i=1./(2.*v);
    double C1=-(k*S0+lcb-1.)*tv2i;
    double S1=(k*C0-lsb)*tv2i;
    const double reltol2=1.0E-4;
    double ai=1.;//a;
    double Csum=(C0+C1/ai)/ai;
   // ai*=a*a;
    double reli=1.;
    double Cn=C1;
    double Cnm=C0;
    double Sn=S1;
    double Snm=S0;;
    double bs=b*lsb;
    double bc=b*lcb;
    double i=2.;
    bool sofar=false;
    double compa=log(1./(1.-b));
    double lob=log(b);
    int ilim=60;//ceil(-gsl_sf_lambert_W0(lob/(log(1.-b)*reltol2))/lob);
  //  if(Csum>compa){ std::cout<<"fat "<<Csum<<" "<<C0<<" "<<C1<<" "<<S0<<" "<<S1<<" "<<compa<<" "<<b<<" "<<k<<" "<<v<<"\n";
  //  exit(1);
  //  }
  //  if(Csum!=0.0){
          double C2=0.;
          double S2=0.;
          double C3=0.;
          double S3=0.;
        while((reli>reltol2 || fabs(Csum)>compa*(1.+reltol2))){
            double Cnsa=Cn;
            double Snsa=Sn;
            Cn=tv2i*((i-1.)*Cnm-k*Snsa-bc);
            Sn=tv2i*((i-1.)*Snm+k*Cnsa-bs);
            Cnm=Cnsa;
            Snm=Snsa;
            double Cnad=Cn/ai;
            reli=fabs(Cnad/Csum);
            Csum+=Cnad;
           // if(Csum==0.0 && Cn!=0.) reli=1.;
        //    else if(Csum==0.) break;
       //     if(reli<=reltol2 && !sofar) sofar=true;
        //    else sofar=false;
            bc*=b;
            bs*=b;
          // bc=lcb*pow(b,i);
          // bs=lsb*pow(b,i);

          if(i==16){
            C2=Cn;
            S2=Sn;
          }
          if(i==17){
            C3=Cn;
            S3=Sn;
     //       std::cout<<"valout "<<C0<<" "<<S0<<" "<<C1<<" "<<S1<<" "<<C2<<" "<<S2<<" "<<C3<<" "<<S3<<" "<<
       //     Cnm<<" "<<Snm<<" "<<Cn<<" "<<Sn<<"\n";
          //  exit(1);
          }
            i+=1.;
          //  ai*=a;

        //  if(fabs(Csum)<tomin) return 0.0;
           /* if(fabs(Csum)>compa && i>1000.){// || fabs(Sn)>compa ){
                std::cout<<"too much csum "<<Csum<<" "<<Cnad<<" "<<C0<<" "<<S0<<" "<<C1<<" "<<S1<<" "<<v<<" "<<b<<"\n "
                <<k<<" "<<i<<" "<<k/sqrt(v)<<" "<<lsb<<" "<<lcb<<" "<<bs<<" "<<bc<<" "<<tv2i<<" "<<compa<<"\n"
                <<" "<<C2<<" "<<S2<<" "<<C3<<" "<<S3<<"\n";
                exit(1);
            }*/
            if(i>ilim) break;
        }
  //  }
   // std::cout<<"finito\n";
    if((Csum>log(1./(1.-b))||std::isnan(Csum)) && k<1.0E-2){
        std::cout<<"Csum "<<Csum<<" C0 "<<C0<<" S0 "<<S0<<" C1 "<<C1<<
        " S1 "<<S1<<" C2 "<<C2<<" S2 "<<S2<<" C3 "<<C3<<" S3 "<<S3<<
        " tv2i "<<tv2i<<" a  "<<a<<" b  "<<log(1./(1.-b))<<" Hsor  "<<Hsor<<" ai  "<<ai<<" b "<<b<<" v "<<v<<" k "<<k<<
        "kk/v"<<k*k/v<<" bb/v "<<b*b/v<<" "<<k*S0<<" "<<lcb/b-1.<<" "<<expf<<" "<<cb<<" "<<sb<<"\n";
    //    if(Csum<0.) std::cout<<"neg val\n";
    exit(1);
     //   std::cout<<"k v b "<<k<<" "<<v<<" "<<b<<"\n";
    }
    return Csum;
}



















std::vector<double> gridinb(){
    std::vector<double> q(Npsi,0.0);
    q[0]=ls;
    q[Npsi-1]=us;
    double h=exp(log(us/ls)/((double) Npsi-1));
    for(int i=1;i<Npsi-1;i++){
        q[i]=q[i-1]*h;
    }
    return q;
}
std::vector<double> gwgrid(){
    std::vector<double> q(Ngw,0.0);
    q[0]=lgw;
    q[Ngw-1]=ugw;
    double h=exp(log(ugw/lgw)/((double) Ngw-1));
    for(int i=1;i<Ngw-1;i++){
        q[i]=q[i-1]*h;
    }
    return q;
}
const double gpf=1./6.;
double gaSff(double a){
    return g0/g0;
}


double Omartevoc(double a){
  //  std::cout<<taubuild<<" "<<(a-aini)/Hsor<<" "<<tauLini<<"\n";
    if(a<taubuild*Hsor+aini){
        return pow(1.-(taubuild-(a-aini)/Hsor)/taubuild,gammab);
    } else if(a>(tauLini)*Hsor+aini){
        if(hrate){
            return pow(tauLini/((a-aini)/Hsor),gamma2h);
        } else {
            return pow(tauLini/((a-aini)/Hsor),gamma2);
        }
    } else {
        return 1.;
    }
}

double Levoa(double a){
//    if(a==aini) return 1.
   // if((a-aini)/(Hsor*tauLini)<1) return 1.;
    if(a<taubuild*Hsor+aini){
        return 1.;
    } else if(a<(tauLini)*Hsor+aini){
        return 1.;
    }
    if(hrate){
        return pow((a-aini)/(Hsor*tauLini),gamma1h);
    } else {
        return pow((a-aini)/(Hsor*tauLini),gamma1);
    }
}


/*double Omartevoc(double a){
  //  std::cout<<taubuild<<" "<<(a-aini)/Hsor<<" "<<tauLini<<"\n";
    if(a<taubuild*Hsor+aini){
        return pow(1.-(taubuild-(a-aini)/Hsor)/taubuild,gammab);
    } else if(){
        if(hrate){
            return pow(tauLini/((a-aini)/Hsor-taubuild+tauLini),gamma2h);
        } else {
            return pow(tauLini/((a-aini)/Hsor-taubuild+tauLini),gamma2);
        }
    }
}

double Levoa(double a){
//    if(a==aini) return 1.
   // if((a-aini)/(Hsor*tauLini)<1) return 1.;
    if(a<taubuild*Hsor+aini){
        return 1.;
    }
    if(hrate){
        return pow((a-aini-(taubuild-tauLini)*Hsor)/(Hsor*tauLini),gamma1h);
    } else {
        return pow((a-aini-(taubuild-tauLini)*Hsor)/(Hsor*tauLini),gamma1);
    }
}*/
/*double sfevo(double a){
//    if(a==aini) return 1.
   // if((a-aini)/(Hsor*tauLini)<1) return 1.;
   switch(sfmodel){
        case 0:
        return 0.;
        break;
        case 1:
        if(a<taubuild*Hsor+aini){
            return initsf*(1.-pow(1.-(taubuild-(a-aini)/Hsor)/taubuild,gammab));
        } else {
            return 0.;
        }
        break;
        case 2:
        return initsf;
        break;
        case 3:
        if(a<taubuild*Hsor+aini){
            return initsf*(1.-pow(1.-(taubuild-(a-aini)/Hsor)/taubuild,gammab));
        } else if(a<(taubuild+sfsf*tauLini)*Hsor+aini){
            return initsf*(1.-pow(1.-(sfsf*tauLini-(a-aini-taubuild*Hsor)/Hsor)/(sfsf*tauLini),gammab));
        } else {
            return 0.;
        }
        default:
        return 0.;
        break;
   }
}*/
double sfevo(double a){
//    if(a==aini) return 1.
   // if((a-aini)/(Hsor*tauLini)<1) return 1.;
   switch(sfmodel){
        case 0:
        return 0.;
        break;
        case 1:
        if(a<tauLini*Hsor+aini){
            return initsf*(1.-pow(1.-(tauLini-(a-aini)/Hsor)/tauLini,gammab));
        } else {
            return 0.;
        }
        break;
        case 2:
        return initsf;
        break;
        case 3:
        if(a<taubuild*Hsor+aini){
            return initsf;//initsf*(1.-pow(1.-(taubuild-(a-aini)/Hsor)/taubuild,gammab));
        } else if(a<(sfsf*tauLini)*Hsor+aini){
            return initsf*(1.-pow(((a-aini-abuild)/Hsor)/(sfsf*tauLini),gammab));
        } else {
            return 0.;
        }
        default:
        return 0.;
        break;
   }
}
















class enspec{
    private:
        double amp;
        double maxamp;
        double t; //here t represents the scalefactor
        double maxL;
        double L;
        double cl;
        double tp;
        double Kc;
        void timevo();
    public:
        enspec();
        enspec(double, double, double, double);
        double get_enk(double);
        double get_enk2(double);
        double get_t();
        void set_t(double);
        double get_L();
        void set_L(double);
        double get_cl();
        void set_cl(double);
        double getamp();
};

enspec::enspec(){}

enspec::enspec(double amp0, double t0, double L0, double cl0){
    cl=cl0;
    Kc=tpi*L/cl;
    maxL=L0;
    maxamp=amp0;
    amp=maxamp*Omartevoc(t0);
    L=maxL*Levoa(t0);
    t=t0;
}

// here k defined as (k*L/2*pi)

double enspec::get_enk(double k){
    k*=L/tpi;
    double k2=k*k;
    if(k<Kc){
        double vala=0.;
        vala=enorm*amp*L*k2*k2/pow(Lc+k2,Kolpot);
   //     if(k>0.9 && k<1.1 && vala>0.1) std::cout<<vala<<"\n";
        return vala;
    } else {
        return 0.0;
    }
}
double enspec::get_enk2(double k){
    k*=L/tpi;
    double k2=k*k;
    if(k<Kc){
        double vala=0.;
        vala=enorm2*amp*L*k2*k2/pow(Lc2+k2,Kolpot2);
   //     if(k>0.9 && k<1.1 && vala>0.1) std::cout<<vala<<"\n";
        return vala;
    } else {
        return 0.0;
    }
}

double enspec::getamp(){
    return amp;
}
double enspec::get_t(){
    return t;
}

void enspec::set_t(double ti){
    t=ti;
    amp=maxamp*Omartevoc(t);
    L=maxL*Levoa(t);
    Kc=tpi*L/cl;
   // std::cout<<"amp : "<<amp<<" "<<L<<" "<<ti<<" "<<maxamp<<" "<<maxL<<"\n";
}

double enspec::get_L(){
    return L;
}

void enspec::set_L(double Li){
    L=Li;
}

double enspec::get_cl(){
    return cl;
}

void enspec::set_cl(double clo){
    cl=clo;
}

double a_to_conft(double a){
    return a/(Hsor*gaSff(a));
}

double conft_to_a(double t){
    return t*Hsor*gaSff(t);
}

double rate0(double dt, double vscaleq,double vscalep){
    double mulf=1.;
    if(addcon) mulf=cos(dt*vscaleq*sqrt(3.))*cos(dt*vscalep*sqrt(3.));
    return exp(-0.5*dt*dt*(vscaleq*vscaleq+vscalep*vscalep))*mulf;
}

const double epox=4./3.;
double rate5(double dt, double vscaleq,double vscalep){
    double efac=pow(dt*vscaleq,epox)+pow(dt*vscalep,epox);
    return exp(-efac*0.5)*cos(0.5*efac*sqrt(3.));
}



























//const double dsnorm=sqrt(3.*3./(2.*4.*kap)*tpi*tpi);
const double dsnorm=sqrt(2./(3.*kap)*tpi*tpi);
const double dscalepref=co1*dsnorm; //assumes K^4 spec
const double funcintpref=8./3.;//2./(3.*pi*pi*pi);//8.*pow(pi,-3)/3.;//8.*4.*2./(tpi*tpi*4.*3.);//pow(pi,-5);//*Omr;//16.0*Omr*rhoc*GN/(9.0*pi*pi*pi*H0*H0);//prefactor for Omgw
class gwenint{
    private:
    enspec vspec;
    enspec bspec;
    enspec hvspec;
    enspec hbspec;
    double vep;
    double veq;
    double bep;
    double beq;
    double vep2;
    double veq2;
    double dscale(double);
    double dscalelag(double);
    std::vector<double> A;
    std::vector<double> plaw;
    std::vector<double> biggrid;
    std::vector<double> yongrid;
    std::vector<double> subintval;
    std::vector<double> Alag;
    std::vector<double> plawlag;
    std::vector<double> yongridlag;
    std::vector<double> subintvallag;
    double gfac;
    double gfacp;
    double soundf;
    double nssoundf;
    double nsoundf;
    double kinfrac;
    double magfrac;
    double dscaleprefl; //L*Omratio
    double dscaleq2;
    double dscaleq;
    double dscalepg;
    double dscalek2;
    double Aplawosc;
    double Aplawosclag;
    void setup_Aplaw();
    void setup_Aplawlag();
    double t;
    double a;
    double k;
    double q;
    double p;
    double k2;
    double p2;
    double L;
    double q2;
    double dsca;
    double tana;
    double tanb;
    double cas2;
    double prsca;
    double prbo;
    double prao;
    double qrsca;
    double intfark;
    double coefp;
    double coefm;
    double coefmq;
    double coefmp;
    double Ctf;
    gsl_integration_workspace *wsp;
    gsl_integration_workspace *wspp;
    gsl_integration_workspace *wspq;
    gsl_integration_qawo_table *qawoit;
    gsl_integration_qawo_table *qawoit2;
    public:
    double bin;
    double ain;
    double iabin;
    gwenint();
    gwenint(double, double,double,enspec,enspec,enspec,enspec);
    void setkat(double, double);
    void setat(double);
    void setq(double);
    void setp(double);
    double getq();
    double getp();
    double drate(double);
    double drates2(double);
    double drateint();
    double intfunc();
    double intintfunc();
    double finalint();
    double getiabin();
    double getbin();
    double getmaxamp();
    double gettana();
    double gettanb();
    double getqrsca();
    double getprsca();
    double getprao();
    double getprbo();
    double calCI(double,double,double,double);
    double calCItrr(double,double,double,double);
    double calCItrr5(double,double,double,double);
    double calCItrial2(double,double,double,double);
    void calC0S0(double,double,double,double*,double*);
    double calCIiter(double,double,double,double);
    double magsvel;
    double cb; //ctab
    double sb; //stab
    double as;
    double ac;
    double Cia;
    double Sia;
    double Ciai;
    double Siai;
    double el0;
    double el1;
    double el2;
    double el3;
    double el4;
    double ilimCI;
    double compa;
    double ampiz;

};

gwenint::gwenint(){}

double gwenint::getprao(){
    return prao;
}
double gwenint::getprbo(){
    return prbo;
}
double gwenint::getprsca(){
    return prsca;
}
double gwenint::getqrsca(){
    return qrsca;
}
double gwenint::gettana(){
    return tana;
}
double gwenint::gettanb(){
    return tanb;
}

double gwenint::getmaxamp(){
    return std::max(vspec.getamp(),bspec.getamp());
}

double gwenint::getiabin(){
    return iabin;
}
double gwenint::getbin(){
    return bin;
}




void gwenint::calC0S0(double vv, double kv, double bv, double *C0, double *S0){
   /* if(bv<numtol){
        *C0=0.0;
        *S0=0.0;
        return;
    }*/
 //   v=16.4095;
 //   k=0.404267;
   // b=0.005552845;
    const double reltol=1.0E-6;
    const double lreltol=log(reltol);
    double sv=sqrt(vv);
    double x=bv*sv;
    double y=0.5*kv/sv;
    double ta=2.0*x;
    double tab=ta*y;
    double x2=x*x;
    double k2=kv*kv;
    double k2i4v=k2*0.25/vv;
    double ctab=cos(tab);//cb;//cos(tab);
    double stab=sin(tab);//sb;//sin(tab);
    double ctabm=1.-ctab;
    if(tab<1.0E-2){
        ctab=1.-tab*tab*(0.5-tab*tab/24.);
        stab=tab*(1.-tab*tab/6.);
        ctabm=tab*tab*(0.5-tab*tab/24.);
    }
    double pref=0.5*sqrt(pi)/sv;//*exp(-k2i4v);
    double ek2i4v=exp(-k2i4v);
    double emx2=exp(-x2)/pi;
    double k2i4vx2=-k2i4v-x2;
    double eext=exp(k2i4vx2);
    double Cs=pref*(ek2i4v*gsl_sf_erf(x)+eext/(2.*pi*x)*ctabm);
    int iopt=ceil(kv/sv);
    int ilim=iopt+10;
    int istart=iopt-10;
    if(istart<1) istart=1;
  //  istart=1;
 //   double check=-4.*b*b*v-4*maxp;
 //   double check=sqrt(4.*(y*y-lreltol));

   // int istart=1;
    /*if(check>=0.){
        double ioptc;
        if(iopt>10) ioptc=iopt-10;
        else ioptc=1;
        istart=floor(std::min(std::max(1.,iopt+floor(sqrt(check))),ioptc));
        // choice of istart
        
    } else istart=1.;*/

    double rel0=1.;
    double bac=2.*y;
    //istart=floor(std::max(bac-check-2.,2.))-1;
//    ilim=ceil(bac+check)+2;
 //   istart=std::max(1,ilim-20);
   // ilim=10000;
    //istart=1.;

   // istart=1;
    int i=istart;
    double tole=1.+reltol;
    while(((rel0>reltol || fabs(Cs)>=bv*(1.+numtol))*tole) ||i<=ilim){
        double i2=(double)i*i;
        double eiy=exp(-0.25*i2+(double)i*y+k2i4vx2);
        double eimy=exp(-0.25*i2-(double)i*y+k2i4vx2);
        double csh=0.5*(eiy+eimy);
        double ssh=0.5*(eiy-eimy);
        double eextl=exp(k2i4vx2-i2*0.25);
        double Csc=pref*(2./pi)*(2.*x*(eextl-ctab*csh)+(double)i*stab*ssh)/(i2+4.*x2);
     //   if(std::isnan(Csc)) std::cout<<i<<" "<<csh<<" "<<ssh<<" "<<emx2<<" "<<pref<<" "<<y<<"\n";
        rel0=fabs(Csc/Cs);
        Cs+=Csc;
        i+=1;
    }
    *C0=Cs;
    //double Ss=-pref/(2.*pi*x)*(eext*stab-ek2i4v*tab);
    double Ss=-pref*(1./(2.*pi*x)*(eext*stab)-ek2i4v*y/pi);
    double Ssi=Ss;
    rel0=1.;
 //   istart=1;
    i=istart;
   // ilim=10000;
    while(((rel0>reltol || fabs(Ss)>=bv*(1.+numtol))*tole) || i<=ilim){
        double i2=(double)i*i;
        double eiy=exp(-0.25*i2+(double)i*y+k2i4vx2);
        double eimy=exp(-0.25*i2-(double)i*y+k2i4vx2);
        double eiy2=exp(-0.25*i2+(double)i*y-k2i4v);
        double eimy2=exp(-0.25*i2-(double)i*y-k2i4v);
        double csh=0.5*(eiy+eimy);
        double ssh=0.5*(eiy-eimy);
        double ssh2=0.5*(eiy2-eimy2);
        double Ssc=-pref*(2./pi)*((2.*x*(stab*csh)+(double)i*ctab*ssh)/(i2+4.*x2)-ssh2/((double)i));
        rel0=fabs(Ssc/Ss);
     //   if(i>1900){
    //    std::cout<<Ss<<" "<<Ssc<<" "<<x<<" "<<y<<" "<<k<<" "<<v<<" "<<b<<" "<<eext<<" "<<stab<<" "
     //   <<iopt<<"\n";
     //   exit(1);
     //   }
        Ss+=Ssc;
        i+=1;
    }
    *S0=Ss;
   /* if(*C0!=0. && *S0==0.){
        std::cout<<"more 1295 S0 "<<*S0<<" i "<<i<<" C0 "<<*C0<<" k "<<kv<<" v "<<vv<<" b "<<bv<<" Ssi "<<Ssi<<
        " ilim "<<ilim<<" istart "<<istart<<"\n";
    }*/
 /*   if(std::isnan(*C0) || std::isnan(*S0)){
        std::cout<<"aua\n";
        exit(1);
    }*/
}


double gwenint::calCIiter(double vv, double kv, double bv, double av){
    if(bv<av*numtol) return 0.;
    vv*=av*av;
  //  b=b/a;
  //  k*=a;
   // a=1.;
    //bool sw=false;
    //if(k<0.1) sw=true;
  //  double as=sin(k*b);
  //  double ac=cos(k*b);
  /*  double ai=a-b;
    double as=sin(k);
    double ac=cos(k);
    double Cia=gsl_sf_Ci(k);
    double Sia=gsl_sf_Si(k);
    double Ciai=gsl_sf_Ci(k*ai);
    double Siai=gsl_sf_Si(k*ai);
    double Cne=ac*(Cia-Ciai)+as*(Sia-Siai);
    double cb=cos(k*b);
    double sb=sin(k*b);

    double el0=Cne;
    double el1=0.;
    double el2=0.;
    double b2=b*b;
    double b3=b*b2;
    double k2=k*k;
    if(k*b<1.0E-4){
        el1=el0-b*(1.+0.5*b)+b2*k2/24.*(4.*b+3.*b2);
        el2=el0-b*(1.+0.5*b+b2/3.+0.25*b3)+k2*b3/120.*(20.+15.*b+12.*b2+10.*b3);
    } else {
        el1=el0+((1.-cb)/k-(1.+b)*sb)/k;
        el2=el0+(-6+k2-cb*(6.-k2*(1.+2.*b+3.*b2))-sb*k*(-2.+k2*(1.+b2*(1.+b)+b)-6.*b))/(k2*k2);
    }*/
//    std::cout<<el0<<" "<<v<<" "<<el1<<" "<<el2<<"\n"; 
    double v2=vv*vv;
    double v3=vv*v2;
    double v4=vv*v3;
    double resu= el0-vv*el1+v2*0.5*el2-v3/6.*el3+v4/24.*el4;
    if(fabs(resu)>compa && fabs(resu)<compa*1.01){
        resu=compa;
    } else if(fabs(resu)>compa*1.01){
     std::cout<<"not good "<<resu<<" v "<<vv<<" k "<<kv*av<<" b "<<bv/av<<" comp "<<compa<<" "<<el0<<" "<<el1<<
     " "<<el2<<" "<<el3<<" "<<el4<<"\n";   
    }
//    if(fabs(resu)>compa) std::cout<<"1337 gotta prob "<<resu<<" comp "<<compa<<" v "<<vv<<" k "<<kv*av<<" b "<<bv/av<<"\n";
    return resu;

    
}



    const double reltol2=5.0E-3;
    const int ilims=12;
double gwenint::calCI(double vv, double kv, double bv, double av){
  // cb=cos(kv*bv);
  // sb=sin(kv*bv);
 /*  if(vv*av*av<=0.9){
       double testir= calCIiter(vv,kv,bv,av);
//       if(kv*bv<1.5 && testir<0.0) std::cout<<" "<<el0<<" "<<el1<<" "<<vv*av*av<<" "<<kv*av<<" "<<testir<<"\n";
       return testir;
   }*/
   double ilimCIst=ilimCI;
   double compast=compa;
   double cbst=cos(kv*bv);//cb;
   double sbst=sin(kv*bv);//sb;
/*   double bmax=sltomin/sqrt(vv);
   if(bv>bmax){ 
       bv=bmax;
       cb=cos(kv*bv);
       sb=sin(kv*bv);
       double dbr=bv/av;
       ilimCI=ceil(-gsl_sf_lambert_W0(log(dbr)/(log(1.-dbr)*reltol2))/log(dbr));
       compa=log(1./(1.-dbr));
       av=bv+aini;
   }*/
   if(bv<av*numtol || bv==av){
//        std::cout<<"overt\n";   
       return 0.;

   }
    double C0=0.;
    double S0=0.;
    bv=bv/av;
    vv*=av*av;
    kv*=av;
 //   b=0.0218789;
 //   k=8.06858;
 //   v=5.4225;
    calC0S0(vv,kv,bv,&C0,&S0);
    if(C0==0.){
  //   std::cout<<"over\n";
     return 0.;
    }
 //   std::cout<<"shitte "<<C0<<" "<<S0<<"\n";
    double b2=bv*bv;
  //  double v2=v*v;
    double expf=exp(-vv*b2);
    //double sb=sin(k*b);
    //double cb=cos(k*b);
    double lsb=sbst*expf;
    double lcb=cbst*expf;
    double tv2i=1./(2.*vv);
    double C1=-(kv*S0+lcb-1.)*tv2i;
    double S1=(kv*C0-lsb)*tv2i;

    double ai=1.;//a;
    double Csum=(C0+C1/ai)/ai;
    double k2=kv*kv;
   // ai*=a*a;
    double reli=1.;
    double Cn=C1;
    double Cnm=C0;
    double Sn=S1;
    double Snm=S0;;
    double bs=bv*lsb;
    double bc=bv*lcb;
    double i=2.;
    bool sofar=false;
    //double compa=log(1./(1.-b));
  //  double lob=log(b);
    //int ilim=ceil(-gsl_sf_lambert_W0(lob/(log(1.-b)*reltol2))/lob);
    double C2=tv2i*(C0-kv*S1-bc);
    double S2=tv2i*(S0+kv*C1-bs);


    bc*=bv;
    bs*=bv;
    double C3=tv2i*(2.*C1-kv*S2-bc);
    double S3=tv2i*(2.*S1+kv*C2-bs);
//    if(C0+C1+C2+C3==0.) std::cout<<"what is wrong :\n";
  /*  if(C0+C1<0.){
        std::cout<<"first "<<C0<<" "<<C1<<" "<<C2<<" "<<C3<<" "<<vv<<" "<<kv<<" "<<bv<<"\n";
        exit(1);
    }*/
    double erel=C0+C1+C2+C3;
    if(fabs(erel)>log(1./(1.-bv))){
        std::cout<<"shitte "<<erel<<" "<<vv<<" "<<kv<<" "<<bv<<" "<<av<<"\n";
    }
    return erel;
    /*
    bc*=bv;
    bs*=bv;
    double C4=tv2i*(3.*C2-kv*S3-bc);
    double S4=tv2i*(3.*S2+kv*C3-bs);
    bc*=bv;
    bs*=bv;
    double C5=tv2i*(4.*C3-kv*S4-bc);
    double S5=tv2i*(4.*S3+kv*C4-bs);
    */


   /* if(S1<.0 && kv*bv<1.5){
        std::cout<<" "<<S1<<" "<<" "<<S0*kv<<" "<<kv*C0<<" "<<lsb<<" "<<bv<<"\n";
        exit(1);
    }*/
    double Ce0=C0;
    double Ce1=C2;
 //   double Ce2=C4;
    double Cu0=C1;
    double Cu1=C3;
//  double Cu2=C5;

    double C4=0.;
    double C5=0.;
    Csum+=C2+C3;
    // Csum+=C4+C5;
    double tv2is=-tv2i*tv2i;
    double Ee=expf*(cbst*(2.*vv*b2)-bv*sbst*kv);
    double Ee2=-expf*cbst;
//    double b2v=bv*bv;
    i=4.;
    // i=6;
//    ilimCI=3;
    while((reli>reltol2 || fabs(Csum)>compa*(1.+reltol2)) && i<=ilimCI){
        //double Ch=0.;
        double Cbck=0.;
        double Cadd=0.;
     //   Ee=expf*(cb*(2.*vv*b2-(i-2.))-sb*bv*kv);
        if((int)i%2==0){
           Cbck=Ce1;
           Ce1=tv2is*((k2-2.*vv*(2.*i-3.))*Ce1+(i-2.)*(i-3.)*Ce0+pow(bv,i-3)*(Ee+Ee2*(i-2.)));
           //Cbck=Ce2;
      //     Ce2=tv2is*(k2*(Ce2-Ce1)-2.*vv*((2.*i-3.)*Ce2-(2.*i-7.)*Ce1)+(i-2.)*(i-3.)*Ce1-(i-4.)*(i-5.)*Ce0)+b2v*Ce2; //homogenous
           //  Ce0=Ce1;
           // Ce1=Ce2;

    //       Ce1=tv2is*((2.*vv*(2.*i-4.))*Ce1-(i-2.)*(i-3.)*Ce0+pow(bv,i-3)*Ee);
           Ce0=Cbck;
           Cadd=Ce1;
        } else {
           Cbck=Cu1;
           Cu1=tv2is*((k2-2.*vv*(2.*i-3.))*Cu1+(i-2.)*(i-3.)*Cu0+pow(bv,i-3)*(Ee+Ee2*(i-2.)));
     //      Cu1=tv2is*((2.*vv*(2.*i-4.))*Cu1-(i-2.)*(i-3.)*Cu0+pow(bv,i-3)*Ee);
           Cu0=Cbck;
           Cadd=Cu1;
        }
        if(i==5) C5=Cadd;
        if(i==4) C4=Cadd;
        reli=fabs(Cadd/Csum);
        Csum+=Cadd;
  //      if(C3>1.0 && kv*bv<1.5) std::cout<<" "<<C3<<"\n";
        i+=1.;
     //   if(i>ilimCI) break;
    }
    if(fabs(Csum)<compa*1.01 && fabs(Csum)>compa){
        Csum=compa;
    }
/*    if(((fabs(Csum)>compa)||std::isnan(Csum))){
        std::cout<<"Csump "<<Csum<<" C0 "<<C0<<" S0 "<<S0<<" C1 "<<C1<<
        " S1 "<<S1<<" C2 "<<C2<<" S2 "<<S2<<" C3 "<<C3<<" S3 "<<S3<<" C4 "<<C4<<" C5 "<<C5<<
        " tv2i "<<tv2i<<" a  "<<av<<" b  "<<compa<<" Hsor  "<<Hsor<<" ai  "<<ai<<" b "<<bv<<" v "<<vv<<" k "<<kv<<
        "\n";
    }*/
   ilimCI=ilimCIst;
   compa=compast;
   cb=cbst;
   sb=sbst;
   /* if(Csum<0.0){ std::cout<<"C0 "<<C0<<" C1 "<<C1<<" C2 "<<C2<<" C3 "<<C3<<" S0 "<<S0<<" S1 "<<S1<<" v "<<vv<<" k "<<kv<<" b "<<bv<<" S "<<Csum<<"\n";
    exit(1);
    }*/
    return Csum;
}

double gwenint::calCItrr(double vv, double kv, double bv, double av){
  // cb=cos(kv*bv);
  // sb=sin(kv*bv);
   if(vv*av*av<=0.9){
       double testir= calCIiter(vv,kv,bv,av);
//       if(kv*bv<1.5 && testir<0.0) std::cout<<" "<<el0<<" "<<el1<<" "<<vv*av*av<<" "<<kv*av<<" "<<testir<<"\n";
       return testir;
   }
   double ilimCIst=ilimCI;
   double compast=compa;
   double cbst=cb;
   double sbst=sb;
   double bmax=sltomin/sqrt(vv);
   if(bv>bmax){ 
       bv=bmax;
       cb=cos(kv*bv);
       sb=sin(kv*bv);
       double dbr=bv/av;
       ilimCI=ceil(-gsl_sf_lambert_W0(log(dbr)/(log(1.-dbr)*reltol2))/log(dbr));
       compa=log(1./(1.-dbr));
       av=bv+aini;
   }
   if(bv<av*numtol || bv==av) return 0.;
    double C0=0.;
    double S0=0.;
    bv=bv/av;
    vv*=av*av;
    kv*=av;
 //   b=0.0218789;
 //   k=8.06858;
 //   v=5.4225;
    calC0S0(vv,kv,bv,&C0,&S0);
    if(C0==0.) return 0.;
 //   std::cout<<"shitte "<<C0<<" "<<S0<<"\n";
    double b2=bv*bv;
  //  double v2=v*v;
    double expf=exp(-vv*b2);
    //double sb=sin(k*b);
    //double cb=cos(k*b);
    double lsb=sb*expf;
    double lcb=cb*expf;
    double tv2i=1./(2.*vv);
    double C1=-(kv*S0+lcb-1.)*tv2i;
    double S1=(kv*C0-lsb)*tv2i;

    double ai=1.;//a;
    double Csum=(C0+C1/ai)/ai;
    double k2=kv*kv;
   // ai*=a*a;
    double reli=1.;
    //double Cn=C1;
   // double Cnm=C0;
   // double Sn=S1;
   // double Snm=S0;;
    double bs=bv*lsb;
    double bc=bv*lcb;
    double i=2.;
    bool sofar=false;
    double tv2is=-tv2i*tv2i;
    //double compa=log(1./(1.-b));
  //  double lob=log(b);
    //int ilim=ceil(-gsl_sf_lambert_W0(lob/(log(1.-b)*reltol2))/lob);
    double C2=tv2i*(C0-kv*S1-bc);
    double corrf=bv*(sb*bv*kv+cb*(1.-2.*vv*b2))/(sb*bv-2.*cb*vv*b2);
    double C3=C2*corrf+tv2is*(k2*(C1-C0*corrf)-(3.*C1-C0*corrf)/tv2i);
    corrf=bv*(sb*bv*kv+cb*(2.-2.*vv*b2))/(sb*bv+cb*(1.-2.*vv*b2));
    double C4=C3*corrf+tv2is*(k2*(C2-C1*corrf)-(5.*C2-3.*C1*corrf)/tv2i+2.*C0);
    double C5=0.;
    double Cn=C4;
    double Cnm=C3;
    double Cnm2=C2;
    double Cnm3=C1;
    double Cnm4=C0;
    Csum+=C3+C4;
    i=5.;
//    ilimCI=3;
    while((reli>reltol2 || fabs(Csum)>compa*(1.+reltol2)) && i<=ilimCI){
        double Cadd=Cn;
        corrf=bv*(sb*bv*kv+cb*((i-2.)-2.*vv*b2))/(sb*bv+cb*((i-3.)-2.*vv*b2));
        Cn=Cn*corrf+tv2is*((k2-2.*i)*(Cnm-Cnm2*corrf)+(3.*Cnm-5.*Cnm2*corrf)/tv2i+i*i*(Cnm3-Cnm4*corrf)-i*(5.*Cnm3-corrf*7.*Cnm4)+6.*Cnm3-12.*corrf*Cnm4);
        Cnm4=Cnm3;
        Cnm3=Cnm2;
        Cnm2=Cnm;
        Cnm=Cadd;
        if(i==5) C5=Cadd;
        //if(i==4) C4=Cadd;
        reli=fabs(Cadd/Csum);
        Csum+=Cadd;
  //      if(C3>1.0 && kv*bv<1.5) std::cout<<" "<<C3<<"\n";
        i+=1.;
     //   if(i>ilimCI) break;
    }
    if(fabs(Csum)<compa*1.01 && fabs(Csum)>compa){
        Csum=compa;
    }
    if(((fabs(Csum)>compa)||std::isnan(Csum))){
        std::cout<<"Csump "<<Csum<<" C0 "<<C0<<" S0 "<<S0<<" C1 "<<C1<<
        " S1 "<<S1<<" C2 "<<C2<<" C3 "<<C3<<" C4 "<<C4<<" C5 "<<C5<<
        " tv2i "<<tv2i<<" a  "<<av<<" b  "<<compa<<" Hsor  "<<Hsor<<" ai  "<<ai<<" b "<<bv<<" v "<<vv<<" k "<<kv<<
        "\n";
    }
   ilimCI=ilimCIst;
   compa=compast;
   cb=cbst;
   sb=sbst;
   /* if(Csum<0.0){ std::cout<<"C0 "<<C0<<" C1 "<<C1<<" C2 "<<C2<<" C3 "<<C3<<" S0 "<<S0<<" S1 "<<S1<<" v "<<vv<<" k "<<kv<<" b "<<bv<<" S "<<Csum<<"\n";
    exit(1);
    }*/
    return Csum;
}

//rework with olver solver
const double erriter=1.0E-4;
const int maxr=8;
double gwenint::calCItrial2(double vv, double kv, double bv, double av){
  // cb=cos(kv*bv);
  // sb=sin(kv*bv);
   if(vv*av*av<=0.8){
       double testir= calCIiter(vv,kv,bv,av);
//       if(kv*bv<1.5 && testir<0.0) std::cout<<" "<<el0<<" "<<el1<<" "<<vv*av*av<<" "<<kv*av<<" "<<testir<<"\n";
       return testir;
   }
   double ilimCIst=ilimCI;
   double compast=compa;
   double cbst=cb;
   double sbst=sb;
   double bmax=sltomin/sqrt(vv);
   if(bv>bmax){ 
       bv=bmax;
       cb=cos(kv*bv);
       sb=sin(kv*bv);
       double dbr=bv/av;
       ilimCI=ceil(-gsl_sf_lambert_W0(log(dbr)/(log(1.-dbr)*reltol2))/log(dbr));
       compa=log(1./(1.-dbr));
       av=bv+aini;
   }
   if(bv<av*numtol || bv==av) return 0.;
    double C0=0.;
    double S0=0.;
    bv=bv/av;
    vv*=av*av;
    kv*=av;
    calC0S0(vv,kv,bv,&C0,&S0);
    if(C0==0.) return 0.;
 //   std::cout<<"shitte "<<C0<<" "<<S0<<"\n";
    double b2=bv*bv;
  //  double v2=v*v;
    double expf=exp(-vv*b2);
    //double sb=sin(k*b);
    //double cb=cos(k*b);
    double lsb=sb*expf;
    double lcb=cb*expf;
    double tv2i=1./(2.*vv);
    double C1=-(kv*S0+lcb-1.)*tv2i;
//    double S1=(kv*C0-lsb)*tv2i;

    double ai=1.;//a;
    double Csum=(C0+C1/ai)/ai;
    double k2=kv*kv;
   // ai*=a*a;
    double reli=1.;
    double Cn=C1;
    double Cnm=C0;
    double bs=bv*lsb;
    double bc=bv*lcb;


    double tv2is=-tv2i*tv2i;
    double Ee=expf*(cb*(2.*vv*b2)-bv*sb*kv);
    double Ee2=-expf*cb;
    int iniL=ceil(0.5*(double)ilimCI)-1;
    int iniN=2*ceil((double)ilimCI*0.5+0.5);

    std::vector<double> ae(iniN,0.);
    std::vector<double> be(iniN,0.);
    std::vector<double> ce(iniN,1./tv2is);
    std::vector<double> de(iniN,0.);
    std::vector<double> ee(iniN,0.);
    std::vector<double> pe(iniN+1,0.);
    std::vector<double> au(iniN,0.);
    std::vector<double> bu(iniN,0.);
    std::vector<double> cu(iniN,1./tv2is);
    std::vector<double> du(iniN,0.);
    std::vector<double> eu(iniN,0.);
    std::vector<double> pu(iniN+1,0.);
    de[0]=Ee/bv;
    du[0]=Ee+Ee2;
    pe[0]=0.;
    pe[1]=1.;
    pu[0]=0.;
    pu[1]=1.;
    ee[0]=C0;
    eu[0]=C1;
    double Csume=C0;
    double Csumu=C1;
    double bve=bv;
    double bvu=b2;
    for(int m=1;m<iniN;m++){

        double mr=(double)m;
        ae[m]=-2.*mr*(2.*mr-1.);
        au[m]=-2.*mr*(2.*mr+1.);
        be[m]=k2-2.*vv*(4.*mr+1.);
        bu[m]=k2-2.*vv*(4.*mr+3.);
        de[m]=bve*Ee+bve*2.*mr*Ee2;//b2*(de[m-1]-Ee2*(m-1));
        du[m]=bvu*Ee+bvu*(2.*mr+1.)*Ee2;
        pe[m+1]=(be[m]*pe[m]-ae[m]*pe[m-1])/ce[m];
        pu[m+1]=(bu[m]*pu[m]-au[m]*pu[m-1])/cu[m];
        ee[m]=(ae[m]*ee[m-1]-de[m]*pe[m])/ce[m];
        eu[m]=(au[m]*eu[m-1]-du[m]*pu[m])/cu[m];
        bve*=b2;
        bvu*=b2;
    }
    double coeffe=fabs(0.5*erriter*ee[iniL-1]/(pe[iniL-1]*pe[iniL]));
    double coeffu=fabs(0.5*erriter*eu[iniL-1]/(pu[iniL-1]*pu[iniL]));
    double tce=ee[iniN-1]/(pe[iniN-1]*pe[iniN]);
    double tcu=eu[iniN-1]/(pu[iniN-1]*pu[iniN]);
    int iniNe=iniN;
    int iniNu=iniN;
    double mce=(double)iniN;
    int mre=(int)mce;
    while(fabs(tce)>coeffe && mre<maxr*iniL){
        ae.push_back(-2.*mce*(2.*mce-1.));
        be.push_back(k2-2.*vv*(4.*mce+1.));
        ce.push_back(1./tv2is);
       // de.push_back(de[mre-1]*b2);
        de.push_back(bve*Ee+bve*2.*mce*Ee2);//b2*(de[m-1]-Ee2*(m-1));
        pe.push_back((be[mre]*pe[mre]-ae[mre]*pe[mre-1])/ce[mre]);
        ee.push_back((ae[mre]*ee[mre-1]-de[mre]*pe[mre])/ce[mre]);
        tce=ee[mre]/(pe[mre]*pe[mre+1]);
        bve*=b2;
        mce+=1.;
        mre++;
    }
    iniNe=mre;
    double mcu=(double)iniN;
    int mru=(int)mcu;
    while(fabs(tcu)>coeffu && mru<maxr*iniL){
        au.push_back(-2.*mcu*(2.*mcu+1.));
        bu.push_back(k2-2.*vv*(4.*mcu+3.));
        cu.push_back(1./tv2is);
        //du.push_back(de[mru]*b2);
        du.push_back(bvu*Ee+bvu*(2.*mcu+1.)*Ee2);
        pu.push_back((bu[mru]*pu[mru]-au[mru]*pu[mru-1])/cu[mru]);
        eu.push_back((au[mru]*eu[mru-1]-du[mru]*pu[mru])/cu[mru]);
        tcu=eu[mru]/(pu[mru]*pu[mru+1]);
        mcu+=1.;
        mru++;
        bvu*=b2;
    }
    iniNu=mru;
    std::vector<double> Ce(iniNe,0.0);
    std::vector<double> Cu(iniNu,0.0);
    Ce[iniNe-1]=0.;//pow(bv, 2*(iniNe)+1)*C0/(3.*(double)(iniNe-1));
    Cu[iniNu-1]=0.;//pow(bv, 2*(iniNu)+1)*C1/(2*(double)iniNu);
  //  Ce[iniNe-2]=ee[iniNe-2]/pe[iniNe-1];
  //  Cu[iniNu-2]=eu[iniNu-2]/pu[iniNu-1];
    for(int m=iniNe-2;m>=0;m--){
        Ce[m]=(ee[m]+pe[m]*Ce[m+1])/pe[m+1];
    }
 //   Ce[0]=C0;
    for(int m=iniNu-3;m>=0;m--){
        Cu[m]=(eu[m]+pu[m]*Cu[m+1])/pu[m+1];
    }
 //   Cu[0]=C1;


    for(int m=1;m<iniL;m++){
        Csum+=Ce[m];
        Csume+=Ce[m];
        Csum+=Cu[m];
        Csumu+=Cu[m];
    }


    if(fabs(Csum)<compa*1.01 && fabs(Csum)>compa){
        Csum=compa;
    }
    if(((fabs(Csum)>compa)||std::isnan(Csum))){
        std::cout<<"Csump "<<Csum<<" C0 "<<C0<<" S0 "<<S0<<" C1 "<<C1<<" iniNe "<<iniNe <<" iniNu "<< iniNu <<" Csu "
        <<Csumu<<" Cse "<<Csume<<" iniL "<<iniL<<" C2 "<<Ce[1]<<" C3 "<<Cu[1]<< 
        " tv2i "<<tv2i<<" a  "<<av<<" compa  "<<compa<<" Hsor  "<<Hsor<<" ai  "<<ai<<" b "<<bv<<" v "<<vv<<" k "<<kv<<
        "\n";
        if(iniNe>6) std::cout<<"test e5 "<<ee[5]<<" e6 "<<ee[6]<<" p2 "<<pe[2]<<" p3 "<<pe[3]<<"\n";
    }
   ilimCI=ilimCIst;
   compa=compast;
   cb=cbst;
   sb=sbst;
   /* if(Csum<0.0){ std::cout<<"C0 "<<C0<<" C1 "<<C1<<" C2 "<<C2<<" C3 "<<C3<<" S0 "<<S0<<" S1 "<<S1<<" v "<<vv<<" k "<<kv<<" b "<<bv<<" S "<<Csum<<"\n";
    exit(1);
    }*/
    return Csum;
}
// here is Millers backward iteration for nonhomogenous
const int largeNM=10;
double gwenint::calCItrr5(double vv, double kv, double bv, double av){
  // cb=cos(kv*bv);
  // sb=sin(kv*bv);
   if(vv*av*av<=0.8){
       double testir= calCIiter(vv,kv,bv,av);
//       if(kv*bv<1.5 && testir<0.0) std::cout<<" "<<el0<<" "<<el1<<" "<<vv*av*av<<" "<<kv*av<<" "<<testir<<"\n";
       return testir;
   }
   double ilimCIst=ilimCI;
   double compast=compa;
   double cbst=cb;
   double sbst=sb;
   double bmax=sltomin/sqrt(vv);
   if(bv>bmax){ 
       bv=bmax;
       cb=cos(kv*bv);
       sb=sin(kv*bv);
       double dbr=bv/av;
       ilimCI=ceil(-gsl_sf_lambert_W0(log(dbr)/(log(1.-dbr)*reltol2))/log(dbr));
       compa=log(1./(1.-dbr));
       av=bv+aini;
   }
   if(bv<av*numtol || bv==av) return 0.;
    double C0=0.;
    double S0=0.;
    bv=bv/av;
    vv*=av*av;
    kv*=av;
    calC0S0(vv,kv,bv,&C0,&S0);
    if(C0==0.) return 0.;
 //   std::cout<<"shitte "<<C0<<" "<<S0<<"\n";
    double b2=bv*bv;
  //  double v2=v*v;
    double expf=exp(-vv*b2);
    //double sb=sin(k*b);
    //double cb=cos(k*b);
    double lsb=sb*expf;
    double lcb=cb*expf;
    double tv2i=1./(2.*vv);
    double C1=-(kv*S0+lcb-1.)*tv2i;
//    double S1=(kv*C0-lsb)*tv2i;

    double ai=1.;//a;
    double Csum=(C0+C1/ai)/ai;
    double k2=kv*kv;
   // ai*=a*a;
    double reli=1.;
    double Cn=C1;
    double Cnm=C0;
    double bs=bv*lsb;
    double bc=bv*lcb;


    double tv2is=-tv2i*tv2i;
    double Ee=expf*(cb*(2.*vv*b2)-bv*sb*kv);
    double Ee2=-expf*cb;
    int iniL=ceil(0.5*(double)ilimCI)-1;
    int iniN=std::max(150,iniL*2);//largeNM*ceil(ilimCI+1);

    std::vector<double> ae(iniN,0.);
    std::vector<double> be(iniN,0.);
    std::vector<double> ce(iniN,1./tv2is);
    std::vector<double> de(iniN,0.);
//    std::vector<double> ee(iniN,0.);
//    std::vector<double> pe(iniN+1,0.);
    std::vector<double> au(iniN,0.);
    std::vector<double> bu(iniN,0.);
    std::vector<double> cu(iniN,1./tv2is);
    std::vector<double> du(iniN,0.);
//    std::vector<double> eu(iniN,0.);
//    std::vector<double> pu(iniN+1,0.);
    de[0]=Ee/bv;
    du[0]=Ee+Ee2;
 //   pe[0]=0.;
  //  pe[1]=1.;
  //  pu[0]=0.;
  //  pu[1]=1.;
  //  ee[0]=C0;
   // eu[0]=C1;
    std::vector<double> ue(iniN,0.0);
    std::vector<double> uu(iniN,0.0);
    std::vector<double> ve(iniN,0.0);
    std::vector<double> vu(iniN,0.0);
    double Ceno=0;//C0*pow(bv,iniN)/((double) iniN);//0.;
    double Cenm=1;//C0*pow(bv,iniN-2)/((double) iniN);//1.; 
    double Cuno=0;//C1*pow(bv,iniN)/((double) iniN);//0.;
    double Cunm=1;//C1*pow(bv,iniN-2)/((double) iniN);//1.;
    uu[iniN-1]=Cuno;
    ue[iniN-1]=Ceno;
    uu[iniN-2]=Cunm;
    ue[iniN-2]=Cenm;

    double bve=bv;
    double bvu=b2;
    double Csume=C0;
    double Csumu=C1;
    for(int m=1;m<iniN;m++){
        double mr=(double)m;
        ae[m]=-2.*mr*(2.*mr-1.);
        au[m]=-2.*mr*(2.*mr+1.);
        be[m]=k2-2.*vv*(4.*mr+1.);
        bu[m]=k2-2.*vv*(4.*mr+3.);
        //de[m]=b2*de[m-1];
        //du[m]=b2*du[m-1];
        de[m]=bve*Ee+bve*2.*mr*Ee2;//b2*(de[m-1]-Ee2*(m-1));
        du[m]=bvu*Ee+bvu*(2.*mr+1.)*Ee2;
     //   pe[m+1]=(be[m]*pe[m]-ae[m]*pe[m-1])/ce[m];
     //   pu[m+1]=(bu[m]*pu[m]-au[m]*pu[m-1])/cu[m];
     //   ee[m]=(ae[m]*ee[m-1]-de[m]*pe[m])/ce[m];
     //   eu[m]=(au[m]*eu[m-1]-du[m]*pu[m])/cu[m];
     bve*=b2;
     bvu*=b2;
    }
    uu[iniN-3]=bu[iniN-2]/au[iniN-2];
    ue[iniN-3]=be[iniN-2]/ae[iniN-2];
    vu[iniN-3]=du[iniN-2]/au[iniN-2];
    ve[iniN-3]=de[iniN-2]/ae[iniN-2];
    int mm=iniN-3;
    if(std::isnan(ue[mm])||std::isnan(uu[mm])||std::isnan(ve[mm])||std::isnan(vu[mm])){
        std::cout<<ue[mm]<<" "<<uu[mm]<<" "<<ve[mm]<<" "<<vu[mm]<<" "<<mm<<"\n";
        exit(1);
    }
    for(int m=iniN-4;m>=0;m--){
        double mr=(double)m;
        ue[m]=(be[m+1]*ue[m+1]-ce[m+1]*ue[m+2])/ae[m+1];
        uu[m]=(bu[m+1]*uu[m+1]-cu[m+1]*uu[m+2])/au[m+1];
        vu[m]=(cu[m+1]*uu[m+1]+du[m+1]*uu[m+1])/au[m+1];
        ve[m]=(ce[m+1]*ue[m+1]+de[m+1]*ue[m+1])/ae[m+1];
        if(std::isnan(ue[m])||std::isnan(uu[m])||std::isnan(ve[m])||std::isnan(vu[m])){
            std::cout<<ue[m+2]<<" "<<uu[m+2]<<" "<<ve[m+2]<<" "<<vu[m+2]<<" "<<m<<"\n";
            exit(1);
        }

    }
    std::vector<double> Ce(iniN,0.0);
    std::vector<double> Cu(iniN,0.0);
    Ce[0]=C0;//pow(bv, 2*(iniNe)+1)*C0/(3.*(double)(iniNe-1));
    Cu[0]=C1;//pow(bv, 2*(iniNu)+1)*C1/(2*(double)iniNu);
  //  Ce[iniNe-2]=ee[iniNe-2]/pe[iniNe-1];
  //  Cu[iniNu-2]=eu[iniNu-2]/pu[iniNu-1];
    for(int m=1;m<iniN-1;m++){
        Ce[m]=(ue[m]*Ce[m-1]-ve[m-1])/ue[m-1];
        Cu[m]=(uu[m]*Cu[m-1]-vu[m-1])/uu[m-1];
        if(std::isnan(ue[m])||std::isnan(uu[m])||std::isnan(ve[m])||std::isnan(vu[m])){
            std::cout<<ue[m]<<" "<<uu[m]<<" "<<ve[m]<<" "<<vu[m]<<"\n";
            exit(1);
        }

    }


    for(int m=1;m<iniL;m++){
        Csum+=Ce[m];
        Csume+=Ce[m];
        Csum+=Cu[m];
        Csumu+=Cu[m];
    }


    if(fabs(Csum)<compa*1.01 && fabs(Csum)>compa){
        Csum=compa;
    }
    if(((fabs(Csum)>compa)||std::isnan(Csum))){
        std::cout<<"Csump "<<Csum<<" C0 "<<C0<<" S0 "<<S0<<" C1 "<<C1<<" Csu "
        <<Csumu<<" Cse "<<Csume<<" iniL "<<iniL<<" C2 "<<Ce[1]<<" C3 "<<Cu[1]<< 
        " tv2i "<<tv2i<<" a  "<<av<<" compa  "<<compa<<" Hsor  "<<Hsor<<" ai  "<<ai<<" b "<<bv<<" v "<<vv<<" k "<<kv<<
        "\n";
    }
   ilimCI=ilimCIst;
   compa=compast;
   cb=cbst;
   sb=sbst;
   /* if(Csum<0.0){ std::cout<<"C0 "<<C0<<" C1 "<<C1<<" C2 "<<C2<<" C3 "<<C3<<" S0 "<<S0<<" S1 "<<S1<<" v "<<vv<<" k "<<kv<<" b "<<bv<<" S "<<Csum<<"\n";
    exit(1);
    }*/
    return Csum;
}


const int intsubi=10000;
void gwenint::setq(double qo){
    q=qo;
    q2=q*q;
    ain=log(fabs(q-k));
    bin=log(q+k);
    tana=atan(ain);
    tanb=atan(bin)-tana;
    iabin=1./(bin-ain);
    dscaleq=dscale(q*L/tpi)/Hsor;
    dscaleq2=dscaleq*dscaleq;
    double ppp=k+q;
    double ppm=fabs(k-q);
    prsca=ppp+ppm;
    prbo=exp(-ppm/prsca);
    prao=exp(-ppp/prsca);
 //   gsl_integration_qawo_table_set(qawoit2,q/Hsor,a-aini,GSL_INTEG_COSINE);
   // if(dscaleq>10.){
  //  std::cout<<"ds "<<dscaleq<<" "<<q<<" "<<Hsor<<"\n";
  //  exit(1);
  //  }
}
void gwenint::setp(double po){
    p=po;
    p2=p*p;
}
double gwenint::getp(){
    return p;
}
double gwenint::getq(){
    return q;
}

const int intsubip=10000;
const int intsubiq=10000;
const double sonAlf=1./sqrt(3.);
const double sonAlf2=sonAlf*sonAlf;
gwenint::gwenint(double kn, double an, double Ln, enspec vsi, enspec bsi, enspec hvsi, enspec hbsi){
    a=an;
    k=kn;
    k2=k*k;
    L=Ln;
    wsp=gsl_integration_workspace_alloc(intsubi);
    wspq=gsl_integration_workspace_alloc(intsubiq);
    wspp=gsl_integration_workspace_alloc(intsubip);
    if(forcing)
    qawoit=gsl_integration_qawo_table_alloc(k/Hsor,a-aini-abuild,GSL_INTEG_COSINE,intsubi);
    else 
    qawoit=gsl_integration_qawo_table_alloc(k/Hsor,a-aini,GSL_INTEG_COSINE,intsubi);
  //  qawoit2=gsl_integration_qawo_table_alloc(k/Hsor,a-aini,GSL_INTEG_COSINE,intsubi);

    setup_Aplaw();
    setup_Aplawlag();
    vspec=vsi;
    bspec=bsi;
    hvspec=hvsi;
    hbspec=hbsi;
    L=vspec.get_L();
    kinfrac=vspec.getamp()/(vspec.getamp()+bspec.getamp());
    if(vspec.getamp()==0.) kinfrac=0.;
    magfrac=1.-kinfrac;
    if(bspec.getamp()==0.) magfrac=0.;
    soundf=kinfrac*sfevo(a);
    nsoundf=kinfrac-soundf+magfrac;
    nssoundf=sqrt(nsoundf);
    dscaleprefl=dscalepref/L*sqrt(std::max(vspec.getamp(),bspec.getamp()));
  //  if(a<aini+Hsor*tauLini) dscaleprefl=dscalepref/L*sqrt(Omvir+Ombir);
   // magsvel=std::min(sqrt((1./3.+2./3.*bspec.getamp())/(1.+2./3.*bspec.getamp())),1./sqrt(3.))/Hsor;
    magsvel=std::min(sqrt(1./3.+sonAlf2/(1.+sonAlf2)),1./sqrt(3.))/Hsor;
    dsca=initds;
    cas2=cos(k/Hsor*a)*(gsl_sf_Ci(k/Hsor*a)-gsl_sf_Ci(k/Hsor*(aini)))+sin(k/Hsor*a)*(gsl_sf_Si(k/Hsor*a)-gsl_sf_Si(k/Hsor*(aini)));
    qrsca=2.*pi/L;
    if(k>qrsca){
        qrsca=k;
    }

    double db=a-aini;
    double dbr=db/a;
    double kar=k*a/Hsor;
    cb=cos(k/Hsor*db);
    sb=sin(k/Hsor*db);
    as=sin(kar);
    ac=cos(kar);
    Cia=gsl_sf_Ci(kar);
    Sia=gsl_sf_Si(kar);
    Ciai=gsl_sf_Ci(k/Hsor*aini);
    Siai=gsl_sf_Si(k/Hsor*aini);

    double b2r=dbr*dbr;
    double b3r=b2r*dbr;
    double k2r=kar*kar;
    if(kar<1.0E-1){
        double lib=log(1.-dbr);
        double b4r=b3r*dbr;
        double b5r=b4r*dbr;
        double b6r=b5r*dbr;
        double b7r=b6r*dbr;
        double b8r=b7r*dbr;
        el0=-lib+0.25*(2.*dbr+b2r+2.*lib)*k2r;
        el1=-dbr-0.5*b2r-lib+k2r/24.*(12.*dbr+6.*b2r+4.*b3r+3.*b4r+12.*lib);
        el2=(-12.*dbr-6.*b2r-4.*b3r-3.*b4r-12.*lib)/12.+k2r/120.*(60.*dbr+30.*b2r+20.*b3r+15.*b4r+12.*b5r+10.*b6r+60.*lib);
        el3=(-60.*dbr-30.*b2r-20.*b3r-15.*b4r-12*b5r-10.*b6r-60.*lib)/60.+k2r/1680.*(840.*dbr+420.*b2r+280.*b3r+210.*b4r+168.*b5r+140.*b6r+120.*b7r+105.*b8r+840.*lib);
        el4=(-840.*dbr-420.*b2r-280.*b3r-210.*b4r-168*b5r-140.*b6r-120.*b7r-105.*b8r-840.*lib)/840.+k2r/5040.*(2520.*dbr+1620.*b2r+840.*b3r+630.*b4r+504.*b5r+420.*b6r+360.*b7r+315.*b8r+280.*dbr*b8r+252.*b8r*b2r+2520.*lib);
       // el1=el0-dbr*(1.+0.5*dbr)+b2r*k2r/24.*(4.*dbr+3.*b2r);
       // el2=el0-dbr*(1.+0.5*dbr+b2r/3.+0.25*b3r)+k2r*b3r/120.*(20.+15.*dbr+12.*b2r+10.*b3r);
    } else {
        double k4r=k2r*k2r;
        double k6r=k4r*k2r;
        double k8r=k6r*k2r;
        double b4r=b3r*dbr;
        double b5r=b4r*dbr;
        double b6r=b5r*dbr;
        double b7r=b6r*dbr;
        el0=ac*(Cia-Ciai)+as*(Sia-Siai);
        el1=el0+((1.-cb)/kar-(1.+dbr)*sb)/kar;
        el2=el0+((1.-cb)*(k2r-6.)-cb*k2r*(2.*dbr+3.*(b2r))-sb*kar*(k2r*(1.+b2r+b3r+dbr)-2.-6.*dbr))/k4r;
        el3=el0+((1.-cb)*(120-6.*k2r+k4r)+cb*(6.*k2r*(4.*dbr+10.*b2r)-k4r*(2.*dbr+3.*b2r+4.*b3r+5.*b4r))-kar*sb*(24.+120*dbr+k2r*(-2.-20.*b3r-12.*b2r-6.*dbr)+k4r*(b4r+b5r+b3r+b2r+dbr)))/k6r;
        el4=el0+((1.-cb)*(-5040.+120.*k2r-6*k4r+k6r)-cb*(120.*k2r*(6.*dbr+21.*b2r)-6.*k4r*(4.*dbr+10.*b2r+20.*b3r+35*b4r)+k6r*(2.*dbr+3.*b2r+4.*b3r+5.*b4r+6.*b5r+7.*b6r))
        -kar*sb*(-720.-5040.*dbr+k2r*(24.+840*b3r+360.*b2r+120*dbr)+k4r*(-2.-42.*b5r-30.*b4r-20.*b3r-12.*b2r-6.*dbr)+k6r*(1.+b6r+b7r+b5r+b4r+b3r+b2r+dbr)))/k8r;
       // el2=el0+(-6+k2r-cb*(6.-k2r*(1.+2.*dbr+3.*b2r))-sb*kar*(-2.+k2r*(1.+b2r*(1.+dbr)+dbr)-6.*dbr))/(k2r*k2r);
    }
    ilimCI=2*ceil(-gsl_sf_lambert_W0(log(dbr)/(log(1.-dbr)*reltol2))/log(dbr));
    compa=log(1./(1.-dbr));






    //std::cout<<dscaleprefl<<" "<<dscalepref<<" "<<L<<" "<<vspec.getamp()<<" "<<bspec.getamp()<<" "<<Ln<<
    //" "<<Omartevoc(a)<<" "<<Levoa(a)<<"\n";
}

void gwenint::setkat(double kn, double an){
    a=an;
    k=kn;
    k2=k*k;
    gsl_integration_qawo_table_free(qawoit);

    vspec.set_t(an);
    bspec.set_t(an);
    hvspec.set_t(an);
    hbspec.set_t(an);
    L=vspec.get_L();
    kinfrac=vspec.getamp()/(vspec.getamp()+bspec.getamp());
    magfrac=1.-kinfrac;
    if(vspec.getamp()==0.) kinfrac=0.;
    if(bspec.getamp()==0.) magfrac=0.;
    soundf=kinfrac*sfevo(a);
    nsoundf=kinfrac-soundf+magfrac;
    nssoundf=sqrt(nsoundf);
    dscaleprefl=dscalepref/L*sqrt(std::max(vspec.getamp(),bspec.getamp()));
  //  if(a<aini+Hsor*tauLini) dscaleprefl=dscalepref/L*sqrt(Omvir+Ombir);
    double dscalek=dscale(k*L/tpi);
    double dscalekp=dscale(0.5*k*L/tpi)/Hsor;
    double rarg=4.;
    double delra=a-aini;
    if(rarg/dscalekp<a-aini) delra=rarg/dscalekp;
    if(forcing){
        if(rarg/dscalekp<a-aini-abuild) delra=rarg/dscalekp;
    }

    if(forcing)
    qawoit=gsl_integration_qawo_table_alloc(k/Hsor,delra,GSL_INTEG_COSINE,intsubi);
    else 
    qawoit=gsl_integration_qawo_table_alloc(k/Hsor,delra,GSL_INTEG_COSINE,intsubi);
    dscalek2=dscalek*dscalek/(Hsor*Hsor);
    cas2=cos(k/Hsor*a)*(gsl_sf_Ci(k/Hsor*a)-gsl_sf_Ci(k/Hsor*(aini)))+sin(k/Hsor*a)*(gsl_sf_Si(k/Hsor*a)-gsl_sf_Si(k/Hsor*(aini)));
    qrsca=2.*pi/L;
    if(k>qrsca){
        qrsca=k;
    }
    double db=a-aini;
    double dbr=db/a;
    double kar=k*a/Hsor;
    cb=cos(k/Hsor*db);
    sb=sin(k/Hsor*db);
    as=sin(kar);
    ac=cos(kar);
    Cia=gsl_sf_Ci(kar);
    Sia=gsl_sf_Si(kar);
    Ciai=gsl_sf_Ci(k/Hsor*aini);
    Siai=gsl_sf_Si(k/Hsor*aini);
    el0=ac*(Cia-Ciai)+as*(Sia-Siai);
    double b2r=dbr*dbr;
    double b3r=b2r*dbr;
    double k2r=kar*kar;
    if(kar<1.0E-1){
        double lib=log(1.-dbr);
        double b4r=b3r*dbr;
        double b5r=b4r*dbr;
        double b6r=b5r*dbr;
        double b7r=b6r*dbr;
        double b8r=b7r*dbr;
        el0=-lib+0.25*(2.*dbr+b2r+2.*lib)*k2r;
        el1=-dbr-0.5*b2r-lib+k2r/24.*(12.*dbr+6.*b2r+4.*b3r+3.*b4r+12.*lib);
        el2=(-12.*dbr-6.*b2r-4.*b3r-3.*b4r-12.*lib)/12.+k2r/120.*(60.*dbr+30.*b2r+20.*b3r+15.*b4r+12.*b5r+10.*b6r+60.*lib);
        el3=(-60.*dbr-30.*b2r-20.*b3r-15.*b4r-12*b5r-10.*b6r-60.*lib)/60.+k2r/1680.*(840.*dbr+420.*b2r+280.*b3r+210.*b4r+168.*b5r+140.*b6r+120.*b7r+105.*b8r+840.*lib);
        el4=(-840.*dbr-420.*b2r-280.*b3r-210.*b4r-168*b5r-140.*b6r-120.*b7r-105.*b8r-840.*lib)/840.+k2r/5040.*(2520.*dbr+1620.*b2r+840.*b3r+630.*b4r+504.*b5r+420.*b6r+360.*b7r+315.*b8r+280.*dbr*b8r+252.*b8r*b2r+2520.*lib);
       // el1=el0-dbr*(1.+0.5*dbr)+b2r*k2r/24.*(4.*dbr+3.*b2r);
       // el2=el0-dbr*(1.+0.5*dbr+b2r/3.+0.25*b3r)+k2r*b3r/120.*(20.+15.*dbr+12.*b2r+10.*b3r);
    } else {
        double k4r=k2r*k2r;
        double k6r=k4r*k2r;
        double k8r=k6r*k2r;
        double b4r=b3r*dbr;
        double b5r=b4r*dbr;
        double b6r=b5r*dbr;
        double b7r=b6r*dbr;
        el0=ac*(Cia-Ciai)+as*(Sia-Siai);
        el1=el0+((1.-cb)/kar-(1.+dbr)*sb)/kar;
        el2=el0+((1.-cb)*(k2r-6.)-cb*k2r*(2.*dbr+3.*(b2r))-sb*kar*(k2r*(1.+b2r+b3r+dbr)-2.-6.*dbr))/k4r;
        el3=el0+((1.-cb)*(120-6.*k2r+k4r)+cb*(6.*k2r*(4.*dbr+10.*b2r)-k4r*(2.*dbr+3.*b2r+4.*b3r+5.*b4r))-kar*sb*(24.+120*dbr+k2r*(-2.-20.*b3r-12.*b2r-6.*dbr)+k4r*(b4r+b5r+b3r+b2r+dbr)))/k6r;
        el4=el0+((1.-cb)*(-5040.+120.*k2r-6*k4r+k6r)-cb*(120.*k2r*(6.*dbr+21.*b2r)-6.*k4r*(4.*dbr+10.*b2r+20.*b3r+35*b4r)+k6r*(2.*dbr+3.*b2r+4.*b3r+5.*b4r+6.*b5r+7.*b6r))
        -kar*sb*(-720.-5040.*dbr+k2r*(24.+840*b3r+360.*b2r+120*dbr)+k4r*(-2.-42.*b5r-30.*b4r-20.*b3r-12.*b2r-6.*dbr)+k6r*(1.+b6r+b7r+b5r+b4r+b3r+b2r+dbr)))/k8r;
       // el2=el0+(-6+k2r-cb*(6.-k2r*(1.+2.*dbr+3.*b2r))-sb*kar*(-2.+k2r*(1.+b2r*(1.+dbr)+dbr)-6.*dbr))/(k2r*k2r);
    }
    ilimCI=2*ceil(-gsl_sf_lambert_W0(log(dbr)/(log(1.-dbr)*reltol2))/log(dbr));
    compa=log(1./(1.-dbr));
}

void gwenint::setat(double an){
    a=an;
    gsl_integration_qawo_table_free(qawoit);
    //if(forcing)
    //qawoit=gsl_integration_qawo_table_alloc(k/Hsor,a-aini-abuild,GSL_INTEG_COSINE,intsubi);
    //else 
    //qawoit=gsl_integration_qawo_table_alloc(k/Hsor,a-aini,GSL_INTEG_COSINE,intsubi);
    vspec.set_t(an);
    bspec.set_t(an);
    hvspec.set_t(an);
    hbspec.set_t(an);
    L=vspec.get_L();
    kinfrac=vspec.getamp()/(vspec.getamp()+bspec.getamp());
    magfrac=1.-kinfrac;
    if(vspec.getamp()==0.) kinfrac=0.;
    if(bspec.getamp()==0.) magfrac=0.;
    soundf=kinfrac*sfevo(a);
    nsoundf=kinfrac-soundf+magfrac;
    nssoundf=sqrt(nsoundf);
    dscaleprefl=dscalepref/L*sqrt(std::max(vspec.getamp(),bspec.getamp()));
  //  if(a<aini+Hsor*tauLini) dscaleprefl=dscalepref/L*sqrt(Omvir+Ombir);
    cas2=cos(k/Hsor*a)*(gsl_sf_Ci(k/Hsor*a)-gsl_sf_Ci(k/Hsor*(aini)))+sin(k/Hsor*a)*(gsl_sf_Si(k/Hsor*a)-gsl_sf_Si(k/Hsor*(aini)));
//    std::cout<<dscaleprefl<<" "<<vspec.getamp()<<" "<<an<<" "<<tauLini*Hsor<<"\n";
    qrsca=2.*pi/L;
    double dscalekp=dscale(0.5*k*L/tpi)/Hsor;
    double rarg=4.;
    double delra=a-aini;
    if(rarg/dscalekp<a-aini) delra=rarg/dscalekp;
    if(forcing){
        if(rarg/dscalekp<a-aini-abuild) delra=rarg/dscalekp;
    }

    if(forcing)
    qawoit=gsl_integration_qawo_table_alloc(k/Hsor,delra,GSL_INTEG_COSINE,intsubi);
    else 
    qawoit=gsl_integration_qawo_table_alloc(k/Hsor,delra,GSL_INTEG_COSINE,intsubi);
    if(k>qrsca){
        qrsca=k;
    }
    double db=a-aini;
    double dbr=db/a;
    double kar=k*a/Hsor;
    cb=cos(k/Hsor*db);
    sb=sin(k/Hsor*db);
    as=sin(kar);
    ac=cos(kar);
    Cia=gsl_sf_Ci(kar);
    Sia=gsl_sf_Si(kar);
    Ciai=gsl_sf_Ci(k/Hsor*aini);
    Siai=gsl_sf_Si(k/Hsor*aini);
    el0=ac*(Cia-Ciai)+as*(Sia-Siai);
    double b2r=dbr*dbr;
    double b3r=b2r*dbr;
    double k2r=kar*kar;
    if(kar<2.0E-1){
        double lib=log(1.-dbr);
        double b4r=b3r*dbr;
        double b5r=b4r*dbr;
        double b6r=b5r*dbr;
        double b7r=b6r*dbr;
        double b8r=b7r*dbr;
        el0=-lib+0.25*(2.*dbr+b2r+2.*lib)*k2r;
        el1=-dbr-0.5*b2r-lib+k2r/24.*(12.*dbr+6.*b2r+4.*b3r+3.*b4r+12.*lib);
        el2=(-12.*dbr-6.*b2r-4.*b3r-3.*b4r-12.*lib)/12.+k2r/120.*(60.*dbr+30.*b2r+20.*b3r+15.*b4r+12.*b5r+10.*b6r+60.*lib);
        el3=(-60.*dbr-30.*b2r-20.*b3r-15.*b4r-12*b5r-10.*b6r-60.*lib)/60.+k2r/1680.*(840.*dbr+420.*b2r+280.*b3r+210.*b4r+168.*b5r+140.*b6r+120.*b7r+105.*b8r+840.*lib);
        el4=(-840.*dbr-420.*b2r-280.*b3r-210.*b4r-168*b5r-140.*b6r-120.*b7r-105.*b8r-840.*lib)/840.+k2r/5040.*(2520.*dbr+1620.*b2r+840.*b3r+630.*b4r+504.*b5r+420.*b6r+360.*b7r+315.*b8r+280.*dbr*b8r+252.*b8r*b2r+2520.*lib);
       // el1=el0-dbr*(1.+0.5*dbr)+b2r*k2r/24.*(4.*dbr+3.*b2r);
       // el2=el0-dbr*(1.+0.5*dbr+b2r/3.+0.25*b3r)+k2r*b3r/120.*(20.+15.*dbr+12.*b2r+10.*b3r);
    } else {
        double k4r=k2r*k2r;
        double k6r=k4r*k2r;
        double k8r=k6r*k2r;
        double b4r=b3r*dbr;
        double b5r=b4r*dbr;
        double b6r=b5r*dbr;
        double b7r=b6r*dbr;
        el0=ac*(Cia-Ciai)+as*(Sia-Siai);
        el1=el0+((1.-cb)/kar-(1.+dbr)*sb)/kar;
        el2=el0+((1.-cb)*(k2r-6.)-cb*k2r*(2.*dbr+3.*(b2r))-sb*kar*(k2r*(1.+b2r+b3r+dbr)-2.-6.*dbr))/k4r;
        el3=el0+((1.-cb)*(120-6.*k2r+k4r)+cb*(6.*k2r*(4.*dbr+10.*b2r)-k4r*(2.*dbr+3.*b2r+4.*b3r+5.*b4r))-kar*sb*(24.+120*dbr+k2r*(-2.-20.*b3r-12.*b2r-6.*dbr)+k4r*(b4r+b5r+b3r+b2r+dbr)))/k6r;
        el4=el0+((1.-cb)*(-5040.+120.*k2r-6*k4r+k6r)-cb*(120.*k2r*(6.*dbr+21.*b2r)-6.*k4r*(4.*dbr+10.*b2r+20.*b3r+35*b4r)+k6r*(2.*dbr+3.*b2r+4.*b3r+5.*b4r+6.*b5r+7.*b6r))
        -kar*sb*(-720.-5040.*dbr+k2r*(24.+840*b3r+360.*b2r+120*dbr)+k4r*(-2.-42.*b5r-30.*b4r-20.*b3r-12.*b2r-6.*dbr)+k6r*(1.+b6r+b7r+b5r+b4r+b3r+b2r+dbr)))/k8r;
       // el2=el0+(-6+k2r-cb*(6.-k2r*(1.+2.*dbr+3.*b2r))-sb*kar*(-2.+k2r*(1.+b2r*(1.+dbr)+dbr)-6.*dbr))/(k2r*k2r);
    }
    ilimCI=2*ceil(-gsl_sf_lambert_W0(log(dbr)/(log(1.-dbr)*reltol2))/log(dbr));
    compa=log(1./(1.-dbr));
//    std::cout<<dscaleprefl<<" "<<dscalepref<<" "<<L<<" "<<vspec.getamp()<<" "<<bspec.getamp()<<" "<<
  //  " "<<Omartevoc(a)<<" "<<Levoa(a)<<"\n";
}
double functap(double K){
    double K2=K*K;
    if(lagrangian){
        return K2*K2*K2/pow(Lc+K2,Kolpot);
    }else {//corresponds to int k^2 E(k) dk
        return K2*K2/pow(Lc+K2,Kolpot); //corresponds to int E(k) dk
    }  
}
double functaplag(double K){
    double K2=K*K;
    return K2*K2*K2/pow(Lc+K2,Kolpot);
}

void gwenint::setup_Aplaw(){
    biggrid=gridinb();

    yongrid.resize(Npsi);
    A.resize(Npsi-1);
    plaw.resize(Npsi-1);
    subintval.resize(Npsi);
    for(int i=0;i<Npsi;i++){
        yongrid[i]=functap(biggrid[i]);
    }
    for(int i=0;i<Npsi-1;i++){
        plaw[i]=log(yongrid[i+1]/yongrid[i])/log(biggrid[i+1]/biggrid[i]);
        A[i]=yongrid[i]*pow(biggrid[i],-plaw[i]);
    }
    Aplawosc=yongrid[Npsi-1]*pow(biggrid[Npsi-1],-infind);
    subintval[0]=A[0]*pow(biggrid[0],plaw[0]+1.)/(plaw[0]+1.);//lowbcoeff*pow(biggrid[0],initpow);
    for(int i=1;i<Npsi;i++){
        subintval[i]=subintval[i-1]+A[i-1]/(plaw[i-1]+1.0)*(pow(biggrid[i],plaw[i-1]+1.0)-pow(biggrid[i-1],plaw[i-1]+1.0));
      //  std::cout<<subintval[i]<<"\n";
    }
 //   std::cout<<subintval[Npsi-1]<<"\n";
}

void gwenint::setup_Aplawlag(){
   // biggrid=gridinb();

    yongridlag.resize(Npsi);
    Alag.resize(Npsi-1);
    plawlag.resize(Npsi-1);
    subintvallag.resize(Npsi);
    for(int i=0;i<Npsi;i++){
        yongridlag[i]=functaplag(biggrid[i]);
    }
    for(int i=0;i<Npsi-1;i++){
        plawlag[i]=log(yongridlag[i+1]/yongridlag[i])/log(biggrid[i+1]/biggrid[i]);
        Alag[i]=yongridlag[i]*pow(biggrid[i],-plawlag[i]);
    }
    Aplawosclag=yongridlag[Npsi-1]*pow(biggrid[Npsi-1],-infindlag);
    subintvallag[0]=Alag[0]*pow(biggrid[0],plawlag[0]+1.)/(plawlag[0]+1.);//lowbcoeff*pow(biggrid[0],initpow);
    for(int i=1;i<Npsi;i++){
        subintvallag[i]=subintvallag[i-1]+Alag[i-1]/(plawlag[i-1]+1.0)*(pow(biggrid[i],plawlag[i-1]+1.0)-pow(biggrid[i-1],plawlag[i-1]+1.0));
      //  std::cout<<subintval[i]<<"\n";
    }
 //   std::cout<<subintval[Npsi-1]<<"\n";
}
double gwenint::dscale(double Q){
    double val=0.0;
    double Qini=Q;
    if(Q>L/dsca) Q=L/dsca;
    // dscaleprefl at the end
    if(Q<=ls){
        val= (subintval[0]); //sqrt(lowbcoeff)*pow(biggrid[0],initpow*0.5);
    //    if(val<1.0E-100)
         //   std::cout<<"dsca "<<val<<" "<<Q<<" "<<lowbcoeff<<" "<<biggrid[0]<<" "<<initpow<<" "<<val<<"\n";
    } else if(Q>=us){
        val= (subintval[Npsi-1]+Aplawosc*(pow(Q,infind)-pow(biggrid[Npsi-1],infind)));
      //  if(val<1.0E-100)
      //      std::cout<<"dsce "<<val<<" "<<Q<<" "<<subintval[Npsi-1]<<" "<<Aplawosc<<" "<<infind<<"\n";
    } else {
        std::vector<double>::const_iterator it=std::upper_bound(biggrid.begin(),biggrid.end(),Q);
        size_t i=it-biggrid.begin()-1;
        val= (subintval[i]+A[i]/(plaw[i]+1.0)*(pow(Q,plaw[i]+1.)-pow(biggrid[i],plaw[i]+1.)));
      //  if(val>1000.){
        //    std::cout<<"dsci "<<val<<" "<<Q<<" "<<L<<" "<<dsca<<" "<<subintval[i]<<" "<<A[i]<<" "<<plaw[i]<<" "<<biggrid[i]<<
          //  " "<<i<<" "<<biggrid[Npsi-1]<<" "<<biggrid[0]<<"\n";
         //   exit(1);
       // }
    }
  //  std::cout<<val<<"\n";
  if(!lagrangian){
      if(puresweep){
          //factor sqrt(2./3.) kept out
          val=Qini*dscaleprefl*sqrt(subintval[Npsi-1]);
      }else {
          val=Qini*dscaleprefl*(sqrt((val*(1.5+Qini*Qini*0.2)+subintval[Npsi-1])/(2.5+0.2*Qini*Qini)));
      }
  } else {
   val*=dscaleprefl*sqrt(val);   
  }// corresponds to  Q^2* int E(k) dk
  //  std::cout<<"val "<<val/Qini<<" "<<dscaleprefl<<" "<<L<<"\n";
   // exit(1);
    return val;
}

double gwenint::dscalelag(double Q){
    double val=0.0;
    double Qini=Q;
    double preval=0.;
    if(Q>L/dsca) Q=L/dsca;
    // dscaleprefl at the end
    if(Q<=ls){
        preval= (subintvallag[0]); //sqrt(lowbcoeff)*pow(biggrid[0],initpow*0.5);
    //    if(val<1.0E-100)
         //   std::cout<<"dsca "<<val<<" "<<Q<<" "<<lowbcoeff<<" "<<biggrid[0]<<" "<<initpow<<" "<<val<<"\n";
    } else if(Q>=us){
        preval= (subintvallag[Npsi-1]+Aplawosclag*(pow(Q,infindlag)-pow(biggrid[Npsi-1],infindlag)));
      //  if(val<1.0E-100)
      //      std::cout<<"dsce "<<val<<" "<<Q<<" "<<subintval[Npsi-1]<<" "<<Aplawosc<<" "<<infind<<"\n";
    } else {
        std::vector<double>::const_iterator it=std::upper_bound(biggrid.begin(),biggrid.end(),Q);
        size_t i=it-biggrid.begin()-1;
        preval= (subintvallag[i]+Alag[i]/(plawlag[i]+1.0)*(pow(Q,plawlag[i]+1.)-pow(biggrid[i],plawlag[i]+1.)));
      //  if(val>1000.){
        //    std::cout<<"dsci "<<val<<" "<<Q<<" "<<L<<" "<<dsca<<" "<<subintval[i]<<" "<<A[i]<<" "<<plaw[i]<<" "<<biggrid[i]<<
          //  " "<<i<<" "<<biggrid[Npsi-1]<<" "<<biggrid[0]<<"\n";
         //   exit(1);
       // }
    }
  //  std::cout<<val<<"\n";
   val=sqrt(preval)*dscaleprefl;  
  // if(val<=0.0) std::cout<<"shit "<<preval<<"\n";
  // val=sqrt(val)*dscaleprefl*sqrt(2/3);  
  //  std::cout<<"val "<<val/Qini<<" "<<dscaleprefl<<" "<<L<<"\n";
   // exit(1);
    return val;
}
const double relvel=1./sqrt(3.)/Hsor;//should be 1./sqrt(3.)
double gwenint::drate(double ap){
   // double tr=0.0;
   // if(a-ap)
    switch(dratech){
        case 0:
         //   return rate0(ap,dscale(q*L/(tpi)),dscale(p*L/(tpi)))/(a-ap);//*cos(ap*k/Hsor);
           { double adap=a-ap;
            double corrfea= 1.;
            double comprf=1.;
            double forcer=1.;
            double comprfq=soundf;
            double comprfp=soundf;
            if(compr) {
                comprfq*=cos(q*relvel*ap);
                comprfp*=cos(p*relvel*ap);
                comprf=(coefm*4.*comprfq*comprfp*veq2*vep2+3.*nsoundf*(coefmq*comprfq*veq2*(vep)+coefmp*comprfp*vep2*(veq)))+(nsoundf*nsoundf*(veq*vep+beq*bep)+Ctf);
            }
            
            if(magson){
                comprfq*=cos(q*magsvel*ap);
                comprfp*=cos(p*magsvel*ap);
                comprf=(coefm*4.*comprfq*comprfp*veq2*vep2+6.*nsoundf*(coefmq*comprfq*veq2*(vep)+coefmp*comprfp*vep2*(veq)))+(nsoundf*nsoundf*(veq*vep+beq+bep)+Ctf);
            }
            if(forcing2 && adap<aini+abuild){
                if(a>=aini+abuild) forcer=(adap-aini)*(adap-aini)/(abuild*abuild);
                else forcer=(adap-aini)*(adap-aini)/((a-aini)*(a-aini));           
           //     if(a>=aini+abuild) forcer=(adap-aini)/(abuild);
            //    else forcer=(adap-aini)/((a-aini));           
            }

            if(simplyi) corrfea=Omartevoc(adap)*Levoa(adap)/(Omartevoc(a)*Levoa(a));//*cos(ap*k/Hsor);
            if(altrate){
                return corrfea*rate5(ap,nssoundf*dscaleq,nssoundf*dscalepg)/(adap)*comprf*forcer;//*cos(ap*k/Hsor);
            } else{
                return corrfea*rate0(ap,nssoundf*dscaleq,nssoundf*dscalepg)/(adap)*comprf*forcer;//*cos(ap*k/Hsor);
            }
           }
          //  std::cout<<tr<<" "<<a-ap<<" "<<q<<" "<<p<<"\n";
         //   return tr;
            break;
        default:

            return 0.0;
            break;
    }
}
double gwenint::drates2(double ap){
   // double tr=0.0;
   // if(a-ap)
    switch(dratech){
        case 3:
         //   return rate0(ap,dscale(q*L/(tpi)),dscale(p*L/(tpi)))/(a-ap);//*cos(ap*k/Hsor);
         {   double adap=a-ap;
            return Omartevoc(adap)*Levoa(adap)/(Omartevoc(a)*Levoa(a))/(adap);//*cos(ap*k/Hsor);
          //  std::cout<<tr<<" "<<a-ap<<" "<<q<<" "<<p<<"\n";
         }
         //   return tr;
            break;
        default:
            return 0.0;
            break;
    }
}

double dratei(double ap, void *data){
    gwenint *gwenintfi=static_cast<gwenint *>(data);
    return gwenintfi->drate(ap);
}


double dratei2(double ap, void *data){
    gwenint *gwenintfi=static_cast<gwenint *>(data);
    return gwenintfi->drates2(ap);
}
const double paral=1.;
const double paral2=paral*paral;
    const double p1=1./sqrt(3.)+1.;
    const double p2=1.-1./sqrt(3.);
double gwenint::drateint(){
   /* 
    double result=0.0;
    double abserr=0.0;
    double epsabi1=errafa*1.0E-2;
    const double epsreli1=1.0E-3;
    gsl_function F;
    F.function=&dratei;
    F.params=this;
    gsl_integration_qawo(&F,0,epsabi1,epsreli1,intsubi,wsp,qawoit,&result,&abserr);
    //std::cout<<result<<"\n";
    return result;
    */
  //  double dscaleq=dscale(q*L/tpi)/Hsor;
  if(a<=abuild+aini && forcing){
    double adif=a-aini;
    double kH=k/Hsor;
    double kadif=kH*adif;
    if(kadif<0.3) return adif/3.;

    double kadif1=kadif*p1;
    double kadif2=kadif*p2;

        //return 0.5*adif/(aini+0.5*abuild);
 //   return (1.-cos(k*adif))/(k*k*adif*(aini+0.5*abuild));
    return 1./(kH*p1*kadif1*kadif1)*(kadif1-sin(kadif1))+1./(kH*p2*kadif2*kadif2)*(kadif2-sin(kadif2));
  }
    switch(dratech){
        case 0:
            {
                double dscalep=dscale(p*L/tpi)/Hsor;
                double dscalelagp=dscalelag(p*L/tpi)/Hsor;
                double dscalelagq=dscalelag(q*L/tpi)/Hsor;
                
                double dscaleint=dscale(1.)/Hsor;
                double dscalelagint=dscalelag(1.)/Hsor;
          //      if(a>aini+Hsor*tauLini){
                //    std::cout<<dscaleint*a<<"\n";
                //std::cout<<sqrt(vspec.getamp())/L*a<<"\n";
         //       std::cout<<sqrt(Omartevoc(a))/Levoa(a)*(a-aini)<<"\n";
              //  std::cout<<gamma2*0.5+gamma1<<"\n";
               
            //    }
                dscalepg=dscalep;
          //      double dscale;

                double val=0.;
               // double ql2=q*q;
               // double pl2=p*p;
                //double nfact=ql2*pl2/(ql2+pl2);
                //if(Hsor*Hsor*(a-aini)*(a-aini)<nfact) return val;
              //  double nfact=q*p*k/(q*k+p*k+q*p);
              //  if(Hsor*(a-aini)<nfact) return val;
              //  if(std::min(p,q)<(a-aini)/Hsor) return val;
           //   if(a-aini<Hsor/k) return val;
          // if(paral*std::max(sqrt(dscaleq2),dscalep)<(a-aini)/Hsor) return val;
          
          if(thetcut && (paral*dscaleq<dscalefac/(a) &&  paral*dscalep<dscalefac/(a)/* ||
              paral*dscaleint<dscalefac/(a) || a<aini+Hsor*tauLini*/)) return val;
   //        if(thetcut && (paral*dscalelagq<dscalefac/(a) || paral*dscalelagp<dscalefac/(a) ||
    //           paral*dscalelagint<dscalefac/(a) || a<aini+Hsor*tauLini)) return val;
       //    if(thetcut && (paral*(dscalelagq+dscalelagp)<dscalefac/(a)/* || paral*dscalelagp<dscalefac/(a) ||
         //      paral*dscalelagint<dscalefac/(a) || a<aini+Hsor*tauLini*/)) return val;
                double vstot=paral2*(dscaleq2/*+dscalek2*/+dscalep*dscalep);

        //  if(thetcut && (paral*dscaleq<dscalefac/(a) ||  paral*dscalep<dscalefac/(a)/* ||
            //  paral*dscaleint<dscalefac/(a) || a<aini+Hsor*tauLini*/)){
          //  if(paral*dscaleq<dscalefac/a) dscaleq=dscalefac/a;
          //  if(paral*dscalep<dscalefac/a) dscalep=dscalefac/a;
         // };
                //vstot+=(Hsor*Hsor)*(1./(q*q)+1./(p*p));
            //    vstot+=k*k/(Hsor*Hsor);
       //         val=calCI(vstot,k/Hsor,a-aini,a);





            double result=0.0;
            double abserr=arr;
            double result2=0.0;
            double abserr2=arr;
            double epsabi1=errafa;//*1.0E-10;
            const double epsreli1=5.0E-3;
            gsl_function F;
            F.function=&dratei;
            F.params=this;
            double sepa=3.5;
            double cala=sepa/sqrt(vstot);

            int cher=2;
            double val2=0.0;
            if(a-aini>cala && cala/a<=0.15 && !lagrangian && appint){
                val=calCI(0.5*vstot,k/Hsor,cala,a);
            //    val2=calCI(vstot,q/Hsor,cala,a);
                cher=0;
                val=(val2+val);
            } else if(1.-aini/a<=0.15 && !lagrangian && appint){
                val=calCI(0.5*vstot,k/Hsor,a-aini,a);
             //   val2=calCI(vstot,q/Hsor,cala,a);
                cher=1;
                val=(val2+val);
            } else {
             /*   if((a-aini)>cala){
                   gsl_integration_qawo_table_set_length(qawoit,cala);
                }*/
                gsl_integration_qawo(&F,0,epsabi1,epsreli1,intsubi,wsp,qawoit,&result,&abserr);
            //    gsl_integration_qawo(&F,0,epsabi1,epsreli1,intsubi,wsp,qawoit2,&result2,&abserr2);
if(forcing)    {
    double kH=k/Hsor;
    result+=0.5*drate(a-(aini+abuild))/(kH*kH*kH*p1*p1*p1*abuild*abuild)*((kH*kH*p1*p1*abuild*abuild-2.)*sin(kH*p1*(aini+abuild-a))+2.*kH*p1*abuild*cos(kH*p1*(aini+abuild-a))+2.*sin(kH*p1*(aini-a)));
    result+=0.5*drate(a-(aini+abuild))/(kH*kH*kH*p2*p2*p2*abuild*abuild)*((kH*kH*p2*p2*abuild*abuild-2.)*sin(kH*p2*(aini+abuild-a))+2.*kH*p2*abuild*cos(kH*p2*(aini+abuild-a))+2.*sin(kH*p2*(aini-a)));
}
    //drate(a-(aini+abuild))/(k*k*abuild)*(cos(k*(a-(aini+abuild)))-cos(k*(a-aini))-k*abuild*sin(k*(a-(aini+abuild))))/(aini+0.5*abuild);
                val=result;
              //  val=(val+result2);
//                val=calCI(vstot,k/Hsor,a-aini,a);
  //              cher=2;
              /*  if((a-aini)>cala){
                   gsl_integration_qawo_table_set_length(qawoit,a-aini);
                    
                }*/

            }
          //  if(fabs(val)>log(1./(1.-(a-aini)/a))) std::cout<<"gotcha \n";
          /*  if(val<0){
                std::cout<<"problem here : "<<val<<" "<<vstot*(a*a)<<" "<<k*(a/Hsor)<<" "<<1.-aini/a<<
                " "<<cher<<" ds "<<dscaleprefl<<" Li  "<<tpi/L<<" amp "<<vspec.getamp()<<" Hsor "<<Hsor<<" a "<<a<<" edd "<<tauLini*Hsor<<
		" k "<<k<<" q "<<q<<ts enthaltene p "<<p<<"\n";
                exit(1);
            }*/





            
      /*          if(val<0.){
                    std::cout<<"vpr "<<vstot*a*a<<" k "<<k/Hsor*a<<" b "<<(a-aini)/a<<" a "<<a<<" val "<<val<<"\n";
                    exit(1);
                }*/
             //   if(fabs(val)>log(a/aini))
           //  std::cout<<val<<" "<<log(a/aini)<<"\n";
                return val;
            }
            break;
        case 1:
            {

                double ali=a-Hsor/k;
                if(aini>ali) return cas2;
                else if(a>ali){
                //std::cout<<"ini\n";
                    double alno=a-(ali-aini);
                    double bno=aini;
                    return
                    cos(k/Hsor*a)*(gsl_sf_Ci(k/Hsor*a)-gsl_sf_Ci(k*ali/Hsor))+sin(k/Hsor*a)*(gsl_sf_Si(k*a/Hsor)-gsl_sf_Si(k*ali/Hsor));
                }
                else return 0.;

             //   return cos(k/Hsor*a)*(gsl_sf_Ci(k/Hsor*alno)-gsl_sf_Ci(k*(bno)/Hsor))+sin(k/Hsor)*(gsl_sf_Si(k*alno/Hsor)-gsl_sf_Si(k*bno/Hsor));

            }
            break;
        case 2:
            return tauLini*Hsor/a;
            break;
        case 3:
        {
            double result=0.0;
            double abserr=arr;
            double result2=0.0;
            double abserr2=arr;
            double epsabi1=errafa;//*1.0E-10;
            const double epsreli1=5.0E-2;
            gsl_function F;
            F.function=&dratei2;
            F.params=this;
            gsl_integration_qawo(&F,0,epsabi1,epsreli1,intsubi,wsp,qawoit,&result,&abserr);
            return result;
          //  return 0.;
        }
            break;
        default:
            return 0.;
            break;
    }
}

double gwenint::intfunc(){
    //double p2=q*(q-k*co)+k*(k-q*co);
    //p=sqrt(p2);
    //double co2=co*co;
    //double q2=q*q;
    //double p2=p*p;
    //double k2=k*k;
    double kq=(q2+k2-p2)/(2.0*q*k);
    double kp=(p2+k2-q2)/(2.0*k*p);
    double pq=(k2-p2-q2)/(2.0*q*p);
    veq=vspec.get_enk(q);
    vep=vspec.get_enk(p);
    veq2=vspec.get_enk2(q);
    vep2=vspec.get_enk2(p);
    beq=bspec.get_enk(q);
    bep=bspec.get_enk(p);
    //5oefp=1.+kq*kq*(2.+kp*kp);
    coefp=1.+kq*kq*(1.+kp*kp)+1.*kp*kp;
    coefm=(1.-kq*kq*(1.-kp*kp)-kp*kp)/coefp;
    //coefmq=kp*kp*(1.-kq*kq)/coefp;
    //coefmp=kq*kq*(1.-kp*kp)/coefp;
    coefmq=(1.+kq*kq)*(1.-pq*pq)/coefp;
    coefmp=(1.+kp*kp)*(1.-pq*pq)/coefp;
   // double prekqp=;
    double Et=(vspec.get_enk(q)*vspec.get_enk(p)+bspec.get_enk(q)*bspec.get_enk(p));
    double Ct=(hvspec.get_enk(q)*hvspec.get_enk(p)+hbspec.get_enk(q)*hbspec.get_enk(p));
    Ctf=4.*Ct*kq*kp/coefp;
    if(k<min || p<min || q<min){
        return 0.0;
    } else {
        double val=0.;
        if(!compr && !magson){
        val=k2*funcintpref*drateint()*(Et*(coefp)+4.0*Ct*kq*kp);
      //  std::cout<<veq<<"\n";
        }
        else 
        val=k2*funcintpref*drateint()*coefp;//*((1.+kq*kq*(2.0+kp*kp)));
       // std::cout<<val<<" "<<k2<<" "<<kq<<" "<<kp<<" "<<Et<<" "<<Ct<<" "<<vspec.get_enk(1.0E0)<<"\n";
        return val;
    }
}









const bool tanvarp=false;
const bool supervarp=true;
double integrateop(double t, void *data){
    gwenint *gwenintfi=static_cast<gwenint *>(data);
	double tn=t;
	double dtn=1.;
//	sinpqt(0,0,t,&tn,&dtn);
	t=tn;	
    if(!tanvarp && !supervarp){
        double tiab=t+gwenintfi->getiabin();
        if(tiab>0.0){
            double p=exp(gwenintfi->getbin()-(1.-t)/tiab);
            if(p<min || p>max) return 0.0;
        //    std::cout<<p<<" "<<gwenintfi->getbin()<<"\n";
            gwenintfi->setp(p);
            return gwenintfi->intfunc()*(1.+gwenintfi->getiabin())/(tiab*tiab)*dtn;
        } else 
            return 0.0;
    } else if(tanvarp && !supervarp){
        double ptan=tan(gwenintfi->gettanb()*t+gwenintfi->gettana());
        double p=exp(ptan);
        if(p<min || p>max) return 0.0;
        gwenintfi->setp(p);
        return gwenintfi->intfunc()*(1.+ptan*ptan)*gwenintfi->gettanb()*dtn;
    } else {
        double apa=gwenintfi->getprsca();
        if(t==0.) return 0.;
        double p=-apa*log(t);
        gwenintfi->setp(p);
        return gwenintfi->intfunc()*dtn/t/p*apa;
        //return gwenintfi->intfunc()*dtn/t*apa;
    }
}

double gwenint::intintfunc(){
    double ab=0.;
    double bb=1.;
    if(supervarp){
        ab=prao;
        bb=prbo;
    }
    gsl_function F;
    const int keyp=GSL_INTEG_GAUSS15;
    F.function=&integrateop;
    F.params=this;
    double result=0.0;
    double abserr=arr;
    double epsabsp=iabserr;//errafa*std::max(vspec.get_enk(q),bspec.get_enk(q))*getmaxamp();
    const double epsrelp=2.0E-2;
  //  std::cout<<"prog "<<result<<" "<<q<<" "<<vspec.get_enk(q)<<"\n";
    gsl_integration_qag(&F,ab,bb,epsabsp,epsrelp,intsubip,keyp,wspp,&result,&abserr);
 //   std::cout<<"zero step\n";
  //  std::cout<<"prog "<<result<<" "<<q<<"\n";
    return result;
}

const bool tanvar=false;
const bool supervarq=true;
double integrateoq(double t,void *data){
    gwenint *gwenintfi=static_cast<gwenint *>(data);
	double tn=t;
	double dtn=1.;
//	sinpqt(0,0,t,&tn,&dtn);
    if(!supervarq)
        t=-1.+2.*tn;	
    if(!tanvar && !supervarq){
        double t2=t*t;
        if(t2<=1.-fp){
            double tpi2=1./(1.-t2);
            double q=exp(t*tpi2);
        //    std::cout<<"t "<<t<<" "<<q<<"\n";
            if(q<min || q>max) return 0.0;
            gwenintfi->setq(q);
            return gwenintfi->intintfunc()*(1.+t2)*tpi2*tpi2*dtn;
        } else 
            return 0.0;
    } else if(tanvar && !supervarq){
        double qtap=tan(0.5*pi*t);
        double q=exp(qtap);
        if(q<min || q>max) return 0.0;
        gwenintfi->setq(q);
        return gwenintfi->intintfunc()*(1.+qtap*qtap)*pi*0.5*dtn*2.;
    } else {
        double apa=gwenintfi->getqrsca();
        double q=-apa*log(t);
        if(t==0.) return 0.;
        if(q<min || q>max) return 0.0;
        gwenintfi->setq(q);
        return gwenintfi->intintfunc()*dtn/t/q*apa;
     //   return gwenintfi->intintfunc()*dtn/t*apa;
    }
}

double gwenint::finalint(){
    const double ab=0.;//other -1.
    const double bb=1.;
    const int keyq=GSL_INTEG_GAUSS15;
    gsl_function F;
    F.function=&integrateoq;
    F.params=this;
    double result=0.0;
    double abserr=0.0;
    const double epsabsq=iabserr;//errafa*pow(std::max(vspec.getamp(),bspec.getamp()),2);
    const double epsrelq=3.0E-2;
    gsl_integration_qag(&F,ab,bb,epsabsq,epsrelq,intsubiq,keyq,wspq,&result,&abserr);
   // if(true) std::cout<<"one step "<<result<<" "<<k<<" "<<a/aini<<"\n";
    return result;
}

double integrateotf(double av,void *data){
    gwenint *gwenintfi=static_cast<gwenint *>(data);
    double a=exp(av);
    if(a>aini*(1.+fp)){
        gwenintfi->setat(a);
        return gwenintfi->finalint();//*a; //factor a gets cancelled by 1/a in formula
    } else 
        return 0.0;
}
// adapt conf. time
double whollyintegrated(double ato,double abe,gwenint &tgwenint,gsl_integration_cquad_workspace *gwwsp){
   double aa=log(ato);
   double ba=log(abe);
   gsl_function F;
   F.function=&integrateotf;
   F.params=&tgwenint;
   double result=0.0;
   double abserr=0;
   size_t nevals=0;
   const double epsabst=errafa*pow(tgwenint.getmaxamp(),2);
   const double epsrelt=5.0E-2;
   gsl_integration_cquad(&F,aa,ba,epsabst,epsrelt,gwwsp,&result,&abserr,&nevals);
   return result;
}



  //  odeiv2_driver *d=gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,mi3,intatol,presicio);
  //     status=gsl_odeiv2_driver_apply(d, &pi1,pim,f);
  //      gsl_odeiv2_driver_free(d);






const int subitq=10000;

void calcspeclogk(){

    std::vector<double> gwgp=gwgrid();
    std::vector<gwenint> gwenvec(Ngw);
    std::vector<gsl_integration_cquad_workspace *> cquadwsp(Ngw);
    std::vector<double> gwenergog(Ngw,0.0);
    std::cout<<"starting "<<aini<<" "<<abuild<<"\n";
    gsl_set_error_handler_off();
  // std::cout<<aini<<" "<<amax<<"\n";
 //  exit(1);
    #pragma omp parallel for schedule(dynamic,1)
    for(int i=0;i<Ngw;i++){
        double gr2=pow(gaSff(amax),2);
        enspec vspe=enspec(Omvir,aini,Lini,initds);
        enspec bspe=enspec(Ombir,aini,Lini,initds);
        enspec hvspe=enspec(helvr*Omvir,aini,Lini,initds);
        enspec hbspe=enspec(helbr*Ombir,aini,Lini,initds);
        gwenvec[i]=gwenint(gwgp[i],aini,Lini,vspe,bspe,hvspe,hbspe);
        cquadwsp[i]=gsl_integration_cquad_workspace_alloc(subitq);
   //     if(gwgp[i]>0.3 && gwgp[i]<0.8) gwenergog[i]=0.;

     gwenergog[i]=gr2*whollyintegrated(aini,amax,gwenvec[i],cquadwsp[i]);
       // std::cout<<"one resolved sc: "<<gwgp[i]<<" i "<<i<<" ama "<<amax<<" ami "<<aini<<" ener "<<gwenergog[i]<<"\n";
        std::cout<<gwgp[i]<<" "<<gwenergog[i]<<"\n";
    }
    std::ofstream output;
    output.open(direc+"/"+fname,std::ios_base::trunc);
    for(int i=0;i<Ngw;i++){
     //   if((gwgp[i]<0.3 || gwgp[i]>0.8) 
        output<<gwgp[i]<<" "<<gwenergog[i]<<"\n";
    }
    output.close();
}
