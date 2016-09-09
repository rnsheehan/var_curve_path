#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the special functions declared in the namespace
// R. Sheehan 23 - 11 - 2015

void special::gser(double *gamser,double a,double x,double *gln)
{
	// This function returns the incomplete gamma function P(a,x), calculated by its series representation
	// Also returns ln[gamma(a)] as gln

	int n;
	double sum,del,ap;

	static const int ITMAX=(100);
	static const double EPS=(3.0e-7);

	*gln=gammln(a);
	if(x<=0.0){
		if(x<0.0)
			std::cerr<<"x less than 0 in routine GSER"<<"\n";
		*gamser=0.0;
		return;
	}else{
		ap=a;
		del=sum=1.0/a;
		for(n=1;n<=ITMAX;n++){
			ap+=1.0;
			del*=x/ap;
			sum+=del;
			if(fabs(del)<fabs(sum)*EPS){
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		std::cerr<<"a too large, ITMAX too small in routine GSER"<<"\n";
		return;
	}
}

void special::gcf(double *gammcf, double a, double x, double *gln)
{
	// This function returns the incomplete gamma function Q(a,x), calculated by its continued fraction representation
	// Also returns ln[gamma(a)] as gln

	int n;
	double gold=0.0,g,fac=1.0,b1=1.0;
	double b0=0.0,anf,ana,an,a1,a0=1.0;

	static const int ITMAX=(100);
	static const double EPS=(3.0e-7);

	*gln=gammln(a);
	a1=x;
	for(n=1;n<=ITMAX;n++){
		an=static_cast<double>(n);
		ana=an-a;
		a0=(a1+a0*ana)*fac;
		b0=(b1+b0*ana)*fac;
		anf=an*fac;
		a1=x*a0+anf*a1;
		b1=x*b0+anf*b1;
		if(a1){
			fac=1.0/a1;
			g=b1*fac;
			if(fabs((g-gold)/g)<EPS){
				*gammcf=exp(-x+a*log(x)-(*gln))*g;
				return;
			}
			gold=g;
		}
	}
	std::cerr<<"a too large, ITMAX too small in routine GCF"<<"\n";
}

double special::gammp(double a,double x)
{
	// Computes the incomplete gamma function P(a,x) from the functions gser and gcf

	double gamser,gammcf,gln;

	if(x<0.0||a<=0.0) 
		std::cerr<<"Invalid arguments in routine GAMMP"<<"\n";
	if(x<(a+1.0)){
		gser(&gamser,a,x,&gln);
		return gamser;
	} 
	else{
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}

double special::gammq(double a,double x)
{
	// Computes the incomplete gamma function Q(a,x)=1-P(a,x) from the functions gser and gcf

	double gamser,gammcf,gln;

	if(x < 0.0 || a <= 0.0){
		std::cerr<<"Invalid arguments in routine GAMMQ"<<"\n";
	}
	if(x<(a+1.0)){
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} 
	else{
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}

double special::gammln(double xx)
{
	// Return the value of ln[gamma(xx)] for xx>0

	double x,tmp,ser;
	static double cof[6]={76.18009173,-86.50532033,24.01409822,-1.231739516,0.120858003e-2,-0.536382e-5};
	int j;

	x=xx-1.0;
	tmp=x+5.5;
	tmp-=(x+0.5)*log(tmp);
	ser=1.0;
	for (j=0;j<=5;j++){
		x+=1.0;
		ser+=(cof[j]/x);
	}
	return -tmp+log(2.50662827465*ser);
}

double special::erff(double x)
{
	// Return the error function erf(x)

	return ( x < 0.0 ? -gammp(0.5, template_funcs::DSQR(x) ) : gammp(0.5, template_funcs::DSQR(x) ) ) ; 
}

double special::erffc(double x)
{
	// Return the complementary error function erfc(x)

	return ( x < 0.0 ? 1.0 + gammp(0.5, template_funcs::DSQR(x) ) : gammq(0.5, template_funcs::DSQR(x) ) ) ; 
}

void special::two_F_one(double a, double b, double c, double x, double &F, double &dF)
{
	// Implementation of the hypergeometric function, {2}_F_{1}(a, b, c; z)
	// This series is valid for |x| < 1, with rapid convergence when |x| <= 1/2
	// Series is divergent when (c - a - b) <= -1
	// Series is absolutely convergent when (c - a - b) > 0
	// Series is conditionally convergent when -1 < (c - a - b) < 0 and x = 1 is excluded
	// See Abramowitz and Stegun, Ch. 15

	// This implementation is based on that given in NRinC, sect. 6.13
	// The function returns the value of F and its derivative at some point x
	// An implementation involving std::complex arguments is also possible
	// R. Sheehan 15 - 4 - 2014

	if(fabs(x) < 1.0){

		int i, nterms; 
		double fac, temp, aa, bb, cc;

		F = dF = 0.0; 
		fac = 1.0; temp = fac;
		aa = a; bb = b; cc = c; 

		nterms = 1000; // compute the partial sum out to nterms

		for(i=1; i<=nterms; i++){
		
			fac *= ((aa*bb)/cc); 
		
			dF += fac; 
		
			fac *= ( ( 1.0 / ( static_cast<double>(i) ) )*x);
		
			F = temp + fac; 
		
			if(F == temp){
				return; 
			}
		
			temp = F; 

			aa += 1.0; bb += 1.0; cc += 1.0; 
		}

	}
	else{
		
		F = dF = 0.0; 

	}
}

double special::factorial(int n)
{
	// return the factorial of n by recursion
	// It is possible to represent 170! as a floating point number
	// Attempting to compute 171! will cause numerical overflow
	// Factorials up to 22! are exact, 23! and higher are numerically approximate
	// See NRinC, sect 6.1. 
	// R. Sheehan 15 - 4 - 2014

	if(n >-1 && n < 170){

		if(n == 0 || n == 1){
			return 1.0; 
		}
		else if(n == 2){
			return 2.0; 
		}
		else if(n == 3){
			return 6.0; 
		}
		else if(n == 4){
			return 24.0; 
		}
		else if(n == 5){
			return 120.0; 
		}
		else if(n == 6){
			return 720.0; 
		}
		else{
			// recursively compute n!
			return ( static_cast<double>(n)*factorial(n-1) ); 
		}

	}
	else{
		std::cerr<<"Cannot compute "<<n<<"! by this method\nc.f. gammln routine for factorials of large numbers"; 
		return 0.0; 
	}
}

double special::bessj0(double x)
{
	//Return the Bessel Function J0(x) for all real x
	double z,ax,xx,y,ans,ans1,ans2;

	if((ax=fabs(x))<8.0){
		y=template_funcs::DSQR(x);
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7+
			y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029533985.0+y*(9494680.718+
			y*(59272.64853+y*(267.8532712+y))));
		ans=ans1/ans2;
	}
	else{
		z=8.0/ax;
		y=template_funcs::DSQR(z);
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4+
			y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2=-0.1562499995e-1+y*(0.1430488765e-3+
			y*(-0.6911147651e-5+y*(0.7621095161e-6-y*0.934945152e-7)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}

double special::bessj1(double x)
{
	//Returns the Bessel function J1(x) for all real x
	double ax,z,xx,y,ans,ans1,ans2;

	if((ax=fabs(x))<8.0){
		y=template_funcs::DSQR(x);
		ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
			+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
		ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
			+y*(99447.43394+y*(376.9991397+y))));
		ans=ans1/ans2;
	}
	else{
		z=8.0/ax;
		y=template_funcs::DSQR(z);
		xx=ax-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
			+y*(0.245752017e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
			+y*(0.8449199096e-5+y*(-0.88228987e-6
			+y*0.105787412e-6)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
		if(x<0.0) ans=-ans;
	}
	return ans;
}

double special::bessj(int n,double x)
{
	//Returns te Bessel function Jn(x) for all real x and n >= 2
	int j,jsum,m;
	double ax,bj,bjm,bjp,sum,tox,ans;

	static const double ACC=40.0;
	static const double BIGNO=1e10;
	static const double BIGNI=1.0e-10;

	if(n<2){
		std::cerr<<"Index n less than 2 in bessj\n";	
	}
	ax=fabs(x);
	if(ax==0.0){
		return 0.0;	
	}
	else if(ax>static_cast<double>(n)){
		//Upwards recurrence from J0 and J1
		tox=2.0/ax;
		bjm=bessj0(ax);
		bj=bessj1(ax);
		for(j=1;j<n;j++){
			bjp=j*tox*bj-bjm;
			bjm=bj;
			bj=bjp;
		}
		ans=bj;
	}
	else{
		//Downwards recurrence from an even m here computed
		tox=2.0/ax;
		m=2*((n+static_cast<int>(sqrt(ACC*n)))/2);
		jsum=0;
		bjp=ans=sum=0.0;
		bj=1.0;
		for(j=m;j>0;j--){
			bjm=j*tox*bj-bjp;
			bjp=bj;
			bj=bjm;
			if(fabs(bj)>BIGNO){
				bj*=BIGNI;
				bjp*=BIGNI;
				ans*=BIGNI;
				sum*=BIGNI;
			}
			if(jsum) sum+=bj;
			jsum=!jsum;
			if(j==n) ans=bjp;
		}
		sum=2.0*sum-bj;
		ans/=sum;
	}
	return x<0.0 && (n&1) ? -ans:ans; 
}

double special::bessel_J(int n,double x)
{
	if(n==0){
		return bessj0(x);	
	}
	else if(n==1){
		return bessj1(x);	
	}
	else{
		return bessj(n,x);	
	}
}

double special::bessy0(double x)
{
	//Returns the Bessel Function Y0(x) for real positive values of x
	double z,xx,y,ans,ans1,ans2;

	if(x<8.0){
		y=template_funcs::DSQR(x);
		ans1=-2957821389.0+y*(7062834065.0+y*(-512359803.6+
			y*(10879881.29+y*(-86327.92757+y*228.4622733))));
		ans2=40076544269.0+y*(745249964.8+y*(7189466.438+
			y*(47447.26470+y*(226.1030244+y))));
		ans=(ans1/ans2)+0.636619772*bessj0(x)*log(x);
	}
	else{
		z=8.0/x;
		y=template_funcs::DSQR(z);
		xx=x-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4+
			y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2=-0.1562499995e-1+y*(0.1430488765e-3+
			y*(-0.6911147651e-5+y*(0.7621095161e-6+y*(-0.934945152e-7))));
		ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
	}
	return ans;
}

double special::bessy1(double x)
{
	//Returns the Bessel Function Y1(x) for real positive x
	double z,xx,y,ans,ans1,ans2;

	if(x<8.0){
		y=template_funcs::DSQR(x);
		ans1=x*(-0.4900604943e13+y*(0.1275274390e13
			+y*(-0.5153438139e11+y*(0.7349264551e9
			+y*(-0.4237922726e7+y*0.8511937935e4)))));
		ans2=0.2499580570e14+y*(0.424441966e12
			+y*(0.3733650367e10+y*(0.2245904002e8
			+y*(0.1020426050e6+y*(0.3549632885e3+y)))));
		ans=(ans1/ans2)+0.636619772*(bessj1(x)*log(x)-1.0/x);
	}
	else{
		z=8.0/x;
		y=template_funcs::DSQR(z);
		xx=x-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
			+y*(0.245752017e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
			+y*(0.8449199096e-5+y*(-0.88228987e-6
			+y*0.105787412e-6)));
		ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
	}
	return ans;
}

double special::bessy(int n,double x)
{
	//Returns the Bessel Function Yn(x) for n >= 2
	int j;
	double by,bym,byp,tox;

	if(n<2){
		std::cerr<<"Index n less than 2 in bessy\n";
	}
	tox=2.0/x;
	by=bessy1(x);
	bym=bessy0(x);
	for(j=1;j<n;j++){
		byp=j*tox*by-bym;
		bym=by;
		by=byp;
	}
	return by;
}

double special::bessel_Y(int n,double x)
{
	if(n==0){
		return bessy0(x);
	}
	else if(n==1){
		return bessy1(x);
	}
	else{
		return bessy(n,x);
	}
}

double special::bessi0(double x)
{
	//Returns the modified Bessel Function I0(x) for any real x
	double ax,ans,y;

	if((ax=fabs(x))<3.75){
		y=x/3.75;
		y*=y;
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	}
	else{
		y=3.75/ax;
		ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.132859e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2))))))));
	}
	return ans;
}

double special::bessi1(double x)
{
	//Returns the modified Bessel Function I1(x) for any real x
	double ax,ans,y;
	if((ax=fabs(x))<3.75){
		y=x/3.75;
		y*=y;
		ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
			+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
	}
	else{
		y=3.75/ax;
		ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
			-y*0.420059e-2));
		ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
			+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
		ans*=(exp(ax)/sqrt(ax));
	}
	return x<0.0?-ans:ans;
}

double special::bessi(int n,double x)
{
	//Returns the modified Bessel Function In(x) for n>=2 for any real x
	int j;
	double bi,bim,bip,tox,ans;

	static const double ACC=40.0;
	static const double BIGNO=1e10;
	static const double BIGNI=1.0e-10;

	if(n<2) std::cerr<<"Index n less than 2 in bessi\n";

	if(x==0.0){
		return 0.0;
	}
	else{
		tox=2.0/fabs(x);
		bip=ans=0.0;
		bi=1.0;
		for(j=2*(n+static_cast<int>(sqrt(ACC*n)));j>0;j--){
			bim=bip+j*tox*bi;
			bip=bi;
			bi=bim;
			if(fabs(bi)>BIGNO){
				ans*=BIGNI;
				bi *=BIGNI;
				bip*=BIGNI;
			}
			if(j==n) ans=bip;
		}
		ans*=bessi0(x)/bi;
		return x<0.0 && (n&1) ? -ans:ans;
	}
}

double special::bessel_I(int n,double x)
{
	if(n==0){
		return bessi0(x);
	}
	else if(n==1){
		return bessi1(x);
	}
	else{
		return bessi(n,x);
	}
}

double special::bessk0(double x)
{
	//Returns the modified Bessel function K0(x) for positive real x
	double y,ans;
	if(x<=2.0){
		y=template_funcs::DSQR(x)/4.0;
		ans=(-log(x/2.0)*bessi0(x))+(-0.57721566+y*(0.42278420
			+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
			+y*(0.10750e-3+y*0.74e-5))))));
	}
	else{
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
			+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
			+y*(-0.251540e-2+y*0.53208e-3))))));
	}
	return ans;
}

double special::bessk1(double x)
{
	//Returns the modified Bessel function K1(x) for real positive x	
	double y,ans;
	if(x<=2.0){
		y=x*x/4.0;
		ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144
			+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
			+y*(-0.110404e-2+y*(-0.4686e-4)))))));
	}
	else{
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
			+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
			+y*(0.325614e-2+y*(-0.68245e-3)))))));
	}
	return ans;
}

double special::bessk(int n,double x)
{
	//Returns the modified Bessel function Kn(x) for n >=2
	int j;
	double bk,bkm,bkp,tox;

	if(n<2) std::cerr<<"Index n less than 2 in bessk\n";

	tox=2.0/x;
	bkm=bessk0(x);
	bk=bessk1(x);
	for(j=1;j<n;j++){
		bkp=bkm+j*tox*bk;
		bkm=bk;
		bk=bkp;
	}
	return bk;
}

double special::bessel_K(int n,double x)
{
	if(n==0){
		return bessk0(x);
	}
	else if(n==1){
		return bessk1(x);
	}
	else{
		return bessk(n,x);
	}
}

void special::fresnel(double x, double *s, double *c)
{
	// Computes the Fresnel integrals S(x) and C(x) for all real x.
	// Taken from NRinC
	// R. Sheehan 26 - 3 - 2009

	int k, n, odd;

	double a,ax,fact,pix2,sign,sum,sumc,sums,term,test;

	std::complex<double> b,cc,d,h,del,cs;

	static const int MAXIT=100;
	static const int TRUE=1;

	static const double FPMIN=1.0e-30;
	static const double XMIN=1.5;
	static const double EPS=(3.0e-12);

	ax=fabs(x);

	if(ax<sqrt(FPMIN)){					// Special case: avoid failure of convergence
		*s=0.0;							// test because of underflow.
		*c=ax;
	}else if(ax<=XMIN){					// Evaluate both series simultaneously.
		sum=sums=0.0;
		sumc=ax;
		sign=1.0;
		fact=PI_2*ax*ax;
		odd=TRUE;
		term=ax;
		n=3;
		for(k=1;k<=MAXIT;k++){
			term *= fact/k;
			sum += sign*term/n;
			test=fabs(sum)*EPS;
			if(odd){
				sign = -sign;
				sums=sum;
				sum=sumc;
			}else{
				sumc=sum;
				sum=sums;
			}
			if(term<test)break;
			odd=!odd;
			n+=2;
		}
//		if (k > MAXIT) nrerror("series failed in frenel");
		*s=sums;							
		*c=sumc;
	}else{								// Evaluate Integrals by use of Lentz's Continued Fraction Expansion Method
		pix2=PI*ax*ax;
		b=std::complex<double>(1.0,-pix2);
		cc=std::complex<double>(1.0/FPMIN,0.0);
		d=h=one/b;
		n = -1;
		for(k=2;k<=MAXIT;k++){
			n+=2;
			a=-n*(n+1);
			b=b+4.0;
			d=one/((a*d)+b);
			cc=(b+(a/cc));
			del=(cc*d);
			h=(h*del);
			if(fabs(del.real()-1.0)+fabs(del.imag())<EPS)break;
		}
//		if (k > MAXIT) nrerror("cf failed in frenel");
		h=(std::complex<double>(ax,-ax)*h);
		cs=(std::complex<double>(0.5,0.5)*(one-(std::complex<double>(cos(0.5*pix2),sin(0.5*pix2))*h)));
		*c=cs.real();
		*s=cs.imag();
	}
	if(x<0.0){						// Use antisymmetry to obtain values for negative x
		*c=-(*c);
		*s=-(*s);
	}
}

double special::Ell_K(double x, bool conjugate)
{
	// Complete elliptic integral of the first kind defined by Hypergeometric function
	// Abramowitz and Stegun, Ch 17, sect 17.3
	// function computes K'(k) = K(1-k) when conjugate = true, otherwise returns K(k)
	// R. Sheehan 23 - 11 - 2015

	try{

		if(x >= 0.0 && x < 1){

			double F, dF; 

			special::two_F_one( 0.5, 0.5, 1, ( conjugate ? 1.0 - x :  x ), F, dF); 
		
			return (F*PI_2);

		}
		else{
			std::string reason; 

			reason = "Error in special::Ell_K(double x, bool conjugate)\n"; 
			reason += "input argument x = " + template_funcs::toString(x, 3) + "\n"; 

			throw std::invalid_argument(reason); 
		}

	}
	catch(std::invalid_argument &e){		
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double special::Ell_E(double x, bool conjugate)
{
	// Complete elliptic integral of the second kind defined by Hypergeometric function
	// Abramowitz and Stegun, Ch 17, sect 17.3
	// function computes E'(k) = E(1-k) when conjugate = true, otherwise returns E(k)
	// R. Sheehan 23 - 11 - 2015

	try{

		if(x >= 0.0 && x < 1){

			double F, dF; 

			special::two_F_one( -0.5, 0.5, 1, ( conjugate ? 1.0 - x :  x ), F, dF); 
		
			return (F*PI_2);

		}
		else{
			std::string reason; 

			reason = "Error in special::Ell_E(double x, bool conjugate)\n"; 
			reason += "input argument x = " + template_funcs::toString(x, 3) + "\n"; 

			throw std::invalid_argument(reason); 
		}

	}
	catch(std::invalid_argument &e){		
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void special::Ell_K_E(double k, double &Kval, double &Eval, bool conjugate)
{
	// complete elliptic integrals of the first and second kind
	// polynomial approximation accurate to within |eps| <= 2e-8 when 0 <= k < 1
	// Abramowitz and Stegun, Ch 17, sect 17.3
	// functions compute K'(k) = K(1-k) when conjugate = true, otherwise returns K(k), or E(k)
	// R. Sheehan 23 - 11 - 2015

	//double arg = (conjugate ? k: 1.0 - k); // decide on the argument, i.e. return K(k) or K'(k) = K(1-k)

	try{

		if(k >= 0.0 && k < 1){

			// polynomial coefficients for K(k)
			double AK[5] = {1.38629436112, 0.09666344259, 0.03590092383, 0.03742563713, 0.01451196212}; 
			double BK[5] = {0.5, 0.12498593597, 0.06880248576, 0.03328355346, 0.00441787012};

			// polynomial coefficients for E(k)
			double AE[5] = {1.0, 0.44325141463, 0.06260601220, 0.04757383546, 0.01736506451}; 
			double BE[5] = {0.0, 0.24998368310, 0.09200180037, 0.04069697526, 0.00526449639};

			double m1 = (1.0 - ( conjugate ? 1.0 - k :  k ) ); 
			double lnm1 = log(1.0/m1); 
			double akval, bkval, aeval, beval; 

			// evaluate all polynomials simultaneously
			akval = AK[4]; bkval = BK[4]; // K(k)
			aeval = AE[4]; beval = BE[4]; // E(k)
			for(int i=3; i>=0; i--){
				// K(k)
				akval = akval*m1 + AK[i]; 
				bkval = bkval*m1 + BK[i]; 

				// B(k)
				aeval = aeval*m1 + AE[i]; 
				beval = beval*m1 + BE[i]; 
			}

			// multiply b polynomials by lnm1
			bkval *= lnm1; beval *= lnm1; 

			// Assign values for K(k) and E(k)
			Kval = (akval + bkval); 
			Eval = (aeval + beval);

		}
		else{
			std::string reason; 

			reason = "Error in special::Ell_K_E(double k, double &Kval, double &Eval, bool conjugate)\n"; 
			reason += "input argument k = " + template_funcs::toString(k, 3) + "\n"; 

			throw std::invalid_argument(reason); 
		}

	}
	catch(std::invalid_argument &e){		
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}