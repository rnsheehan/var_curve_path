#ifndef SPECIAL_FUNCTIONS_H
#define SPECIAL_FUNCTIONS_H

// Declaration of namespace special
// Contains implementations of mathematical objects known as special functions
// Implementation is done for real variables for most functions will also work with type double replace
// with type complex<double> 
// List of functions is not exhaustive, functions are added to namespace as required
// R. Sheehan 23 - 11 - 2015

// For connections between special functions see http://www.johndcook.com/blog/special_function_diagram/ 

namespace special{

	// Gamma Function
	double gammln(double x); // Computes the value of ln[gamma(xx)] for xx>0

	// Incomplete Gamma Function
	double gammp(double a, double x); // Computes the incomplete gamma function P(a,x) from the functions gser and gcf
	double gammq(double a, double x); // Computes the incomplete gamma function Q(a,x)=1-P(a,x) from the functions gser and gcf

	void gser(double *gamser,double a,double x,double *gln); // Computes the incomplete gamma function P(a,x), calculated by its series representation
	void gcf(double *gammcf,double a,double x,double *gln); // Computes the incomplete gamma function Q(a,x), calculated by its continued fraction representation
	
	// Error Function
	double erff(double x); // erf(x)
	double erffc(double x); // 1 - erf(x)

	// Factorial Function
	double factorial(int n); 

	// Hypergeometric Function
	void two_F_one(double a, double b, double c, double x, double &F, double &dF); 

	// Bessel Functions
	double bessel_J(int n, double x); // Bessel Function of the 1st kind Jnu(x)

	double bessel_Y(int n,double x); // Bessel Function of the 2nd kind Ynu(x)

	double bessel_I(int n,double x); // Modified Bessel Function Inu(x)

	double bessel_K(int n,double x); // Modified Bessel Function Knu(x)

	//Bessel Function of the 1st kind
	double bessj0(double x);
	double bessj1(double x);
	double bessj(int n,double x);

	//Bessel Function of the 2nd kind
	double bessy0(double x);
	double bessy1(double x);
	double bessy(int n,double x);

	//Modified Bessel Function Inu(x)
	double bessi0(double x);
	double bessi1(double x);
	double bessi(int n,double x);

	//Modified Bessel Function Knu(x)
	double bessk0(double x);
	double bessk1(double x);
	double bessk(int n,double x);
	
	// Fresnel Integrals
	void fresnel(double x, double *s, double *c); // Fresnel Integrals 

	// Complete Elliptic Integrals of the First and Second Kinds

	double Ell_K(double x, bool conjugate = false); // Elliptic integral of first kind defined by Hypergeometric Function
	double Ell_E(double x, bool conjugate = false); // Elliptic integral of second kind defined by Hypergeometric Function

	void Ell_K_E(double k, double &Kval, double &Eval, bool conjugate = false); // polynomial approximation to K(k) and E(k) and their conjugates
}

#endif