#ifndef TEMPLATES_H
#define TEMPLATES_H

// Template fucntions for commonly performed mathematical operations
// R. Sheehan 25 - 3 - 2013

namespace template_funcs{

	template <class T> void Quad_Solve(T a,T b,T c,T &rt1,T &rt2)
	{
		// Solve the quadratic equation a x^{2} + b x + c = 0
	
		T q;

		q = -0.5*( b + Signum(b) * sqrt( ( template_funcs::DSQR(b)-4.0*a*c ) ) );
	
		rt1=(q/a);
	
		rt2=(c/q); 
	}
	
	template <class T> T Signum(T a)
	{
		// The sign operator
		T darg;
		//return ((darg=(a))==0.0?0.0:(darg=(a))>=0.0?1.0:-1.0);
		return ( (darg=(a)) >= (T)(0) ? (T)(1) : -(T)(1) ); // Setting the Sign of zero to be 1
	}

	template <class T> T DSQR(T a)
	{
		// Efficient squaring operator
		// Write injuries in dust, benefits in marble
		T darg;
		return ( (darg=(a)) == (T)(0) ? (T)(0) : darg*darg );
	}

	template <class T> T SIGN(T a,T b)
	{
		return ((b)>(T)(0)?fabs(a):-fabs(a));
	}

	template <class T> std::string toString(const T & t)
	{
		// This is just too convenient not to use
		// Is there a version that can include something similar to %0.5d ? 
		// There has to be, look into the setw method for strings
		// Requires the string-stream (sstream) standard library
		// R. Sheehan 16 - 5 - 2011
    
		std::ostringstream oss; // create a stream
		oss << t;				// insert value to stream
		return oss.str();		// return as a string
	}

	template <class T> std::string toString(const T &t,int places)
	{
		// toString function that allows for the
		// number of decimal places to be specified 
		// far too convenient
		// R. Sheehan 17 - 5 - 2011

		std::ostringstream oss; // create a stream

		oss<<std::fixed<<std::setprecision(places)<<t; // insert value to stream

		return oss.str(); // return as a string
	}

}

#endif