#ifndef ATTACH_H
#include "Attach.h"
#endif

// Source code for the XC_Bend class
// R. Sheehan 7 - 8 - 2012

// Constructors
XC_Bend::XC_Bend()
{
	// Default Constructor
	// R. Sheehan 7 - 8 - 2012

	coeffs_defined = false; 
	
	bend_type=Npts=signT=0; 
	
	width=t2=Reqc=Teqc=Leqc=Lbend=Lbend_2=Loverlap=Loverlayer=Lovertotal=0.0;
	xs=ys=sxc=syc=kt=keqc=kmax=0.0;
	r=alpha=sigma=nu=g=f=c1=c2=c3=c4=c5=c6=FCc3=FSc3=0.0;

	// populate the vector with the various bend types
	bend_types.push_back(CC); bend_types.push_back(LC); bend_types.push_back(TC); bend_types.push_back(QC);
}

XC_Bend::XC_Bend(int type, double W, double R, double T)
{
	// Constructor for CC, LC, TC bends
	// R. Sheehan 7 - 8 - 2012

	// populate the vector with the various bend types
	bend_types.push_back(CC); bend_types.push_back(LC); bend_types.push_back(TC); bend_types.push_back(QC);
	
	set_params(type,W,R,T);
}

XC_Bend::~XC_Bend()
{
	// Deconstructor

	clear_memory(); 
}

// Methods

// Getters
int XC_Bend::get_bend_type()
{
	// return the wg bend type
	// R. Sheehan 7 - 8 - 2012
	
	return bend_type; 
}

int XC_Bend::get_Npts()
{
	// return the number of points used to represent the bend
	// R. Sheehan 7 - 8 - 2012
	
	return Npts; 
}

std::string XC_Bend::get_bend_name()
{
	// return the name of the bend as a std::string

	std::string name; 

	try{

		if(bend_type == CC){
			return "CC"; 
		}
		else if(bend_type == LC){
			return "LC"; 
		}
		else if(bend_type == TC){
			return "TC"; 
		}
		else if(bend_type == QC){
			return "QC"; 
		}
		else{
			std::string reason; 
			reason = "Error: Bend Type not defined in XC_Bend::get_bend_name()\n"; 
			throw std::logic_error(reason); 
		}

	}
	catch(std::logic_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double XC_Bend::get_bend_length()
{
	// return the bend length
	// R. Sheehan 8 - 8 - 2012
	
	return Lbend; 
}

double XC_Bend::get_max_curvature()
{
	// return the maximum curvature
	// R. Sheehan 25 - 5 - 2016

	return kmax; 
}

// Setters
void XC_Bend::set_bend_type(int val)
{
	// set the bend type
	// R. Sheehan 7 - 8 - 2012
	
	bend_type = val;
}

//void XC_Bend::set_n_pts(int val)
//{
//	// set the number of points along the bend
//	// R. Sheehan 7 - 8 - 2012
//	
//	Npts = val;
//}

void XC_Bend::set_bend_rad(double val)
{
	// assign a value to the bend radius

	Reqc = val; 
}

void XC_Bend::set_params(int type, double W, double R, double T)
{
	// define parameters for CC, LC, TC bends
	// R. Sheehan 7 - 8 - 2012

	try{

		// check whether the bend type that was input is a valid value
		bool c1 = ( std::find( bend_types.begin(), bend_types.end(), type) != bend_types.end() ); 
		bool c2 = W > 0.0 ? true : false; // waveguide width > 0
		bool c3 = R > 0.0 ? true : false; // bend radius > 0
		bool c4 = T > 0.0 && T <= PI  ? true : false; // bend angle 0 < T <= \pi

		if(c1 && c2 && c3 && c4){
		
			bend_type = type;
	
			Npts = 1001; // this is a default value, no need to change it really
	
			signT = 1; 
	
			width = W; // waveguide width
			t2 = 0.5*width; // half waveguide width
			Reqc = R; // equivalent circle radius
			keqc = 1.0 / Reqc; // equivalent circle curvature
			Teqc = T; // bend angle 
			Leqc = R*T; // equivalent circle length
			Lbend = Leqc; // actual length of the bend
			Lbend_2 = 0.5*Lbend; // half the bend length
	
			/*Loverlap=0.01; //Size of overlap is 10 nm, represented here in um
			Loverlayer=1; // Default length of overlayer is 1 um extra
			Lovertotal = Loverlap+Loverlayer;*/
	
			// No longer using this, any errors encountered using it are yours to fix
			// R. Sheehan 8 - 9 - 2016
			Loverlap = Loverlayer = Lovertotal = 0.0; 
	
			xs=ys=sxc=syc=ks=0.0; 
	
			XC_coeffs(); 	
		}
		else{
			std::string reason = "Error: void XC_Bend::set_params(int type, double W, double R, double T)\n"; 
			if(!c1) reason += "bend type: " + template_funcs::toString(type) + " is not a valid choice\n"; 
			if(!c2) reason += "path width: " + template_funcs::toString(W, 3) + " cannot be negative\n"; 
			if(!c3) reason += "bend radius: " + template_funcs::toString(R, 3) + " cannot be negative\n"; 
			if(!c4) reason += "bend angle: " + template_funcs::toString(T, 3) + " must be in range 0 < T <= pi\n"; 

			throw std::invalid_argument(reason); 
		}	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void XC_Bend::XC_coeffs()
{
	// Assign the bend coefficients depending on the bend type
	// R. Sheehan 7 - 8 - 2012
	
	if(bend_type == LC){
		LC_coeffs();
	}
	else if(bend_type == TC){
		TC_coeffs();
	}	
	else if(bend_type == QC){
		QC_coeffs(); 
	}
}

void XC_Bend::LC_coeffs()
{
	// Assign values to the coefficients necessary to define a LC bend
	// R. Sheehan 7 - 8 - 2012
	
	double S, C; 

	alpha=(4.0/(Reqc*Lbend));
	c1=0.5*sqrt(PI*Reqc*Lbend);
	c2=(1.0/c1);
	//c3=sqrt(Teqc/PI); // Very subtle error here
	c3=sqrt(Lbend/(Reqc*PI)); 
	
	special::fresnel(c3, &S, &C); 
	FCc3 = C; 
	FSc3 = S; 

	coeffs_defined = true; 
}

void XC_Bend::TC_coeffs()
{
	// Assign values to the coefficients necessary to define a TC bend
	// R. Sheehan 7 - 8 - 2012

	double S, C; 
	
	f=0.25; // This is the fraction of the total length assigned to the first linear portion

	// Case 1: k_{max} > k_{CC}
	// standard case
	TC_type = 1; 
	g=(1.0/(1.0-f)); 
	kt=(g/Reqc); 		 
	
	// Case 2: k_{max} == k_{CC}
	/*TC_type = 2; 
	g=(1.0/(1.0-f)); 
	kt=(1.0/Reqc); 
	Lbend *= g;*/ 

	// Case 3: k_{max} < k_{CC}
	/*TC_type = 3; 
	r = 0.8; 
	g=(1.0/(1.0-f)); 
	g *= 1.0/r; 
	kt=(r/Reqc); 
	Lbend *= g;*/

	// alpha=(template_funcs::DSQR(kt)/(kt*Lbend-(Lbend/Reqc))); 
	// Very subtle error here
	// Error occurs because you're effectively re-defining theta with every iteration, 
	// when theta is most definitely fixed
	
	alpha=(template_funcs::DSQR(kt)/(kt*Lbend-Teqc));
	sigma=(kt/alpha); 
	nu=(Lbend-sigma);

	c1=sqrt(PI/alpha);
	c2=(1.0/c1);
	c3=(kt/(sqrt((alpha*PI))));
	c4=(kt*(Lbend-sigma));
	c5=(kt*(Lbend-1.5*sigma));
	c6=(0.5*kt*sigma);

	special::fresnel(c3, &S, &C); 
	FCc3 = C; 
	FSc3 = S; 

	coeffs_defined = true; 
}

void XC_Bend::QC_coeffs()
{
	// Assign values to the parameters associated with the QC_Bend

	alpha = ( (6.0*Teqc)/(Lbend*template_funcs::DSQR(Lbend)) ); // max curvature 
	nu = ( 2.0/(3.0*Lbend) ); // path-length scaling for x_coord, y_coord
	g = (1.0/nu); // 3 L / 2, used in the calculation of the x_coord

	coeffs_defined = true; 
}

void XC_Bend::define_eqc_coords()
{
	// Define the coordinates that define the equivalent circle
	// R. Sheehan 7 - 8 - 2012

	try{

		bool c1 = Npts > 1 ? true : false; 
		bool c2 = Teqc > 0.0 && Teqc <= PI  ? true : false; // bend angle 0 < T <= \pi

		if(c1 && c2){

			double ds=(Teqc/(Npts-1));
			double s;

			wg_Xeqc.resize(Npts+1);
			wg_Yeqc.resize(Npts+1);

			for(int i = 1; i <= Npts ; i++){
				s = (i-1)*ds; 
				wg_Xeqc[i] = CC_x_coord(s);
				wg_Yeqc[i] = CC_y_coord(s);
			}
		
		}
		else{
			std::string reason = "Error: void XC_Bend::define_eqc_coords()\n"; 
			if(!c1) reason += "Npts: " + template_funcs::toString(Npts) + " is not defined\n"; 
			if(!c2) reason += "bend angle: " + template_funcs::toString(Teqc, 3) + " must be in range 0 < T <= pi\n"; 

			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void XC_Bend::define_curvature()
{
	// define the bend curvature as a function of path-length
	// R. Sheehan 7 - 3 - 2009
	
	try{

		bool c1 = Npts > 1 ? true : false; 
		bool c2 = Lbend > 0.0 ? true : false; // L = R theta > 0

		if(c1 && c2){

			double ds = (Lbend/(Npts-1));
			double s;

			kmax = 0.0; 

			wg_curv.resize(Npts+1);

			for(int i=1;i<=Npts;i++){
				s = (i-1)*ds; 
				wg_curv[i] = curvature(s);
				if(kmax < wg_curv[i]) kmax = wg_curv[i]; // determine the max value of the curvature
			}
		
		}
		else{
			std::string reason = "Error: void XC_Bend::define_curvature()\n"; 
			if(!c1) reason += "Npts: " + template_funcs::toString(Npts) + " is not defined\n"; 
			if(!c2) reason += "Lbend: " + template_funcs::toString(Lbend, 3) + " cannot be negative\n"; 

			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}	
}

void XC_Bend::define_bend_angle()
{
	// Define the bending angle as a function of path-length
	// R. Sheehan 7 - 3 - 2009

	try{

		bool c1 = Npts > 1 ? true : false; 
		bool c2 = Lbend > 0.0 ? true : false; // L = R theta > 0

		if(c1 && c2){

			double ds = (Lbend/(Npts-1));
			double s;

			wg_bend_ang.resize(Npts+1);

			for(int i=1;i<=Npts;i++){
				s = (i-1)*ds; 
				wg_bend_ang[i] = bend_angle(s);
			}	
	
			re_scale_bend_angle();
		
		}
		else{
			std::string reason = "Error: void XC_Bend::define_bend_angle()\n"; 
			if(!c1) reason += "Npts: " + template_funcs::toString(Npts) + " is not defined\n"; 
			if(!c2) reason += "Lbend: " + template_funcs::toString(Lbend, 3) + " cannot be negative\n"; 

			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	} 
}

void XC_Bend::define_bend(bool loud)
{
	// define the bend coordinates based on the bend type
	// Repeat bend definition calculation until scale factors are identical to one
	// This will give the correct values for the bend whose parameters match those of the equivalent circle
	// R. Sheehan 8 - 8 - 2012
	
	if(bend_type == CC){
		
		define_bend_coords();
		
	}
	else{
		// Only need to re-scale in the LC and TC cases
		
		int n_iter = 1; 
		int max_iter = 30;
		
		double Lbend_old = 0.0; 

		//define_bend_coords(); 
		
		while(n_iter < max_iter){		
			Lbend_old = Lbend; 
				
			define_bend_coords();		
			re_scale_coords();
			
			// Convergence to the correct bend parameters is defined by correct length
			// Notice also that the scaling factors xs, ys get closer and closer to one
			// Would probably take more iterations for those to strictly converge to one because of division of floats. 
			if(fabs(Lbend - Lbend_old) < 1.0E-9){
				break; 
			}
			
			n_iter++;
			
			if(loud){
				std::cout<<"Iteration "<<n_iter<<", Lbend = "<<std::setprecision(15)<<Lbend<<"\n";
			}
		}
	
	}
}

void XC_Bend::define_bend_coords()
{
	// Define the coordinates of the waveguide centre
	// R. Sheehan 31 - 3 - 2009

	try{

		bool c1 = Npts > 1 ? true : false; 
		bool c2 = Lbend > 0.0 ? true : false; // L = R theta > 0

		if(c1 && c2){

			double ds = (Lbend/(Npts-1));
			double s;

			wg_X.resize(Npts+1);
			wg_Y.resize(Npts+1);

			for(int i=1;i<=Npts;i++){
				s = (i-1)*ds; 
				wg_X[i]=x_coord(s);
				wg_Y[i]=y_coord(s);
			}
		
		}
		else{
			std::string reason = "Error: void XC_Bend::define_bend_coords()\n"; 
			if(!c1) reason += "Npts: " + template_funcs::toString(Npts) + " is not defined\n"; 
			if(!c2) reason += "Lbend: " + template_funcs::toString(Lbend, 3) + " cannot be negative\n"; 

			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void XC_Bend::re_scale_coords(bool loud)
{
	// Scale the bend coordinates so that the end-points meet at the correct location
	// R. Sheehan 31 - 3 - 2009

	try{

		bool c1 = Npts > 1 ? true : false; 
		bool c2 = Teqc > 0.0 && Teqc <= PI  ? true : false; // bend angle 0 < T < \pi/2
		bool c3 = Reqc > 0.0 ? true : false; // bend radius > 0
		bool c4 = wg_X.size() > 0 ? true : false; 
		bool c5 = wg_Y.size() > 0 ? true : false; 
		bool c6 = wg_Xeqc.size() > 0 ? true : false; 
		bool c7 = wg_Yeqc.size() > 0 ? true : false; 
		bool c8 = (c1 && c2 && c3 && c4 && c5 && c6 && c7) ? true : false; 

		if( c8 ){

			if(bend_type == CC){
				// No need to re-scale for CC bends
				// that's the bend we're trying to match
		
				xs = ys = 1.0;
				Lbend = Leqc; 
			}
			else{
				if(Teqc==PI){
					xs=1.0;
					ys=(2.0*Reqc)/(wg_Y[Npts]);
				}
				else{
					xs=(wg_Xeqc[Npts]/(wg_X[Npts]));
					ys=(wg_Yeqc[Npts]/(wg_Y[Npts]));
				}
		
				if(loud){
					std::cout<<std::setprecision(15)<<"xs = "<<xs<<", ys = "<<ys<<"\n";
					std::cout<<std::setprecision(15)<<"______: x_last = "<<wg_Xeqc[Npts]<<", y_last = "<<wg_Yeqc[Npts]<<"\n";
					std::cout<<std::setprecision(15)<<"Before: x_last = "<<wg_X[Npts]<<", y_last = "<<wg_Y[Npts]<<"\n"<<"\n";
				}
		
				for(int i=1;i<=Npts;i++){
					wg_X[i]*=xs;
					wg_Y[i]*=ys;
				}
		
				//std::cout<<std::setprecision(15)<<"_After: x_last = "<<wg_X.last()<<", y_last = "<<wg_Y.last()<<"\n";

				//Need to determine the real-length of the curve here
				Lbend = 0.0; 
				for(int i=2;i<=Npts;i++){
					Lbend += sqrt( template_funcs::DSQR( wg_X[i] - wg_X[i-1] ) + template_funcs::DSQR( wg_Y[i] - wg_Y[i-1] ) );
				}
		
				Lbend_2 = 0.5*Lbend; 
				//std::cout<<"Lbend = "<<Lbend<<", Lbend_2 = "<<Lbend_2<<"\n";
		
				XC_coeffs(); // re-define bend coefficients for the corrected bend length
			}		
		}
		else{
			std::string reason = "Error: void XC_Bend::re_scale_coords(bool loud)\n"; 
			if(!c1) reason += "Npts: " + template_funcs::toString(Npts) + " is not defined\n"; 
			if(!c2) reason += "bend angle: " + template_funcs::toString(Teqc, 3) + " must be in range 0 < T <= pi\n"; 
			if(!c3) reason += "radius: " + template_funcs::toString(Reqc, 3) + " must be positive\n"; 
			if(!c4) reason += "wg_X contains no data\n"; 
			if(!c5) reason += "wg_Y contains no data\n"; 
			if(!c6) reason += "wg_Xeqc contains no data\n"; 
			if(!c7) reason += "wg_Yeqc contains no data\n"; 

			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void XC_Bend::re_scale_curvature()
{
	// Calculate the curvature profile of the bend from the scaled bend coordinates

	try{

		bool c1 = Npts > 0 ? true : false;
		bool c4 = wg_X.size() > 0 ? true : false; 
		bool c5 = wg_Y.size() > 0 ? true : false; 
		bool c6 = wg_curv.size() > 0 ? true : false; 

		if(c1 && c4 && c5 && c6){

			double Rc=0.0;
			double d1,d2,d3,A;	
			 
			point P1, P2, P3;
	
			wg_R.resize(Npts+1);

			for(int j=2;j<=Npts-1;j++){
		
				// Define the points along the bend
				P1.set_X(wg_X[j-1]); P1.set_Y(wg_Y[j-1]);	
				P2.set_X(wg_X[j]); P2.set_Y(wg_Y[j]);		
				P3.set_X(wg_X[j+1]); P3.set_Y(wg_Y[j+1]);

				// compute the distances between the points
				d1 = coord_geom::dist( P1, P2 );
				d2 = coord_geom::dist( P2, P3 );
				d3 = coord_geom::dist( P1, P3 );

				// Compute the area of the triangle formed by the points on the bend
				A = coord_geom::area( P1, P2, P3 );

				// Compute the radius of curvature at a given point along the bend
				Rc = ( d1 * d2 * d3 ) / ( 4.0 * A ); // Determines the radius of curvature at the jth point along the curve, excluding the endpoints
		
				wg_R[j] = Rc; 
			}

			//Rescale the curvature and the bending angle
			std::vector<double> tmpvec(Npts+1);

			for(int i=2;i<=Npts-1;i++){
				tmpvec[i] = (1.0/wg_R[i]);
			}

			double kmax,knmax;

			kmax = *std::max_element( wg_curv.begin(), wg_curv.end() ); 

			knmax = *std::max_element( tmpvec.begin(), tmpvec.end() ); 

			ks = (knmax/kmax); // set curvature scale factor

			// scale the curvature profile
			for(int i=1;i<=Npts;i++){
				wg_curv[i]*=ks;
			}
	
			//std::cout<<"ks = "<<ks<<"\n";
	
			// Add in the radius of curvature values at the start and end of the bend
			// Convention is that if kappa = 0 then R is set to 0, not infinity

			(fabs(wg_curv[1]) > EPS ? wg_R[1] = 1.0/wg_curv[1] : wg_R[1] = 0.0);
	
			(fabs(wg_curv[Npts]) > EPS ? wg_R[Npts] = 1.0/wg_curv[Npts] : wg_R[Npts] = 0.0);
	
			tmpvec.clear();		
		}
		else{
			std::string reason = "Error: void XC_Bend::define_eqc_coords()\n"; 
			if(!c1) reason += "Npts: " + template_funcs::toString(Npts) + " is not defined\n"; 
			if(!c4) reason += "wg_X contains no data\n"; 
			if(!c5) reason += "wg_Y contains no data\n"; 
			if(!c6) reason += "wg_curv contains no data\n"; 

			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void XC_Bend::re_scale_bend_angle()
{
	// Re-scale the bend angle ordinate values in case of mis-alignment
	// Ensure that the bend angle is scaled to the correct value at the end point
	// R. Sheehan 26 - 7 - 2013
	
	try{

		bool c1 = Npts > 0 ? true : false; 
		bool c2 = Teqc > 0.0 && Teqc <= PI_2  ? true : false; // bend angle 0 < T < \pi/2
		bool c3 = wg_bend_ang.size() > 0 ? true : false ; 

		if(c1 && c2 && c3){

			// Does this mean that bend angle is scaled only for LC bends? 
			// Does it not need to be scaled for other bend types? 
			// R. Sheehan 9 - 9 - 2016
			if(bend_type == LC){		
				double angle_scale = (Teqc/(wg_bend_ang[Npts]));
		
				for(int i=1; i<=Npts; i++){
					wg_bend_ang[i]*=angle_scale; 
				}
			}
		
		}
		else{
			std::string reason = "Error: void XC_Bend::re_scale_bend_angle()\n"; 
			if(!c1) reason += "Npts: " + template_funcs::toString(Npts) + " is not defined\n"; 
			if(!c2) reason += "bend angle: " + template_funcs::toString(Teqc, 3) + " must be in range 0 < T < pi/2\n"; 
			if(!c3) reason += "wg_bend_ang contains no data\n"; 

			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void XC_Bend::define_waveguide_edges()
{
	// Define the points that make the upper and lower edges of the waveguide
	// input is the set of points that make up the waveguide centre
	// R. Sheehan 8 - 8 - 2012

	try{

		bool c1 = Npts > 0 ? true : false; 
		bool c4 = wg_X.size() > 0 ? true : false; 
		bool c5 = wg_Y.size() > 0 ? true : false;
		bool c6 = wg_bend_ang.size() > 0 ? true : false ; 

		if(c1 && c4 && c5 && c6){

			double m2,cc,dk2,x1,x2,tf,y1,y2,d1,d2;
			double A,B,C;
			point p1, p2, pc;

			pc.set_X(0.0); pc.set_Y(Reqc); // Define the coordinates of the centre point of the equivalent circle

			// Firstly define the coordinates for the upper and lower edges
			for(int i=1;i<=Npts;i++){
				if(wg_Y[i]<EPS){
					wg_Xl.push_back(wg_X[i]);
					wg_Xu.push_back(wg_X[i]);
					wg_Yl.push_back(-1.0*t2);
					wg_Yu.push_back(t2);
				}
				else{
					//Calculate the parameters
			
					tf = tan(wg_bend_ang[i]);
			
					if(fabs(tf)<EPS){
						wg_Xl.push_back(wg_X[i]);
						wg_Xu.push_back(wg_X[i]);
						wg_Yl.push_back(-1.0*t2);
						wg_Yu.push_back(t2);
					}
					else if(fabs(wg_bend_ang[i]-PI_2)<EPS){
						wg_Xl.push_back(wg_X[i]+t2);
						wg_Xu.push_back(wg_X[i]-t2);
						wg_Yl.push_back(wg_Y[i]);
						wg_Yu.push_back(wg_Y[i]);
					}
					else{
						m2=-1.0/tf;
						dk2=(template_funcs::DSQR(wg_X[i])+template_funcs::DSQR(wg_Y[i]));
						cc=(wg_Y[i]-wg_X[i]*m2);
				
						//Form the quadratic
						A=(1.0+template_funcs::DSQR(m2));
						B=(2.0*m2*cc-2.0*wg_X[i]-2.0*m2*wg_Y[i]);
						C=(dk2+template_funcs::DSQR(cc)-template_funcs::DSQR(t2)-2.0*cc*wg_Y[i]);

						//Find x upper and x lower
						//x increases for lower => x1
						//x decreases for upper => x2
						template_funcs::Quad_Solve(A,B,C,x1,x2);

						y1=m2*x1+cc;
						y2=m2*x2+cc;

						p1.set_X(x1); p1.set_Y(y1);
						p2.set_X(x2); p2.set_Y(y2);

						d1=coord_geom::dist(p1,pc);
						d2=coord_geom::dist(p2,pc);

						if(d1>d2){
							wg_Xl.push_back(x1);
							wg_Xu.push_back(x2);
							wg_Yl.push_back(y1);
							wg_Yu.push_back(y2);
						}
						else{
							wg_Xl.push_back(x2);
							wg_Xu.push_back(x1);
							wg_Yl.push_back(y2);
							wg_Yu.push_back(y1);
						}
					}
				}
			}	
		
		}
		else{
			std::string reason = "Error: void XC_Bend::define_waveguide_edges()\n"; 
			if(!c1) reason += "Npts: " + template_funcs::toString(Npts) + " is not defined\n"; 
			if(!c4) reason += "wg_X contains no data\n"; 
			if(!c5) reason += "wg_Y contains no data\n"; 
			if(!c6) reason += "wg_bend_ang contains no data\n"; 

			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
	
	

}

void XC_Bend::define_waveguide_outline()
{
	// Define the outline of the waveguide based on the upper and lower edges
	// R. Sheehan 8 - 8 - 2012

	try{

		bool c1 = Npts > 0 ? true : false; 
		bool c4 = wg_Xl.size() > 0 ? true : false; 
		bool c5 = wg_Yl.size() > 0 ? true : false;
		bool c6 = wg_Xu.size() > 0 ? true : false; 
		bool c7 = wg_Yu.size() > 0 ? true : false;

		if(c1 && c4 && c5 && c6 && c7){

			//To account for the sense of the angle change the sign on the y-values
			//Similarly for when the radius is negative
			if(signT == 1){
				// Add overlayer section if necessary
				if(Lovertotal > 0.0){
					waveguide.push_back(point(-Lovertotal,t2));
				}

				// Define upper portion of curve
				for(int i=0;i<Npts;i++){
					waveguide.push_back(point(wg_Xu[i],wg_Yu[i]));
				}

				if(Lovertotal > 0.0){
					// Add overlayer section if necessary
					double dy=fabs(wg_Yu[Npts-1]-wg_Yu[Npts-2]);
					double dx=fabs(wg_Xu[Npts-1]-wg_Xu[Npts-2]);
					double th=atan(dy/dx);

					waveguide.push_back(point(wg_Xu[Npts-1]+Lovertotal*cos(th),wg_Yu[Npts-1]+Lovertotal*sin(th)));

					waveguide.push_back(point(wg_Xl[Npts-1]+Lovertotal*cos(th),wg_Yl[Npts-1]+Lovertotal*sin(th)));
				}
		
				// Define lower portion of curve
				for(int i=(Npts-1);i>=0;i--){
					waveguide.push_back(point(wg_Xl[i],wg_Yl[i]));
				}

				if(Lovertotal > 0.0){
					// Close curve
					waveguide.push_back(point(-Lovertotal,-t2));
					waveguide.push_back(point(-Lovertotal,t2));		
				}
				else{
					// Close curve
					waveguide.push_back(point(wg_Xu[0],wg_Yu[0]));
				}
			}
			else{
				// Change the sign of the upper and lower curves
				// x values stay the same
				for(int i=0;i<Npts;i++){
					wg_Yl[i]*=-1.0;
					wg_Yu[i]*=-1.0;
				}

				// Add overlayer section if necessary
				if(Lovertotal > 0.0){
					waveguide.push_back(point(-Lovertotal,t2));
				}
        
				// Define upper portion of curve
				for(int i=0;i<Npts;i++){
					waveguide.push_back(point(wg_Xl[i],wg_Yl[i]));
				}

				// Add overlayer section if necessary
				if(Lovertotal > 0.0){
					double dy=fabs(wg_Yl[Npts-1]-wg_Yl[Npts-2]);
					double dx=fabs(wg_Xl[Npts-1]-wg_Xl[Npts-2]);
					double th=atan(dy/dx);

					waveguide.push_back(point(wg_Xl[Npts-1]+Lovertotal*cos(th),wg_Yl[Npts-1]-Lovertotal*sin(th)));
					waveguide.push_back(point(wg_Xu[Npts-1]+Lovertotal*cos(th),wg_Yu[Npts-1]-Lovertotal*sin(th)));
				}

				// Add lower portion of curve
				for(int i=(Npts-1);i>=0;i--){
					waveguide.push_back(point(wg_Xu[i],wg_Yu[i]));
				}

				if(Lovertotal > 0.0){
					// Close curve
					waveguide.push_back(point(-Lovertotal,-t2));
					waveguide.push_back(point(-Lovertotal,t2));		
				}
				else{
					// Close curve
					waveguide.push_back(point(wg_Xl[0],wg_Yl[0]));
				}
			}    
		}
		else{
			std::string reason = "Error: void XC_Bend::define_waveguide_outline()\n"; 
			if(!c1) reason += "Npts: " + template_funcs::toString(Npts) + " is not defined\n"; 
			if(!c4) reason += "wg_Xl contains no data\n"; 
			if(!c5) reason += "wg_Yl contains no data\n"; 
			if(!c6) reason += "wg_Xu contains no data\n"; 
			if(!c7) reason += "wg_Yu contains no data\n"; 

			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void XC_Bend::output_curves(bool loud)
{
	// Send the data to files
	// R. Sheehan 8 - 8 - 2012

	try{

		bool c1 = Npts > 0 ? true : false; 
		bool c2 = Teqc > 0.0 && Teqc <= PI_2  ? true : false; // bend angle 0 < T < \pi/2
		bool c3 = Lbend > 0.0 ? true : false; // bend angle 0 < T < \pi/2
		bool c4 = wg_Xeqc.size() > 0 ? true : false; 
		bool c5 = wg_Yeqc.size() > 0 ? true : false; 
		bool c6 = wg_curv.size() > 0 ? true : false; 
		bool c7 = wg_bend_ang.size() > 0 ? true : false; 
		bool c8 = wg_X.size() > 0 ? true : false; 
		bool c9 = wg_Y.size() > 0 ? true : false; 
		bool c10 = wg_Xl.size() > 0 ? true : false; 
		bool c11 = wg_Yl.size() > 0 ? true : false; 
		bool c12 = wg_Xu.size() > 0 ? true : false; 
		bool c13 = wg_Yu.size() > 0 ? true : false; 
		bool c14 = waveguide.size() > 0 ? true : false; 
		bool c15 = (c1 && c2 && c3 && c4 && c5 && c6 && c7 && c8 && c9 && c10 && c11 && c12 && c13 && c14) ? true : false; 

		if(c15){

			double s,ds;
			ds=(Lbend/(Npts-1));

			if(loud){
				std::cout<<"\n";
				std::cout<<"The length of the bend is "<<Lbend<<" microns\n";
				std::cout<<"\n";
			}

			std::ofstream outp;
			std::string info; 
			std::string filename; 

			// info std::string for filenames
			info = "_R_"+template_funcs::toString(Reqc,2)+"_T_"+template_funcs::toString(Teqc,4)+"_W_"+template_funcs::toString(width,2)+dottxt; 

			// Equivalent circle data
			filename = "Equivalent_Circle_Centre"+info; 

			outp.open(filename.c_str(),std::ios_base::out|std::ios_base::trunc);

			if(outp.is_open()){
				for(int i=1;i<=Npts;i++){
					outp<<std::setprecision(15)<<wg_Xeqc[i]<<" , "<<wg_Yeqc[i]<<"\n";
				}

				outp.close();
			}

			// curvature versus path-length data
			if(bend_type == TC){
				filename = "TC_Curvature_"+template_funcs::toString(TC_type)+info; 
			}
			else if(bend_type == LC){
				filename = "LC_Curvature"+info;
			}
			else if(bend_type == QC){
				filename = "QC_Curvature"+info;
			}
			else{
				filename = "CC_Curvature"+info;
			}

			outp.open(filename.c_str(),std::ios_base::out|std::ios_base::trunc);

			if(outp.is_open()){
				for(int i=1;i<=Npts;i++){
					s = (i-1)*ds; 
					//outp<<std::setprecision(15)<<s<<" , "<<wg_curv[i]<<" , "<<curvature(s)<<"\n";
					outp<<std::setprecision(15)<<s<<" , "<<wg_curv[i]<<"\n";
				}

				outp.close();
			}

			// Bend Angle Data
			if(bend_type == TC){
				filename = "TC_Bend_Angle_"+template_funcs::toString(TC_type)+info; 
			}
			else if(bend_type == LC){
				filename = "LC_Bend_Angle"+info; 
			}
			else if(bend_type == QC){
				filename = "QC_Bend_Angle"+info; 
			}
			else{
				filename = "CC_Bend_Angle"+info; 
			}

			outp.open(filename.c_str(),std::ios_base::out|std::ios_base::trunc);

			if(outp.is_open()){
				for(int i=1;i<=Npts;i++){
					s = (i-1)*ds; 
					//outp<<std::setprecision(15)<<s<<" , "<<wg_bend_ang[i]<<" , "<<bend_angle(s)<<"\n";
					outp<<std::setprecision(15)<<s<<" , "<<wg_bend_ang[i]<<"\n";
				}

				outp.close();
			}

			// Bend Centre Data
			if(bend_type == TC){
				filename = "TC_Bend_Centre_"+template_funcs::toString(TC_type)+info; 
			}
			else if(bend_type == LC){
				filename = "LC_Bend_Centre"+info; 
			}
			else if(bend_type == QC){
				filename = "QC_Bend_Centre"+info; 
			}
			else{
				filename = "CC_Bend_Centre"+info; 
			}

			outp.open(filename.c_str(),std::ios_base::out|std::ios_base::trunc);

			if(outp.is_open()){
				for(int i=1;i<=Npts;i++){
					outp<<std::setprecision(15)<<wg_X[i]<<" , "<<wg_Y[i]<<"\n";
				}

				outp.close();
			}

			// Lower Edge Data
			if(bend_type == TC){
				filename = "TC_Lower_Edge_"+template_funcs::toString(TC_type)+info; 
			}
			else if(bend_type == LC){
				filename = "LC_Lower_Edge"+info; 
			}
			else if(bend_type == QC){
				filename = "QC_Lower_Edge"+info; 
			}
			else{
				filename = "CC_Lower_Edge"+info; 
			}

			outp.open(filename.c_str(),std::ios_base::out|std::ios_base::trunc);

			if(outp.is_open()){
				for(int i=0;i<static_cast<int>(wg_Xl.size());i++){
					outp<<std::setprecision(15)<<wg_Xl[i]<<" , "<<wg_Yl[i]<<"\n";
				}

				outp.close();
			}

			// Upper Edge Data
			if(bend_type == TC){
				filename = "TC_Upper_Edge_"+template_funcs::toString(TC_type)+info; 
			}
			else if(bend_type == LC){
				filename = "LC_Upper_Edge"+info; 
			}
			else if(bend_type == QC){
				filename = "QC_Upper_Edge"+info; 
			}
			else{
				filename = "CC_Upper_Edge"+info; 
			}

			outp.open(filename.c_str(),std::ios_base::out|std::ios_base::trunc);

			if(outp.is_open()){
				for(int i=0;i<static_cast<int>(wg_Xu.size());i++){
					outp<<std::setprecision(15)<<wg_Xu[i]<<" , "<<wg_Yu[i]<<"\n";
				}

				outp.close();
			}

			// Waveguide Outline Data
			if(bend_type == TC){
				filename = "TC_Waveguide_"+template_funcs::toString(TC_type)+info; 
			}
			else if(bend_type == LC){
				filename = "LC_Waveguide"+info; 
			}
			else if(bend_type == QC){
				filename = "QC_Waveguide"+info; 
			}
			else{
				filename = "CC_Waveguide"+info; 
			}

			outp.open(filename.c_str(),std::ios_base::out|std::ios_base::trunc);

			if(outp.is_open()){
				for(int i=0;i<static_cast<int>(waveguide.size());i++){
					outp<<std::setprecision(15)<<waveguide[i].get_X()<<" , "<<waveguide[i].get_Y()<<"\n";
				}

				outp.close();
			}
		
		}
		else{
			std::string reason = "Error: void XC_Bend::output_curves(bool loud)\n"; 
			if(!c1) reason += "Npts: " + template_funcs::toString(Npts) + " is not defined\n"; 
			if(!c2) reason += "bend angle: " + template_funcs::toString(Teqc, 3) + " must be in range 0 < T < pi/2\n"; 
			if(!c3) reason += "Lbend: " + template_funcs::toString(Lbend, 3) + " must be positive\n"; 
			if(!c4) reason += "wg_Xeqc contains no data\n"; 
			if(!c5) reason += "wg_Yeqc contains no data\n"; 
			if(!c6) reason += "wg_curv contains no data\n"; 
			if(!c7) reason += "wg_bend_ang contains no data\n"; 
			if(!c8) reason += "wg_X contains no data\n"; 
			if(!c9) reason += "wg_Y contains no data\n"; 
			if(!c10) reason += "wg_Xl contains no data\n"; 
			if(!c11) reason += "wg_Yl contains no data\n"; 
			if(!c12) reason += "wg_Xu contains no data\n"; 
			if(!c13) reason += "wg_Yu contains no data\n"; 
			if(!c14) reason += "waveguide contains no data\n"; 

			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void XC_Bend::clear_memory()
{
	// Remove the memory used by the object
	
	wg_curv.clear(); 
	wg_R.clear(); 
	wg_bend_ang.clear(); 
	wg_Xeqc.clear(); 
	wg_Yeqc.clear(); 
	wg_X.clear();
	wg_Y.clear();

	wg_Xl.clear();
	wg_Yl.clear();
	wg_Xu.clear();
	wg_Yu.clear();

	waveguide.clear();
}

// Curvature versus pathlength profiles
double XC_Bend::curvature(double s)
{
	// Return the bend curvature
	// R. Sheehan 7 - 8 - 2012

	try{

		if(coeffs_defined){

			if(s < 0.0 || s > Lbend){
				return 0.0;
			}
			else{
				if(bend_type == LC){
					return LC_curvature(s);
				}
				else if(bend_type == TC){
					return TC_curvature(s);
				}
				else if(bend_type == QC){
					return QC_curvature(s); 
				}
				else{
					return CC_curvature(s);
				}
			}
		}
		else{
			std::string reason = "Error: double XC_Bend::curvature_coord(double s)\n"; 
			reason += "Calculation coefficients have not been defined\n"; 
			throw std::invalid_argument(reason); 
			return 0.0; 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double XC_Bend::CC_curvature(double x)
{
	// CC bends have constant curvature
	// R. Sheehan 26 - 3 - 2009
	
	return keqc; 
}

double XC_Bend::TC_curvature(double s)
{
	// Define the curvature at a point along a TC bend
	// R. Sheehan 26 - 3 - 2009
	
	if(s>=0.0 && s<=sigma){
		return alpha*s;
	}
	else if(s>sigma && s<=nu){
		return kt;
	}
	else{
		return alpha*(Lbend-s);
	}
}

double XC_Bend::LC_curvature(double s)
{
	// Define the curvature at a point along a LC bend
	// R. Sheehan 26 - 3 - 2009
	
	if(s>=0.0 && s<=Lbend_2){
		return (alpha*s);
	}
	else{
		return (alpha*(Lbend-s));
	}
}

double XC_Bend::QC_curvature(double s)
{
	// Define the curvature at a point along a QC bend
	// R. Sheehan 15 - 4 - 2014

	if(s>=0.0 && s<=Lbend){
		return alpha*s*(Lbend-s); 
	}
	else{
		return 0.0; 
	}
}

// Bend-angle versus pathlength profiles
double XC_Bend::bend_angle(double s)
{
	// Return the bend angle
	// R. Sheehan 7 - 8 - 2012

	try{

		if(coeffs_defined){
			if(bend_type == LC){
				return LC_bend_angle(s);
			}
			else if(bend_type == TC){
				return TC_bend_angle(s); 
			}
			else if(bend_type == QC){
				return QC_bend_angle(s); 
			}
			else{
				return CC_bend_angle(s); 
			}		
		}
		else{
			std::string reason = "Error: double XC_Bend::x_coord(double s)\n"; 
			reason += "Calculation coefficients have not been defined\n"; 
			throw std::invalid_argument(reason); 
			return 0.0; 
		}	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double XC_Bend::CC_bend_angle(double s)
{
	// Define the bend-angle at a point along a CC bend
	// R. Sheehan 26 - 3 - 2009
	
	return keqc*s;
}

double XC_Bend::TC_bend_angle(double s)
{
	// Define the bend-angle at a point along a TC bend
	// R. Sheehan 26 - 3 - 2009
	
	if(s>=0.0 && s<=sigma){
		return 0.5*alpha*template_funcs::DSQR(s);
	}
	else if(s>sigma && s<=nu){
		return kt*(s-0.5*sigma);
	}
	else{
		return -0.5*alpha*(template_funcs::DSQR(s)-2.0*Lbend*s+template_funcs::DSQR(Lbend))+kt*(Lbend-sigma);
	}
}

double XC_Bend::LC_bend_angle(double s)
{
	// Define the bend-angle at a point along a LC bend
	// R. Sheehan 26 - 3 - 2009
	
	if(s>=0.0 && s<=Lbend_2){
		return 0.5*alpha*template_funcs::DSQR(s);
	}
	else{
		//return alpha*(Lbend*s-0.5*template_funcs::DSQR(s)-(3.0/8.0)*template_funcs::DSQR(Lbend))+(Teqc/2.0);
		return alpha*(Lbend*s-0.5*template_funcs::DSQR(s)-(3.0/8.0)*template_funcs::DSQR(Lbend))+(Lbend/(2.0*Reqc));
	}
}

double XC_Bend::QC_bend_angle(double s)
{
	// Define the bend-angle at a point along a QC bend
	// R. Sheehan 15 - 4 - 2014

	if(s>=0.0 && s<=Lbend){
		return (alpha * s * ( (Lbend_2*s) - (template_funcs::DSQR(s)/3.0) ) ); 
	}
	else{
		return 0.0; 
	}
}

// Waveguide centre coordinates
double XC_Bend::x_coord(double s)
{
	// Return the bend coordinate in x
	// R. Sheehan 7 - 8 - 2012

	try{

		if(coeffs_defined){
			if(bend_type == LC){
				return LC_x_coord(s);
			}
			else if(bend_type == TC){
				return TC_x_coord(s); 
			}
			else if(bend_type == QC){
				return QC_x_coord(s); 
			}
			else{
				// For the CC bend you need to pass an angle value not a path-length value
				// R. Sheehan 26 - 7 - 2013
				s/=Reqc; 
				return CC_x_coord(s); 
			}
		}
		else{
			std::string reason = "Error: double XC_Bend::x_coord(double s)\n"; 
			reason += "Calculation coefficients have not been defined\n"; 
			throw std::invalid_argument(reason); 
			return 0.0; 
		}	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double XC_Bend::CC_x_coord(double s)
{
	// x-coordinate of a point on a CC bend
	// R. Sheehan 26 - 3 - 2009
	
	// Updated R. Sheehan 26 - 7 - 2013
	// Remember that you're supposed to be passing angluar value to this function
	
	return Reqc*sin(s); 
}

double XC_Bend::TC_x_coord(double s)
{
	// x-coordinate of a point on a TC bend
	// R. Sheehan 26 - 3 - 2009
	
	double pos;
	double C,S;

	if(s>=0.0 && s<=sigma){
		pos=c2*s;
		special::fresnel(pos, &S, &C); 
		return (c1*C);
	}
	else if(s>sigma && s<=nu){
		return ((2.0/kt)*sin(0.5*kt*(s-sigma))*cos(0.5*kt*s)+c1*FCc3);
	}
	else{
		pos=c2*(s-Lbend);
		special::fresnel(pos,&S,&C);
		return (c1*(cos(c4)*(C+FCc3)+sin(c4)*(S+FSc3)+FCc3)+(1.0/kt)*(sin(c5)-sin(c6)));
	}
}

double XC_Bend::LC_x_coord(double s)
{
	// x-coordinate of a point on a LC bend
	// R. Sheehan 26 - 3 - 2009
	
	double pos;
	double C,S;

	if(s>=0.0 && s<=Lbend_2){
		pos=c2*s;
		special::fresnel(pos, &S, &C); 
		return (c1*C);
	}
	else{
		pos=c2*(s-Lbend);
		special::fresnel(pos,&S,&C);//Take advantage of the fact that this function computes both integrals at once
		return (c1*(cos(Teqc)*(C+FCc3)+sin(Teqc)*(S+FSc3)+FCc3));
	}	
}

double XC_Bend::QC_x_coord(double s)
{
	// x-coordinate on a QC bend
	// Need to compute a sum over a series of hypergeometric functions
	// R. Sheehan 15 - 4 - 2014

	int k, nterms,p1,p2,p3; 
	double a, b, c, z, F, dF; // arguments to the hypergeometric function
	double term, series, t1, t2, t3, t4, t5; 

	z = nu*s; 

	nterms = 10; // this has been shown to provide sufficient accuracy by A. Walsh
	
	series = 0.0; 
	
	for(k = 0; k<=nterms; k++){
		
		term = 0.0; 

		// powers
		p1 = 2*k; p2 = 4*k+1; p3 = 6*k+1; 

		// value of the hypergeometric function stored in F
		a = static_cast<double>(p2); // 1+4k
		b = static_cast<double>(-p1); // -2k
		c = static_cast<double>(2+4*k); // 4k+2
		special::two_F_one(a, b, c, z, F, dF); // compute {2}_F_{1}

		// terms in the sum
		t1 = ( 1.0 / special::factorial(p1) ); // 1/(2k)! 
		t2 = pow(alpha,p1); // alpha^{2k}
		t3 = 1.0 / ( a );  // 1 / 4k+1
		t4 = pow(z,p2); // z^{4k+1}
		t5 = pow(g,p3) / pow(3.0,p1); // g^{6k+1} / 3^{2k} , g = 1 / nu

		term = ( t1 * t2 * t3 * t4 * t5 * F ); 

		if( k % 2 == 1){
			term *= -1.0; // multiply by (-1)^{k}
		}

		series += term; 
	}

	return series; 
}

double XC_Bend::y_coord(double s)
{
	// Return the bend coordinate in x
	// R. Sheehan 7 - 8 - 2012

	try{

		if(coeffs_defined){
			if(bend_type == LC){
				return LC_y_coord(s);
			}
			else if(bend_type == TC){
				return TC_y_coord(s); 
			}
			else if(bend_type == QC){
				return QC_y_coord(s); 
			}
			else{
				// For the CC bend you need to pass an angle value not a path-length value
				// R. Sheehan 26 - 7 - 2013
				s/=Reqc; 
				return CC_y_coord(s); 
			}		
		}
		else{
			std::string reason = "Error: double XC_Bend::y_coord(double s)\n"; 
			reason += "Calculation coefficients have not been defined\n"; 
			throw std::invalid_argument(reason); 
			return 0.0; 
		}	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double XC_Bend::CC_y_coord(double s)
{
	// y-coordinate of a point on a CC bend
	// R. Sheehan 26 - 3 - 2009
	
	// Updated R. Sheehan 26 - 7 - 2013
	// Remember that you're supposed to be passing angluar value to this function
	
	return Reqc*(1.0-cos(s)); 
}

double XC_Bend::TC_y_coord(double s)
{
	// y-coordinate of a point on a TC bend
	// R. Sheehan 26 - 3 - 2009
	
	double pos;
	double C,S;

	if(s>=0.0 && s<=sigma){
		pos=c2*s;
		special::fresnel(pos, &S, &C); 
		return (c1*S);
	}
	else if(s>sigma && s<=nu){
		return ((2.0/kt)*sin(0.5*kt*(s-sigma))*sin(0.5*kt*s)+c1*FSc3);
	}
	else{
		pos=c2*(s-Lbend);
		special::fresnel(pos,&S,&C);
		return (c1*(sin(c4)*(C+FCc3)-cos(c4)*(S+FSc3)+FSc3)+(1.0/kt)*(cos(c6)-cos(c5)));
	}
}

double XC_Bend::LC_y_coord(double s)
{
	// y-coordinate of a point on a LC bend
	// R. Sheehan 26 - 3 - 2009
	
	double pos;
	double C,S;

	if(s>=0.0 && s<=Lbend_2){
		pos=c2*s;
		special::fresnel(pos, &S, &C); 
		return (c1*S);
	}
	else{
		pos=c2*(s-Lbend);
		special::fresnel(pos,&S,&C);//Take advantage of the fact that this function computes both integrals at once
		return (c1*(sin(Teqc)*(C+FCc3)-cos(Teqc)*(S+FSc3)+FSc3));
	}
}

double XC_Bend::QC_y_coord(double s)
{
	// y-coordinate on a QC bend
	// Need to compute a sum over a series of hypergeometric functions
	// R. Sheehan 15 - 4 - 2014

	int k, nterms,p1,p2; 
	double a, b, c, z, F, dF; // arguments to the hypergeometric function
	double term, series, t1, t2, t3, t4, t5; 

	z = nu*s; 

	nterms = 10; // this has been shown to provide sufficient accuracy by A. Walsh
	
	series = 0.0; 
	
	for(k = 0; k<=nterms; k++){
		
		term = 0.0; 

		// powers
		p1 = 2*k+1; p2 = 4*k+3; 

		// value of the hypergeometric function stored in F
		a = static_cast<double>(-p1); // -2k-1
		b = static_cast<double>(p2); // 4k+3
		c = 1.0+b; // 4k+4
		special::two_F_one(a, b, c, z, F, dF); // compute {2}_F_{1}

		// terms in the sum
		t1 = ( 1.0 / special::factorial(p1) ); // 1/(2k+1)! 
		t2 = pow(alpha,p1); // alpha^{2k+1}
		t3 = 1.0 / ( static_cast<double>(8*k+6) );  // 1 / 8k+6
		t4 = ( pow(Lbend,p1) / pow(4.0,k) ); // L^{2k+1} / 4^{k}
		t5 = pow(s, p2); // s^{4k+3}

		term = ( t1 * t2 * t3 * t4 * t5 * F ); 

		if( k % 2 == 1){
			term *= -1.0; // multiply by (-1)^{k}
		}

		series += term; 
	}

	return series; 
}