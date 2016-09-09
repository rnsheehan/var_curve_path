#ifndef BEND_PATHS_H
#define BEND_PATHS_H

// Bend description class incorporating bends having different curvature types
// First implemented separately for TC, LC and NLC bends March 2009
// Mathematics behind curve definitions to be found in Notebook 984, and probably PhD Thesis
// Units of bend radius assumed to be microns
// => units of curvature = inverse microns
// => no factor of 1000 floating about the place generally confusing the hell out of everybody
// => kappa[mm^{-1}] = 1000 / R[um] for your erudition
// Updated to a single object 
// R. Sheehan 7 - 8 - 2012

// Updated to include a quadratic curvature bend
// R. Sheehan 15 - 4 - 2014

class XC_Bend{
public:
	// Constructors
	XC_Bend();
	XC_Bend(int type, double W, double R, double T);
	~XC_Bend(); // Deconstructor

	// Methods
	
	// Getters
	int get_bend_type();
	int get_Npts();

	std::string get_bend_name(); 
	
	double get_bend_length();
	double get_max_curvature(); 
	
	// Setters
	void set_bend_type(int val);
	void set_bend_rad(double val); 
	//void set_n_pts(int val);
	
	void set_params(int type, double W, double R, double T);
	
	// Bend Definition Methods
	void define_eqc_coords(); // define the coordinates of the equivalent circle
	void define_curvature(); // define the bend curvature based on the bend type
	void define_bend_angle(); // define the bend angle based on the bend type
	void define_bend(bool loud = false); // iteratively define the wg bend centre coordinates
	void define_bend_coords(); // define the bend coordinates based on the bend type
	//void define_waveguide_edges(std::vector<double> &X_c, std::vector<double> &Y_c); // Define the upper and lower edges of the waveguide
	void define_waveguide_edges(); // Define the upper and lower edges of the waveguide
	void define_waveguide_outline(); // Define the outline of the waveguide based on the bend type
	void output_curves(bool loud = false); // send the data to the files
	void clear_memory(); // de-allocate memory used

	// The user does not need access to these methods
private:
	void XC_coeffs(); // assign values to the XC bend coeffs
	void LC_coeffs(); // assign values to the LC bend coeffs
	void TC_coeffs(); // assign values to the TC bend coeffs	
	void QC_coeffs(); // assign values to the QC bend coeffs

	void re_scale_coords(bool loud = false); // Scale the bend coordinates so that the end-points meet at the correct location
	void re_scale_curvature(); // Scale the bend curvature so that the end-points meet at the correct location
	void re_scale_bend_angle(); // Scale the bend angle so that the end-points have the correct value

	// Curvature versus pathlength profiles
	double curvature(double s);
	double CC_curvature(double s);
	double TC_curvature(double s);
	double LC_curvature(double s); 
	double QC_curvature(double s);

	// Bend-angle versus pathlength profiles
	double bend_angle(double s);
	double CC_bend_angle(double s);
	double TC_bend_angle(double s);
	double LC_bend_angle(double s);
	double QC_bend_angle(double s);

	// Waveguide centre coordinates
	double x_coord(double s);	
	double y_coord(double s);

	double CC_x_coord(double s);
	double TC_x_coord(double s);
	double LC_x_coord(double s);
	double QC_x_coord(double s);

	double CC_y_coord(double s);
	double TC_y_coord(double s);
	double LC_y_coord(double s);
	double QC_y_coord(double s);

private:
	// Members

	bool coeffs_defined; 

	int bend_type; // parameter that determines bend curvature profile
	int TC_type; // this distinguishes between three different cases of TC bend
	int Npts; // num point used to plot bend bath
	int signT; // direction of the bend 
	
	// Bend parameters common to all bends
	double width; // waveguide width
	double t2; // half waveguide width
	double Reqc; // equivalent circle radius
	double keqc; // equivalent circle curvature
	double kmax; // maximum curvature of given bend
	double Teqc; // bend angle
	double Leqc; // equivalent circle length
	double Lbend; // actual length of the bend
	double Lbend_2; // half the bend length

	double Loverlap; // This is the amount by which the two concatenated waveguides must overlap in PICDraw.cpp
	double Loverlayer; // This is the length that the overlayer protrudes beyond the waveguide in PICDraw.cpp
	double Lovertotal; // Loverlap + Loverlayer
	
	// Bend scaling factors
	double xs; 
	double ys;
	double sxc;
	double syc;
	double kt;
	double ks; 
	
	// Parameters for LC & TC & QC bends
	double alpha;
	double sigma;
	double nu;
	double r; 
	double g;
	double f;
	double c1;
	double c2;
	double c3;
	double c4;
	double c5;
	double c6;
	double FCc3;
	double FSc3;	
	
	// Various containers

	std::vector<int> bend_types; // vector holding different bend type identifiers

	std::vector<double> wg_curv; // curvature along the bend
	std::vector<double> wg_R; // Radius of curvature along the curve	
	std::vector<double> wg_bend_ang; // bending angle along the bend
	std::vector<double> wg_Xeqc; // x coordinates of an equivalent waveguide based on a circle
	std::vector<double> wg_Yeqc; //y coordinates of an equivalent waveguide based on a circle
	std::vector<double> wg_X; // x coordinates of the waveguide centre
	std::vector<double> wg_Y; // y coordinates of the waveguide centre

	std::vector<double> wg_Xl; // x coordinates of lower branch of waveguide
	std::vector<double> wg_Yl; // y coordinates of lower branch of waveguide
	std::vector<double> wg_Xu; // x coordinates of upper branch of waveguide
	std::vector<double> wg_Yu; // y coordinates of upper branch of waveguide

	std::vector<point> waveguide; // points describing the outline of the waveguide bend
};

#endif