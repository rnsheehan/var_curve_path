#ifndef ATTACH_H
#include "Attach.h"
#endif

// sample call: Variable_Curvature_Paths 1003 2.5 500 PI/3

int main(int argc, char *argv[])
{
	try{
		
		if(argc >= 4){
			
			// List off the input parameters
			// Program needs 5 or more parameters to run, remember that the name of the program is also considered a parameter
			// argv[0] = program name
			// argv[1] = bend type
			// argv[2] = waveguide width
			// argv[3] = bend radius
			// argv[4] = bend angle

			std::cout<<argc-1<<" parameters were input into the program\n"; 
			for(int count = 1; count < argc; count++){
				std::cout<<"argv["<<count<<"] = "<<argv[count]<<"\n"; 
			}
			std::cout<<"\n";

			std::cout<<"Calculation underway\n"; 
			
			// Define input parameters
			int bnd_type = atoi( argv[1] ); 
			double wg_width = atof( argv[2] );
			double wg_rad = atof( argv[3] ); 
			double wg_ang = atof( argv[4] ); 

			XC_Bend the_bend(bnd_type, wg_width, wg_rad, wg_ang);

			the_bend.define_eqc_coords(); // strictly necessary for computing the correct curvature profile
			the_bend.define_curvature(); // strictly necessary for computing the correct curvature profile
			the_bend.define_bend_angle(); // technically not needed to compute bend path shape
			the_bend.define_bend(); // strictly necessary for computing the correct curvature profile
	
			// These are necessary for defining the waveguide outline
			the_bend.define_waveguide_edges(); 

			the_bend.define_waveguide_outline();

			the_bend.output_curves();

			the_bend.clear_memory();

			std::cout<<"\nCalculation complete\n"; 
		}
		else{
			std::string reason = "Error: Slab_WG_Slv\n";
			reason += "Insufficient number of input arguments\n"; 

			throw std::invalid_argument(reason); 
		}

	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

	return 0; 
}