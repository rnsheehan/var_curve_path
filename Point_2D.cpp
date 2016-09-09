#ifndef ATTACH_H
#include "Attach.h"
#endif

point::point()
{
	// Default constructor

	quadrant = 0; 

	x = y = 0.0; 
}

point::point(double xpos, double ypos)
{
	// Constructor

	set_XY(xpos, ypos); 
}

void point::set_XY(double xpos, double ypos)
{
	x = xpos; 
	y = ypos; 

	set_quad(); 
}

void point::set_quad()
{
	//Defines the quadrant of the element
	//Follows the trigonemtric convention of quadrant definition based on the unit circle
	//R. Sheehan 14 - 1 - 2010

	if(x>=0.0){
		if(y>=0.0){
			quadrant=1; // x>=0 && y>=0
		}
		else{
			quadrant=4; // x>=0 && y<0
		}
	}
	else{
		if(y>=0.0){
			quadrant=2; // x<0 && y>=0
		}
		else{
			quadrant=3; // x<0 && y<0
		}
	}
}