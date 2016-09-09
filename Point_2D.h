#ifndef POINT_2D_H
#define POINT_2D_H

// Class that implements a point in the plane, i.e. the point (x, y) in cartesian coordinates
// R. Sheehan 21 - 5 - 2009

class point{
public:
	// Constructors
	point(); 
	point(double xpos, double ypos); 

	//Members

	// Getters
	inline int get_quad(){return quadrant;} // return the quadrant in which the node lives
	
	inline double get_X(){return x;} // get the value of the node position along the x dirn
	inline double get_Y(){return y;} // get the value of the node position along the y dirn
	
	inline void set_X(double xpos){x = xpos; }
	inline void set_Y(double ypos){y = ypos; }

	void set_XY(double xpos, double ypos); // assign the position of the node
	void set_quad(); // define the quadrant in which the node lives, standard trigonometric interpretation of quadrants

private:
	int quadrant; // in which quadrant does the point live

	double x; // value of x-coordinate
	double y; // value of y-coordinate
}; 

#endif