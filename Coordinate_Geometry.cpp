#ifndef ATTACH_H
#include "Attach.h"
#endif

double coord_geom::dist(point p,point q)
{
	// Compute the distance between the two points p and q
	
	double dx, dy;

	dx = ( p.get_X() - q.get_X() );
	
	dy = ( p.get_Y() - q.get_Y() );

	return sqrt( template_funcs::DSQR(dx) + template_funcs::DSQR(dy) );
}

double coord_geom::area(point p,point q,point r)
{
	// Compute the are bounded by the points p, q, r
	
	double d1, d2, d3, u;

	d1 = ( q.get_Y() - r.get_Y() );
	d2 = ( r.get_Y() - p.get_Y() );
	d3 = ( p.get_Y() - q.get_Y() );
	
	u = p.get_X() * d1 + q.get_X() * d2 + r.get_X() * d3;

	return 0.5*fabs(u);
}