#ifndef COORDINATE_GEOMETRY_H
#define COORDINATE_GEOMETRY_H

// namespace that implements various formulae that are useful in the field of coordinate geometry
// R. Sheehan 7 - 8 - 2012

namespace coord_geom{

	double dist(point p, point q); // distance between two points in the plane

	double area(point p, point q, point r); // area contained by three points in the plane
}

#endif