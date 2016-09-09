#ifndef ATTACH_H
#define ATTACH_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <complex>
#include <vector>
#include <cmath>

// Constants
static const double EPS=(3.0e-12);

static const double p=(atan(1.0));
static const double Two_PI=(8.0*p);
static const double PI=(4.0*p);
static const double PI_2=(2.0*p);
static const double PI_3=((4.0/3.0)*p);
static const double PI_4=(p);
static const double PI_5=((4.0/5.0)*p);
static const double PI_6=((2.0/3.0)*p);

// Bend types
static const int STR=1000; // straight waveguide
static const int CC=1001; // constant curvature - i.e. circular
static const int LC=1002; // linear curvature
static const int TC=1003; // trapezoidal curvature
static const int QC=1004; // quadratic curvature

static const std::complex<double> one(1.0, 0.0); 

static const std::string dottxt=".txt";

#include "Useful.h"
#include "Templates.h"
#include "Special_Functions.h" // namespace that contains implementations of various mathematical special functions
#include "Point_2D.h"
#include "Coordinate_Geometry.h"
#include "Bend_Paths.h"

#endif