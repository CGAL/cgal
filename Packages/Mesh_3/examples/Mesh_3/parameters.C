#include "parameters.h"

const int number_of_initial_points = 50;

const double curvature_bound = 0.001;

// bound on squared radius of surface Delaunay balls
const double size_bound =  0.5;

// bound on the squared radius of circumspheres
const double tets_size_bound = 0.0;

const double aspect_ratio_bound = 30.; // in degrees

const double tets_aspect_ratio_bound = 2.5;

const double enclosing_sphere_radius = 6.;

const double precision = 1e-6;

const bool bipolar_oracle = true;
