#include "parameters.h"

int number_of_initial_points;

bool bipolar_oracle;

Double_options double_options;

void init_parameters()
{
  double_options["curvature_bound"] = 1000000000000000.;
  
  // bound on radius of surface Delaunay balls
  double_options["size_bound"] = 10000000000000.;

  // bound on the squared radius of circumspheres
  double_options["tets_size_bound"] = 10000000000.;
  double_options["facets_aspect_ratio_bound"] = 0.; // in degrees
  double_options["tets_aspect_ratio_bound"] = 0.;
  double_options["enclosing_sphere_radius"] = 6.;
  double_options["precision"] = 1e-6;
  double_options["center_x"] = 0.0;
  double_options["center_y"] = 0.0;
  double_options["center_z"] = 0.0;
  double_options["number_of_initial_points"] = 20;
  double_options["number_of_pump"] = 1;
  double_options["pumping_bound"] = 0.45;
  double_options["sliver_test"] = 0.15;
  bipolar_oracle = true;
}
