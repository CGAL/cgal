#include <CGAL/leda_integer.h>
typedef leda_integer RT;

#include <CGAL/Homogeneous.h>
typedef CGAL::Homogeneous<RT> Kernel;

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Timer.h>
#include <sstream>
#include <fstream>
#include "grid_generator.h"

typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Aff_transformation_3 Aff_transformation_3;
typedef Nef_polyhedron::Vector_3 Vector_3;
typedef CGAL::grid_generator<Nef_polyhedron> ggen;

bool cgal_nef3_timer_on = false;

int main(int argc, char* argv[]) {

  assert(argc>1 && argc < 5);
  
  int nx = argc>2 ? std::atoi(argv[2]) : 2;
  int ny = argc>3 ? std::atoi(argv[3]) : nx;

  std::ifstream in(argv[1]);
  Nef_polyhedron Nin;
  in >> Nin;

  std::ostringstream out1;
  ggen g(out1, Nin,false);
  g.print(nx,1,1);
  std::istringstream in1(out1.str());
  Nef_polyhedron N1;
  in1 >> N1;
  CGAL_assertion(N1.is_valid());

  Nin.transform(Aff_transformation_3(0,-1,0,
				     1,0,0,
				     0,0,1,1));

  std::ostringstream out2;
  ggen g2(out2, Nin,false);
  g2.print(1,ny,1);
  std::istringstream in2(out2.str());
  Nef_polyhedron N2;
  in2 >> N2;
  CGAL_assertion(N2.is_valid());

  char* dir="nef3/quadratic_";
  char* horizontal="horizontal";
  char* vertical="vertical";
  char* translated="vertical_translated";
  char* suffix=".nef3";
  char* us="_";
  
  char* full_suffix = new char[strlen(argv[2])+strlen(suffix)+2];
  strcpy(full_suffix, us);
  strcat(full_suffix, argv[2]);
  strcat(full_suffix, suffix);

  char* full_horizontal = new char[strlen(dir)+strlen(horizontal)+strlen(full_suffix)+1];
  strcpy(full_horizontal, dir);
  strcat(full_horizontal, horizontal);
  strcat(full_horizontal, full_suffix);
  std::ofstream out_horizontal(full_horizontal);
  out_horizontal << N2;

  char* full_vertical = new char[strlen(dir)+strlen(vertical)+strlen(full_suffix)+1];
  strcpy(full_vertical, dir);
  strcat(full_vertical, vertical);
  strcat(full_vertical, full_suffix);
  std::ofstream out_vertical(full_vertical);
  out_vertical << N1;

  N1.transform(Aff_transformation_3(CGAL::TRANSLATION, Vector_3(0,0,1)));

  char* full_translated = new char[strlen(dir)+strlen(translated)+strlen(full_suffix)+1];
  strcpy(full_translated, dir);
  strcat(full_translated, translated);
  strcat(full_translated, full_suffix);
  std::ofstream out_translated(full_translated);
  out_translated << N1;
}
