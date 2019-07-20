#include <CGAL/leda_integer.h>
typedef leda_integer RT;

#include <CGAL/Homogeneous.h>
typedef CGAL::Homogeneous<RT> Kernel;

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Timer.h>
#include <sstream>
#include <fstream>
#include "tetrahedron_generator.h"
#include "grid_generator.h"

typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Vector_3 Vector_3;
typedef Nef_polyhedron::Aff_transformation_3 Aff_transformation_3;
typedef tetrahedron_generator<Kernel> tgen;
typedef CGAL::grid_generator<Nef_polyhedron> ggen;

bool cgal_nef3_timer_on = false;

int main(int argc, char* argv[]) {

  assert(argc>1 && argc < 7);
  
  int nx = argc>2 ? std::atoi(argv[2]) : 2;
  int ny = argc>3 ? std::atoi(argv[3]) : 2;
  int nz = argc>4 ? std::atoi(argv[4]) : 2;

  std::ifstream in(argv[1]);
  Nef_polyhedron Nin;
  in >> Nin;
  Nin.transform(Aff_transformation_3(CGAL::SCALING,2,1));
  std::ostringstream out1;
  ggen g(out1, Nin);
  g.print(nx,ny,nz);
  std::istringstream in1(out1.str());
  Nef_polyhedron N1;
  in1 >> N1;
  RT s = g.size_x();
  N1.transform(Aff_transformation_3(CGAL::TRANSLATION,Vector_3(s,s,s,2)));
  CGAL_assertion(N1.is_valid());

  std::ostringstream out2;
  CGAL::Random r;
  if(argc>5) {
    tgen t2(out2,s,std::atoi(argv[5]));
    t2.create_tetrahedra(nx+1,ny+1,nz+1);
  } else {
    tgen t2(out2,s);
    t2.create_tetrahedra(nx+1,ny+1,nz+1);    
  }
  std::istringstream in2(out2.str());
  Nef_polyhedron N2;
  in2 >> N2;
  CGAL_assertion(N2.is_valid());

  char* dir="nef3/";
  char* tetrahedra="tetrahedra";
  char* grid="grid";
  char* suffix=".nef3";
  char* us="_";
  
  char* full_suffix = new char[strlen(argv[2])+strlen(argv[3])+strlen(argv[4])+strlen(argv[5])+strlen(suffix)+5];
  strcpy(full_suffix, us);
  strcat(full_suffix, argv[2]);
  strcat(full_suffix, us);
  strcat(full_suffix, argv[3]);
  strcat(full_suffix, us);
  strcat(full_suffix, argv[4]);
  strcat(full_suffix, us);
  strcat(full_suffix, argv[5]);
  strcat(full_suffix, suffix);

  char* full_tetrahedra = new char[strlen(dir)+strlen(tetrahedra)+strlen(full_suffix)+1];
  strcpy(full_tetrahedra, dir);
  strcat(full_tetrahedra, tetrahedra);
  strcat(full_tetrahedra, full_suffix);
  std::ofstream out_tetrahedra(full_tetrahedra);
  out_tetrahedra << N2;

  char* full_grid = new char[strlen(dir)+strlen(grid)+strlen(full_suffix)+1];
  strcpy(full_grid, dir);
  strcat(full_grid, grid);
  strcat(full_grid, full_suffix);
  std::ofstream out_grid(full_grid);
  out_grid << N1;  
}
