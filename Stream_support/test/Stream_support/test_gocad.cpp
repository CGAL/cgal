
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/helpers.h>

#include <CGAL/IO/GOCAD.h>
#include <iostream>
#include <sstream>

template<class FaceGraph, class Point>
bool test_io()
{
  FaceGraph fg;
  /*const char* tet = "GOCAD TSurf 1                   \n"
                    "HEADER {                        \n"
                    "name:Tetrahedron                \n"
                    "*border:on                      \n"
                    "*border*bstone:on               \n"
                    "}                               \n"
                    "GOCAD_ORIGINAL_COORDINATE_SYSTEM\n"
                    "NAME Default                    \n"
                    "AXIS_NAME \"X\" \"Y\" \"Z\"     \n"
                    "AXIS_UNIT \"m\" \"m\" \"m\"     \n"
                    "ZPOSITIVE Elevation             \n"
                    "END_ORIGINAL_COORDINATE_SYSTEM  \n"
                    "TFACE                           \n"
                    "VRTX 0 0 0 0                    \n"
                    "VRTX 1 1 0 0                    \n"
                    "VRTX 2 0 1 0                    \n"
                    "VRTX 3 0 0 1                    \n"
                    "TRGL 0 2 1                      \n"
                    "TRGL 2 0 3                      \n"
                    "TRGL 1 2 3                      \n"
                    "TRGL 0 1 3                      \n"
                    "END                             \n";*/
  CGAL::make_tetrahedron(Point(0, 0, 0), Point(1, 1, 0),
                         Point(2, 0, 1), Point(3, 0, 0), fg);
  std::ostringstream out;
  out << fg;
  if(out.fail())
  {
    std::cerr<<"Tetrahedron writing failed."<<std::endl;
    return false;
  }
  FaceGraph fg2;
  std::istringstream in( out.str());
  in >> fg2;

  if(in.fail()){
    std::cerr<<"Tetrahedron reading failed."<<std::endl;
    return false;
  }

  if(num_vertices(fg2) != 4
     || num_faces(fg2) != 4)
  {
    std::cerr<<"Facegraph construction failed."<<std::endl;
    return false;
  }

  return true;
}

int main()
{
  typedef CGAL::Exact_predicates_exact_constructions_kernel Epeck;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;

  if(!test_io<CGAL::Surface_mesh<Epick::Point_3>, Epick::Point_3>())
  {
    return 1;
  }
  if(!test_io<CGAL::Surface_mesh<Epeck::Point_3>,Epeck::Point_3>())
  {
    return 1;
  }
  if(!test_io<CGAL::Polyhedron_3<Epick>, Epick::Point_3>())
  {
    return 1;
  }
  if(!test_io<CGAL::Polyhedron_3<Epeck>, Epeck::Point_3>())
  {
    return 1;
  }
  return 0;
}
