
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/helpers.h>

#include <CGAL/IO/GOCAD.h>
#include <iostream>
#include <sstream>
#include <string>

template<class FaceGraph, class Point>
bool test_io()
{
  FaceGraph fg;
  CGAL::make_tetrahedron(Point(0, 0, 0), Point(1, 1, 0),
                         Point(2, 0, 1), Point(3, 0, 0), fg);
  std::ostringstream out;
  CGAL::write_gocad(fg, out, "tetrahedron");
  if(out.fail())
  {
    std::cerr<<"Tetrahedron writing failed."<<std::endl;
    return false;
  }
  FaceGraph fg2;
  std::istringstream in( out.str());
  std::string name, color;
  CGAL::read_gocad(fg2, in, name, color);
  if(name != "tetrahedron"){
    std::cerr<<"reading error: tetrahedron != "<<name<<std::endl;
    return 1;
  }
  if( !color.empty()){
    std::cerr<<"reading error: there should be no color."<<std::endl;
    return 1;
  }

  if(in.fail()){
    std::cerr<<"Tetrahedron reading failed."<<std::endl;
    return false;
  }

  if(num_vertices(fg2) != 4){
    std::cerr<<"Wrong number of vertices: 4 != "<<num_vertices(fg2)<<std::endl;
    return false;
  }

   if(num_faces(fg2) != 4)
  {
    std::cerr<<"Wrong number of faces: 4 != "<<num_faces(fg2)<<std::endl;
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
