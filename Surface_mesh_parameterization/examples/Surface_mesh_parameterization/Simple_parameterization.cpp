#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <CGAL/IO/Surface_mesh_parameterization/File_off.h>
#include <CGAL/parameterize.h>

#include <CGAL/Polygon_mesh_processing/measure.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Polyhedron_3<Kernel> Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

int main(int argc, char * argv[])
{
  std::ifstream in((argc>1)?argv[1]:"../data/mushroom_big_hole.off");
  if(!in){
    std::cerr << "Problem loading the input data" << std::endl;
    return 1;
  }

  Mesh sm;
  in >> sm;

  halfedge_descriptor hd = CGAL::Polygon_mesh_processing::longest_border(sm).first;

  CGAL::Unique_hash_map<vertex_descriptor,Point_2> uvhm;
  boost::associative_property_map<CGAL::Unique_hash_map<vertex_descriptor,Point_2> > uvmap(uvhm);

  CGAL::parameterize(sm, hd, uvmap);

  std::ofstream out("result.off");
  CGAL::Parameterization::output_uvmap_to_off(sm, hd, uvmap, out);

  return 0;
}
