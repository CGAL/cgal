#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> SMesh;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

template<typename MAP>
struct Bot
{
  Bot(MAP map):map(map){}
  template<typename VD, typename T>
  void operator()(const T&,VD vd) const
  {
    put(map, vd, get(map, vd)+Kernel::Vector_3(-2.0,0.0,1.0));
  }
  MAP map;

};

template<typename MAP>
struct Top
{
  Top(MAP map):map(map){}

  template<typename VD, typename T>
  void operator()(const T&, VD vd) const
  {
    put(map, vd, get(map, vd)+Kernel::Vector_3(0.0,2.0,-1.0));
  }

  MAP map;
};

template <class Mesh>
void test_mesh(const char* filename)
{
  Mesh in, out;
  std::ifstream input(filename);

  if (!input || !(input >> in))
  {
    std::cerr << "Error: cannot read Surface Mesh : " << filename << "\n";
    assert(!CGAL::is_empty(in));
    assert(false);
    return ;
  }
  CGAL::Polygon_mesh_processing::extrude_mesh(in, out, Kernel::Vector_3(0.0, 0.0, -1.0));
  std::ofstream extruded_off("extruded.off");
  extruded_off << out;
  extruded_off.close();
  out.clear();

  typedef typename boost::property_map<Mesh, CGAL::vertex_point_t>::type VPMap;
  Bot<VPMap> bot(get(CGAL::vertex_point, out));
  Top<VPMap> top(get(CGAL::vertex_point, out));
  CGAL::Polygon_mesh_processing::extrude_mesh(in, out, bot, top);
  std::ofstream gen_extruded_off("gen_extruded.off");
  gen_extruded_off << out;
  gen_extruded_off.close();
  std::cerr << "All done." << std::endl;
}

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/quad.off";
  test_mesh<SMesh>(filename);
  test_mesh<Polyhedron>(filename);
  return 0;
}
