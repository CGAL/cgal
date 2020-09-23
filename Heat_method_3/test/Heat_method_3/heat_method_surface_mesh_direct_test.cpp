#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <fstream>
#include <iostream>


typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point_3;
typedef CGAL::Surface_mesh<Point_3>                          Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
typedef Surface_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;
typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Surface_mesh, CGAL::Heat_method_3::Direct> Heat_method;



int main(int argc, char* argv[])
{
  //read in mesh
  Surface_mesh sm;
  const char* filename = (argc > 1) ? argv[1] : "./data/sphere.off";
  std::ifstream in(filename);
  in >> sm;
  //property map for the distance values to the source set
  Vertex_distance_map vertex_distance = sm.add_property_map<vertex_descriptor, double>("v:distance", 0).first;

  Heat_method hm(sm);

  //add the first vertex as the source set
  vertex_descriptor source = *(vertices(sm).first);
  hm.add_source(source);
  assert(hm.sources().size() == 1);
  hm.estimate_geodesic_distances(vertex_distance);

  Point_3 sp = sm.point(source);

  std::cout << "source: " << sp  << " " << source << std::endl;
  vertex_descriptor vfar;
  double sdistance = 0;

  for(vertex_descriptor vd : vertices(sm)){
    if(get(vertex_distance,vd) > sdistance){
      vfar = vd;
      sdistance = get(vertex_distance,vd);
    }
  }
  assert(sdistance > 2.9);
  assert(sdistance < CGAL_PI);

  hm.add_source(vfar);
    assert(hm.sources().size() == 2);
  hm.estimate_geodesic_distances(vertex_distance);

  sdistance = 0;
  for(vertex_descriptor vd : vertices(sm)){
    if(get(vertex_distance,vd) > sdistance){
      sdistance = get(vertex_distance,vd);
    }
  }

  assert(sdistance > 1.4);
  assert(sdistance < CGAL_PI/2.0);

  hm.remove_source(vfar);
  assert(hm.sources().size() == 1);

  hm.clear_sources();

  // add range of sources
  std::vector<vertex_descriptor> vrange;
  vrange.push_back(source);
    vrange.push_back(vfar);
  hm.add_sources(vrange);
  assert(hm.sources().size() == 2);
  hm.estimate_geodesic_distances(vertex_distance);
  sdistance = 0;
  for(vertex_descriptor vd : vertices(sm)){
    if(get(vertex_distance,vd) > sdistance){
      sdistance = get(vertex_distance,vd);
    }
  }


  assert(sdistance > 1.4);
  assert(sdistance < CGAL_PI/2.0);

  // do it again for one source
  hm.clear_sources();
  assert(hm.sources().size() == 0);
  hm.add_source(source);
  hm.estimate_geodesic_distances(vertex_distance);
   for(vertex_descriptor vd : vertices(sm)){
    if(get(vertex_distance,vd) > sdistance){
      sdistance = get(vertex_distance,vd);
    }
  }
  assert(sdistance > 2.9);
  assert(sdistance < CGAL_PI);


  CGAL::Heat_method_3::estimate_geodesic_distances(sm, vertex_distance, source, CGAL::Heat_method_3::Direct());
  sdistance = 0;
  for(vertex_descriptor vd : vertices(sm)){
    if(get(vertex_distance,vd) > sdistance){
      sdistance = get(vertex_distance,vd);
    }
  }

  assert(sdistance > 2.9);
  assert(sdistance < CGAL_PI);


  std::cout << "done" << std::endl;
  return 0;
}
