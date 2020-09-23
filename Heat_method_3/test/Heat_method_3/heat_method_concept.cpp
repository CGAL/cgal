#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

#include "Heat_method_traits_3.h"

#include <fstream>
#include <iostream>


typedef Heat_method_traits_3                                 Kernel;
typedef Kernel::Point_3                                      Point_3;
typedef CGAL::Surface_mesh<Point_3>                          Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
typedef Surface_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;
typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Surface_mesh,
                                                               CGAL::Heat_method_3::Intrinsic_Delaunay,
                                                               boost::property_map< Surface_mesh, CGAL::vertex_point_t>::const_type,
                                                               CGAL::Eigen_solver_traits<Eigen::SimplicialLDLT<CGAL::Eigen_sparse_matrix<double>::EigenType > >,
                                                               Kernel> Heat_method;



int main()
{
  //read in mesh
  Surface_mesh sm;

  if(sm.is_empty()){
    return 0;
  }
  //the heat intensity will hold the distance values from the source set
  Vertex_distance_map heat_intensity = sm.add_property_map<vertex_descriptor, double>("v:heat_intensity", 0).first;

  Heat_method hm(sm);

  //add the first vertex as the source set
  vertex_descriptor source = *(vertices(sm).first);
  hm.add_source(source);
  hm.estimate_geodesic_distances(heat_intensity);

  //Point_3 sp = sm.point(source);

  vertex_descriptor vfar;
  // double sdistance = 0;

  for(vertex_descriptor vd : vertices(sm)){
    std::cout << vd << "  is at distance " << get(heat_intensity, vd) << " from " << source << std::endl;
    /*
    if(squared_distance(sp,sm.point(vd)) > sdistance){
      vfar = vd;
      sdistance = squared_distance(sp,sm.point(vd));
    }
    */
  }
  hm.add_source(vfar);
  hm.estimate_geodesic_distances(heat_intensity);

  for(vertex_descriptor vd : vertices(sm)){
    std::cout << vd << "  is at distance " << get(heat_intensity, vd) << " " << std::endl;
  }

  return 0;
}
