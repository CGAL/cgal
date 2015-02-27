#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/mesh_segmentation.h>

#include <fstream>
#include <map>

typedef CGAL::Simple_cartesian<double>                               Kernel;
typedef Kernel::Point_3                                              Point;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef Polyhedron::Facet_iterator                                   Facet_iterator;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor           vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator             vertex_iterator;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor         halfedge_descriptor;

typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron>        Mean_curvature_skeleton;
typedef Mean_curvature_skeleton::Skeleton                            Skeleton;

typedef boost::graph_traits<Skeleton>::vertex_descriptor             vertex_desc;
typedef boost::graph_traits<Skeleton>::vertex_iterator               vertex_iter;

int main()
{
  Polyhedron mesh;
  std::ifstream input("data/161.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Cannot open data/161.off" << std::endl;
    return 1;
  }

  Skeleton skeleton;
 
  CGAL::extract_mean_curvature_flow_skeleton(mesh, skeleton);

  // create a property-map for sdf values (it is an adaptor for this case)
  typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
  Facet_double_map internal_sdf_map;
  boost::associative_property_map<Facet_double_map> sdf_property_map(internal_sdf_map);

  std::vector<double> distances;
  distances.resize(boost::num_vertices(mesh));

  vertex_iter gvb, gve;
  for (boost::tie(gvb, gve) = boost::vertices(skeleton); gvb != gve; ++gvb)
  {
    vertex_desc gv = *gvb;
    Point skel = skeleton[gv].point;
    for (size_t i = 0; i < skeleton[gv].vertices.size(); ++i)
    {
      Point surf = skeleton[gv].vertices[i]->point();
      distances[skeleton[gv].vertices[i]] = sqrt(squared_distance(skel, surf));
    }
  }

  // compute sdf values with skeleton
  for (Facet_iterator f = mesh.facets_begin(); f != mesh.facets_end(); ++f)
  {
    Polyhedron::Halfedge_const_handle he = f->facet_begin();
    int vid1 = he->vertex()->id();
    int vid2 = he->next()->vertex()->id();
    int vid3 = he->next()->next()->vertex()->id();
    double dis1 = distances[vid1];
    double dis2 = distances[vid2];
    double dis3 = distances[vid3];
    double avg_dis = (dis1 + dis2 + dis3) / 3.0;
    sdf_property_map[f] = avg_dis;
  }

  CGAL::sdf_values_postprocessing(mesh, sdf_property_map);

  // create a property-map for segment-ids (it is an adaptor for this case)
  typedef std::map<Polyhedron::Facet_const_handle, int> Facet_int_map;
  Facet_int_map internal_segment_map;
  boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);

  // segment the mesh using default parameters for number of levels, and smoothing lambda
  // Note that you can use your own scalar value, instead of using SDF calculation computed using the CGAL function
  int number_of_segments = CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map);
  std::cout << "Number of segments: " << number_of_segments << std::endl;

  return 0;
}

