// WARNING : non-functioning example. This exmaple is to understand how to create and use the feature graph for tetrahedral mesh creation. 

#include <CGAL/Image_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Meshing definition
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>

// #include <CGAL/Mesh_3/polylines_to_protect.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/Feature_graph/detect_sharp_features.h>

/// [Meshing definition]
#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Image_domain;
typedef CGAL::Mesh_domain_with_polyline_features_3<Image_domain> Mesh_domain;

typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
/// [Meshing definition]

/// [Detector definition]
typedef CGAL::Feature_graph::Detect_sharp_features_on_labeled_image<K> Feature_graph_detector;
/// [Detector definition]

// To avoid verbose function and named parameters call
namespace params = CGAL::parameters;

template<typename P, typename G>
struct Polyline_visitor
{
  std::vector<std::vector<P> >& polylines;
  G& graph;
  Polyline_visitor(typename std::vector<std::vector<P> >& lines, G& p_graph)
    : polylines(lines), graph(p_graph)
  {}

  void start_new_polyline()
  {
    std::vector<P> V;
    polylines.push_back(V);
  }

  void add_node(typename boost::graph_traits<G>::vertex_descriptor vd)
  {
    std::vector<P>& polyline = polylines.back();
    polyline.push_back(graph[vd]);
  }

  void end_polyline()
  {
    // ignore degenerated polylines
    if(polylines.back().size() < 2)
      polylines.resize(polylines.size() - 1);
  }
};

int main(int argc, char* argv[])
{
  const std::string fname = (argc>1)?argv[1]:CGAL::data_file_path("images/quadDomainCube.inr");

  /// [Loads image]
  CGAL::Image_3 image;
  if(!image.read(fname)){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }
  /// [Loads image]

  /// [Feature detection]
  Feature_graph_detector feature_graph_detector;
  const auto& feature_graph = feature_graph_detector(image);
  /// [Feature detection]

  /// [Domain creation]
  std::vector<std::vector<K::Point_3>> feature_graph_polylines;
  Polyline_visitor visitor(feature_graph_polylines, feature_graph);
  split_graph_into_polylines(feature_graph, visitor);
  Mesh_domain domain
    = Mesh_domain::create_labeled_image_mesh_domain(image,
         params::input_features = feature_graph_polylines);
  /// [Domain creation]

  CGAL::Bbox_3 bbox = domain.bbox();
  double diag = CGAL::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                           CGAL::square(bbox.ymax() - bbox.ymin()) +
                           CGAL::square(bbox.zmax() - bbox.zmin()));
  double sizing_default = diag * 0.05;

  /// [Mesh criteria]
  /// Note that `edge_size` is needed with 1D-features
  Mesh_criteria criteria(params::edge_size = sizing_default,
    params::facet_angle = 30,
    params::facet_size = sizing_default,
    params::facet_distance = sizing_default / 10,
    params::cell_radius_edge_ratio = 0,
    params::cell_size = 0
  );
  /// [Mesh criteria]

  /// [Meshing]
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      params::no_exude(),
                                      params::no_perturb());
  /// [Meshing]

  // Output
  CGAL::dump_c3t3(c3t3, "out");

  return 0;
}
