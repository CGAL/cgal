//#define POLY 

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#ifdef POLY
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#else
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/properties_Surface_mesh.h>
#endif
#include <CGAL/boost/graph/helpers.h>

#include <CGAL/triangulate_hole.h>
#include <CGAL/internal/Hole_filling/Triangulate_hole_Polyhedron_3.h>

#include <CGAL/assertions.h>

#include <cassert>
#include <vector>
#include <set>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
#ifdef POLY
typedef CGAL::Polyhedron_3<Kernel>       Polyhedron;
#else
typedef CGAL::Surface_mesh<Kernel::Point_3> Polyhedron;
#endif
typedef boost::graph_traits<Polyhedron>::face_descriptor Facet_handle;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor        Vertex_handle;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor    Halfedge_handle;
typedef boost::graph_traits<Polyhedron>::halfedge_iterator    Halfedge_iterator;
typedef CGAL::Halfedge_around_face_circulator<Polyhedron> Halfedge_around_facet_circulator;
typedef boost::property_map<Polyhedron,CGAL::vertex_point_t>::type Point_property_map;

void read_poly_with_borders(const char* file_name, Polyhedron& poly, std::vector<Halfedge_handle>& border_reps) {
  border_reps.clear();
  poly.clear();

  std::ifstream input(file_name);
  if ( !input || !(input >> poly)  || (num_vertices(poly) == 0)){
    std::cerr << "  Error: can not read file." << std::endl;
    assert(false);
  }

  std::set<Halfedge_handle> border_map;
  BOOST_FOREACH(Halfedge_handle h,  halfedges(poly)){
    if(face(h,poly)== boost::graph_traits<Polyhedron>::null_face() && border_map.find(h) == border_map.end()){
      border_reps.push_back(h);
      Halfedge_around_facet_circulator hf_around_facet(h,poly), done(hf_around_facet);
      do {
        bool insertion_ok = border_map.insert(*hf_around_facet).second;
        CGAL_assertion(insertion_ok);
      } while(++hf_around_facet != done);
    }
  }
}

/******************************************************************/
// This part test internal functions with Weight_min_max_dihedral_and_area
template<class Polyhedron, class Iterator>
CGAL::internal::Weight_min_max_dihedral_and_area
           calculate_weight_for_patch(Polyhedron& poly,Iterator begin, Iterator end) 
{
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor Halfedge_handle;
  typedef CGAL::internal::Weight_min_max_dihedral_and_area Weight;

  Point_property_map ppmap = get(CGAL::vertex_point, poly);
  Weight res(0,0);
  for(; begin!=end; ++begin) {
    Halfedge_handle edge_it = halfedge(*begin, poly);
    double ang_max = 0;
    for(int i = 0; i < 3; ++i) {
      double angle = 180 - CGAL::abs( 
                                     CGAL::Mesh_3::dihedral_angle(ppmap[target(edge_it,poly)],
                                                                  ppmap[source(edge_it,poly)],
                                                                  ppmap[target(next(edge_it,poly),poly)],
                                                                  ppmap[target(next(opposite(edge_it,poly),poly),poly)]) );
      edge_it = next(edge_it,poly);
      ang_max = (std::max)(angle, ang_max);
    }

    double area = std::sqrt(CGAL::squared_area(ppmap[target(edge_it,poly)],
                                               ppmap[target(next(edge_it,poly),poly)],
                                               ppmap[target(prev(edge_it,poly),poly)]));
    res = res + Weight(ang_max,area);
  }
  return res;
}


void test_triangulate_hole_weight(const char* file_name, bool use_DT) {
  typedef CGAL::internal::Weight_min_max_dihedral_and_area Weight;

  std::cerr << "test_triangulate_hole_weight + useDT: " << use_DT << std::endl;
  std::cerr << "  File: "<< file_name  << std::endl;
  Polyhedron poly;
  std::vector<Halfedge_handle> border_reps;
  read_poly_with_borders(file_name, poly, border_reps);

  for(std::vector<Halfedge_handle>::iterator it = border_reps.begin(); it != border_reps.end(); ++it) {
    std::vector<Facet_handle> patch;
    Weight w_algo = CGAL::internal::triangulate_hole_Polyhedron(poly, *it, back_inserter(patch), use_DT).second;
    if(patch.empty()) { continue; }
    Weight w_test = calculate_weight_for_patch(poly, patch.begin(), patch.end());

    std::cerr << "  Weight returned by algo   : " << w_algo << std::endl;
    std::cerr << "  Weight calculated by test : " << w_test << std::endl;

    const double epsilon = 1e-10;
    if(std::abs(w_algo.w.first - w_test.w.first) > epsilon) {
      assert(false);
    }
    if(std::abs(w_algo.w.second - w_test.w.second) > epsilon) {
      assert(false);
    }
  }

  std::cerr << "  Done!" << std::endl;
}
/******************************************************************/

void test_triangulate_hole(const char* file_name) {
  std::cerr << "test_triangulate_hole:" << std::endl;
  std::cerr << "  File: "<< file_name  << std::endl;
  Polyhedron poly;
  std::vector<Halfedge_handle> border_reps;
  read_poly_with_borders(file_name, poly, border_reps);

  for(std::vector<Halfedge_handle>::iterator it = border_reps.begin(); it != border_reps.end(); ++it) {
    std::vector<Facet_handle> patch;
    CGAL::triangulate_hole(poly, *it, back_inserter(patch));
    if(patch.empty()) {
      std::cerr << "  Error: empty patch created." << std::endl;
      assert(false);
    }
  }

  if(!poly.is_valid() || ! is_closed(poly)) {
    std::cerr << "  Error: patched polyhedron is not valid or closed." << std::endl;
    assert(false);
  }
  
  std::cerr << "  Done!" << std::endl;
}

void test_triangulate_hole_should_be_no_output(const char* file_name) {
  std::cerr << "test_triangulate_hole_should_be_no_output:" << std::endl;
  std::cerr << "  File: "<< file_name  << std::endl;
  Polyhedron poly;
  std::vector<Halfedge_handle> border_reps;
  read_poly_with_borders(file_name, poly, border_reps);

  for(std::vector<Halfedge_handle>::iterator it = border_reps.begin(); it != border_reps.end(); ++it) {
    std::vector<Facet_handle> patch;
    CGAL::triangulate_hole(poly, *it, back_inserter(patch), false);
    if(!patch.empty()) {
      std::cerr << "  Error: patch should be empty" << std::endl;
      assert(false);
    }

    CGAL::triangulate_hole(poly, *it, back_inserter(patch), true);
    if(!patch.empty()) {
      std::cerr << "  Error: patch should be empty" << std::endl;
      assert(false);
    }
  }

  std::cerr << "  Done!" << std::endl;
}

void test_triangulate_and_refine_hole(const char* file_name) {
  std::cerr << "test_triangulate_and_refine_hole:" << std::endl;
  std::cerr << "  File: "<< file_name  << std::endl;
  Polyhedron poly;
  std::vector<Halfedge_handle> border_reps;
  read_poly_with_borders(file_name, poly, border_reps);

  for(std::vector<Halfedge_handle>::iterator it = border_reps.begin(); it != border_reps.end(); ++it) {
    std::vector<Facet_handle> patch_facets;
    std::vector<Vertex_handle> patch_vertices;
    CGAL::triangulate_and_refine_hole(poly, *it, 
      back_inserter(patch_facets), back_inserter(patch_vertices));

    if(patch_facets.empty()) {
      std::cerr << "  Error: empty patch created." << std::endl;
      assert(false);
    }
  }

  if(!poly.is_valid() || ! is_closed(poly)) {
    std::cerr << "  Error: patched polyhedron is not valid or closed." << std::endl;
    assert(false);
  }

  std::cerr << "  Done!" << std::endl;
}

void test_triangulate_refine_and_fair_hole(const char* file_name) {
  std::cerr << "test_triangulate_refine_and_fair_hole:" << std::endl;
  std::cerr << "  File: "<< file_name  << std::endl;
  Polyhedron poly;
  std::vector<Halfedge_handle> border_reps;
  read_poly_with_borders(file_name, poly, border_reps);

  for(std::vector<Halfedge_handle>::iterator it = border_reps.begin(); it != border_reps.end(); ++it) {
    std::vector<Facet_handle> patch_facets;
    std::vector<Vertex_handle> patch_vertices;
    CGAL::triangulate_refine_and_fair_hole(poly, *it, back_inserter(patch_facets), back_inserter(patch_vertices));

    if(patch_facets.empty()) {
      std::cerr << "  Error: empty patch created." << std::endl;
      assert(false);
    }
  }

  if(!poly.is_valid() || ! is_closed(poly)) {
    std::cerr << "  Error: patched polyhedron is not valid or closed." << std::endl;
    assert(false);
  }

  std::cerr << "  Done!" << std::endl;
}

void test_ouput_iterators_triangulate_hole(const char* file_name) {
  std::cerr << "test_ouput_iterators_triangulate_hole:" << std::endl;
  std::cerr << "  File: "<< file_name  << std::endl;

  Polyhedron poly, poly_2;
  std::vector<Halfedge_handle> border_reps, border_reps_2;

  read_poly_with_borders(file_name, poly, border_reps);
  read_poly_with_borders(file_name, poly_2, border_reps_2);

  std::vector<Halfedge_handle>::iterator it_2 = border_reps_2.begin();
  for(std::vector<Halfedge_handle>::iterator it = border_reps.begin(); it != border_reps.end(); ++it, ++it_2) {
    std::vector<Facet_handle> patch;
    CGAL::triangulate_hole(poly, *it, back_inserter(patch));

    std::vector<Facet_handle> patch_2 = patch;
    Facet_handle* output_it = CGAL::triangulate_hole(poly_2, *it_2, &*patch_2.begin());

    if(patch.size() != (output_it - &*patch_2.begin())) {
      std::cerr << "  Error: returned facet output iterator is not valid!" << std::endl;
      std::cerr << "  " << patch.size() << " vs " << (output_it - &*patch_2.begin()) << std::endl;
      assert(false);
    }
  }
  std::cerr << "  Done!" << std::endl;
}

void test_ouput_iterators_triangulate_and_refine_hole(const char* file_name) {
  std::cerr << "test_ouput_iterators_triangulate_and_refine_hole:" << std::endl;
  std::cerr << "  File: "<< file_name  << std::endl;

  Polyhedron poly, poly_2;
  std::vector<Halfedge_handle> border_reps, border_reps_2;;

  read_poly_with_borders(file_name, poly, border_reps);
  read_poly_with_borders(file_name, poly_2, border_reps_2);

  std::vector<Halfedge_handle>::iterator it_2 = border_reps_2.begin();
  for(std::vector<Halfedge_handle>::iterator it = border_reps.begin(); it != border_reps.end(); ++it, ++it_2) {
    std::vector<Facet_handle> patch_facets;
    std::vector<Vertex_handle> patch_vertices;
    CGAL::triangulate_and_refine_hole(poly, *it, back_inserter(patch_facets), back_inserter(patch_vertices));
    // create enough space to hold outputs
    std::vector<Facet_handle> patch_facets_2 = patch_facets;
    std::vector<Vertex_handle> patch_vertices_2 = patch_vertices;
    if(patch_vertices_2.empty()) { patch_vertices_2.push_back(Vertex_handle()); } //just allocate space for dereferencing

    std::pair<Facet_handle*, Vertex_handle*> output_its = 
      CGAL::triangulate_and_refine_hole(poly_2, *it_2, &*patch_facets_2.begin(), &*patch_vertices_2.begin());

    if(patch_facets.size() != (output_its.first - &*patch_facets_2.begin())) {
      std::cerr << "  Error: returned facet output iterator is not valid!" << std::endl;
      std::cerr << "  " << patch_facets.size() << " vs " << (output_its.first - &*patch_facets_2.begin()) << std::endl;
      assert(false);
    }

    if(patch_vertices.size() != (output_its.second - &*patch_vertices_2.begin())) {
      std::cerr << "  Error: returned vertex output iterator is not valid!" << std::endl;
      std::cerr << "  " << patch_vertices.size() << " vs " << (output_its.second - &*patch_vertices_2.begin()) << std::endl;
      assert(false);
    }
  }
  std::cerr << "  Done!" << std::endl;
}


void test_triangulate_refine_and_fair_hole_compile() {
  typedef CGAL::Eigen_solver_traits<
    Eigen::SparseLU<
      CGAL::Eigen_sparse_matrix<double>::EigenType,
      Eigen::COLAMDOrdering<int> > >
  Default_solver;

  Polyhedron poly;
  std::vector<Halfedge_handle> border_reps;
  
  std::vector<Facet_handle> patch_facets;
  std::vector<Vertex_handle> patch_vertices;

  // use all param
  read_poly_with_borders("data/elephant_quad_hole.off", poly, border_reps);
  CGAL::triangulate_refine_and_fair_hole
  (poly, border_reps[0], back_inserter(patch_facets), back_inserter(patch_vertices),
    CGAL::internal::Uniform_weight_fairing<Polyhedron>(poly), Default_solver());

  // default solver
  read_poly_with_borders("data/elephant_quad_hole.off", poly, border_reps);
  CGAL::triangulate_refine_and_fair_hole
    (poly, border_reps[0], back_inserter(patch_facets), back_inserter(patch_vertices),
    CGAL::internal::Uniform_weight_fairing<Polyhedron>(poly), CGAL::Default());

  // default solver and weight
  read_poly_with_borders("data/elephant_quad_hole.off", poly, border_reps);
  CGAL::triangulate_refine_and_fair_hole
    (poly, border_reps[0], back_inserter(patch_facets), back_inserter(patch_vertices));
}

int main() {

  std::vector<std::string> input_files;
  input_files.push_back("data/elephant_triangle_hole.off");
  input_files.push_back("data/elephant_quad_hole.off");
  input_files.push_back("data/mech-holes-shark.off");
  // std::cerr.precision(15);
  for(std::vector<std::string>::iterator it = input_files.begin(); it != input_files.end(); ++it) {
    test_triangulate_hole(it->c_str());
    test_triangulate_and_refine_hole(it->c_str());
    test_triangulate_refine_and_fair_hole(it->c_str());
    test_ouput_iterators_triangulate_and_refine_hole(it->c_str());
    test_ouput_iterators_triangulate_hole(it->c_str());
    test_triangulate_hole_weight(it->c_str(), true);
    test_triangulate_hole_weight(it->c_str(), false);
    std::cerr << "------------------------------------------------" << std::endl;
  }
  test_triangulate_hole_weight("data/RedCircleBox.off", true);
  test_triangulate_hole_weight("data/RedCircleBox.off", false);

  test_triangulate_hole_should_be_no_output("data/non_manifold_vertex.off");
  test_triangulate_hole_should_be_no_output("data/two_tris_collinear.off");

  test_triangulate_refine_and_fair_hole_compile();
  std::cerr << "All Done!" << std::endl;
}
