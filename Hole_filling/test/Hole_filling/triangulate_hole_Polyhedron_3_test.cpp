#include <CGAL/Hole_filling.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/assertions.h>

#include <cassert>
#include <vector>
#include <set>

typedef CGAL::Simple_cartesian<double>   Kernel;
typedef CGAL::Polyhedron_3<Kernel>       Polyhedron;
typedef Polyhedron::Facet_handle         Facet_handle;
typedef Polyhedron::Vertex_handle        Vertex_handle;
typedef Polyhedron::Halfedge_handle      Halfedge_handle;
typedef Polyhedron::Halfedge_iterator    Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;

void read_poly_with_borders(const char* file_name, Polyhedron& poly, std::vector<Halfedge_handle>& border_reps) {
  border_reps.clear();
  poly.clear();

  std::ifstream input(file_name);
  if ( !input || !(input >> poly) || poly.empty() ){
    std::cerr << "  Error: can not read file." << std::endl;
    assert(false);
  }

  std::set<Halfedge_handle> border_map;
  for(Halfedge_iterator it = poly.halfedges_begin(); it != poly.halfedges_end(); ++it){
    if(it->is_border() && border_map.find(it) == border_map.end()){
      border_reps.push_back(it);
      Halfedge_around_facet_circulator hf_around_facet = it->facet_begin();
      do {
        bool insertion_ok = border_map.insert(hf_around_facet).second;
        CGAL_assertion(insertion_ok);
      } while(++hf_around_facet != it->facet_begin());
    }
  }
}

void test_triangulate_hole(const char* file_name) {
  std::cerr << "test_triangulate_hole:" << std::endl;
  std::cerr << "  File: "<< file_name  << std::endl;
  Polyhedron poly;
  std::vector<Halfedge_iterator> border_reps;
  read_poly_with_borders(file_name, poly, border_reps);

  for(std::vector<Halfedge_iterator>::iterator it = border_reps.begin(); it != border_reps.end(); ++it) {
    std::vector<Facet_handle> patch;
    CGAL::triangulate_hole(poly, *it, back_inserter(patch));
    if(patch.empty()) {
      std::cerr << "  Error: empty patch created." << std::endl;
      assert(false);
    }
  }

  if(!poly.is_valid() || !poly.is_closed()) {
    std::cerr << "  Error: patched polyhedron is not valid or closed." << std::endl;
    assert(false);
  }
  
  std::cerr << "  Done!" << std::endl;
}

void test_triangulate_and_refine_hole(const char* file_name) {
  std::cerr << "test_triangulate_and_refine_hole:" << std::endl;
  std::cerr << "  File: "<< file_name  << std::endl;
  Polyhedron poly;
  std::vector<Halfedge_iterator> border_reps;
  read_poly_with_borders(file_name, poly, border_reps);

  for(std::vector<Halfedge_iterator>::iterator it = border_reps.begin(); it != border_reps.end(); ++it) {
    std::vector<Facet_handle> patch_facets;
    std::vector<Vertex_handle> patch_vertices;
    CGAL::triangulate_and_refine_hole(poly, *it, 
      back_inserter(patch_facets), back_inserter(patch_vertices));

    if(patch_facets.empty()) {
      std::cerr << "  Error: empty patch created." << std::endl;
      assert(false);
    }
  }

  if(!poly.is_valid() || !poly.is_closed()) {
    std::cerr << "  Error: patched polyhedron is not valid or closed." << std::endl;
    assert(false);
  }

  std::cerr << "  Done!" << std::endl;
}


void test_triangulate_refine_and_fair_hole(const char* file_name) {
  std::cerr << "test_triangulate_refine_and_fair_hole:" << std::endl;
  std::cerr << "  File: "<< file_name  << std::endl;
  Polyhedron poly;
  std::vector<Halfedge_iterator> border_reps;
  read_poly_with_borders(file_name, poly, border_reps);

  for(std::vector<Halfedge_iterator>::iterator it = border_reps.begin(); it != border_reps.end(); ++it) {
    std::vector<Facet_handle> patch_facets;
    std::vector<Vertex_handle> patch_vertices;
    CGAL::triangulate_refine_and_fair_hole(poly, *it, back_inserter(patch_facets), back_inserter(patch_vertices));

    if(patch_facets.empty()) {
      std::cerr << "  Error: empty patch created." << std::endl;
      assert(false);
    }
  }

  if(!poly.is_valid() || !poly.is_closed()) {
    std::cerr << "  Error: patched polyhedron is not valid or closed." << std::endl;
    assert(false);
  }

  std::cerr << "  Done!" << std::endl;
}

void test_ouput_iterators_triangulate_hole(const char* file_name) {
  std::cerr << "test_ouput_iterators_triangulate_hole:" << std::endl;
  std::cerr << "  File: "<< file_name  << std::endl;

  Polyhedron poly, poly_2;
  std::vector<Halfedge_iterator> border_reps, border_reps_2;

  read_poly_with_borders(file_name, poly, border_reps);
  read_poly_with_borders(file_name, poly_2, border_reps_2);

  std::vector<Halfedge_iterator>::iterator it_2 = border_reps_2.begin();
  for(std::vector<Halfedge_iterator>::iterator it = border_reps.begin(); it != border_reps.end(); ++it, ++it_2) {
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
  std::vector<Halfedge_iterator> border_reps, border_reps_2;;

  read_poly_with_borders(file_name, poly, border_reps);
  read_poly_with_borders(file_name, poly_2, border_reps_2);

  std::vector<Halfedge_iterator>::iterator it_2 = border_reps_2.begin();
  for(std::vector<Halfedge_iterator>::iterator it = border_reps.begin(); it != border_reps.end(); ++it, ++it_2) {
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
  std::vector<Halfedge_iterator> border_reps;
  
  std::vector<Facet_handle> patch_facets;
  std::vector<Vertex_handle> patch_vertices;

  // use all param
  read_poly_with_borders("data/elephant_quad_hole.off", poly, border_reps);
  CGAL::triangulate_refine_and_fair_hole<Default_solver>
  (poly, border_reps[0], back_inserter(patch_facets), back_inserter(patch_vertices), 
    CGAL::internal::Uniform_weight_fairing<Polyhedron>());
  // default weight
  read_poly_with_borders("data/elephant_quad_hole.off", poly, border_reps);
  CGAL::triangulate_refine_and_fair_hole<Default_solver>
  (poly, border_reps[0], back_inserter(patch_facets), back_inserter(patch_vertices));
  // default solver
  read_poly_with_borders("data/elephant_quad_hole.off", poly, border_reps);
  CGAL::triangulate_refine_and_fair_hole
    (poly, border_reps[0], back_inserter(patch_facets), back_inserter(patch_vertices),
    CGAL::internal::Uniform_weight_fairing<Polyhedron>());
  // default solver and weight
  read_poly_with_borders("data/elephant_quad_hole.off", poly, border_reps);
  CGAL::triangulate_refine_and_fair_hole<Default_solver>
    (poly, border_reps[0], back_inserter(patch_facets), back_inserter(patch_vertices));
}

int main() {
  std::vector<std::string> input_files;
  input_files.push_back("data/elephant_triangle_hole.off");
  input_files.push_back("data/elephant_quad_hole.off");
  input_files.push_back("data/mech-holes-shark.off");

  for(std::vector<std::string>::iterator it = input_files.begin(); it != input_files.end(); ++it) {
    test_triangulate_hole(it->c_str());
    test_triangulate_and_refine_hole(it->c_str());
    test_triangulate_refine_and_fair_hole(it->c_str());
    test_ouput_iterators_triangulate_and_refine_hole(it->c_str());
    test_ouput_iterators_triangulate_hole(it->c_str());
    std::cerr << "------------------------------------------------" << std::endl;
  }

  test_triangulate_refine_and_fair_hole_compile();
  std::cerr << "All Done!" << std::endl;
}