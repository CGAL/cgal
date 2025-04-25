
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <fstream>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

namespace PMP = CGAL::Polygon_mesh_processing;

struct My_visitor : public PMP::Autorefinement::Default_visitor
{
  My_visitor(std::size_t nb_input, std::size_t expected_nb_output)
    : nb_input(nb_input)
    , expected_nb_output(expected_nb_output)
  {}

  ~My_visitor()
  {
    for(std::size_t i=0; i<tgt_check.size(); ++i)
    {
      assert( tgt_check[i]==1 );
    }
  }
  void number_of_output_triangles(std::size_t nbt)
  {
    tgt_check.assign(expected_nb_output, 0);
    std::cout << nbt << " " << expected_nb_output << std::endl;
    assert(nbt==expected_nb_output);
  }

  void verbatim_triangle_copy(std::size_t tgt_id, std::size_t src_id)
  {
    assert(src_id<nb_input);
    assert(tgt_id<expected_nb_output);
    assert(tgt_check.size()==expected_nb_output && tgt_check[tgt_id]==0);
    tgt_check[tgt_id]=1;
  }

  void new_subtriangle(std::size_t tgt_id, std::size_t src_id)
  {
    assert(src_id<nb_input);
    assert(tgt_id<expected_nb_output);
    assert(tgt_check.size()==expected_nb_output && tgt_check[tgt_id]==0);
    tgt_check[tgt_id]=1;
  }

  void delete_triangle(std::size_t src_id)
  {
    assert(src_id<nb_input);
  }

  std::size_t nb_input;
  std::size_t expected_nb_output;
  std::vector<int> tgt_check;
};

template <class Tag>
void test(const char* fname, std::size_t nb_vertices_after_autorefine, std::size_t expected_nb_output, Tag tag)
{
  std::cout << "Running tests on " << fname;
  if (std::is_same_v<Tag, CGAL::Sequential_tag>)
    std::cout << " (Sequential)\n";
  else
    std::cout << " (Parallel)\n";

  std::vector<K::Point_3> points;
  std::vector< std::vector<std::size_t> > triangles;
  if (!CGAL::IO::read_polygon_soup(fname, points, triangles))
  {
    std::cerr << "  Input mesh is not a valid file." << std::endl;
    exit(EXIT_FAILURE);
  }

// Testing autorefine()
  My_visitor visitor(triangles.size(), expected_nb_output);
  PMP::autorefine_triangle_soup(points, triangles, CGAL::parameters::visitor(visitor).concurrency_tag(tag).apply_iterative_snap_rounding(true));
  std::cout << points.size() << " " << nb_vertices_after_autorefine << std::endl;
  std::cout << triangles.size() << " " << expected_nb_output << std::endl;
  assert( nb_vertices_after_autorefine==points.size());
  assert( expected_nb_output==triangles.size());
  assert( !PMP::does_triangle_soup_self_intersect(points, triangles) );
//  CGAL::IO::write_polygon_soup("/tmp/debug.off", points, triangles);
}

int main(int argc, const char** argv)
{
  // file expected_nb_of_vertices expected_nb_of_faces (after repair)
  for (int i=0;i<(argc-1)/3; ++i)
  {
    test(argv[1+3*i], atoi(argv[1+3*i+1]), atoi(argv[1+3*i+2]), CGAL::Sequential_tag());
#ifdef CGAL_LINKED_WITH_TBB
    test(argv[1+3*i], atoi(argv[1+3*i+1]), atoi(argv[1+3*i+2]), CGAL::Parallel_tag());
#endif
  }
}
