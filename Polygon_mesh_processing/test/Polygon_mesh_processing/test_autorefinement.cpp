
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <CGAL/IO/polygon_soup_io.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <fstream>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>             Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

template <class TriangleMesh>
struct My_exp_visitor :
  public CGAL::Polygon_mesh_processing::Corefinement::Default_visitor<TriangleMesh>
{
  void after_subface_creations(TriangleMesh&){++(*i);}

  My_exp_visitor()
    : i (new int(0) )
  {}

  std::shared_ptr<int> i;
};

struct My_visitor
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

  std::size_t nb_input;
  std::size_t expected_nb_output;
  std::vector<int> tgt_check;
};

void test_coref_based(const char* fname, std::size_t nb_polylines, std::size_t total_nb_points,
                      std::size_t nb_vertices_after_autorefine, bool all_fixed, std::size_t nb_vertices_after_fix,
                      bool triple_intersection)
{
  std::cout << "Running tests (coref based) on " << fname << "\n";
  std::ifstream input(fname);

  Mesh mesh;
  if (!input || !(input >> mesh))
  {
    std::cerr << "  Input mesh is not a valid off file." << std::endl;
    exit(EXIT_FAILURE);
  }
  input.close();
  std::size_t nb_vertices_before_autorefine = num_vertices(mesh);

// Testing PMP::experimental::surface_self_intersection()
  try{
    std::vector< std::vector<K::Point_3> >polylines;
    PMP::experimental::surface_self_intersection(mesh, std::back_inserter(polylines));
    assert(polylines.size() == nb_polylines);
    std::size_t total_nb_pt=0;
    for(const std::vector<K::Point_3>& polyline : polylines)
      total_nb_pt+=polyline.size();
    assert(total_nb_points == total_nb_pt);
    assert( !triple_intersection );
  }
  catch(const PMP::Corefinement::Triple_intersection_exception&)
  {
    assert( triple_intersection );
  }

// Testing PMP::experimental::autorefine()
  try{
    My_exp_visitor<Mesh> visitor;
    PMP::experimental::autorefine(mesh,
      CGAL::parameters::visitor(visitor));
    mesh.collect_garbage();
    assert( nb_vertices_after_autorefine==num_vertices(mesh));
    assert( (nb_vertices_before_autorefine!=nb_vertices_after_autorefine)== (*(visitor.i) != 0) );
    assert( !triple_intersection );
  }
  catch(const PMP::Corefinement::Triple_intersection_exception&)
  {
    assert( triple_intersection );
  }

// Testing PMP::experimental::autorefine_and_remove_self_intersections()
  try{
    input.open(fname);
    mesh.clear();
    input >> mesh;
    bool res=PMP::experimental::autorefine_and_remove_self_intersections(mesh);
    assert(res==all_fixed);
    mesh.collect_garbage();
    assert( nb_vertices_after_fix==num_vertices(mesh));
    assert( !triple_intersection );
  }
  catch(const PMP::Corefinement::Triple_intersection_exception&)
  {
    assert( triple_intersection );
  }
}

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
  PMP::autorefine_triangle_soup(points, triangles, CGAL::parameters::visitor(visitor).concurrency_tag(tag));
  assert( nb_vertices_after_autorefine==points.size());
  assert( expected_nb_output==triangles.size());
  assert( !PMP::does_triangle_soup_self_intersect(points, triangles) );
//  CGAL::IO::write_polygon_soup("/tmp/debug.off", points, triangles);
}

int main(int argc, const char** argv)
{
  // file nb_polylines total_nb_points nb_vertices_after_autorefine all_fixed nb_vertices_after_fix triple_intersection
  for (int i=0;i<(argc-1)/9; ++i)
  {
    test_coref_based(argv[1+9*i], atoi(argv[1+9*i+1]), atoi(argv[1+9*i+2]),
                     atoi(argv[1+9*i+3]), atoi(argv[1+9*i+4])==0?false:true, atoi(argv[1+9*i+5]), atoi(argv[1+9*i+6])==0?false:true);
    test(argv[1+9*i], atoi(argv[1+9*i+7]), atoi(argv[1+9*i+8]), CGAL::Sequential_tag());
#ifdef CGAL_LINKED_WITH_TBB
    test(argv[1+9*i], atoi(argv[1+9*i+7]), atoi(argv[1+9*i+8]), CGAL::Parallel_tag());
#endif
  }
}
