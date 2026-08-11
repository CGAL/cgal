//#define POLY

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>


#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Weights/uniform_weights.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

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

void read_poly(const std::string file_name, Polyhedron& poly) {
  poly.clear();

  std::ifstream input(file_name);
  if ( !input || !(input >> poly)  || (num_vertices(poly) == 0)){
    std::cerr << "  Error: can not read file." << std::endl;
    assert(false);
  }
}

void detect_borders(Polyhedron& poly, std::vector<Halfedge_handle>& border_reps)
{
  border_reps.clear();
  std::set<Halfedge_handle> border_map;
  for(Halfedge_handle h :  halfedges(poly)){
    if(face(h,poly)== boost::graph_traits<Polyhedron>::null_face() && border_map.find(h) == border_map.end()){
      border_reps.push_back(h);
      Halfedge_around_facet_circulator hf_around_facet(h,poly), done(hf_around_facet);
      do {
        bool insertion_ok = border_map.insert(*hf_around_facet).second;
        assert(insertion_ok);
      } while(++hf_around_facet != done);
    }
  }
}

void read_poly_with_borders(const std::string file_name, Polyhedron& poly, std::vector<Halfedge_handle>& border_reps)
{
  read_poly(file_name, poly);
  detect_borders(poly, border_reps);
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
                                     CGAL::approximate_dihedral_angle(ppmap[target(edge_it,poly)],
                                                                      ppmap[source(edge_it,poly)],
                                                                      ppmap[target(next(edge_it,poly),poly)],
                                                                      ppmap[target(next(opposite(edge_it,poly),poly),poly)]) );
      edge_it = next(edge_it,poly);
      ang_max = (std::max)(angle, ang_max);
    }

    double area = CGAL::sqrt(CGAL::squared_area(ppmap[target(edge_it,poly)],
                                               ppmap[target(next(edge_it,poly),poly)],
                                               ppmap[target(prev(edge_it,poly),poly)]));
    res = res + Weight(ang_max,area);
  }
  return res;
}


void test_triangulate_hole_weight(const std::string file_name, bool use_DT, std::size_t nb_remaining_holes) { //don't test with cdt as we are testing the weights and there is no weight in the cdt version
  typedef CGAL::internal::Weight_min_max_dihedral_and_area Weight;

  std::cout << "test_triangulate_hole_weight + useDT: " << use_DT << std::endl;
  std::cout << "  File: "<< file_name  << std::endl;
  Polyhedron poly;
  std::vector<Halfedge_handle> border_reps;
  read_poly_with_borders(file_name, poly, border_reps);

  for(std::vector<Halfedge_handle>::iterator it = border_reps.begin(); it != border_reps.end(); ++it) {
    std::vector<Facet_handle> patch;
    CGAL::Polygon_mesh_processing::Hole_filling::Default_visitor visitor;
    Weight w_algo = CGAL::Polygon_mesh_processing::internal::triangulate_hole_polygon_mesh(
      poly, *it, back_inserter(patch), get(CGAL::vertex_point, poly), use_DT, Kernel(), false/*use_cdt*/, false/*skip_cubic*/, visitor, 0).second;
    if(patch.empty()) { continue; }
    Weight w_test = calculate_weight_for_patch(poly, patch.begin(), patch.end());

    const double epsilon = 1e-10;
    if( std::abs(w_algo.w.first - w_test.w.first) > epsilon ||
        std::abs(w_algo.w.second - w_test.w.second) > epsilon )
    {
      std::cerr << "  Weight returned by algo   : " << w_algo << std::endl;
      std::cerr << "  Weight calculated by test : " << w_test << std::endl;
      assert(false);
    }
  }

  detect_borders(poly, border_reps);
  assert(border_reps.size()==nb_remaining_holes);
  std::cout << "  Done!" << std::endl;
}
/******************************************************************/

void test_triangulate_hole(const std::string file_name, bool use_cdt) {
  std::cout << "test_triangulate_hole:" << std::endl;
  std::cout << "  File: "<< file_name  << std::endl;
  Polyhedron poly;
  std::vector<Halfedge_handle> border_reps;
  read_poly_with_borders(file_name, poly, border_reps);

  for(std::vector<Halfedge_handle>::iterator it = border_reps.begin(); it != border_reps.end(); ++it) {
    std::vector<Facet_handle> patch;
    CGAL::Polygon_mesh_processing::triangulate_hole(poly, *it,
                                                    CGAL::parameters::
                                                      face_output_iterator(std::back_inserter(patch)).
                                                      use_2d_constrained_delaunay_triangulation(use_cdt));
    if(patch.empty()) {
      std::cerr << "  Error: empty patch created." << std::endl;
      assert(false);
    }
  }

  if(!poly.is_valid() || ! is_closed(poly)) {
    std::cerr << "  Error: patched polyhedron is not valid or closed." << std::endl;
    assert(false);
  }

  std::cout << "  Done!" << std::endl;
}

void test_triangulate_hole_should_be_no_output(const std::string file_name, bool use_cdt) {
  std::cout << "test_triangulate_hole_should_be_no_output:" << std::endl;
  std::cout << "  File: "<< file_name  << std::endl;
  Polyhedron poly;
  std::vector<Halfedge_handle> border_reps;
  read_poly_with_borders(file_name, poly, border_reps);

  for(std::vector<Halfedge_handle>::iterator it = border_reps.begin(); it != border_reps.end(); ++it) {
    std::vector<Facet_handle> patch;
    CGAL::Polygon_mesh_processing::triangulate_hole(poly, *it,
      CGAL::parameters::use_delaunay_triangulation(false)
        .face_output_iterator(back_inserter(patch))
        .use_2d_constrained_delaunay_triangulation(use_cdt));
    if(!patch.empty()) {
      std::cerr << "  Error: patch should be empty" << std::endl;
      assert(false);
    }

    CGAL::Polygon_mesh_processing::triangulate_hole(poly, *it,
      CGAL::parameters::use_delaunay_triangulation(true).face_output_iterator(back_inserter(patch)));
    if(!patch.empty()) {
      std::cerr << "  Error: patch should be empty" << std::endl;
      assert(false);
    }
  }

  std::cout << "  Done!" << std::endl;
}

void test_triangulate_and_refine_hole(const std::string file_name, bool use_cdt) {
  std::cout << "test_triangulate_and_refine_hole:" << std::endl;
  std::cout << "  File: "<< file_name  << std::endl;
  Polyhedron poly;
  std::vector<Halfedge_handle> border_reps;
  read_poly_with_borders(file_name, poly, border_reps);

  for(std::vector<Halfedge_handle>::iterator it = border_reps.begin(); it != border_reps.end(); ++it) {
    std::vector<Facet_handle> patch_facets;
    std::vector<Vertex_handle> patch_vertices;
    CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(poly, *it,
      CGAL::parameters::
        face_output_iterator(std::back_inserter(patch_facets)).
        vertex_output_iterator(std::back_inserter(patch_vertices)).
        use_2d_constrained_delaunay_triangulation(use_cdt));

    if(patch_facets.empty()) {
      std::cerr << "  Error: empty patch created." << std::endl;
      assert(false);
    }
  }

  if(!poly.is_valid() || ! is_closed(poly)) {
    std::cerr << "  Error: patched polyhedron is not valid or closed." << std::endl;
    assert(false);
  }

  std::cout << "  Done!" << std::endl;
}

void test_triangulate_refine_and_fair_hole(const std::string file_name, bool use_cdt) {
  std::cout << "test_triangulate_refine_and_fair_hole:" << std::endl;
  std::cout << "  File: "<< file_name  << std::endl;
  Polyhedron poly;
  std::vector<Halfedge_handle> border_reps;
  read_poly_with_borders(file_name, poly, border_reps);

  for(std::vector<Halfedge_handle>::iterator it = border_reps.begin(); it != border_reps.end(); ++it) {
    std::vector<Facet_handle> patch_facets;
    std::vector<Vertex_handle> patch_vertices;
    CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(poly, *it,
      CGAL::parameters::
        face_output_iterator(std::back_inserter(patch_facets)).
        vertex_output_iterator(std::back_inserter(patch_vertices)).
        use_2d_constrained_delaunay_triangulation(use_cdt));

    if(patch_facets.empty()) {
      std::cerr << "  Error: empty patch created." << std::endl;
      assert(false);
    }
  }

  if(!poly.is_valid() || ! is_closed(poly)) {
    std::cerr << "  Error: patched polyhedron is not valid or closed." << std::endl;
    assert(false);
  }

  std::cout << "  Done!" << std::endl;
}

void test_output_iterators_triangulate_hole(const std::string file_name, bool use_cdt) {
  std::cout << "test_output_iterators_triangulate_hole:" << std::endl;
  std::cout << "  File: "<< file_name  << std::endl;

  Polyhedron poly, poly_2;
  std::vector<Halfedge_handle> border_reps, border_reps_2;

  read_poly_with_borders(file_name, poly, border_reps);
  read_poly_with_borders(file_name, poly_2, border_reps_2);

  std::vector<Halfedge_handle>::iterator it_2 = border_reps_2.begin();
  for(std::vector<Halfedge_handle>::iterator it = border_reps.begin(); it != border_reps.end(); ++it, ++it_2) {
    std::vector<Facet_handle> patch;
    CGAL::Polygon_mesh_processing::triangulate_hole(poly, *it,
      CGAL::parameters::
        face_output_iterator(std::back_inserter(patch)).
        use_2d_constrained_delaunay_triangulation(use_cdt));

    std::vector<Facet_handle> patch_2 = patch;
    Facet_handle* output_it =
      CGAL::Polygon_mesh_processing::triangulate_hole(poly_2, *it_2, CGAL::parameters::face_output_iterator(& *patch_2.begin()));

    if(patch.size() != (std::size_t)(output_it - &*patch_2.begin())) {
      std::cerr << "  Error: returned facet output iterator is not valid!" << std::endl;
      std::cerr << "  " << patch.size() << " vs " << (output_it - &*patch_2.begin()) << std::endl;

      assert(false);
    }
  }
  std::cout << "  Done!" << std::endl;
}

void test_output_iterators_triangulate_and_refine_hole(const std::string file_name, bool use_cdt) {
  std::cout << "test_output_iterators_triangulate_and_refine_hole:" << std::endl;
  std::cout << "  File: "<< file_name  << std::endl;

  Polyhedron poly, poly_2;
  std::vector<Halfedge_handle> border_reps, border_reps_2;;

  read_poly_with_borders(file_name, poly, border_reps);
  read_poly_with_borders(file_name, poly_2, border_reps_2);

  std::vector<Halfedge_handle>::iterator it_2 = border_reps_2.begin();
  for(std::vector<Halfedge_handle>::iterator it = border_reps.begin(); it != border_reps.end(); ++it, ++it_2) {
    std::vector<Facet_handle> patch_facets;
    std::vector<Vertex_handle> patch_vertices;
    CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(poly, *it,
      CGAL::parameters::
        face_output_iterator(std::back_inserter(patch_facets)).
        vertex_output_iterator(std::back_inserter(patch_vertices)).
        use_2d_constrained_delaunay_triangulation(use_cdt));
    // create enough space to hold outputs
    std::vector<Facet_handle> patch_facets_2 = patch_facets;
    std::vector<Vertex_handle> patch_vertices_2 = patch_vertices;
    if(patch_vertices_2.empty()) { patch_vertices_2.push_back(Vertex_handle()); } //just allocate space for dereferencing

    std::pair<Facet_handle*, Vertex_handle*> output_its =
      CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(poly_2, *it_2,
        CGAL::parameters::
          face_output_iterator(&*patch_facets_2.begin()).
          vertex_output_iterator(&*patch_vertices_2.begin()).
          use_2d_constrained_delaunay_triangulation(use_cdt));

    if(patch_facets.size() != (std::size_t) (output_its.first - &*patch_facets_2.begin())) {
      std::cout << "  Error: returned facet output iterator is not valid!" << std::endl;
      std::cout << "  " << patch_facets.size() << " vs " << (output_its.first - &*patch_facets_2.begin()) << std::endl;
      assert(false);
    }

    if(patch_vertices.size() != (std::size_t) (output_its.second - &*patch_vertices_2.begin())) {
      std::cerr << "  Error: returned vertex output iterator is not valid!" << std::endl;
      std::cerr << "  " << patch_vertices.size() << " vs " << (output_its.second - &*patch_vertices_2.begin()) << std::endl;
      assert(false);
    }
  }
  std::cout << "  Done!" << std::endl;
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
  read_poly_with_borders("elephant_quad_hole.off", poly, border_reps);
  CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole
  (poly, border_reps[0],
   CGAL::parameters::
     face_output_iterator(back_inserter(patch_facets)).
     vertex_output_iterator(back_inserter(patch_vertices)).
     weight_calculator(CGAL::Weights::Uniform_weight<Polyhedron>()).
     sparse_linear_solver(Default_solver()).
     use_2d_constrained_delaunay_triangulation(false));

  // default solver
  read_poly_with_borders("elephant_quad_hole.off", poly, border_reps);
  CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole
    (poly, border_reps[0],
     CGAL::parameters::
       face_output_iterator(back_inserter(patch_facets)).
       vertex_output_iterator(back_inserter(patch_vertices)).
       weight_calculator(CGAL::Weights::Uniform_weight<Polyhedron>()));

  // default solver and weight
  read_poly_with_borders("elephant_quad_hole.off", poly, border_reps);
  CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole
    (poly, border_reps[0],
     CGAL::parameters::
       face_output_iterator(back_inserter(patch_facets)).
       vertex_output_iterator(back_inserter(patch_vertices)));
}

void generate_elephant_with_hole()
{
  Polyhedron poly;
  read_poly(CGAL::data_file_path("meshes/elephant.off"), poly);
  int i=0;
  for(Facet_handle fd : faces(poly))
    if (++i==229)
    {
      Halfedge_handle nh=opposite(halfedge(fd,poly), poly);
      CGAL::Euler::remove_face(halfedge(fd, poly), poly);
      std::ofstream output("elephant_triangle_hole.off");
      output << poly;
      output.close();
      CGAL::Euler::remove_face(nh, poly);
      output.open("elephant_quad_hole.off");
      output << poly;
      return;
    }
}
// visitors for triangulate_faces
struct Triangulate_face_visitor
  : public CGAL::Polygon_mesh_processing::Triangulate_faces::Default_visitor<Polyhedron>
{
  using Mesh = CGAL::Surface_mesh<Kernel::Point_3>;
  using Vertex_index = Mesh::Vertex_index;
  using Face_index = Mesh::Face_index;

  std::array<Vertex_handle,3> forbidden_face;
  Triangulate_face_visitor(Vertex_index v0, Vertex_index v1, Vertex_index v2)
  {
    forbidden_face=CGAL::make_array(v0,v1,v2);
    std::sort(forbidden_face.begin(), forbidden_face.end());
  }

  bool accept_face(Face_index,
                   Vertex_index v0, Vertex_index v1, Vertex_index v2) const
  {
    auto a = CGAL::make_array(v0, v1, v2);
    std::sort(a.begin(), a.end());

    return a!=forbidden_face;
  }
};

//visitors for triangulate_hole_polyline
struct Triangulate_face_visitor_reject_all
  : public CGAL::Polygon_mesh_processing::Triangulate_faces::Default_visitor<Polyhedron>
{
  using Mesh = CGAL::Surface_mesh<Kernel::Point_3>;
  using Vertex_index = Mesh::Vertex_index;
  using Face_index = Mesh::Face_index;

  bool accept_face(Face_index,
                   Vertex_index , Vertex_index , Vertex_index) const
  {
    return false;
  }
};

struct Hole_filling_visitor
  : public CGAL::Polygon_mesh_processing::Hole_filling::Default_visitor
{
  std::array<int, 3> forbidden_face;
  Hole_filling_visitor(int i0, int i1, int i2)
  {
    forbidden_face=CGAL::make_array(i0,i1,i2);
    std::sort(forbidden_face.begin(), forbidden_face.end());
  }


  bool accept_triangle(int i0, int i1, int i2) const
  {
    auto a =CGAL::make_array(i0,i1,i2);
    std::sort(a.begin(), a.end());
    return a!=forbidden_face;
  }
};

struct Hole_filling_visitor_reject_all
  : public CGAL::Polygon_mesh_processing::Hole_filling::Default_visitor
{
  bool accept_triangle(int , int , int ) const
  {
    return false;
  }
};

// visitors for triangulate_polygons()
struct Triangulate_polygon_visitor
  : public CGAL::Polygon_mesh_processing::Triangulate_polygons::Default_visitor
{
  std::array<std::size_t,3> forbidden_face;
  Triangulate_polygon_visitor(std::size_t v0, std::size_t v1, std::size_t v2)
  {
    forbidden_face=CGAL::make_array(v0,v1,v2);
    std::sort(forbidden_face.begin(), forbidden_face.end());
  }

  bool accept_face(std::size_t,
                   std::size_t v0, std::size_t v1, std::size_t v2) const
  {
    auto a = CGAL::make_array(v0, v1, v2);
    std::sort(a.begin(), a.end());

    return a!=forbidden_face;
  }
};

struct Triangulate_polygon_visitor_reject_all
  : public CGAL::Polygon_mesh_processing::Triangulate_polygons::Default_visitor
{
  constexpr
  bool accept_face(std::size_t,std::size_t,std::size_t,std::size_t) const
  {
    return false;
  }
};

void test_with_forbidden_triangles()
{
  using Mesh = CGAL::Surface_mesh<Kernel::Point_3>;
  namespace PMP = CGAL::Polygon_mesh_processing;
  namespace params = CGAL::parameters;

  Mesh mesh;
  auto vpm = get(CGAL::vertex_point, mesh);
  CGAL::make_hexahedron(CGAL::Bbox_3(-1,-1,-1,1,1,1), mesh, params::do_not_triangulate_faces(true));
  auto fit = faces(mesh).begin();
  while (fit!=faces(mesh).end())
  {
    auto h = halfedge(*fit, mesh);
    if (get(vpm, source(h,mesh)).z()==1 && get(vpm, target(h,mesh)).z()==1 && get(vpm, target(next(h,mesh),mesh)).z()==1)
      break;
    ++fit;
  }

  auto hv = CGAL::Euler::add_center_vertex(halfedge(*fit, mesh), mesh);
  auto v = target(hv, mesh);
  put(vpm, v, Kernel::Point_3(0.5,0.5,1));
  std::vector<Halfedge_handle> hedges;
  for (auto h : CGAL::halfedges_around_target(hv, mesh))
    hedges.push_back(h);
  for (Halfedge_handle h : hedges)
  {
    Kernel::Point_3 mp = CGAL::midpoint(get(vpm, source(h, mesh)), get(vpm, target(h, mesh)));
    auto hnew = CGAL::Euler::split_edge(h, mesh);
    put(vpm, target(hnew, mesh), mp);
  }
  PMP::triangulate_faces(mesh);
  auto h = CGAL::Euler::remove_center_vertex(halfedge(v, mesh), mesh);

  // testing triangulate_faces
  {
    auto m1=mesh;
    m1.collect_garbage();
    auto m2=m1;
    auto m3=m1;
    auto m4=m1;

    //look for the quad face
    Mesh::Halfedge_index hq;
    for (Mesh::Halfedge_index h : halfedges(m1))
    {
      if (!CGAL::is_triangle(h, m1))
      {
        hq=h;
      }
    }

    std::pair diag1(source(hq,m1), target(next(hq,m1), m1));
    std::pair diag2(target(hq,m2), source(prev(hq,m2), m2));
    Triangulate_face_visitor vis1(source(hq,m1), target(hq,m1), target(next(hq,m1), m1));
    Triangulate_face_visitor vis2(source(hq,m2), target(hq,m2), source(prev(hq,m2), m2));

    PMP::triangulate_faces(m1, params::visitor(vis1));
    PMP::triangulate_faces(m2, params::visitor(vis2));
    PMP::triangulate_faces(m3, params::visitor(Triangulate_face_visitor_reject_all()));

    assert(CGAL::is_triangle_mesh(m1));
    assert(CGAL::is_triangle_mesh(m2));
    assert(!CGAL::is_triangle_mesh(m3));

    assert(!halfedge(diag1.first, diag1.second, m1).second);
    assert(!halfedge(diag2.first, diag2.second, m2).second);

    auto mp = CGAL::midpoint(m4.point(source(hq,m4)), m4.point(target(hq, m4)));
    m4.point(target(CGAL::Euler::split_edge(hq, m4), m4))=mp;

    auto m5=m4;
    auto m6=m4;

    PMP::triangulate_faces(m4);

    Triangulate_face_visitor vis5(source(hq,m4), target(hq,m4), target(next(hq,m4), m4)); // m4 on purpose
    std::pair diag5(source(hq,m4), source(prev(hq,m4), m4));
    PMP::triangulate_faces(m5, params::visitor(vis5));
    PMP::triangulate_faces(m6, params::visitor(Triangulate_face_visitor_reject_all()));
    assert(CGAL::is_triangle_mesh(m4));
    assert(CGAL::is_triangle_mesh(m5));
    assert(!CGAL::is_triangle_mesh(m6));
    assert(!halfedge(diag5.first, diag5.second, m5).second);
  }

  CGAL::Euler::remove_face(h, mesh);
  mesh.collect_garbage();

  for (auto h2 : halfedges(mesh))
  {
    if (is_border(h2, mesh))
    {
      h=h2;
      break;
    }
  }

  std::array< std::pair<bool, bool>, 4 > configs = {{ {true, true}, {true, false}, {false, true}, {false, false}}};

  for (auto c : configs)
  {
    auto mesh0 = mesh;
    auto mesh1 = mesh;
    auto mesh2 = mesh;
    auto mesh3 = mesh;
    std::pair diag(target(h, mesh), target(next(next(h,mesh),mesh),mesh));
    std::pair diag1(source(h, mesh), target(next(h,mesh),mesh));
  // TODO: test that the visitor is using the vertices as documented!
    PMP::triangulate_hole(mesh0, h, params::visitor(Hole_filling_visitor(0,1,2)).use_delaunay_triangulation(c.first).use_2d_constrained_delaunay_triangulation(c.second));
    PMP::triangulate_hole(mesh1, h, params::visitor(Hole_filling_visitor(0,1,3)).use_delaunay_triangulation(c.first).use_2d_constrained_delaunay_triangulation(c.second));
    PMP::triangulate_hole(mesh2, h, params::visitor(Hole_filling_visitor_reject_all()).use_delaunay_triangulation(c.first).use_2d_constrained_delaunay_triangulation(c.second));

    assert(CGAL::is_closed(mesh0));
    assert(CGAL::is_closed(mesh1));
    assert(!CGAL::is_closed(mesh2));
    assert(!halfedge(diag.first, diag.second, mesh0).second);
    assert(!halfedge(diag1.first, diag1.second, mesh1).second);

    auto mp = CGAL::midpoint(mesh3.point(source(h,mesh3)), mesh3.point(target(h, mesh3)));
    mesh3.point(target(CGAL::Euler::split_edge(h, mesh3), mesh3))=mp;

    auto mesh4=mesh3;
    auto mesh5=mesh3;

    std::map<Mesh::Vertex_index, int> vmap;
    int i=0;
    for (auto haf : CGAL::halfedges_around_face(h, mesh3))
      vmap[target(haf,mesh3)]=i++;

    PMP::triangulate_hole(mesh3, h, params::use_delaunay_triangulation(c.first).use_2d_constrained_delaunay_triangulation(c.second));
    Hole_filling_visitor vis4(vmap.at(source(h,mesh3)), vmap.at(target(h,mesh3)), vmap.at(target(next(h,mesh3), mesh3))); // mesh3 on purpose
    std::pair diag4(source(h,mesh3), source(prev(h,mesh3), mesh3));

    PMP::triangulate_hole(mesh4, h, params::visitor(vis4).use_delaunay_triangulation(c.first).use_2d_constrained_delaunay_triangulation(c.second));
    PMP::triangulate_hole(mesh5, h, params::visitor(Hole_filling_visitor_reject_all()).use_delaunay_triangulation(c.first).use_2d_constrained_delaunay_triangulation(c.second));

    assert(CGAL::is_closed(mesh3));
    assert(CGAL::is_closed(mesh4));
    assert(!CGAL::is_closed(mesh5));
    assert(!halfedge(diag4.first, diag4.second, mesh4).second);
  }

  // now test triangulate_hole_polyline
  using P = Kernel::Point_3;
  std::vector<P> quad = {P(0,0,0), P(1,0,0), P(1,1,0), P(0,1,0)};
  std::vector<std::tuple<int,int,int> > triangles1, triangles2, triangles3;
  PMP::triangulate_hole_polyline(quad, std::back_inserter(triangles1), params::visitor(Hole_filling_visitor(0,1,2)));
  PMP::triangulate_hole_polyline(quad, std::back_inserter(triangles2), params::visitor(Hole_filling_visitor(0,1,3)));
  PMP::triangulate_hole_polyline(quad, std::back_inserter(triangles3), params::visitor(Hole_filling_visitor_reject_all()));

  auto triangle_not_present = [](std::tuple<int,int,int> t, int i0, int i1, int i2)
  {
    auto a = CGAL::make_array(get<0>(t), get<1>(t), get<2>(t));
    std::sort(a.begin(), a.end());
    return a[0]!=i0 || a[1]!=i1 || a[2]!=i2;
  };

  assert(triangles1.size()==2);
  assert(triangles2.size()==2);
  assert(triangles3.empty());

  assert(triangle_not_present(triangles1[0], 0, 1, 2) && triangle_not_present(triangles1[1], 0, 1, 2));
  assert(triangle_not_present(triangles2[0], 0, 1, 3) && triangle_not_present(triangles2[1], 0, 1, 3));

  std::vector<P> octo = {P(0,0,0), P(0.5,-0.5,0), P(1,0,0), P(1.5,0.5,0), P(1,1,0), P(0.5, 1.5, 0), P(0,1,0), P(-0.5,0.5,0)};
  std::vector<std::tuple<int,int,int> > triangles4, triangles5, triangles6;
  PMP::triangulate_hole_polyline(octo, std::back_inserter(triangles4));
  assert(triangles4.size()==6);
  auto [i0,i1,i2] = triangles4[0];
  PMP::triangulate_hole_polyline(octo, std::back_inserter(triangles5), params::visitor(Hole_filling_visitor(i0,i1,i2)));
  PMP::triangulate_hole_polyline(octo, std::back_inserter(triangles6), params::visitor(Hole_filling_visitor_reject_all()));
  assert(triangles5.size()==6);
  assert(triangles6.empty());
  if (i2<i0) std::swap(i0,i2);
  if (i1<i0) std::swap(i0,i1);
  if (i2<i1) std::swap(i2,i1);
  assert(triangle_not_present(triangles1[0], i0, i1, i2) && triangle_not_present(triangles1[1], i0, i1, i2));

  // finally test triangulate_polygons
  std::vector<P> points_quad = { P(-1,0,0), P(0,0,0), P(1,0,0), P(1,1,0), P(0,1,0)};
  std::vector<std::vector<std::size_t>> polygons_quad1 = {{0,1,4}, {1,2,3,4}}, polygons_quad2=polygons_quad1, polygons_quad3=polygons_quad1;
  PMP::triangulate_polygons(points_quad, polygons_quad1, params::visitor(Triangulate_polygon_visitor(1,2,3)));
  PMP::triangulate_polygons(points_quad, polygons_quad2, params::visitor(Triangulate_polygon_visitor(1,2,4)));
  PMP::triangulate_polygons(points_quad, polygons_quad3, params::visitor(Triangulate_polygon_visitor_reject_all()));

  assert(polygons_quad1.size()==3);
  assert(polygons_quad2.size()==3);
  assert(polygons_quad3.size()==2);

  auto polygon_not_present = [](std::vector<std::size_t> v, std::size_t i0, std::size_t i1, std::size_t i2)
  {
    std::sort(v.begin(), v.end());
    return v[0]!=i0 || v[1]!=i1 || v[2]!=i2;
  };

  assert(polygon_not_present(polygons_quad1[0], 1, 2, 3) && polygon_not_present(polygons_quad1[1], 1, 2, 3));
  assert(polygon_not_present(polygons_quad2[0], 1, 2, 4) && polygon_not_present(polygons_quad2[1], 1, 2, 4));

  std::vector<P> points_octo = {P(-1,0,0), P(0,0,0), P(0.5,-0.5,0), P(1,0,0), P(1.5,0.5,0), P(1,1,0), P(0.5, 1.5, 0), P(0,1,0), P(-0.5,0.5,0)};
  std::vector<std::vector<std::size_t>> polygons_octo1 ={{0,1,7}, {1,2,3,4,5,6,7,8}}, polygons_octo2=polygons_octo1, polygons_octo3=polygons_octo1;


  PMP::triangulate_polygons(points_octo, polygons_octo1, params::visitor(Triangulate_polygon_visitor(1,2,3)));
  assert(polygons_octo1.size()==7);
  auto vals = polygons_octo1[1];
  std::sort(vals.begin(), vals.end());
  PMP::triangulate_polygons(points_octo, polygons_octo2, params::visitor(Triangulate_polygon_visitor(vals[0],vals[1],vals[2])));
  PMP::triangulate_polygons(points_octo, polygons_octo3, params::visitor(Triangulate_polygon_visitor_reject_all()));
  assert(polygons_octo2.size()==7);
  assert(polygons_octo3.size()==2);
  for (int i=0;i<7;++i)
  {
    assert(polygon_not_present(polygons_octo2[i], vals[0],vals[1],vals[2]));
  }

}


int main()
{
  std::cout.precision(17);
  std::cerr.precision(17);
  test_with_forbidden_triangles();
  generate_elephant_with_hole();
  std::vector<std::string> input_files;
  input_files.push_back("elephant_triangle_hole.off");
  input_files.push_back("elephant_quad_hole.off");
  input_files.push_back(CGAL::data_file_path("meshes/mech-holes-shark.off"));
  // std::cerr.precision(15);
  for(std::vector<std::string>::iterator it = input_files.begin(); it != input_files.end(); ++it) {
    test_triangulate_hole(it->c_str(), true);
    test_triangulate_hole(it->c_str(), false);
    test_triangulate_and_refine_hole(it->c_str(), true);
    test_triangulate_and_refine_hole(it->c_str(), false);
    test_triangulate_refine_and_fair_hole(it->c_str(), true);
    test_triangulate_refine_and_fair_hole(it->c_str(), false);
    test_output_iterators_triangulate_and_refine_hole(it->c_str(), true);
    test_output_iterators_triangulate_and_refine_hole(it->c_str(), false);
    test_output_iterators_triangulate_hole(it->c_str(), true);
    test_output_iterators_triangulate_hole(it->c_str(), false);
    test_triangulate_hole_weight(it->c_str(), true, 0);
    test_triangulate_hole_weight(it->c_str(), false, 0);
    std::cout << "------------------------------------------------" << std::endl;
  }
  test_triangulate_hole_should_be_no_output("data-repair/non_manifold_vertex.off", true);
  test_triangulate_hole_should_be_no_output("data-repair/non_manifold_vertex.off", false);
  test_triangulate_hole_should_be_no_output("data-repair/two_tris_collinear.off", true);
  test_triangulate_hole_should_be_no_output("data-repair/two_tris_collinear.off", false);

  test_triangulate_refine_and_fair_hole_compile();
  std::cout << "All Done!" << std::endl;
}
