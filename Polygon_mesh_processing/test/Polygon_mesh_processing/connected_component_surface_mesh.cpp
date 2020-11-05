#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/property_map.h>

#include <boost/property_map/function_property_map.hpp>

#include <iostream>
#include <fstream>
#include <cstring>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Simple_cartesian<double>     Kernel;
typedef Kernel::Point_3                    Point;
typedef Kernel::Compare_dihedral_angle_3   Compare_dihedral_angle_3;
typedef CGAL::Surface_mesh<Point>          Mesh;

template <typename G, typename GT>
struct Constraint
  : public boost::put_get_helper<bool, Constraint<G, GT> >
{
  typedef typename GT::FT                                  FT;

  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
  typedef boost::readable_property_map_tag                 category;
  typedef bool                                             value_type;
  typedef bool                                             reference;
  typedef edge_descriptor                                  key_type;

  Constraint() : g(nullptr), gt(nullptr), bound(0) {}
  Constraint(const G& g, const GT& gt, const FT bound) : g(&g), gt(&gt), bound(bound) {}

  bool operator[](const edge_descriptor e) const
  {
    const Mesh& rg = *g;

    return gt->compare_dihedral_angle_3_object()(
             rg.point(source(e, rg)),
             rg.point(target(e, rg)),
             rg.point(target(next(halfedge(e, rg), rg), rg)),
             rg.point(target(next(opposite(halfedge(e, rg), rg), rg), rg)),
          bound) == CGAL::SMALLER;
  }

  const G* g;
  const GT* gt;
  FT bound;
};

template <typename G, typename GT>
struct Face_descriptor_area_functor
{
  typedef typename boost::graph_traits<G>::face_descriptor     face_descriptor;

  Face_descriptor_area_functor(const G& g, const GT& gt) : g(g), gt(gt) { }

  typename GT::FT operator()(const face_descriptor f) const
  {
    const auto& vpm = get(CGAL::vertex_point, g);

    return gt.compute_area_3_object()(get(vpm, source(halfedge(f, g), g)),
                                      get(vpm, target(halfedge(f, g), g)),
                                      get(vpm, target(next(halfedge(f, g), g), g)));
  }

  const G& g;
  const GT& gt;
};

void test_CC_with_default_size_map(Mesh sm,
                                   const Kernel& k)
{
  std::cout << " -- test with default size map -- " << std::endl;

  typedef boost::graph_traits<Mesh>::face_descriptor                      face_descriptor;
  typedef Kernel::FT                                                      FT;

  const FT bound = std::cos(0.7 * CGAL_PI);

  std::vector<face_descriptor> cc;
  face_descriptor fd = *faces(sm).first;
  CGAL::Polygon_mesh_processing::connected_component(fd, sm, std::back_inserter(cc));

  std::cerr << "connected components without edge constraints" << std::endl;
  std::cerr << cc.size() << " faces in the CC of " << fd << std::endl;
  assert(cc.size() == 1452);

  std::cerr << "\nconnected components with edge constraints (dihedral angle < 3/4 pi)" << std::endl;
  Mesh::Property_map<face_descriptor,std::size_t> fccmap;
  fccmap = sm.add_property_map<face_descriptor, std::size_t>("f:CC").first;
  std::size_t num = PMP::connected_components(sm, fccmap);

  std::cerr << "The graph has " << num << " connected components (face connectivity)" << std::endl;
  assert(num == 3);

  std::vector<face_descriptor> one_face_per_cc(num);
  std::vector<std::size_t> cc_size(num,0);

  for(face_descriptor f : faces(sm))
  {
    //    std::cout  << f << " in connected component " << fccmap[f] << std::endl;
    std::size_t ccid=fccmap[f];
    if (++cc_size[ccid]==1)
      one_face_per_cc[ccid]=f;
  }

  std::size_t id_of_cc_to_remove = std::distance(cc_size.begin(),
                                                 std::min_element(cc_size.begin(), cc_size.end()));

  Mesh copy1 = sm;
  Mesh copy2 = sm;

  // remove cc from copy1
  std::vector<face_descriptor> ff;
  for (std::size_t i=0;i<num;++i)
    if (i!=id_of_cc_to_remove)
      ff.push_back(one_face_per_cc[i]);

  // default face size map, but explicitely passed
  PMP::keep_connected_components(copy1, ff,
     PMP::parameters::edge_is_constrained_map(Constraint<Mesh, Kernel>(copy1, k, bound))
                     .face_size_map(CGAL::Constant_property_map<face_descriptor, std::size_t>(1)));

  // remove cc from copy2
  ff.clear();
  ff.push_back(one_face_per_cc[id_of_cc_to_remove]);
  PMP::remove_connected_components(copy2, ff,
     PMP::parameters::edge_is_constrained_map(Constraint<Mesh, Kernel>(copy2, k, bound)));

  std::cerr << "We keep the " << num-1 << " largest components" << std::endl;
  PMP::keep_largest_connected_components(sm, num-1,
    PMP::parameters::edge_is_constrained_map(Constraint<Mesh, Kernel>(sm, k, bound)));

  sm.collect_garbage();
  copy1.collect_garbage();
  copy2.collect_garbage();

  assert( num_vertices(sm)==num_vertices(copy1) && num_vertices(copy1)==num_vertices(copy2) );
  assert( num_edges(sm)==num_edges(copy1) && num_edges(copy1)==num_edges(copy2) );
  assert( num_faces(sm)==num_faces(copy1) && num_faces(copy1)==num_faces(copy2) );

  {
    Mesh m;
    Point p(0,0,0), q(1,0,0), r(0,1,0), s(0,0,1);
    CGAL::make_tetrahedron(p,q,r,s,m);
    CGAL::make_triangle(p,q,r,m);
    CGAL::make_tetrahedron(p,q,r,s,m);
    PMP::keep_large_connected_components(m, 4);

    assert(vertices(m).size() == 8);
  }
}

void test_CC_with_area_size_map(Mesh sm,
                                const Kernel& k)
{
  std::cout << " -- test with area size map -- " << std::endl;

  typedef boost::graph_traits<Mesh>::face_descriptor                      face_descriptor;

  Face_descriptor_area_functor<Mesh, Kernel> f(sm, k);
  std::vector<face_descriptor> faces_to_remove; // faces that would be removed but are not because we're doing a dry run
  std::size_t nv = num_vertices(sm);
  std::size_t num = PMP::internal::number_of_connected_components(sm);

  std::cout << "We keep the " << 2 << " largest components" << std::endl;
  std::size_t res = PMP::keep_largest_connected_components(sm, 2,
                                                           PMP::parameters::face_size_map(boost::make_function_property_map<face_descriptor>(f))
                                                                           .dry_run(true)
                                                                           .output_iterator(std::back_inserter(faces_to_remove)));

  // didn't actually remove anything
  assert(PMP::internal::number_of_connected_components(sm) == num);
  assert(num_vertices(sm) == nv);

  if(num > 2)
  {
    assert(res == num - 2);
    assert(!faces_to_remove.empty());
  }

  PMP::keep_largest_connected_components(sm, 2,
                                         PMP::parameters::face_size_map(
                                           boost::make_function_property_map<face_descriptor>(f)));
  assert(vertices(sm).size() == 1459);

  {
    Mesh m;

    Point p(0,0,0), q(1,0,0), r(0,1,0), s(0,0,1);
    CGAL::make_tetrahedron(p,q,r,s,m);
    CGAL::make_tetrahedron(p,q,r,s,m);

    Point t(100,100,100);
    CGAL::make_triangle(p,q,t,m);

    Face_descriptor_area_functor<Mesh, Kernel> f(m, k);
    PMP::keep_large_connected_components(m, 10,
                                         CGAL::parameters::face_size_map(
                                           boost::make_function_property_map<face_descriptor>(f)));
    assert(vertices(m).size() == 3);

    PMP::keep_largest_connected_components(m, 1);
    assert(PMP::internal::number_of_connected_components(m) == 1);
    PMP::keep_largest_connected_components(m, 0);
    assert(is_empty(m));
    assert(PMP::internal::number_of_connected_components(m) == 0);
  }
}

int main(int /*argc*/, char** /*argv*/)
{
  const char* filename = "data/blobby_3cc.off";
  Mesh sm;
  std::ifstream in(filename);
  assert(in.good());
  in >> sm;

  Kernel k;

  std::cout << "VEF " << num_vertices(sm) << " " << num_edges(sm) << " " << num_faces(sm) << "\n";

  test_CC_with_default_size_map(sm, k);
  test_CC_with_area_size_map(sm, k);

  return EXIT_SUCCESS;
}
