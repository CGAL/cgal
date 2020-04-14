#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/optimal_bounding_box.h>

#include <boost/property_map/function_property_map.hpp>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
typedef K::Point_3                                             Point;
typedef K::Aff_transformation_3                                Aff_transformation;

typedef CGAL::Surface_mesh<Point>                              Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor   vertex_descriptor;

struct Aff_tr_fct
{
  Aff_tr_fct() : m_at(nullptr), m_sm(nullptr) { }
  Aff_tr_fct(const Aff_transformation& at, const Surface_mesh& sm) : m_at(&at), m_sm(&sm) { }

  Point operator()(const vertex_descriptor v) const { return m_at->transform(m_sm->point(v)); }

private:
  const Aff_transformation* m_at;
  const Surface_mesh* m_sm;
};

int main(int argc, char** argv)
{
  std::ifstream input((argc > 1) ? argv[1] : "data/pig.off");
  Surface_mesh sm;
  if (!input || !(input >> sm) || sm.is_empty())
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  // get the transformation that yields the optimal bounding box
  Aff_transformation at;
  CGAL::oriented_bounding_box(sm, at);

  // functor to apply the affine transformation to a vertex of the mesh
  Aff_tr_fct aff_tr_fct(at, sm);
  auto aff_tr_vpm = boost::make_function_property_map<vertex_descriptor>(aff_tr_fct);

  // rotated AABB tree
  typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh, decltype(aff_tr_vpm)> AABB_face_graph_primitive;
  typedef CGAL::AABB_traits<K, AABB_face_graph_primitive>                              AABB_face_graph_traits;

  CGAL::AABB_tree<AABB_face_graph_traits> tree(faces(sm).begin(), faces(sm).end(), sm, aff_tr_vpm);

  return EXIT_SUCCESS;
}
