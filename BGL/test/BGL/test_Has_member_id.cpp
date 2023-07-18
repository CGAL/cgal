#include <CGAL/assertions.h>
#include <CGAL/boost/graph/internal/Has_member_id.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Triangulation_vertex_base_with_id_2.h>
#include <CGAL/Triangulation_face_base_with_id_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3> Polyhedron_with_ids;

struct StructWithId
{
  std::size_t m_id;
  StructWithId() : m_id(12) {}
  //with overload
  std::size_t id() const { return m_id; }
  std::size_t& id() { return m_id; }
};
struct StructNoId
{};

int main()
{
  using namespace CGAL::internal;

  static_assert(!Has_member_id<StructNoId>::value);
  static_assert(Has_member_id<StructWithId>::value);
  static_assert(!Has_member_id<Polyhedron::Face>::value);
  static_assert(Has_member_id<Polyhedron_with_ids::Facet>::value);
  static_assert(Has_member_id<Polyhedron_with_ids::FBase>::value);
  static_assert(
    (Has_member_id<Polyhedron_with_ids::Items::Face_wrapper<Polyhedron_with_ids::HDS, K>::Face>::value));

  static_assert(!Has_member_id<CGAL::Triangulation_face_base_2<K> >::value);
  static_assert(Has_member_id<CGAL::Triangulation_vertex_base_with_id_2<K> >::value);
  static_assert(Has_member_id<CGAL::Triangulation_face_base_with_id_2<K> >::value);

  return 0;
}
