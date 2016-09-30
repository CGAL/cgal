#include <CGAL/assertions.h>
#include <CGAL/boost/graph/internal/Has_member_id.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

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

  CGAL_static_assertion(!Has_member_id<StructNoId>::value);
  CGAL_static_assertion(Has_member_id<StructWithId>::value);
  CGAL_static_assertion(!Has_member_id<Polyhedron::Face>::value);
  CGAL_static_assertion(Has_member_id<Polyhedron_with_ids::Facet>::value);
  CGAL_static_assertion(Has_member_id<Polyhedron_with_ids::FBase>::value);
  CGAL_static_assertion(
    (Has_member_id<Polyhedron_with_ids::Items::Face_wrapper<Polyhedron_with_ids::HDS, K>::Face>::value));

  return 0;
}
