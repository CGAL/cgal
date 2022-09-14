#define CGAL_NO_MPZF_DIVISION_OPERATOR 1

#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>

using SCK = CGAL::Simple_cartesian<double>;

struct Any {
  template <class T> operator const T&();
};

template <typename Predicate>
void check_pred() {
  std::cerr << std::is_invocable_v<Predicate, Any>;
  std::cerr << std::is_invocable_v<Predicate, Any, Any>;
  std::cerr << std::is_invocable_v<Predicate, Any, Any, Any>;
  std::cerr << std::is_invocable_v<Predicate, Any, Any, Any, Any>;
  std::cerr << std::is_invocable_v<Predicate, Any, Any, Any, Any, Any>;
  std::cerr << std::is_invocable_v<Predicate, Any, Any, Any, Any, Any, Any>;
  std::cerr << std::is_invocable_v<Predicate, Any, Any, Any, Any, Any, Any, Any>;
  std::cerr << std::is_invocable_v<Predicate, Any, Any, Any, Any, Any, Any, Any, Any>;

  // The following asserts that no predicate from the kernel has more than
  // 8 arguments (actually the assertions are only from 9 to 12 arguments).
  static_assert(!std::is_invocable_v<Predicate, Any, Any, Any, Any, Any, Any, Any, Any, Any>);
  static_assert(!std::is_invocable_v<Predicate, Any, Any, Any, Any, Any, Any, Any, Any, Any, Any>);
  static_assert(!std::is_invocable_v<Predicate, Any, Any, Any, Any, Any, Any, Any, Any, Any, Any, Any>);
  static_assert(!std::is_invocable_v<Predicate, Any, Any, Any, Any, Any, Any, Any, Any, Any, Any, Any, Any>);
  std::cerr << '\n';
}

int main()
{
#define CGAL_Kernel_pred(P, Pf)  \
  std::cerr << #P << ": "; \
  check_pred<SCK::P>();
#include <CGAL/Kernel/interface_macros.h>

  // Bug with predicates with multiple overload of the call operator with the
  // same number of arguments: the call with `Any` is ambiguous.
  static_assert(std::is_invocable_v<SCK::Angle_3, Any, Any, Any>);
  static_assert(!std::is_invocable_v<SCK::Angle_3, Any, Any, Any, Any>); // AMBIGUOUS CALL
  static_assert(!std::is_invocable_v<SCK::Is_degenerate_2, Any>); // AMBIGUOUS CALL
  return 0;
}


/*

WORK IN PROGRESS:

In the CGAL Kernel:
  - 2D: 49 predicates
  - 3D: 50 predicates


## Try to detect all possible types of arguments of predicates, from the doc

```
[lrineau@fernand]~/Git/cgal-master/build-doc/doc_output/Kernel_23/xml% grep -h '         <type>' classKernel_1_1(^(Construct|Compute|Assign)*).xml | sed 's/<ref[^>]*>//; s|</ref>||; s| &amp;|\&|' | sed 's/Kernel::/K::/'| sort | uniq -c | sort -n  | grep -v _3
```

3D: (14 types of arguments)

const K::Direction_3&
const K::Triangle_3&
const K::Circle_3&
const K::Ray_3&
const K::Segment_3&
const K::Iso_cuboid_3&
const K::Line_3&
const K::Tetrahedron_3&
const K::FT&
const K::Plane_3&
const K::Sphere_3&
const K::Vector_3&
const K::Weighted_point_3&
const K::Point_3&

2D: (10 types arguments)

const K::Vector_2&
const K::Direction_2&
const K::Iso_rectangle_2&
const K::Ray_2&
const K::Circle_2&
const K::Triangle_2&
const K::FT&
const K::Segment_2&
const K::Weighted_point_2&
const K::Line_2&
const K::Point_2&

*/
