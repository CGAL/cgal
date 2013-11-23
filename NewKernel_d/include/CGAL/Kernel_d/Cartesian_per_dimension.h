#ifndef CGAL_KD_CARTESIAN_PER_DIM_H
#define CGAL_KD_CARTESIAN_PER_DIM_H
#include <CGAL/functor_tags.h>
#include <CGAL/Dimension.h>
#include <CGAL/predicates/sign_of_determinant.h>

// Should probably disappear.

namespace CGAL {
template <class Dim_, class R_, class Derived_>
struct Cartesian_per_dimension : public R_ {};
}

#endif
