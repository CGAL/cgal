#ifndef CGAL_AOC_TRAITS_1_H
#define CGAL_AOC_TRAITS_1_H

#include <CGAL/enum.h>

namespace CGAL {

/// \ingroup PkgArrangementOnCurve1Concepts
/// \cgalConcept
class AocTraits_1 {
public:
  /// A point guaranteed to lie on the master curve
  typedef unspecified_type Point_1;

  /// Functor evaluating relative positions along the curve parameterization
  class Compare_x_1 {
  public:
    Comparison_result operator()(const Point_1& p1, const Point_1& p2) const;
  };

  Compare_position_1 compare_x_1_object() const;
};

} // namespace CGAL

#endif
