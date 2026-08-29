#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_IO_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_IO_H

#include <iostream>

#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1op_left_shift
 *
 * Inserts the arrangement object `arr` into the output stream `os`. Only the
 * basic geometric and topological features of the arrangement are
 * inserted. Auxiliary data that may be attached to the vertices and edges is
 * ignored.
 */
  template <typename GeometryTraits_1, typename TopologyTraits, bool BinarySearch>
std::ostream& operator<<(std::ostream& os,
                         const Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>& arr);

/*! \ingroup PkgArrangementOnCurve1op_right_shift
 *
 * Extracts an arrangement from a given input stream `is`. Only the basic
 * geometric and topological features of the arrangement are read and no
 * auxiliary data is attached to the vertices and edges.
 */
  template <typename GeometryTraits_1, typename TopologyTraits, bool BinarySearch>
std::istream& operator>>(std::istream& is,
                         Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>& arr);

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
