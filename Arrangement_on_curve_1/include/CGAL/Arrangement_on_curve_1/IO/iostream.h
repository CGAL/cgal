#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_IO_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_IO_H

#include <iostream>
#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

template <typename GeometryTraits_1, typename TopologyTraits>
std::ostream& operator<<(std::ostream& os,
                         const Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits>& arr)
{
  // Your output serialization logic here...
  return os;
}

template <typename GeometryTraits_1, typename TopologyTraits>
std::istream& operator>>(std::istream& is,
                         Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits>& arr)
{
  // Your input deserialization logic here...
  return is;
}

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
