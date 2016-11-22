#ifndef POLYGON_TRIANGULATION_IS_VALID_TRAITS_2_H
#define POLYGON_TRIANGULATION_IS_VALID_TRAITS_2_H

namespace CGAL {

template <class Kernel_>
class Polygon_triangulation_is_valid_traits_2
{
  typedef  Kernel_ Kernel;
 public:
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Compare_xy_2 Compare_xy_2;
};

}
#endif // POLYGON_TRIANGULATION_IS_VALID_2_H
