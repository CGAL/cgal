#ifndef MONOTONE_POLYGON_TRIANGULATION_TRAITS_2
#define MONOTONE_POLYGON_TRIANGULATION_TRAITS_2

namespace CGAL {

template <class Diagonal_, class Kernel_>
class Monotone_polygon_triangulation_traits_2
{
 private:
  typedef Kernel_ Kernel;

 public:
  typedef Diagonal_ Diagonal;
  typedef typename Kernel::Less_yx_2 Less_yx_2;
  typedef typename Kernel::Orientation_2 Orientation_2;

  Less_yx_2 less_yx_2_object() const {
    return Less_yx_2();
  }
  Orientation_2 orientation_2_object() const {
    return Orientation_2();
  }
};

}

#endif // MONOTONE_POLYGON_TRIANGULATION_TRAITS_2
