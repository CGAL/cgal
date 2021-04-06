#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>

template <typename K>
typename K::Point_2 point_2 (int x, int y)
{ return typename K::Point_2 (typename K::FT(x), typename K::FT(y)); }
template <typename K>
typename K::Point_3 point_3 (int x, int y, int z)
{ return typename K::Point_3 (typename K::FT(x), typename K::FT(y), typename K::FT(z)); }

// Generate default objects so that the test compiles *AND* runs fine
template <typename K, typename O> O default_object(const O&) { return O(); }
template <typename K> typename K::Point_2 default_object(const typename K::Point_2&)
{ return point_2<K>(0,0); }
template <typename K> typename K::Segment_2 default_object(const typename K::Segment_2&)
{ return typename K::Segment_2(point_2<K>(0,0), point_2<K>(0,1)); }
template <typename K> typename K::Circle_2 default_object(const typename K::Circle_2&)
{ return typename K::Circle_2(point_2<K>(0,0), typename K::FT(1.0)); }
template <typename K> typename K::Triangle_2 default_object(const typename K::Triangle_2&)
{ return typename K::Triangle_2(point_2<K>(0,0), point_2<K>(0,1), point_2<K>(1,0)); }
template <typename K> typename K::Iso_rectangle_2 default_object(const typename K::Iso_rectangle_2&)
{ return typename K::Iso_rectangle_2(point_2<K>(0,0), point_2<K>(1,1)); }
template <typename K> typename K::Point_3 default_object(const typename K::Point_3&)
{ return point_3<K>(0,0,0); }
template <typename K> typename K::Segment_3 default_object(const typename K::Segment_3&)
{ return typename K::Segment_3(point_3<K>(0,0,0), point_3<K>(0,0,1)); }
template <typename K> typename K::Sphere_3 default_object(const typename K::Sphere_3&)
{ return typename K::Sphere_3(point_3<K>(0,0,0), typename K::FT(1.0)); }
template <typename K> typename K::Triangle_3 default_object(const typename K::Triangle_3&)
{ return typename K::Triangle_3(point_3<K>(0,0,0), point_3<K>(0,0,1), point_3<K>(0,1,0)); }
template <typename K> typename K::Tetrahedron_3 default_object(const typename K::Tetrahedron_3&)
{ return typename K::Tetrahedron_3(point_3<K>(0,0,0), point_3<K>(0,0,1),
                                   point_3<K>(0,1,0), point_3<K>(1,0,0)); }
template <typename K> typename K::Iso_cuboid_3 default_object(const typename K::Iso_cuboid_3&)
{ return typename K::Iso_cuboid_3(point_3<K>(0,0,0), point_3<K>(1,1,1)); }

template <typename Kernel, typename Object, int dim>
void test_2d()
{
  std::array<Object, 1> dummy = { default_object<Kernel>(Object()) };
  typename Kernel::Line_2 line;
  typename Kernel::Point_2 centroid;
  CGAL::linear_least_squares_fitting_2 (dummy.begin(), dummy.end(), line, centroid,
                                        CGAL::Dimension_tag<dim>(), Kernel(),
                                        CGAL::Eigen_diagonalize_traits<typename Kernel::FT, 2>());
}

template <typename Kernel, typename Object, int dim>
void test_3d()
{
  std::array<Object, 1> dummy = { default_object<Kernel>(Object()) };
  typename Kernel::Line_3 line;
  typename Kernel::Plane_3 plane;
  typename Kernel::Point_3 centroid;
  CGAL::linear_least_squares_fitting_3 (dummy.begin(), dummy.end(), line, centroid,
                                        CGAL::Dimension_tag<dim>(), Kernel(),
                                        CGAL::Eigen_diagonalize_traits<typename Kernel::FT, 3>());
  CGAL::linear_least_squares_fitting_3 (dummy.begin(), dummy.end(), plane, centroid,
                                        CGAL::Dimension_tag<dim>(), Kernel(),
                                        CGAL::Eigen_diagonalize_traits<typename Kernel::FT, 3>());
}

template <typename Kernel>
void test_kernel()
{
  test_2d<Kernel, typename Kernel::Point_2, 0>();
  test_2d<Kernel, typename Kernel::Segment_2, 1>();
  test_2d<Kernel, typename Kernel::Segment_2, 0>();
  test_2d<Kernel, typename Kernel::Circle_2, 2>();
  test_2d<Kernel, typename Kernel::Circle_2, 1>();
  test_2d<Kernel, typename Kernel::Triangle_2, 2>();
  test_2d<Kernel, typename Kernel::Triangle_2, 1>();
  test_2d<Kernel, typename Kernel::Triangle_2, 0>();
  test_2d<Kernel, typename Kernel::Iso_rectangle_2, 2>();
  test_2d<Kernel, typename Kernel::Iso_rectangle_2, 1>();
  test_2d<Kernel, typename Kernel::Iso_rectangle_2, 0>();
  test_3d<Kernel, typename Kernel::Point_3, 0>();
  test_3d<Kernel, typename Kernel::Segment_3, 1>();
  test_3d<Kernel, typename Kernel::Segment_3, 0>();
  test_3d<Kernel, typename Kernel::Sphere_3, 3>();
  test_3d<Kernel, typename Kernel::Sphere_3, 2>();
  test_3d<Kernel, typename Kernel::Triangle_3, 2>();
  test_3d<Kernel, typename Kernel::Triangle_3, 1>();
  test_3d<Kernel, typename Kernel::Triangle_3, 0>();
  test_3d<Kernel, typename Kernel::Tetrahedron_3, 3>();
  test_3d<Kernel, typename Kernel::Tetrahedron_3, 2>();
  test_3d<Kernel, typename Kernel::Tetrahedron_3, 1>();
  test_3d<Kernel, typename Kernel::Tetrahedron_3, 0>();
  test_3d<Kernel, typename Kernel::Iso_cuboid_3, 3>();
  test_3d<Kernel, typename Kernel::Iso_cuboid_3, 2>();
  test_3d<Kernel, typename Kernel::Iso_cuboid_3, 1>();
  test_3d<Kernel, typename Kernel::Iso_cuboid_3, 0>();
}

int main()
{
  test_kernel<CGAL::Simple_cartesian<float> >();
  test_kernel<CGAL::Simple_cartesian<double> >();
  test_kernel<CGAL::Exact_predicates_inexact_constructions_kernel>();
  test_kernel<CGAL::Exact_predicates_exact_constructions_kernel>();

  return EXIT_SUCCESS;
}
