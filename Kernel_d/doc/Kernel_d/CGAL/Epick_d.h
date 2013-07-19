
namespace CGAL {

/*!
\ingroup PkgKernelDKernels

A model for `Kernel_d` that uses %Cartesian coordinates to represent the
geometric objects. The parameter `Dimension` is the dimension of the
ambient Euclidean space. It may be either `Dimension_tag<d>` or
`Dynamic_dimension_tag`. It supports construction of points from `double`
Cartesian coordinates. It provides exact geometric predicates, but
inexact geometric constructions. The geometric predicates are made exact
without sacrificing speed thanks to the use of filters.

This kernel is default constructible and copyable. It does not carry any
state so it is possible to use objects created by one instance with
functors created by another one.

Note that this kernel does not completely conform to the `Kernel_d`
concept: it is missing the constructions `Lift_to_paraboloid_d` and
`Project_along_d_axis_d` which do not make sense with a single fixed
dimension.


\cgalModels `Kernel_d`
\cgalModels `DelaunayTriangulationTraits`

\sa `CGAL::Epick_d<Dimension>::Point_d`
\sa `CGAL::Cartesian_d<FieldNumberType>`
\sa `CGAL::Homogeneous_d<RingNumberType>`

*/
template< typename Dimension >
class Epick_d {
public:
/*!
represents a point in the Euclidean space
\cgalModels `DefaultConstructible`
\cgalModels `Copyable`
\cgalModels `Assignable`
*/
class Point_d {
public:
/*! introduces a point with coordinates (x0, x1, ...). */
Point_d(double x0, double x1, ...);

/*! returns the i'th coordinate of a point.r
    \pre `i` is non-negative and less than the dimension */
double operator[](int i)const;
}
/// @}

}; /* end Epick_d */
} /* end namespace CGAL */
