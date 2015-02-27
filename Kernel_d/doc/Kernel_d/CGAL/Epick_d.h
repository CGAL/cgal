
namespace CGAL {

/*!
\ingroup PkgKernelDKernels

A model for `Kernel_d` that uses %Cartesian coordinates to represent the
geometric objects. 

This kernel is default constructible and copyable. It does not carry any
state so it is possible to use objects created by one instance with
functors created by another one.

This kernel supports construction of points from `double`
%Cartesian coordinates. It provides exact geometric predicates, but
the geometric constructions are not guaranteed to be exact. The geometric
predicates are made exact without sacrificing speed thanks to the use of
filters.

\tparam DimensionTag is a tag representing the dimension of the
ambient Euclidean space. It may be either `Dimension_tag<d>` where `d` is
an integer
or `Dynamic_dimension_tag`. In the latter case, the dimension of the space is specified for each point when it is constructed, so it does not need to be known at compile-time.


\attention Only the interfaces specific to this class are listed below. Refer to the
concepts for the rest.

\attention Known bugs: the functor `Intersect_d` is not yet implemented. `Contained_in_affine_hull` assumes that the iterators refer to an affinely independent family. `Orientation_d` only works for points, not vectors.

\attention Ancient compilers like gcc-4.2 or icc 14 are not supported, but gcc-4.4 and
icc 15 work.

\attention This kernel requires the \ref thirdpartyEigen "Eigen" library.



\cgalModels `Kernel_d`
\cgalModels `DelaunayTriangulationTraits`

\sa `CGAL::Cartesian_d<FieldNumberType>`
\sa `CGAL::Homogeneous_d<RingNumberType>`

*/
template< typename DimensionTag >
class Epick_d {
public:
/*!
represents a point in the Euclidean space
\cgalModels `DefaultConstructible`
\cgalModels `Assignable`
*/
class Point_d {
public:
/*! introduces a point with coordinates (x0, x1, ...) where the number of
    coordinates matches the dimension.
    \pre `DimensionTag` is a fixed dimension, not `Dynamic_dimension_tag`. */
Point_d(double x0, double x1, ...);

/*! introduces a point with coordinate set `[first,end)`.
    \pre If `DimensionTag` is a fixed dimension, it matches `distance(first,end)`.
    \cgalRequires The value type of `InputIterator` is convertible to `double`.
    */
template<typename InputIterator>
Point_d(InputIterator first, InputIterator end);

/*! returns the i'th coordinate of a point.
    \pre `i` is non-negative and less than the dimension. */
double operator[](int i)const;

/*! returns an iterator pointing to the zeroth Cartesian coordinate. */
Cartesian_const_iterator_d cartesian_begin()const;
/*! returns an iterator pointing beyond the last Cartesian coordinate. */
Cartesian_const_iterator_d cartesian_end()const;
};
/// @}

}; /* end Epick_d */
} /* end namespace CGAL */
