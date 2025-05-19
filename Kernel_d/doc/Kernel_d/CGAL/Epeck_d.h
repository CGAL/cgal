
namespace CGAL {

/*!
\ingroup PkgKernelDKernels

A model for `Kernel_d`, minus `Kernel_d::Point_of_sphere_d`, that uses %Cartesian coordinates to represent the
geometric objects.

This kernel is default constructible and copyable. It does not carry any
state so it is possible to use objects created by one instance with
functors created by another one.

This kernel supports construction of points from `double` %Cartesian
coordinates. It provides exact geometric predicates and constructions. The
geometric predicates are made exact without sacrificing speed thanks to the use
of filters. The geometric constructions are made exact without sacrificing
speed thanks to a lazy mechanism, similar to
`Exact_predicates_exact_constructions_kernel`. A construction creates an
approximate object, and stores a directed acyclic graph (DAG) of the operation
and arguments used. When an operation needs more precision on an object than is
currently available, which should be rare, \cgal reconstructs exactly all the
ancestors of the object and replaces this part of the graph with exact objects.
This should be transparent for users, those details do not affect the
functionality, but they can cause surprising running time where the costly part
of an algorithm is not the construction itself, but a seemingly trivial use
afterwards that causes exact reconstruction of a large part of the structure.

`Sphere_d` is represented internally with a center and a squared radius, and the exact type used by the kernel is a rational type. This means that `Kernel_d::Point_of_sphere_d` cannot be implemented exactly in dimensions up to 3 (higher dimensions would require a decomposition of an integer as a sum of `d` squares), so currently it is not provided at all.

\tparam DimensionTag is a tag representing the dimension of the
ambient Euclidean space. It may be either `Dimension_tag<d>` where `d` is
an integer
or `Dynamic_dimension_tag`. In the latter case, the dimension of the space is specified for each point when it is constructed, so it does not need to be known at compile-time.


\attention Only the interfaces specific to this class are listed below. Refer to the
concepts for the rest.

\attention Known bugs: the functor `Kernel_d::Intersect_d` is not yet implemented. `Kernel_d::Contained_in_affine_hull` assumes that the iterators refer to an affinely independent family. `Kernel_d::Orientation_d` only works for points, not vectors. Visual Studio 2015 is not supported.

\attention This kernel requires the \ref thirdpartyEigen "Eigen" library.



\cgalModels{Kernel_d,DelaunayTriangulationTraits,RegularTriangulationTraits,SearchTraits,RangeSearchTraits}

\sa `CGAL::Cartesian_d<FieldNumberType>`
\sa `CGAL::Homogeneous_d<RingNumberType>`
\sa `CGAL::Epick_d<DimensionTag>`

*/
template< typename DimensionTag >
struct Epeck_d {
/*!
represents a point in the Euclidean space
\cgalModels{DefaultConstructible,Assignable}
*/
class Point_d {
public:
/*! introduces a point with coordinates (x0, x1, ...) where the number of
    coordinates matches the dimension.
    \pre `DimensionTag` is a fixed dimension, not `Dynamic_dimension_tag`. */
Point_d(double x0, double x1, ...);

/*! introduces a point with coordinate set `[first,end)`.
    \pre If `DimensionTag` is a fixed dimension, it matches `distance(first,end)`.
    \tparam ForwardIterator has its value type that is convertible to `double`.
    */
template<typename ForwardIterator>
Point_d(ForwardIterator first, ForwardIterator end);

/*! returns the i-th coordinate of a point.
    \pre `i` is non-negative and less than the dimension. */
double operator[](int i)const;

/*! returns an iterator pointing to the zeroth %Cartesian coordinate. */
Cartesian_const_iterator_d cartesian_begin()const;
/*! returns an iterator pointing beyond the last %Cartesian coordinate. */
Cartesian_const_iterator_d cartesian_end()const;
};

/*!
represents a weighted point in the Euclidean space
\cgalModels{DefaultConstructible,Assignable}
*/
class Weighted_point_d {
public:
/*! introduces a weighted point with point p and weight w. */
Weighted_point_d(const Point_d& p, const double& w);
/*! extracts the point of a weighted point. */
Point_d point() const;
/*! extracts the weight of a weighted point. */
double weight() const;
};

/*! \cgalModels{Kernel_d::Center_of_sphere_d}
 */
class Construct_circumcenter_d {
public:
/*! returns the center of the sphere defined by `A=tuple[first,last)`. The sphere is centered in the affine hull of A and passes through all the points of A. The order of the points of A does not matter.
    \pre A is affinely independent.
    \tparam ForwardIterator has `Epeck_d::Point_d` as value type.
    */
template<typename ForwardIterator>
Point_d operator()(ForwardIterator first, ForwardIterator last);
};
class Compute_squared_radius_d {
public:
/*! returns the radius of the sphere defined by `A=tuple[first,last)`. The sphere is centered in the affine hull of A and passes through all the points of A. The order of the points of A does not matter.
    \pre A is affinely independent.
    \tparam ForwardIterator has `Epeck_d::Point_d` as value type.
    */
template<class ForwardIterator>
FT operator()(ForwardIterator first, ForwardIterator last);
};
class Compute_squared_radius_smallest_orthogonal_sphere_d {
public:
/*! returns the radius of the sphere defined by `A=tuple[first,last)`. The sphere is centered in the affine hull of A and orthogonal to all the spheres of A. The order of the points of A does not matter.
    \pre A is affinely independent.
    \tparam ForwardIterator has `Epeck_d::Weighted_point_d` as value type.
    */
template<class ForwardIterator>
FT operator()(ForwardIterator first, ForwardIterator last);
};
/*! \cgalModels{Kernel_d::Side_of_bounded_sphere_d}
 */
class Side_of_bounded_sphere_d {
public:
/*! returns the relative position of point p to the sphere defined by `A=tuple[first,last)`. The sphere is centered in the affine hull of A and passes through all the points of A. The order of the points of A does not matter.
    \pre A is affinely independent.
    \tparam ForwardIterator has `Epeck_d::Point_d` as value type.
    */
template<class ForwardIterator>
Bounded_side operator()(ForwardIterator first, ForwardIterator last, const Point_d&p);
};
class Power_side_of_bounded_power_sphere_d {
public:
/*! returns the relative position of weighted point p to the sphere defined by `A=tuple[first,last)`. The sphere is centered in the affine hull of A and orthogonal to all the spheres of A. The order of the points of A does not matter.
    \pre A is affinely independent.
    \tparam ForwardIterator has `Epeck_d::Weighted_point_d` as value type.
    */
template<class ForwardIterator>
Bounded_side operator()(ForwardIterator first, ForwardIterator last, const Weighted_point_d&p);
};
class Construct_power_sphere_d {
public:
/*! returns a weighted point on the affine hull of the weighted points of `A=tuple[first,last)` at power distance 0 of each of them. In other words, this returns the smallest sphere orthogonal to the spheres of A.
    \pre A is affinely independent.
    \tparam ForwardIterator has `Epeck_d::Weighted_point_d` as value type.
    */
template<class ForwardIterator>
Weighted_point_d operator()(ForwardIterator first, ForwardIterator last);
};
class Compute_power_product_d {
public:
/*! returns the power product (aka power distance) of the weighted points `pw` and `qw`, that is, the squared Euclidean distance between the points minus their weights.
 */
FT operator()(Weighted_point_d pw, Weighted_point_d qw);
};
Construct_circumcenter_d construct_circumcenter_d_object();
Compute_power_product_d compute_power_product_d_object();
Compute_squared_radius_d compute_squared_radius_d_object();
Compute_squared_radius_smallest_orthogonal_sphere_d compute_squared_radius_smallest_orthogonal_sphere_d_object();
Construct_power_sphere_d construct_power_sphere_d_object();
Power_side_of_bounded_power_sphere_d power_side_of_bounded_power_sphere_d_object();
}; /* end Epeck_d */
} /* end namespace CGAL */
