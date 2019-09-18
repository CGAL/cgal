namespace CGAL {

/*!
\defgroup PkgInterpolationSurfaceNeighbors3 3D Surface Neighbors Functions
\ingroup PkgInterpolation2SurfaceNeighbor

Given a set of sample points issued from a surface and a query point
`p`, the functions `surface_neighbors_3()` compute the neighbors of `p` on
the surface within the sample points. If the sampling is sufficiently
dense, the neighbors are provably close to the point `p` on the
surface (cf. the manual pages and
\cgalCite{bf-lcss-02},\cgalCite{cgal:f-csapc-03}). They are defined to
be the neighbors of `p` in the regular triangulation dual
to the power diagram which is equivalent to the intersection of the
Voronoi cell of the query point `p` with the tangent plane to the
surface at `p`.

The functions \c surface_neighbors_certified_3() also return, in
addition, a Boolean value that certifies whether or not, the Voronoi
cell of `p` can be affected by points that lie outside the input
range, i.e.\ outside the ball centered on `p` passing through the
furthest sample point from `p` in the range `[first, beyond)`. If the sample
points are collected by a k-nearest neighbor or a range search
query, this permits to verify that a large enough neighborhood has
been considered.

\cgalHeading{Requirements}

<OL>
<LI>`Dt` is equivalent to the class
`Delaunay_triangulation_3`.
<LI>`OutputIterator::value_type` is equivalent to
`Dt::Point_3`, i.e.\ a point type.
<LI>`ITraits` is equivalent to the class `Voronoi_intersection_2_traits_3<K>`.
</OL>

\sa `CGAL::Voronoi_intersection_2_traits_3<K>`
\sa PkgInterpolationSurfaceNeighborCoordinates3

\cgalHeading{Implementation}

These functions compute the regular triangulation of
the sample points and the point `p` using a traits class
equivalent to `Voronoi_intersection_2_traits_3<K>`. They determine
the neighbors of `p` in this triangulation. The functions which
certify the result need to compute, in addition, the Voronoi vertices
of the cell of `p` in this diagram.
*/
/// @{

/*!
The sample points \f$ \mathcal{P}\f$ are provided in the range
`[first, beyond)`.
`InputIterator::value_type` is the point type `Kernel::Point_3`. The
tangent plane is defined by the point `p` and the vector `normal`. The
parameter `K` determines the kernel type that will instantiate the
template parameter of `Voronoi_intersection_2_traits_3<K>`.

The surface neighbors of `p` are computed which are the
neighbors of `p` in the regular triangulation that is dual to
the intersection of the 3D Voronoi diagram of \f$ \mathcal{P}\f$ with
the tangent plane. The point sequence that is computed by the
function is placed starting at `out`. The function returns an
iterator that is placed past-the-end of the resulting point
sequence.
*/
template <class OutputIterator, class InputIterator, class Kernel>
OutputIterator
surface_neighbors_3(InputIterator first, InputIterator beyond,
                    const typename Kernel::Point_3& p,
                    const typename Kernel::Vector_3& normal,
                    OutputIterator out,
                    const Kernel& K);

/*!
Same as above only that the traits class must be instantiated by
the user. `ITraits` must be equivalent to
`Voronoi_intersection_2_traits_3<K>`.
*/
template <class OutputIterator, class InputIterator, class ITraits>
OutputIterator surface_neighbors_3(InputIterator first, InputIterator beyond,
                                   const typename ITraits::Point_2& p,
                                   OutputIterator out,
                                   const ITraits& traits);

/*!
Similar to the first function. The additional third return
value is `true` if the furthest point in the range
`[first, beyond)` is further
away from `p` than twice the distance from `p` to the
furthest vertex of the intersection of the Voronoi cell of `p`
with the tangent plane defined be `(p,normal)`. It is
`false` otherwise.
*/
template <class OutputIterator, class InputIterator, class Kernel>
std::pair< OutputIterator, bool >
surface_neighbors_certified_3(InputIterator first, InputIterator beyond,
                              const typename Kernel::Point_3& p,
                              const typename Kernel::Vector_3& normal,
                              OutputIterator out,
                              const Kernel& K);

/*!
Same as above except that this function
takes the maximal distance from `p` to the points in the range
`[first, beyond)` as additional parameter.
*/
template <class OutputIterator, class InputIterator, class Kernel>
std::pair< OutputIterator, bool >
surface_neighbors_certified_3(InputIterator first, InputIterator beyond,
                              const typename Kernel::Point_3& p,
                              const typename Kernel::Vector_3& normal,
                              const typename Kernel::FT& max_distance,
                              OutputIterator out,
                              const Kernel& kernel);

/*!
Same as above only that the traits class must be instantiated by the user.
`ITraits` must be equivalent to `Voronoi_intersection_2_traits_3<K>`. There is no
parameter `max_distance`.
*/
template <class OutputIterator, class InputIterator, class ITraits>
std::pair< OutputIterator, bool >
surface_neighbors_certified_3(InputIterator first, InputIterator beyond,
                              const typename ITraits::Point_2& p,
                              OutputIterator out,
                              const ITraits& traits);

/*!
Same as above with the parameter `max_distance`.
*/
template <class OutputIterator, class InputIterator, class ITraits>
std::pair< OutputIterator, bool >
surface_neighbors_certified_3(InputIterator first, InputIterator beyond,
                              const typename ITraits::Point_2& p,
                              const typename ITraits::FT& max_distance,
                              OutputIterator out,
                              const ITraits& traits);

/*!
Computes the surface neighbor coordinates with respect to the points
that are vertices of the Delaunay triangulation `dt`. The type `Dt`
must be equivalent to `Delaunay_triangulation_3<Gt, Tds>`. The
optional parameter `start` is used for the used as a starting place
for the search of the conflict zone. It may be the result of the call
`dt.locate(p)`. This function instantiates the template parameter
`ITraits` to be `Voronoi_intersection_2_traits_3<Dt::Geom_traits>`.

This function allows to filter some potential neighbors of the
query point `p` from \f$ \mathcal{P}\f$ via its three-dimensional
Delaunay triangulation. All surface neighbors of `p` are
necessarily neighbors in the Delaunay triangulation of \f$ \mathcal{P}
\cup \{p\}\f$.
*/
template < class Dt, class OutputIterator >
OutputIterator
surface_neighbors_3(const Dt& dt,
                    const typename Dt::Geom_traits::Point_3& p,
                    const typename Dt::Geom_traits::Vector_3& normal,
                    OutputIterator out,
                    typename Dt::Cell_handle start = typename Dt::Cell_handle());

/*!
Same as above only that the parameter `traits` instantiates
the geometric traits class. Its type `ITraits` must be
equivalent to `Voronoi_intersection_2_traits_3<K>`.
*/
template < class Dt, class OutputIterator, class ITraits>
OutputIterator surface_neighbors_3(const Dt& dt,
                                   const typename ITraits::Point_2& p,
                                   OutputIterator out,
                                   const ITraits& traits,
                                   typename Dt::Cell_handle start = typename Dt::Cell_handle());

/// @}

} /* namespace CGAL */

