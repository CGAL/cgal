namespace CGAL {

/*!
\defgroup PkgInterpolationSurfaceNeighborCoordinates3 3D Surface Neighbor Coordinates Functions
\ingroup PkgInterpolation2SurfaceNeighbor

The functions `surface_neighbor_coordinates_3()` compute natural neighbor coordinates for 
surface points associated to a finite set of sample points issued from 
the surface. The coordinates are computed from the intersection of the 
Voronoi cell of the query point `p` with the tangent plane to the 
surface at `p`. If the sampling is sufficiently dense, the 
coordinate system meets the properties described in the manual pages 
and in \cgalCite{bf-lcss-02},\cgalCite{cgal:f-csapc-03}. The query 
point `p` needs to lie inside the convex hull of the projection of 
the sample points onto the tangent plane at `p`. 

The functions `surface_neighbor_coordinates_certified_3()` return, in
addition, a second Boolean value (the fourth value of the quadruple)
that certifies whether or not, the Voronoi cell of `p` can be affected
by points that lie outside the input range, i.e.\ outside the ball
centered on `p` passing through the furthest sample point from `p` in
the range `[first, beyond)`. If the sample points are collected by a
`k`-nearest neighbor or a range search query, this permits to check
whether the neighborhood which has been considered is large enough.

\cgalHeading{Requirements}

<OL> 
<LI>`Dt` is equivalent to the class 
`Delaunay_triangulation_3`. 
<LI>The value type of `OutputIterator` is equivalent to 
`std::pair<Dt::Point_3, Dt::Geom_traits::FT>`, i.e.\ a pair 
associating a point and its natural neighbor coordinate. 
<LI>`ITraits` is equivalent to the class `Voronoi_intersection_2_traits_3<K>`. 
</OL> 

\sa `CGAL::linear_interpolation()`
\sa `CGAL::sibson_c1_interpolation()` 
\sa `CGAL::farin_c1_interpolation()`
\sa `CGAL::Voronoi_intersection_2_traits_3<K>`
\sa PkgInterpolationSurfaceNeighbors3

\cgalHeading{Implementation}

This functions construct the regular triangulation of the input points 
instantiated with `Voronoi_intersection_2_traits_3<Kernel>` or `ITraits` if provided. 
They return the result of the function call 
`PkgInterpolationRegularNeighborCoordinates2`
with the regular triangulation and `p` as arguments. 

*/
/// @{

/*!
The sample points \f$ \mathcal{P}\f$ are provided in the range
`[first`, beyond)`.
The value type of `InputIterator` is the point type
`Kernel::Point_3`. The tangent plane is defined by the point
`p` and the vector `normal`. The parameter `K`
determines the kernel type that will instantiate
the template parameter of `Voronoi_intersection_2_traits_3<K>`. 

The natural neighbor coordinates for `p` are computed in the
power diagram that results from the intersection of the `3D` Voronoi
diagram of \f$ \mathcal{P}\f$ with the tangent plane. The sequence of
point/coordinate pairs that is computed by the function is placed
starting at `out`. The function returns a triple with an
iterator that is placed past-the-end of the resulting sequence of
point/coordinate pairs, the normalization factor of the coordinates
and a Boolean value which is set to true iff the coordinate
computation was successful, i.e.\ if `p` lies inside the convex
hull of the projection of the points \f$ \mathcal{P}\f$ onto the tangent
plane.
*/
template <class OutputIterator, class InputIterator, class
Kernel> CGAL::Triple< OutputIterator, typename Kernel::FT, bool >
surface_neighbor_coordinates_3(InputIterator first, InputIterator
beyond, const typename Kernel::Point_3& p, const typename
Kernel::Vector_3& normal, OutputIterator out, const Kernel& K);

/*!
the same as above only that the traits class
must be instantiated by the user. `ITraits` must be equivalent
to `Voronoi_intersection_2_traits_3<K>`.
*/
template <class OutputIterator, class InputIterator, class
ITraits> CGAL::Triple< OutputIterator, typename ITraits::FT, bool >
surface_neighbor_coordinates_3(InputIterator first, InputIterator
beyond, const typename ITraits::Point_2& p,OutputIterator out, const
ITraits& traits);

/*!
Similar to the first function. The additional fourth return
value is `true` if the furthest point in the range
`[first, beyond)` is further
away from `p` than twice the distance from `p` to the
furthest vertex of the intersection of the Voronoi cell of `p`
with the tangent plane defined by `(p,normal)`. It is
`false` otherwise.
*/
template <class OutputIterator, class InputIterator, class
Kernel> CGAL::Quadruple< OutputIterator, typename Kernel::FT, bool,
bool > surface_neighbor_coordinates_certified_3(InputIterator first,
InputIterator beyond, const typename Kernel::Point_3& p, const
typename Kernel::Vector_3& normal, OutputIterator out, const Kernel&
K);

/*!
The same as above except that this function takes the
maximal distance from p to the points in the range
`[first, beyond)` as additional parameter.
*/
template <class OutputIterator,
class InputIterator, class Kernel> CGAL::Quadruple< OutputIterator,
typename Kernel::FT, bool, bool >
surface_neighbor_coordinates_certified_3(InputIterator first,
InputIterator beyond, const typename Kernel::Point_3& p, const
typename Kernel::FT& max_distance, OutputIterator out, const Kernel&
kernel);

/*!
The same as above only
that the traits class must be instantiated by the user and without
the parameter `max_distance`. `ITraits` must be equivalent
to `Voronoi_intersection_2_traits_3<K>`.
*/
template <class OutputIterator, class InputIterator, class
ITraits> CGAL::Quadruple< OutputIterator, typename ITraits::FT, bool, bool >
surface_neighbor_coordinates_certified_3(InputIterator first,
InputIterator beyond, const typename ITraits::Point_2& p,
OutputIterator out, const ITraits& traits);

/*!
The same as above with the parameter
`max_distance`.
*/
template <class OutputIterator, class InputIterator, class
ITraits> CGAL::Quadruple< OutputIterator, typename ITraits::FT, bool, bool >
surface_neighbor_coordinates_certified_3(InputIterator first,
InputIterator beyond, const typename ITraits::Point_2& p, const
typename ITraits::FT& max_distance, OutputIterator out, const
ITraits& traits);

/*!
computes the surface neighbor coordinates with respect to the points
that are vertices of the Delaunay triangulation `dt`. The type `Dt`
must be equivalent to `Delaunay_triangulation_3<Gt, Tds>`. The
optional parameter `start` is used as a starting place for the search
of the conflict zone. It may be the result of the call
`dt.locate(p)`. This function instantiates the template parameter
`ITraits` to be `Voronoi_intersection_2_traits_3<Dt::Geom_traits>`.

This function allows to filter some potential neighbors of the 
query point `p` from \f$ \mathcal{P}\f$ via its three-dimensional 
Delaunay triangulation. All surface neighbors of `p` are 
necessarily neighbors in the Delaunay triangulation of \f$ \mathcal{P} 
\cup \{p\}\f$. 
*/
template < class Dt, class OutputIterator >
CGAL::Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
surface_neighbor_coordinates_3(const Dt& dt, const typename
Dt::Geom_traits::Point_3& p, const typename
Dt::Geom_traits::Vector_3& normal, OutputIterator out, typename
Dt::Cell_handle start = typename Dt::Cell_handle());

/*!
The same as above only that the parameter traits instantiates
the geometric traits class. Its type `ITraits` must be
equivalent to `Voronoi_intersection_2_traits_3<K>`.
*/
template < class Dt, class OutputIterator,
class ITraits> CGAL::Triple< OutputIterator, typename
Dt::Geom_traits::FT, bool > surface_neighbor_coordinates_3(const Dt& dt,
const typename Dt::Geom_traits::Point_3& p, OutputIterator out,
const ITraits& traits, typename Dt::Cell_handle start = typename
Dt::Cell_handle());

/// @}

} /* namespace CGAL */

