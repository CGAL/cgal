namespace CGAL {

/*!
\ingroup PkgPolygonPartitioning2

\brief determines if the polygons in the range [`poly_first`, `poly_beyond`)
define a valid convex partition of the polygon defined by the points in the 
range [`point_first`, `point_beyond`). 
A convex partition is valid if the 
polygons do not overlap, the union of the polygons is the same as the original 
polygon given by the sequence of points, and if each partition polygon is 
convex. 
The function returns `true` iff the partition is valid and otherwise
returns `false`.

\pre The points in the range [`point_first`, `point_beyond`) define a simple, counterclockwise-oriented polygon.




\cgalHeading{Requires}

- `Traits` is a model of the concept 
  `ConvexPartitionIsValidTraits_2`. 
- `std::iterator_traits<InputIterator>::%value_type` should be `Traits::Point_2`, 
  which should also be the type of the points stored in an object 
  of type `Traits::Polygon_2`. 
- `std::iterator_traits<ForwardIterator>::%value_type` should be `Traits::Polygon_2`. 


The default traits class `Default_traits` is `Partition_traits_2`, 
with the representation type determined by `std::iterator_traits<InputIterator>::%value_type`. 

\sa `CGAL::approx_convex_partition_2()` 
\sa `CGAL::greene_approx_convex_partition_2()` 
\sa `CGAL::optimal_convex_partition_2()` 
\sa `CGAL::partition_is_valid_2()` 
\sa `CGAL::is_convex_2()` 

\cgalHeading{Implementation}

This function calls `partition_is_valid_2()` using the function object 
`Is_convex_2` to determine the convexity of each partition polygon. 
Thus the time required by this function is \f$ O(n \log n + e \log e)\f$ where 
\f$ n\f$ is the total number of vertices in the partition polygons and \f$ e\f$ the 
total number of edges. 

\cgalHeading{Example}

See the example presented with the function `approx_convex_partition_2()` 
for an illustration of the use of this function. 

*/

template<class InputIterator, class ForwardIterator, class Traits>
bool
convex_partition_is_valid_2 (InputIterator point_first, 
InputIterator point_beyond,
ForwardIterator poly_first, 
ForwardIterator poly_beyond,
const Traits& traits = Default_traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgPolygonPartitioning2

\brief returns `true` iff the polygons in the range [`poly_first`, 
`poly_beyond`) define a valid partition of the polygon defined by the 
points in the range [`point_first`, `point_beyond`) and 
`false` otherwise. A valid partition is one in 
which the polygons are nonoverlapping and the union of the polygons is the 
same as the original polygon. 
Each polygon must also satisfy the property 
tested by `Traits::Is_valid()`. 

\pre Points in the range [`point_first`, `point_beyond`) define a simple, counterclockwise-oriented polygon.

\cgalHeading{Requires}

- `Traits` is a model of the concept 
  `PartitionIsValidTraits_2` and the 
  concept defining the requirements for the validity test 
  implemented by `Traits::Is_valid()`. 
- `std::iterator_traits<InputIterator>::%value_type` should be `Traits::Point_2`, 
  which should also be the type of the points stored in an object 
  of type `Traits::Polygon_2`. 
- `std::iterator_traits<ForwardIterator>::%value_type` should be 
  `Traits::Polygon_2`. 

The default traits class `Default_traits` is `Partition_traits_2`, 
with the representation type determined by `std::iterator_traits<InputIterator>::%value_type`. 

\sa `CGAL::approx_convex_partition_2()` 
\sa `CGAL::greene_approx_convex_partition_2()` 
\sa `CGAL::is_y_monotone_2()` 
\sa `CGAL::optimal_convex_partition_2()` 
\sa `CGAL::Partition_is_valid_traits_2<Traits, PolygonIsValid>` 
\sa `CGAL::y_monotone_partition_2()` 
\sa `CGAL::is_convex_2()` 

\cgalHeading{Implementation}

This function requires \f$ O(n \log n + e \log e + \Sigma_{i=1}^p m_i)\f$ where \f$ n\f$ 
is the total number of vertices of the \f$ p\f$ partition polygons, \f$ e\f$ is the 
total number of edges of the partition polygons and \f$ m_i\f$ is the time required 
by `Traits::Is_valid()` to test if partition polygon \f$ p_i\f$ is valid. 

\cgalHeading{Example}

See the example presented with the function `optimal_convex_partition_2()` 
for an illustration of the use of this function. 
*/

template<class InputIterator, class ForwardIterator, class Traits>
bool
partition_is_valid_2 (InputIterator point_first, InputIterator point_beyond,
ForwardIterator poly_first, ForwardIterator poly_beyond,
const Traits& traits = Default_traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgPolygonPartitioning2

\brief determines if the polygons in the range [`poly_first`, `poly_beyond`)
define a valid \f$ y\f$-monotone partition of the simple, counterclockwise-oriented polygon represented by the points 
in the range [`point_first`, `point_beyond`). 
A valid partition is one in 
which the polygons are nonoverlapping and the union of the polygons is the 
same as the original polygon and each polygon is \f$ y\f$-monotone 

\pre P
The function returns `true` iff the partition is valid and otherwise
returns false.



\cgalHeading{Requires}

- `Traits` is a model of the concept 
  `YMonotonePartitionIsValidTraits_2`. 
- `std::iterator_traits<InputIterator>::%value_type` should be `Traits::Point_2`, 
  which should also be the type of the points stored in an object 
  of type `Traits::Polygon_2`. 
- `std::iterator_traits<ForwardIterator>::%value_type` should be 
  `Traits::Polygon_2`. 

The default traits class `Default_traits` is `Partition_traits_2`, 
with the representation type determined by `std::iterator_traits<InputIterator>::%value_type`. 

\sa `CGAL::y_monotone_partition_2()` 
\sa `CGAL::is_y_monotone_2()` 
\sa `CGAL::partition_is_valid_2()` 
\sa `CGAL::Partition_is_valid_traits_2<Traits, PolygonIsValid>` 

\cgalHeading{Implementation}

This function uses the function `partition_is_valid_2()` together with 
the function object `Is_y_monotone_2` to determine if each polygon 
is \f$ y\f$-monotone or not. Thus the time required is \f$ O(n \log n + e \log e)\f$ 
where \f$ n\f$ is the total number of vertices of the partition polygons and 
\f$ e\f$ is the total number of edges. 

\cgalHeading{Example}

See the example presented with the function `y_monotone_partition_2()` 
for an illustration of the use of this function. 

*/

template<class InputIterator, class ForwardIterator, class Traits>
bool
y_monotone_partition_is_valid_2 (InputIterator point_first, 
InputIterator point_beyond,
ForwardIterator poly_first, 
ForwardIterator poly_beyond,
const Traits& traits = Default_traits);

} /* namespace CGAL */
