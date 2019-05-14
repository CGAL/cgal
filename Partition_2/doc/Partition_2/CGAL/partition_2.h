namespace CGAL {

/*!
\ingroup PkgPolygonPartitioning2

\brief computes a partition of the polygon defined by the points in
the range [`first`, `beyond`) into convex polygons. The
counterclockwise-oriented partition polygons are written to the
sequence starting at position `result`. The past-the-end iterator for
the resulting sequence of polygons is returned.  

The number of convex polygons produced is 
no more than four times the minimal number. 

\pre The points in the range [`first`, `beyond`) define a simple
counterclockwise-oriented polygon.

\cgalHeading{Requirements}

<OL> 
<LI>`Traits` is a model of the concept 
`PartitionTraits_2` 
and, for the purposes of checking the postcondition that the partition 
produced is valid, it should also be a model of 
the concept `ConvexPartitionIsValidTraits_2`. 
<LI>`std::iterator_traits<OutputIterator>::%value_type` should be `Traits::Polygon_2`. 
<LI>`std::iterator_traits<InputIterator>::%value_type` should be `Traits::Point_2`, 
which should also be the type of the points stored in an object 
of type `Traits::Polygon_2`. 
<LI>Points in the range `[first, beyond)` must define a simple polygon 
whose vertices are oriented counterclockwise. 
</OL> 

The default traits class `Default_traits` is `Partition_traits_2`, 
with the representation type determined by `std::iterator_traits<InputIterator1>::%value_type`. 

\sa `CGAL::convex_partition_is_valid_2()` 
\sa `CGAL::greene_approx_convex_partition_2()` 
\sa `CGAL::optimal_convex_partition_2()` 
\sa `CGAL::partition_is_valid_2()` 
\sa `CGAL::Partition_is_valid_traits_2<Traits, PolygonIsValid>` 
\sa `CGAL::y_monotone_partition_2()` 

\cgalHeading{Implementation}

This function implements the algorithm of Hertel and Mehlhorn 
\cgalCite{hm-ftsp-83} and is based on the class 
`Constrained_triangulation_2`. Given a triangulation of 
the polygon, the function requires \f$ O(n)\f$ time and 
space for a polygon with \f$ n\f$ vertices. 

\cgalHeading{Example}

The following program computes an approximately optimal 
convex partitioning of a polygon using the default 
traits class and stores the partition polygons in the list 
`partition_polys`. 

\cgalExample{Partition_2/approx_convex_partition_2.cpp} 

*/
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator approx_convex_partition_2(InputIterator first,
InputIterator beyond,
OutputIterator result,
const Traits& traits = Default_traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgPolygonPartitioning2

\brief computes a partition of the polygon defined 
by the points in the range [`first`, `beyond`) into convex 
polygons. The counterclockwise-oriented partition polygons are written to
the sequence starting at position `result`. 
The number of convex polygons produced is 
no more than four times the minimal number. 
The past-the-end iterator for 
the resulting sequence of polygons is returned.



\pre The points in the range [`first`, `beyond`) define a simple, counterclockwise-oriented polygon.

\cgalHeading{Requirements}

<OL> 
<LI>`Traits` is a model of the concepts `PartitionTraits_2` 
and `YMonotonePartitionTraits_2`. 
For the purpose of 
checking the validity of the \f$ y\f$-monotone partition produced as 
a preprocessing step for the convex partitioning, it must also 
be a model of `YMonotonePartitionIsValidTraits_2`. 
For the purpose of checking 
the postcondition that the convex partition is valid, `Traits` 
must also be a model of `ConvexPartitionIsValidTraits_2`. 
<LI>`std::iterator_traits<OutputIterator>::%value_type` is equivalent to 
`Traits::Polygon_2`. 
<LI>`std::iterator_traits<InputIterator>::%value_type` is equivalent to 
`Traits::Point_2`, 
which should also be equivalent to the type of the points stored in 
an object of type `Traits::Polygon_2`. 
</OL> 

The default traits class `Default_traits` is `Partition_traits_2`, 
with the representation type determined by `std::iterator_traits<InputIterator>::%value_type`. 

\sa `CGAL::approx_convex_partition_2()` 
\sa `CGAL::convex_partition_is_valid_2()` 
\sa `CGAL::optimal_convex_partition_2()` 
\sa `CGAL::partition_is_valid_2()` 
\sa `CGAL::y_monotone_partition_2()` 

\cgalHeading{Implementation}

This function implements the approximation algorithm of 
Greene \cgalCite{g-dpcp-83} and requires \f$ O(n \log n)\f$ time and \f$ O(n)\f$ space 
to produce a convex partitioning given a \f$ y\f$-monotone partitioning of a 
polygon with \f$ n\f$ vertices. The function `y_monotone_partition_2()` 
is used to produce the monotone partition. 

\cgalHeading{Example}

The following program computes an approximately optimal 
convex partitioning of a polygon using the default 
traits class and stores the partition polygons in the list 
`partition_polys`. 

\cgalExample{Partition_2/greene_approx_convex_partition_2.cpp} 

*/

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator greene_approx_convex_partition_2(InputIterator first,
InputIterator beyond,
OutputIterator result,
const Traits& traits = Default_traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgPolygonPartitioning2

\brief computes a partition of the polygon defined
by the points in the range [`first`, `beyond`) into convex
polygons. The counterclockwise-oriented partition polygons are written to
the sequence starting at position `result`. 
The number of convex polygons produced is minimal. 
The past-the-end iterator for
the resulting sequence of polygons is returned.



\pre The points in the range [`first`, `beyond`) define a simple, counterclockwise-oriented polygon.

\cgalHeading{Requirements}

<OL> 
<LI>`Traits` is a model of the concept `OptimalConvexPartitionTraits_2`. 
For the purposes of checking the 
postcondition that the partition is valid, `Traits` should 
also be a model of `ConvexPartitionIsValidTraits_2`. 

<LI>`std::iterator_traits<OutputIterator>::%value_type` should be 
`Traits::Polygon_2`. 
<LI>`std::iterator_traits<InputIterator>::%value_type` should be `Traits::Point_2`, 
which should also be the type of the points stored in an object 
of type `Traits::Polygon_2`. 
</OL> 

The default traits class `Default_traits` is `Partition_traits_2`, 
with the representation type determined by `std::iterator_traits<InputIterator>::%value_type`. 

\sa `CGAL::approx_convex_partition_2()` 
\sa `CGAL::convex_partition_is_valid_2()` 
\sa `CGAL::greene_approx_convex_partition_2()` 
\sa `CGAL::partition_is_valid_2()` 
\sa `CGAL::Partition_is_valid_traits_2<Traits, PolygonIsValid>` 

\cgalHeading{Implementation}

This function implements the dynamic programming algorithm of Greene 
\cgalCite{g-dpcp-83}, which requires \f$ O(n^4)\f$ time and \f$ O(n^3)\f$ space to 
produce a partitioning of a polygon with \f$ n\f$ vertices. 

\cgalHeading{Example}

The following program computes an optimal 
convex partitioning of a polygon using the default 
traits class and stores the partition polygons in the list 
`partition_polys`. 
It then asserts that the partition produced is valid. The 
traits class used for testing the validity is derived from the 
traits class used to produce the partition with the function object 
class `Is_convex_2` used 
to define the required `Is_valid` type. 
(Note that this assertion is superfluous unless the 
postcondition checking for `optimal_convex_partition_2()` has been 
turned off.) 

\cgalExample{Partition_2/optimal_convex_partition_2.cpp} 

*/

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator optimal_convex_partition_2(InputIterator first,
InputIterator beyond,
OutputIterator result,
const Traits& traits = Default_traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgPolygonPartitioning2

\brief computes a partition of the polygon defined 
by the points in the range [`first`, `beyond`) into \f$ y\f$-monotone 
polygons. The counterclockwise-oriented partition polygons are written to
the sequence starting at position `result`. The past-the-end iterator for 
the resulting sequence of polygons is returned.

\pre The points in the range [`first`, `beyond`) define a simple, counterclockwise-oriented polygon.


\cgalHeading{Requirements}

<OL> 
<LI>`Traits` is a model of the concept 
`YMonotonePartitionTraits_2` 
and, for the purposes 
of checking the postcondition that the partition is valid, it should 
also be a model of `YMonotonePartitionIsValidTraits_2`. 
<LI>`std::iterator_traits<OutputIterator>::%value_type` should be 
`Traits::Polygon_2`. 
<LI>`std::iterator_traits<InputIterator>::%value_type` should be `Traits::Point_2`, 
which should also be the type of the points stored in an object 
of type `Traits::Polygon_2`. 
</OL> 

The default traits class `Default_traits` is `Partition_traits_2`, 
with the representation type determined by `std::iterator_traits<InputIterator>::%value_type`. 

\sa `CGAL::approx_convex_partition_2()` 
\sa `CGAL::greene_approx_convex_partition_2()` 
\sa `CGAL::optimal_convex_partition_2()` 
\sa `CGAL::partition_is_valid_2()` 
\sa `CGAL::y_monotone_partition_is_valid_2()` 

\cgalHeading{Implementation}

This function implements the algorithm presented by de Berg <I>et al.</I> 
\cgalCite{bkos-cgaa-97} which requires \f$ O(n \log n)\f$ time 
and \f$ O(n)\f$ space for a polygon with \f$ n\f$ vertices. 

\cgalHeading{Example}

The following program computes a \f$ y\f$-monotone partitioning 
of a polygon using the default 
traits class and stores the partition polygons in the list 
`partition_polys`. It then asserts that each partition polygon 
produced is, in fact, \f$ y\f$-monotone and that the partition is valid. 
(Note that these assertions are superfluous unless the postcondition 
checking for `y_monotone_partition_2()` has been turned off.) 

\cgalExample{Partition_2/y_monotone_partition_2.cpp} 

*/
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator y_monotone_partition_2(InputIterator first, 
InputIterator beyond,
OutputIterator result, 
const Traits& traits = Default_traits);

} /* namespace CGAL */
