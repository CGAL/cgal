namespace CGAL {

/*!
\ingroup PkgPolygonPartitioning2

determines if the sequence of points in the range 
[`first`, `beyond`) defines a \f$ y\f$-monotone 
polygon or not. If so, the function returns `true`, otherwise it
returns `false`. 

\cgalHeading{Requires}

- `Traits` is a model of the concept `IsYMonotoneTraits_2`. 
- `std::iterator_traits<InputIterator>::%value_type` should be `Traits::Point_2`. 

The default traits class `Default_traits` is the kernel in which the 
type `std::iterator_traits<InputIterator>::%value_type` is defined. 

\sa `CGAL::Is_y_monotone_2<Traits>` 
\sa `CGAL::y_monotone_partition_2()` 
\sa `CGAL::y_monotone_partition_is_valid_2()` 

\cgalHeading{Implementation}

This function requires \f$ O(n)\f$ time for a polygon with \f$ n\f$ vertices. 

\cgalHeading{Example}

The following program computes a \f$ y\f$-monotone partitioning 
of a polygon using the default 
traits class and stores the partition polygons in the list 
`partition_polys`. It then asserts that each of the partition 
polygons is, in fact, a \f$ y\f$-monotone polygon and that the partition 
is valid. (Note that the 
assertions are superfluous unless the postcondition checking done 
by `y_monotone_partition_2()` has been turned off during compilation.) 

\cgalExample{Partition_2/y_monotone_partition_2.cpp} 

*/

template<class InputIterator, class Traits>
bool
is_y_monotone_2(InputIterator first, InputIterator beyond,
const Traits& traits);

} /* namespace CGAL */
