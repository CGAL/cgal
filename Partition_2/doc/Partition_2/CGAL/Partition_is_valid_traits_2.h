namespace CGAL {

/*!
\ingroup PkgPolygonPartitioning2

Class that derives a traits class for `partition_is_valid_2()` from 
a given traits class by defining the validity testing function object 
in terms of a supplied template parameter. 

\cgalModels `PartitionIsValidTraits_2`

\sa `CGAL::Is_convex_2<Traits>` 
\sa `CGAL::Is_vacuously_valid<Traits>` 
\sa `CGAL::Is_y_monotone_2<Traits>` 
\sa `CGAL::Partition_traits_2<R>` 

\cgalHeading{Example}

See the example presented with the function `optimal_convex_partition_2()` 
for an illustration of the use of this traits class. 

*/
template< typename Traits, typename PolygonIsValid >
class Partition_is_valid_traits_2 : Traits {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef PolygonIsValid Is_valid; 

/*!

*/ 
typedef Traits::Point_2 Point_2; 

/*!

*/ 
typedef Traits::Polygon_2 Polygon_2; 

/*!

*/ 
typedef Traits::Less_xy_2 Less_xy_2; 

/*!

*/ 
typedef Traits::Left_turn_2 Left_turn_2; 

/*!

*/ 
typedef Traits::Orientation_2 Orientation_2; 

/// @} 

/// \name Operations 
/// The constructors and member functions for creating instances of
/// the above types are inherited from `Traits`. In addition, the
/// following member function is defined:
/// @{

/*!
function returning an instance of `Is_valid` 
*/ 
Is_valid 
is_valid_object(const Traits& traits) const; 

/// @}

}; /* end Partition_is_valid_traits_2 */
} /* end namespace CGAL */
