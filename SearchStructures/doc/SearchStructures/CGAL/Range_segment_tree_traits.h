
namespace CGAL {

/*!
\ingroup PkgRangeSegmentTreesDTraitsClasses

The class `Range_segment_tree_traits_set_2` is a range and segment tree traits class for the 
2-dimensional point class from the \cgal kernel. The class is 
parameterized with a representation class `R`. 

*/
template< typename R >
class Range_segment_tree_traits_set_2 {
public:

/// \name Types 
/// @{

/*!

*/ 
Point_2<R> Key; 

/*!

*/ 
std::pair<Key, Key> Interval; 

/// @}

}; /* end Range_segment_tree_traits_set_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgRangeSegmentTreesDTraitsClasses

The class `Range_segment_tree_traits_set_3` is a range and segment tree traits class for the 3-dimensional 
point class from the \cgal kernel. 
The class is parameterized with a representation class `R`. 

*/
template< typename R >
class Range_segment_tree_traits_set_3 {
public:

/// \name Types 
/// @{

/*!

*/ 
Point_3<R> Key; 

/*!

*/ 
std::pair<Key, Key> Interval; 

/// @}

}; /* end Range_segment_tree_traits_set_3 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgRangeSegmentTreesDTraitsClasses

The class `Range_tree_traits_map_2` is a range tree traits class for the 
2-dimensional point class from the \cgal kernel, where data of 
type `T` is associated to each key. The class is 
parameterized with a representation class `R` and the type of 
the associated data `T`. 

*/
template< typename R, typename T >
class Range_tree_traits_map_2 {
public:

/// \name Types 
/// @{

/*!

*/ 
std::pair<R::Point_2,T> Key; 

/*!

*/ 
std::pair<R::Point_2, R::Point_2 > Interval; 

/// @}

}; /* end Range_tree_traits_map_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgRangeSegmentTreesDTraitsClasses

The class `Range_tree_traits_map_3` is a range and segment tree traits class for the 3-dimensional 
point class from the \cgal kernel, where data of 
type `T` is associated to each key. 
The class is parameterized with a representation class `R` and the type of 
the associated data `T`. 

*/
template< typename R, typename T >
class Range_tree_traits_map_3 {
public:

/// \name Types 
/// @{

/*!

*/ 
std::pair<R::Point_3,T> Key; 

/*!

*/ 
std::pair<R::Point_3, R::Point_3> Interval; 

/// @}

}; /* end Range_tree_traits_map_3 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgRangeSegmentTreesDTraitsClasses

The class `Segment_tree_traits_map_2` is a segment tree traits class for the 
2-dimensional point class from the \cgal kernel, where data of 
type `T` is associated to each interval. The class is 
parameterized with a representation class `R` and the type of 
the associated data `T`. 

*/
template< typename R, typename T >
class Segment_tree_traits_map_2 {
public:

/// \name Types 
/// @{

/*!

*/ 
R::Point_2 Key; 

/*!

*/ 
std::pair<std::pair<Key,Key>,T> Interval; 

/// @}

}; /* end Segment_tree_traits_map_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgRangeSegmentTreesDTraitsClasses

The class `Segment_tree_traits_map_3` is a segment tree traits class for the 3-dimensional 
point class from the \cgal kernel, where data of 
type `T` is associated to each interval. 
The class is parameterized with a representation class `R` and the type of 
the associated data `T`. 

*/
template< typename R, typename T >
class Segment_tree_traits_map_3 {
public:

/// \name Types 
/// @{

/*!

*/ 
std::pair<R::Point_3 > Key; 

/*!

*/ 
std::pair<std::pair<Key, Key>,T> Interval; 

/// @}

}; /* end Segment_tree_traits_map_3 */
} /* end namespace CGAL */
