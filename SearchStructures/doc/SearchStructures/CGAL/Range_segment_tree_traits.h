
namespace CGAL {

/*!
\ingroup PkgSearchStructuresTraitsClasses

The class `Range_segment_tree_set_traits_2` is a range and segment tree traits class for the
2-dimensional point class from the \cgal kernel. The class is
parameterized with a representation class `R`.

\cgalModels `RangeSegmentTreeTraits_k`

*/
template< typename R >
class Range_segment_tree_set_traits_2 {
public:

/// \name Types
/// @{

/*!

*/
  typedef R::Point_2 Key;

/*!

*/
typedef std::pair<Key, Key> Interval;

/// @}

}; /* end Range_segment_tree_set_traits_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSearchStructuresTraitsClasses

The class `Range_segment_tree_set_traits_3` is a range and segment tree traits class for the 3-dimensional
point class from the \cgal kernel.
The class is parameterized with a representation class `R`.

\cgalModels `RangeSegmentTreeTraits_k`

*/
template< typename R >
class Range_segment_tree_set_traits_3 {
public:

/// \name Types
/// @{

/*!

*/
  typedef R::Point_3 Key;

/*!

*/
typedef std::pair<Key, Key> Interval;

/// @}

}; /* end Range_segment_tree_set_traits_3 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSearchStructuresTraitsClasses

The class `Range_tree_map_traits_2` is a range tree traits class for the
2-dimensional point class from the \cgal kernel, where data of
type `T` is associated to each key. The class is
parameterized with a representation class `R` and the type of
the associated data `T`.

\cgalModels `RangeSegmentTreeTraits_k`

*/
template< typename R, typename T >
class Range_tree_map_traits_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::pair<R::Point_2,T> Key;

/*!

*/
typedef std::pair<R::Point_2, R::Point_2 > Interval;

/// @}

}; /* end Range_tree_map_traits_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSearchStructuresTraitsClasses

The class `Range_tree_map_traits_3` is a range and segment tree traits class for the 3-dimensional
point class from the \cgal kernel, where data of
type `T` is associated to each key.
The class is parameterized with a representation class `R` and the type of
the associated data `T`.

\cgalModels `RangeSegmentTreeTraits_k`

*/
template< typename R, typename T >
class Range_tree_map_traits_3 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::pair<R::Point_3,T> Key;

/*!

*/
typedef std::pair<R::Point_3, R::Point_3> Interval;

/// @}

}; /* end Range_tree_map_traits_3 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSearchStructuresTraitsClasses

The class `Segment_tree_map_traits_2` is a segment tree traits class for the
2-dimensional point class from the \cgal kernel, where data of
type `T` is associated to each interval. The class is
parameterized with a representation class `R` and the type of
the associated data `T`.

\cgalModels `RangeSegmentTreeTraits_k`

*/
template< typename R, typename T >
class Segment_tree_map_traits_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef R::Point_2 Key;

/*!

*/
typedef std::pair<std::pair<Key,Key>,T> Interval;

/// @}

}; /* end Segment_tree_map_traits_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSearchStructuresTraitsClasses

The class `Segment_tree_map_traits_3` is a segment tree traits class for the 3-dimensional
point class from the \cgal kernel, where data of
type `T` is associated to each interval.
The class is parameterized with a representation class `R` and the type of
the associated data `T`.

\cgalModels `RangeSegmentTreeTraits_k`

*/
template< typename R, typename T >
class Segment_tree_map_traits_3 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::pair<R::Point_3 > Key;

/*!

*/
typedef std::pair<std::pair<Key, Key>,T> Interval;

/// @}

}; /* end Segment_tree_map_traits_3 */
} /* end namespace CGAL */
