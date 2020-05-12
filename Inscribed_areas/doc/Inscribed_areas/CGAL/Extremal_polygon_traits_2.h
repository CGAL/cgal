namespace CGAL {

/*!
\ingroup PkgInscribedAreasRef

\cgalAdvancedClass
\cgalAdvancedBegin
The class `Extremal_polygon_area_traits_2` provides the types and
operations needed to compute a maximum area \f$ k\f$-gon \f$ P_k\f$ that can
be inscribed into a given convex polygon \f$ P\f$ using the function
`extremal_polygon_2`.
\cgalAdvancedEnd

\tparam K must be a model of `Kernel`.

\cgalModels `ExtremalPolygonTraits_2`

\sa `CGAL::maximum_area_inscribed_k_gon_2()`
\sa `CGAL::maximum_perimeter_inscribed_k_gon_2()`
\sa `CGAL::extremal_polygon_2()`
\sa `CGAL::Extremal_polygon_perimeter_traits_2<K>`
\sa `ExtremalPolygonTraits_2`

*/
template< typename K >
struct Extremal_polygon_area_traits_2 {

/// \name Types
/// @{

/*!
typedef to `K::FT`.
*/
typedef unspecified_type FT;

/*!
typedef to `K::Point_2`.
*/
typedef unspecified_type Point_2;

/*!
typedef to `K::Less_xy_2`.
*/
typedef unspecified_type Less_xy_2;

/*!
typedef to `K::Orientation_2`.
*/
typedef unspecified_type Orientation_2;

/*!
AdaptableBinaryFunction class `op`:
`Point_2` \f$ \times\f$ `Point_2` \f$ \rightarrow\f$ `FT`.
For a fixed `Point_2` \f$ root\f$, `op`\f$ (p,\,q)\f$ returns
twice the area of the triangle \f$ (root,\, q,\, p)\f$.
*/
typedef unspecified_type Operation;

/// @}

/// \name Operations
/// @{

/*!
returns 3.
*/
int min_k() const;

/*!
returns `FT(0)`.
*/
FT init(const Point_2& p, const Point_2& q)
const;

/*!
returns `Operation` where `p` is the fixed
\f$ root\f$ point.
*/
Operation operation( const Point_2& p)
const;

/*!
writes the vertices of
[`points_begin`, `points_end`) forming a maximum area
triangle rooted at `points_begin[0]` to o and returns the
past-the-end iterator for that sequence (== `o + 3`).
*/
template < class RandomAccessIterator, class
OutputIterator > OutputIterator compute_min_k_gon(
RandomAccessIterator points_begin, RandomAccessIterator points_end, FT&
max_area, OutputIterator o) const;

/*!

*/
Less_xy_2 less_xy_2_object();

/*!

*/
Orientation_2 orientation_2_object();

/// @}

}; /* end Extremal_polygon_area_traits_2 */

/*!
\ingroup PkgInscribedAreasRef

\cgalAdvancedClass
\cgalAdvancedBegin
The class `Extremal_polygon_perimeter_traits_2` provides the
types and operations needed to compute a maximum perimeter \f$
k\f$-gon \f$ P_k\f$ that can be inscribed into a given convex polygon
\f$ P\f$ using the function `extremal_polygon_2()`.
\cgalAdvancedEnd

\tparam K must be a model of `Kernel`.

\cgalModels `ExtremalPolygonTraits_2`

\sa `CGAL::maximum_area_inscribed_k_gon_2()`
\sa `CGAL::maximum_perimeter_inscribed_k_gon_2()`
\sa `CGAL::extremal_polygon_2()`
\sa `CGAL::Extremal_polygon_area_traits_2<K>`
\sa `ExtremalPolygonTraits_2`

*/
template< typename K >
struct Extremal_polygon_perimeter_traits_2 {

/// \name Types
/// @{

/*!
typedef to `K::FT`.
*/
typedef unspecified_type FT;

/*!
typedef to `K::Point_2`.
*/
typedef unspecified_type Point_2;

/*!
typedef to `K::Less_xy_2`.
*/
typedef unspecified_type Less_xy_2;

/*!
typedef to `K::Orientation_2`.
*/
typedef unspecified_type Orientation_2;

/*!
AdaptableBinaryFunction class `op`:
`Point_2` \f$ \times\f$ `Point_2` \f$ \rightarrow\f$ `FT`.
For a fixed `Point_2` \f$ root\f$, `op`\f$ (p,\,q)\f$ returns
\f$ d(r,\,p) + d(p,\,q) - d(r,\,q)\f$ where \f$ d\f$ denotes the Euclidean
distance.
*/
typedef unspecified_type Operation;

/// @}

/// \name Operations
/// @{

/*!
returns 2.
*/
int min_k() const;

/*!
returns twice the Euclidean distance between `p` and
`q`.
*/
FT init(const Point_2& p, const Point_2& q)
const;

/*!
returns `Operation` where `p` is the fixed
\f$ root\f$ point.
*/
Operation operation( const Point_2& p)
const;

/*!
writes the pair
(`points_begin[0]`, `p`) where `p` is drawn from
[`points_begin`, `points_end`) such that the Euclidean
distance between both points is maximized (maximum perimeter
2-gon rooted at `points_begin[0]`) to o and returns the
past-the-end iterator for that sequence (== `o + 2`).
*/
template < class RandomAccessIterator, class
OutputIterator > OutputIterator compute_min_k_gon(
RandomAccessIterator points_begin, RandomAccessIterator points_end, FT&
max_area, OutputIterator o) const;

/*!

*/
Less_xy_2 less_xy_2_object();

/*!

*/
Orientation_2 orientation_2_object();

/// @}

}; /* end Extremal_polygon_perimeter_traits_2 */
} /* end namespace CGAL */
