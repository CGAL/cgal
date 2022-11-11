
/*!
\ingroup PkgSegmentDelaunayGraph2Concepts
\cgalConcept

The concept `SegmentDelaunayGraphStorageSite_2` provides the
requirements for the storage sites of a segment Delaunay graph. The
storage sites are sites that are used to store the information of a
site in a more compact form (that uses less storage). This is achieved
by storing handles to points instead of points.

\cgalRefines `DefaultConstructible`
\cgalRefines `CopyConstructible`
\cgalRefines `Assignable`

\cgalHasModel `CGAL::Segment_Delaunay_graph_storage_site_2<Gt>`

\sa `SegmentDelaunayGraphTraits_2`
\sa `CGAL::Segment_Delaunay_graph_site_2<K>`
\sa `CGAL::Segment_Delaunay_graph_storage_site_2<Gt>`
\sa `CGAL::Segment_Delaunay_graph_traits_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>`

*/

class SegmentDelaunayGraphStorageSite_2 {
public:

/// \name Types
/// @{

/*!
The site type.
*/
typedef unspecified_type Site_2;

/*!
The type for a handle to a point.
*/
typedef typename std::set<typename Site_2::Point_2>::iterator
Point_handle;

/// @}

/// \name Creation
/// In addition to the default and copy constructors, the following
/// static methods should be available for constructing sites:
/// @{

/*!
Constructs a storage site from a
point handle. The storage site represents the point associated with
the point handle `hp`.
*/
SegmentDelaunayGraphStorageSite_2
construct_storage_site_2(Point_handle hp);

/*!
Constructs a storage site from two point handles. The storage site
represents the segment the endpoints of which are the points
associated with the point handles `hp1` and `hp2`.
*/
SegmentDelaunayGraphStorageSite_2
construct_storage_site_2(Point_handle hp1, Point_handle hp2);

/*!
Constructs
a storage site from four point handles. The storage site represents
the point of intersection of the segments the endpoints of which are
the points associated with the point handles `hp1`, `hp2` and
`hq1` and `hq2`, respectively.
*/
SegmentDelaunayGraphStorageSite_2
construct_storage_site_2(Point_handle hp1,
Point_handle hp2, Point_handle hq1, Point_handle hq2);

/*!
Constructs
a site from four point handles and a Boolean. The storage site
represents a segment. If `b` is `true`, the first endpoint
of the segment is the point associated with the handle `hp1` and
the second endpoint is the point of intersection of the segments the
endpoints of which are the point associated with the point handles
`hp1`, `hp2` and `hq1`, `hq2`, respectively. If
`b` is `false`, the first endpoint of the represented
segment is the one mentioned above, whereas the second endpoint if
the point associated with the point handle `hp2`.
*/
SegmentDelaunayGraphStorageSite_2
construct_storage_site_2(Point_handle hp1, Point_handle hp2,
Point_handle hq1, Point_handle hq2, bool b);

/*!
Constructs a storage site from six
point handles. The storage site represents of segment the endpoints
of which are points of intersection of two pairs of segments, the
endpoints of which are `hp1`, `hp2`/`hq1`, `hq2` and
`hp1`, `hp2`/`hr1`, `hr2`, respectively.
*/
SegmentDelaunayGraphStorageSite_2
construct_storage_site_2(Point_handle hp1,
Point_handle hp2, Point_handle hq1, Point_handle hq2, Point_handle
hr1, Point_handle hr2);

/// @}

/// \name Predicates
/// @{

/*!
Returns `true` if the storage site
represents a valid point or segment.
*/
bool is_defined();

/*!
Returns `true` if the storage site
represents a point.
*/
bool is_point();

/*!
Returns `true` if the storage site
represents a segment.
*/
bool is_segment();

/*!
Returns `true` if the storage site
represents an input point or a segment defined by two input
points. Returns `false` if it represents a point of intersection
of two segments, or if it represents a segment, at least one
endpoint of which is a point of intersection of two segments.
*/
bool is_input();

/*!
Returns `true` if the
`i`-th endpoint of the corresponding site is an input
point. Returns `false` if the `i`-th endpoint of the
corresponding site is the intersection of two segments.
\pre `i` must be at most \f$ 1\f$, and `ss.is_segment()` must be `true`.
*/
bool is_input(unsigned int i);

/// @}

/// \name Access Functions
/// @{

/*!
Returns a storage site object representing the segment
that supports the segment represented by the storage site.
The returned storage site represents a site, both endpoints
of which are input points.
\pre `ss.is_segment()` must be `true`.
*/
SegmentDelaunayGraphStorageSite_2 supporting_site();

/*!
Returns a storage site that represents the first endpoint of the
represented segment.
\pre `ss.is_segment()` must be `true`.
*/
SegmentDelaunayGraphStorageSite_2 source_site();

/*!
Returns a storage site that represents the second endpoint of the
represented segment.
\pre `ss.is_segment()` must be `true`.
*/
SegmentDelaunayGraphStorageSite_2 target_site();

/*!
Returns a storage site object representing the `i`-th
segment that supports the point of intersection represented
by the storage site.
The returned storage site represents a site, both endpoints
of which are input points.
\pre `i` must be at most \f$ 1\f$, `ss.is_point()` must be `true` and `ss.is_input()` must be `false`.
*/
SegmentDelaunayGraphStorageSite_2 supporting_site(unsigned int i);

/*!
Returns a storage site object representing the `i`-th
segment that supports the \f$ i\f$-th endpoint of the site
which is not the supporting segment of the site.
The returned storage site represents a site, both endpoints
of which are input points.
\pre `i` must be at most \f$ 1\f$, `ss.is_segment()` must be `true` and `ss.is_input(i)` must be `false`.
*/
SegmentDelaunayGraphStorageSite_2 crossing_site(unsigned int i);

/*!
Returns the site represented by the storage
site.
*/
Site_2 site();

/*!
Returns a handle associated with
the represented point. \pre `is_point()` and `is_input()` must both be `true`.
*/
Point_handle point();

/*!
Returns a handle to the source point of the supporting site of the
this site. \pre `is_segment()` must be `true`.
*/
Point_handle source_of_supporting_site();

/*!
Returns a handle to the target point of the supporting site of the
this site. \pre `is_segment()` must be `true`.
*/
Point_handle target_of_supporting_site();

/*!
Returns a handle to the source point of the `i`-th supporting
site of the this site.
\pre `is_point()` must be `true`, `is_input()` must be `false` and `i` must either be `0` or `1`.
*/
Point_handle source_of_supporting_site(unsigned int i);

/*!
Returns a handle to the target point of the `i`-th supporting
site of the this site.
\pre `is_point()` must be `true`, `is_input()` must be `false` and `i` must either be `0` or `1`.
*/
Point_handle target_of_supporting_site(unsigned int i);

/*!
Returns a handle to the source point of the `i`-th crossing site
of the this site.
\pre `is_segment()` must be `true`, `is_input(i)` must be `false` and `i` must either be `0` or `1`.
*/
Point_handle source_of_crossing_site(unsigned int i);

/*!
Returns a handle to the target point of the `i`-th supporting
site of the this site.
\pre `is_segment()` must be `true`, `is_input(i)` must be `false` and `i` must either be `0` or `1`.
*/
Point_handle target_of_crossing_site(unsigned int i);

/// @}

}; /* end SegmentDelaunayGraphStorageSite_2 */

