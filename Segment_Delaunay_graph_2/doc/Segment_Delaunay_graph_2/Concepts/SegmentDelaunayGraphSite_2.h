
/*!
\ingroup PkgSegmentDelaunayGraph2Concepts
\cgalConcept

The concept `SegmentDelaunayGraphSite_2` provides the 
requirements for the sites of a segment Delaunay graph. 

\cgalRefines `DefaultConstructible` 
\cgalRefines `CopyConstructible` 
\cgalRefines `Assignable` 

\cgalHasModel `CGAL::Segment_Delaunay_graph_site_2<K>` 

\sa `SegmentDelaunayGraphTraits_2` 
\sa `CGAL::Segment_Delaunay_graph_site_2<K>` 
\sa `CGAL::Segment_Delaunay_graph_traits_2<K,MTag>` 
\sa `CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K,MTag>` 
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>` 
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>` 

*/

class SegmentDelaunayGraphSite_2 {
public:

/// \name Types 
/// @{

/*!
The point type. 
*/ 
typedef unspecified_type Point_2; 

/*!
The segment type. 
*/ 
typedef unspecified_type Segment_2; 

/*!
The field number type. 
*/ 
typedef unspecified_type FT; 

/*!
The ring number type. 
*/ 
typedef unspecified_type RT; 

/// @} 

/// \name Creation 
/// In addition to the default and copy constructors the following
/// static methods are available for constructing sites:
/// @{

/*!
Constructs a site from a point: the site 
represents the point `p`. 
*/ 
static SegmentDelaunayGraphSite_2 construct_site_2(Point_2 p); 

/*!
Constructs a site from two points: the site represents 
the (open) segment `(p1,p2)`. 
*/ 
static SegmentDelaunayGraphSite_2 construct_site_2(Point_2 p1, Point_2 p2); 

/*!
Constructs a site from four points: the site represents the 
point of intersection of the segments `(p1,p2)` and 
`(q1,q2)`. 
*/ 
static SegmentDelaunayGraphSite_2 construct_site_2(Point_2 p1, Point_2 p2, 
                                                   Point_2 q1, Point_2 q2); 

/*!
Constructs a site from four points and a Boolean: the 
site represents a segment. If `b` is `true` the endpoints 
are `p1` and \f$ p_\times\f$, otherwise \f$ p_\times\f$ and 
`p2`. \f$ p_\times\f$ is the point of intersection of the segments 
`(p1,p2)`,`(q1,q2)`. 
*/ 
static SegmentDelaunayGraphSite_2 construct_site_2(Point_2 p1, Point_2 p2, 
                                                   Point_2 q1, Point_2 q2, bool b); 

/*!
Constructs a site from six 
points: the site represents the segment with endpoints the points of 
intersection of the pairs of segments `(p1,p2)`,`(q1,q2)` 
and `(p1,p2)`,`(r1,r2)`. 
*/ 
static SegmentDelaunayGraphSite_2 construct_site_2(Point_2 p1, Point_2 p2, 
                                                   Point_2 q1, Point_2 q2, 
                                                   Point_2 r1, Point_2 r2); 

/// @} 

/// \name Predicates 
/// @{

/*!
Returns `true` if the site 
represents a valid point or segment. 
*/ 
bool is_defined(); 

/*!
Returns `true` if the site represents 
a point. 
*/ 
bool is_point(); 

/*!
Returns `true` if the site 
represents a segment. 
*/ 
bool is_segment(); 

/*!
Returns `true` if the site 
represents an input point or a segment defined by two input 
points. Returns `false` if it represents a point of intersection 
of two segments, or if it represents a segment, at least one 
endpoint of which is a point of intersection of two segments. 
*/ 
bool is_input(); 

/*!
Returns `true` if the 
`i`-th endpoint of the site is an input point. Returns `false` 
if the `i`-th endpoint of the site is the intersection of two 
segments. 
\pre `i` must be at most \f$ 1\f$, and `s.is_segment()` must be `true`. 
*/ 
bool is_input(unsigned int i); 

/// @} 

/// \name Access Functions 
/// @{

/*!
Returns the point represented by the 
site `s`. 
\pre `s.is_point()` must be `true`. 
*/ 
Point_2 point() const; 

/*!
Returns the segment represented 
by the site `s`. 
\pre `s.is_segment()` must be `true`. 
*/ 
Segment_2 segment() const; 

/*!
Returns the source endpoint of the 
segment. Note that this method can construct an inexact point if the 
number type used is inexact. 
\pre `s.is_segment()` must be `true`. 
*/ 
Point_2 source() const; 

/*!
Returns the target endpoint of the 
segment. Note that this method can construct an inexact point if the 
number type used is inexact. 
\pre `s.is_segment()` must be `true`. 
*/ 
Point_2 target() const; 

/*!
Returns a segment site object representing the segment 
that supports the segment represented by the site. Both 
endpoints of the returned site are input points. 
\pre `s.is_segment()` must be `true`. 
*/ 
SegmentDelaunayGraphSite_2 supporting_site(); 

/*!
Returns a segment site object representing the `i`-th 
segment that supports the point of intersection represented 
by the site. Both endpoints of the returned site are input 
points. 
\pre `i` must be at most \f$ 1\f$, `s.is_point()` must be `true` and `s.is_input()` must be `false`. 
*/ 
SegmentDelaunayGraphSite_2 supporting_site(unsigned int i); 

/*!
Returns a segment site object representing the `i`-th 
segment that supports the \f$ i\f$-th endpoint of the site 
which is not the supporting segment of the site. Both 
endpoints of the returned site are input points. 
\pre `i` must be at most \f$ 1\f$, `s.is_segment()` must be `true` and `s.is_input(i)` must be `false`. 
*/ 
SegmentDelaunayGraphSite_2 crossing_site(unsigned int i); 

/*!
Returns a point site object representing the source point of 
the site. 
\pre `s.is_segment()` must be `true`. 
*/ 
SegmentDelaunayGraphSite_2 source_site(); 

/*!
Returns a point site object representing the target point of 
the site. 
\pre `s.is_segment()` must be `true`. 
*/ 
SegmentDelaunayGraphSite_2 target_site(); 

/*!
Returns the source point of the supporting site of the this site. 
\pre `is_segment()` must be `true`. 
*/ 
Point_2 source_of_supporting_site(); 

/*!
Returns the target point of the supporting site of the this site. 
\pre `is_segment()` must be `true`. 
*/ 
Point_2 target_of_supporting_site(); 

/*!
Returns the source point of the `i`-th supporting site of the 
this site. 
\pre `is_point()` must be `true`, `is_input()` must be `false` and `i` must either be `0` or `1`. 
*/ 
Point_2 source_of_supporting_site(unsigned int i); 

/*!
Returns the target point of the `i`-th supporting site of the 
this site. 
\pre `is_point()` must be `true`, `is_input()` must be `false` and `i` must either be `0` or `1`. 
*/ 
Point_2 target_of_supporting_site(unsigned int i); 

/*!
Returns the source point of the `i`-th crossing site of the 
this site. 
\pre `is_segment()` must be `true`, `is_input(i)` must be `false` and `i` must either be `0` or `1`. 
*/ 
Point_2 source_of_crossing_site(unsigned int i); 

/*!
Returns the target point of the `i`-th supporting site of the 
this site. 
\pre `is_segment()` must be `true`, `is_input(i)` must be `false` and `i` must either be `0` or `1`. 
*/ 
Point_2 target_of_crossing_site(unsigned int i); 

/// @}

}; /* end SegmentDelaunayGraphSite_2 */

