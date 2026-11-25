
namespace CGAL {

/*!
\ingroup PkgTriangulation2IntersectionTags

\deprecated This class is deprecated since \cgal 5.1, as it was ambiguous. It is equivalent to
`No_constraint_intersection_requiring_constructions_tag`. Users should instead use the tags
`No_constraint_intersection_tag` or `No_constraint_intersection_requiring_constructions_tag`,
depending on their needs.
*/
struct No_intersection_tag{};

/*!
\ingroup PkgTriangulation2IntersectionTags

Intersection tag for constrained triangulations, when constraints are not allowed to intersect,
except at a single common extremity, and vertices are not allowed to lie in the interior of a constraint.

More precisely, the following conditions must hold at all times:
- two constraints are either disjoint, or intersect at a single point that is an extremity of both.
  In particular, 'T'-junctions, crossings, and overlapping or identical constraints are forbidden;
- no vertex of the triangulation lies in the interior of a constraint. In particular, inserting
  a point in the interior of a constraint, or a constraint whose interior contains a vertex, is forbidden.

As a consequence, constraints are never split: each constraint corresponds to exactly one
constrained edge of the triangulation. Whether a sequence of insertions satisfies these conditions
depends only on the inserted points and constraints, and not on the order in which they are inserted.

An exception of type `Constrained_triangulation_2::Intersection_of_constraints_exception`
is thrown upon any insertion that would break these conditions.
*/
struct No_constraint_intersection_tag{};

/*!
\ingroup PkgTriangulation2IntersectionTags

Intersection tag for constrained triangulations, when constraints are allowed to intersect,
except if representing the intersection requires the construction of a new point.

More precisely, consider the arrangement of the constraints and of the vertices already in the triangulation.
A new constraint can be inserted if and only if adding it to this arrangement creates no vertex other
than the existing vertices and the extremities of the new constraint. In other words, it may not cross
an existing constraint at a point that lies in the interior of both and that is not already a vertex
of the triangulation. All other configurations are allowed, for example 'T'-junctions,
overlapping or identical constraints, constraints passing through existing vertices,
and points inserted in the interior of constraints. Constraints are split at the vertices lying
in their interior, so that a constraint may correspond to several constrained edges of the triangulation.

An exception of type `Constrained_triangulation_2::Intersection_of_constraints_exception`
is thrown upon insertion of a constraint that would break this condition.

This is the default tag.
*/
struct No_constraint_intersection_requiring_constructions_tag{};

/*!
\ingroup PkgTriangulation2IntersectionTags

Intersection tag for constrained triangulations, when constraints may intersect in any way,
and the geometric traits class provides exact predicates but approximate constructions.

No configuration is forbidden. When two constraints cross (see `No_constraint_intersection_requiring_constructions_tag`),
their intersection point is constructed and inserted as a new vertex, and both constraints are split at this vertex.

All other configurations are handled as with `No_constraint_intersection_requiring_constructions_tag`.
As constructions are approximate, a constructed vertex may not lie exactly on the constraints,
and constrained edges may thus slightly deviate from the input constraints;
the algorithm is designed to remain robust in this setting.

When using this tag, be wary of inserting the same constraint multiple times, as this may
cause problems in case of intersecting constraints. Indeed, when an inserted segment is subdivided, the resulting subsegments
are slightly displaced by numerical error and may therefore intersect the next insertion
of the same segment, creating additional vertices at each new insertion.
Functions are provided to insert unique constraints.

The geometric traits class must be a model of `ConstrainedTriangulationTraits_2`.
No exception of type `Constrained_triangulation_2::Intersection_of_constraints_exception` is thrown.
*/
struct Exact_predicates_tag{};

/*!
\ingroup PkgTriangulation2IntersectionTags

Intersection tag for constrained triangulations, when constraints may intersect in any way,
and the geometric traits class provides exact predicates and exact constructions.

Same as `Exact_predicates_tag`, except that intersection points are computed exactly
and thus lie exactly on the constraints. Using `Constrained_triangulation_plus_2` is recommended
with this tag, as it avoids cascading intersection computations.

The geometric traits class must be a model of `ConstrainedTriangulationTraits_2`.
No exception of type `Constrained_triangulation_2::Intersection_of_constraints_exception` is thrown.
*/
struct Exact_intersections_tag{};

/*!
\ingroup PkgTriangulation2TriangulationClasses

A constrained triangulation is a triangulation of a set of points
which has to include among its edges
a given set of polylines joining the points.

The given polylines are
called *constraints* and the corresponding
edges in the triangulation are called *constrained edges*.


The endpoints of constrained edges are of course vertices of the
triangulation. However the triangulation may include
other vertices as well. Any number of constrained edges may share the same endpoint.

The intersection tag `Itag` determines which configurations of constraints are allowed
(see Section \ref PkgTriangulation2IntersectionTags).

In order to retrieve the constrained edges of a constraint, or
the constraints overlapping with a constrained edge, we provide
the class `Constrained_triangulation_plus_2`. This class maintains
a constraint hierarchy data structure. See
Section \ref Section_2D_Triangulations_Constrained_Plus for details.
This class should also be used when doing exact intersection computations
as it avoids the cascading of intersection computations.

\image html constraints.png
\image latex constraints.png

\tparam Traits is a geometric traits class and must be a model
of the concept `TriangulationTraits_2`.
When intersections of input constraints are supported (i.e., when
the intersection tag is `Exact_predicates_tag` or `Exact_intersections_tag`)
the geometric traits class is required to provide additional function object types
to compute the intersection of two segments.
It has then to be a model of the concept
`ConstrainedTriangulationTraits_2`.

\tparam Tds must be a model of the concept `TriangulationDataStructure_2` or `Default`.
The information about constrained edges is stored in the faces of the triangulation. Thus the nested `Face`
type of a constrained triangulation offers additional functionalities to deal with this information.
These additional functionalities induce additional requirements on the face base class
plugged into the triangulation data structure of a constrained Delaunay triangulation.
The face base of a constrained Delaunay triangulation has to be a model of the concept
`ConstrainedTriangulationFaceBase_2`.

\tparam Itag is the intersection tag which serves to choose between the different
strategies to deal with constraints intersections.
It must be one of the tags described in Section \ref PkgTriangulation2IntersectionTags.

\cgal provides default instantiations for the template parameters `Tds` and `Itag`.
If `Gt` is the geometric traits class parameter, the default triangulation data structure
is the class is the class `Triangulation_data_structure_2<Triangulation_vertex_base_2<Gt>, Constrained_triangulation_face_base_2<Gt> >`.
The default intersection tag is `No_constraint_intersection_requiring_constructions_tag`.

\sa `CGAL::Triangulation_2<Traits,Tds>`
\sa `TriangulationDataStructure_2`
\sa `TriangulationTraits_2`
\sa `ConstrainedTriangulationTraits_2`
\sa `ConstrainedTriangulationFaceBase_2`

\cgalHeading{Implementation}

The insertion of a constrained edge runs in time
proportional to the number of triangles intersected by this edge.

*/
template< typename Traits, typename Tds, typename Itag >
class Constrained_triangulation_2 : public Triangulation_2<Traits,Tds> {
public:

/// \name Types
/// @{


/*!
\deprecated The type of constraints.
*/
typedef std::pair<Point,Point> Constraint;

/*!
A bidirectional iterator to visit
all the edges `e` of the triangulation which are constrained.
The order of visit is undefined.
The value type of this iterator is `Edge`.
*/
typedef unspecified_type Constrained_edges_iterator;

/*!
A range type to iterate over the constrained edges.
*/
typedef Iterator_range<Constrained_edges_iterator> Constrained_edges;

/*!
The intersection tag which decides how
intersections between input constraints are dealt with.
*/
typedef Itag Intersection_tag;

/*!
Exception, derived from `std::exception`, thrown when the intersection tag is
`No_constraint_intersection_tag` or `No_constraint_intersection_requiring_constructions_tag`,
and the insertion of a point or of a constraint would break the conditions associated with the tag
(see Section \ref PkgTriangulation2IntersectionTags).
*/
class Intersection_of_constraints_exception;

/// @}

/// \name Creation
/// @{

/*!
%Default constructor.
*/
Constrained_triangulation_2();

/*!
Copy constructor: All faces and vertices
are duplicated and the constrained status of edges
is copied.
*/
Constrained_triangulation_2(const
Constrained_triangulation_2& ct1);

/// @}

/// \name Queries
/// @{

/*!
returns `true` if edge `e` is a constrained edge.
*/
bool is_constrained(Edge e) const;

/*!
returns `true` if at least one of the edges incident to vertex `v`
is constrained.
*/
bool are_there_incident_constraints(Vertex_handle v) const;

/*!
outputs the constrained edges incident to `v`
into the output iterator `out` and returns the resulting
output iterator.
\tparam OutputItEdges is an `OutputIterator` with `Edge` as value
type.
*/
template<class OutputItEdges>
OutputItEdges incident_constraints(Vertex_handle v,
OutputItEdges out) const;

/*!
returns an iterator that enumerates the constrained edges.
*/
Constrained_edges_iterator constrained_edges_begin() const;

/*!
returns the past-the-end iterator.
*/
Constrained_edges_iterator constrained_edges_end() const;

/*!
returns a range of constrained edges.
*/
Constrained_edges constrained_edges() const;

/// @}

/// \name Insertion and Removal
///
/// @{

/*!
inserts point `p` and restores the status (constrained or not) of all
the touched edges: if `p` lies in the interior of a constrained edge,
this edge is split into two constrained edges. If present, `f` is used as a hint
for the location of `p`.

If the intersection tag is `No_constraint_intersection_tag` and `p` lies in the interior
of a constrained edge, an exception of type `Intersection_of_constraints_exception` is thrown.
*/
Vertex_handle insert(Point p, Face_handle f = Face_handle() );

/*!
inserts point `p` in the triangulation at the location given by `(lt,loc,i)`.
Same as above for constrained edges.
\sa `Triangulation_2::locate()`
*/
Vertex_handle
insert(const Point& p,
Locate_type& lt,
Face_handle loc, int li );

/*!
Equivalent to `insert(p)`.
*/
Vertex_handle push_back(const Point& p);

/*!
inserts points `a` and `b` in this order, and inserts the line segment `ab` as a
constraint. Removes the faces crossed by segment `ab` and creates new
faces instead. If a vertex `c` lies in the interior of segment `ab`, constraint `ab` is
replaced by the two constraints `ac` and `cb`. Apart from the insertion of
`a` and `b`, the algorithm runs in time proportional to the number of
removed triangles.
*/
void insert_constraint(Point a, Point b);

/*!
Equivalent to `insert_constraint(c.first, c.second)`.
*/
void push_back(const std::pair<Point,Point>& c);

/*!
inserts the line segment `s` whose endpoints are the vertices `va` and `vb` as a constraint.
The triangles intersected by `s` are removed and new ones are created.
*/
void insert_constraint(const Vertex_handle & va, const Vertex_handle & vb);

/*!
inserts a polyline defined by the points in the range `[first,last)`.
The polyline is considered as a polygon if the first and last point are equal or if  `close = true`. This enables for example passing the vertex range of a `Polygon_2`.
All the points of the polyline are first inserted, and then each segment of the polyline is inserted as a constraint.
\tparam PointIterator must be an `InputIterator` with the value type `Point`.
*/
template < class PointIterator>
void insert_constraint(PointIterator first, PointIterator last, bool close=false);

/*!
removes a vertex `v`.
\pre Vertex `v` is not incident to a constrained edge.
*/
void remove(Vertex_handle v);

/*!
makes the edges incident to vertex `v` unconstrained edges.
*/
void remove_incident_constraints(Vertex_handle v);

/*!
makes edge `(f,i)` unconstrained.
*/
void remove_constrained_edge(Face_handle f, int i);

/*!
checks the validity of the triangulation and the consistency
of the constrained marks in edges.
*/
bool
is_valid(bool verbose = false, int level = 0) const;

/// @}

}; /* end Constrained_triangulation_2 */

/*!
writes the triangulation as for `Triangulation_2<Traits,Tds>` and, for each face `f`, and integers `i=0,1,2`,
writes "C" or "N" depending whether edge
`(f,i)` is constrained or not.
\relates Constrained_triangulation_2
*/
template <typename  Traits, typename Tds, typename Itag>
std::ostream & operator<<(std::ostream& os, const Constrained_triangulation_2<Traits,Tds,Itag> &ct);

/*!
reads a triangulation from stream `is` and assigns it to c`t`. Data in the stream must have the same format `operator<<` uses.
Note that `ct` is first cleared.
\relates Constrained_triangulation_2
*/
template <typename  Traits, typename Tds, typename Itag>
std::istream& operator>>(std::istream& is,Constrained_triangulation_2<Traits,Tds,Itag> Ct& ct);

} /* end namespace CGAL */
