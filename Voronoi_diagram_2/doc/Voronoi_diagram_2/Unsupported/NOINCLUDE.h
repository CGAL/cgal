CONVINFO Missing include from ../doc_html/Voronoi_diagram_2/Voronoi_diagram_2_ref/Face.h

namespace CGAL {
CONVERROR Additional namespace Voronoi_diagram_2<DG,AT,AP>:: required
/*!
\ingroup PkgVoronoiDiagramAdaptor2

The class `Face` is the class provided by the 
`Voronoi_diagram_2<DG,AT,AP>` class for Voronoi faces. Below we 
present its interface. 

\models `DefaultConstructible`
\models `CopyConstructible`
\models `Assignable`
\models `EqualityComparable`
\models `LessThanComparable` 

\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>` 
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>::Vertex` 
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>::Halfedge` 
\sa `DelaunayGraph_2` 

*/

class Face {
public:

/// \name Types 
/// @{

/*! 
A type for the vertices of the Voronoi diagram. 
*/ 
typedef Hidden_type Vertex; 

/*! 
A type for the halfedges of the Voronoi diagram. 
*/ 
typedef Hidden_type Halfedge; 

/*! 
Handle for the vertices of the Voronoi diagram. 
*/ 
typedef Hidden_type Vertex_handle; 

/*! 
Handle for the faces of the Voronoi diagram. 
*/ 
typedef Hidden_type Face_handle; 

/*! 
Handle for the halfedges of the Voronoi 
diagram. 
*/ 
typedef Hidden_type Halfedge_handle; 

/*! 
A type for a bidirectional 
circulator over the halfedges on the boundary of the face. The value 
type of the circulator is 
`CGAL::Voronoi_diagram_2<DG,AT,AP>::Halfedge`, and is convertible 
to `Halfedge_handle`. 
*/ 
typedef Hidden_type Ccb_halfedge_circulator; 

/*! 
A type for the Delaunay graph. It is a 
model of the `DelaunayGraph_2` concept. 
*/ 
typedef Hidden_type Delaunay_graph; 

/*! 
A type for the handle of the dual vertex. 
*/ 
typedef Delaunay_graph::Vertex_handle 
Delaunay_vertex_handle; 

/// @} 

/// \name Access Methods 
/// @{

/*! 
Returns an incident halfedge on 
the boundary of `f`. 
*/ 
Halfedge_handle halfedge(); 

/*! 
Returns a bidirectional 
circulator for traversing the halfedges on the boundary of 
`f`. The halfedges are traversed in counterclockwise order. 
*/ 
Ccb_halfedge_circulator ccb(); 

/*! 
Returns a handle to the corresponding dual vertex in the Delaunay graph. 
*/ 
Delaunay_vertex_handle dual(); 

/// @} 

/// \name Predicate Methods 
/// @{

/*! 
Returns `true` iff the face is 
an unbounded face in the Voronoi diagram. 
*/ 
bool is_unbounded(); 

/*! 
Returns 
`true` iff `e` is a halfedge of the boundary of 
`f`. 
*/ 
bool is_halfedge_on_ccb(Halfedge e); 

/*! 
Returns `true` iff the following 
conditions are met: the face is not rejected by the chosen 
adaptation policy; 
all its adjacent halfedges do not have zero length; all its adjacent 
halfedges return the face as their adjacent face. 
*/ 
bool is_valid(); 

/// @}

}; /* end Face */
} /* end namespace CGAL */
CONVINFO Missing include from ../doc_html/Voronoi_diagram_2/Voronoi_diagram_2_ref/Halfedge.h

namespace CGAL {
CONVERROR Additional namespace Voronoi_diagram_2<DG,AT,AP>:: required
/*!
\ingroup PkgVoronoiDiagramAdaptor2

The class `Halfedge` is the class provided by the 
`Voronoi_diagram_2<DG,AT,AP>` class for Voronoi halfedges. 
Below we present its interface. 

\models `DefaultConstructible`
\models `CopyConstructible`
\models `Assignable`
\models `EqualityComparable`
\models `LessThanComparable`

\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>` 
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>::Vertex` 
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>::Face` 
\sa `DelaunayGraph_2` 

*/

class Halfedge {
public:

/// \name Types 
/// @{

/*! 
A type for the vertices of the Voronoi diagram. 
*/ 
typedef Hidden_type Vertex; 

/*! 
A type for the faces of the Voronoi diagram. 
*/ 
typedef Hidden_type Face; 

/*! 
Handle for the vertices of the Voronoi diagram. 
*/ 
typedef Hidden_type Vertex_handle; 

/*! 
Handle for the faces of the Voronoi diagram. 
*/ 
typedef Hidden_type Face_handle; 

/*! 
Handle for the halfedges of the Voronoi 
diagram. 
*/ 
typedef Hidden_type Halfedge_handle; 

/*! 
A type for a bidirectional 
circulator over the halfedges of the boundary of a 
Voronoi face. The value type of the circulator is 
`CGAL::Voronoi_diagram_2<DG,AT,AP>::Halfedge` and is convertible to 
`Halfedge_handle`. 
*/ 
typedef Hidden_type Ccb_halfedge_circulator; 

/*! 
A type for the Delaunay graph. It is a 
model of the `DelaunayGraph_2` concept. 
*/ 
typedef Hidden_type Delaunay_graph; 

/*! 
A type for 
the dual edge in the Delaunay graph. 
*/ 
typedef Delaunay_graph::Edge Delaunay_edge; 

/*! 
A type for vertex handles in the Delaunay graph. 
*/ 
typedef Delaunay_graph::Vertex_handle Delaunay_vertex_handle; 

/// @} 

/// \name Access Methods 
/// @{

/*! 
Returns the twin halfedge. 
*/ 
Halfedge_handle twin(); 

/*! 
Same as `e.twin()`. 
*/ 
Halfedge_handle opposite(); 

/*! 
Returns the next halfedge in the 
counterclockwise sense around the boundary of the face that `e` 
is incident to. 
*/ 
Halfedge_handle next(); 

/*! 
Returns the previous halfedge in the 
counterclockwise sense around the boundary of the adjacent face. 
*/ 
Halfedge_handle previous(); 

/*! 
Returns the face that `e` is 
incident to. 
*/ 
Face_handle face(); 

/*! 
Returns the source vertex of 
`e`. 
\pre The source vertex must exist, i.e., `has_source()` must return `true`. 
*/ 
Vertex_handle source(); 

/*! 
Returns the target vertex of 
`e`. 
\pre The target vertex must exist, i.e., `has_target()` must return `true`. 
*/ 
Vertex_handle target(); 

/*! 
Returns a bidirectional 
circulator to traverse the halfedges on the boundary of the Voronoi 
face containing `e`. The circulator is initialized to 
`e`. Applying `operator++` (resp. `operator-`) to this 
circulator returns the next halfedge on the boundary of the face 
containing `e` in the counterclockwise (resp. clockwise) sense. 
*/ 
Ccb_halfedge_circulator ccb(); 

/*! 
Returns the 
corresponding dual edge in the Delaunay graph. 
*/ 
Delaunay_edge dual(); 

/// @}

/// \name
/// In the four methods below we consider Voronoi halfedges to be
/// "parallel" to the \f$ x\f$-axis, oriented from left to right.
/// @{

/*! 
Returns a handle to the vertex in 
the Delaunay graph corresponding to the defining site above 
the Voronoi edge. 
*/ 
Delaunay_vertex_handle up(); 

/*! 
Returns a handle to the vertex 
in the Delaunay graph corresponding to the defining site below 
the Voronoi edge. 
*/ 
Delaunay_vertex_handle down(); 

/*! 
Returns a handle to the vertex in 
the Delaunay graph corresponding to the defining site to the left of 
the Voronoi edge. 
\pre `has_source()` must be `true`. 
*/ 
Delaunay_vertex_handle left(); 

/*! 
Returns a handle to the vertex in 
the Delaunay graph corresponding to the defining site to the right of 
the Voronoi edge. 
\pre `has_target()` must be `true`. 
*/ 
Delaunay_vertex_handle right(); 

/// @} 

/// \name Predicate Methods 
/// @{

/*! 
Returns `true` iff the halfedge 
corresponds to a bisecting segment or a bisecting ray oriented 
appropriately so that its apex is its source. 
*/ 
bool has_source(); 

/*! 
Returns `true` iff the halfedge 
corresponds to a bisecting segment or a bisecting ray oriented 
appropriately so that its apex is its target. 
*/ 
bool has_target(); 

/*! 
Returns `true` iff the source or 
the target of the halfedge does not exist, i.e., if either of 
`has_source()` or `has_target()` return `false`. 
*/ 
bool is_unbounded(); 

/*! 
Returns `true` iff the Voronoi 
edge is an entire bisector. 
*/ 
bool is_bisector(); 

/*! 
Returns `true` iff the Voronoi 
edge has both a source and a target Voronoi vertex. 
*/ 
bool is_segment(); 

/*! 
Returns `true` iff the Voronoi 
edge has either a source or a target Voronoi vertex, but not both; 
in other words it is a bisecting ray. 
*/ 
bool is_ray(); 

/*! 
Returns `true` if the following 
conditions are met: the halfedge is not a rejected 
edge with respect to the chosen adaptation policy; 
the twin edge of its twin edge is itself; its adjacent face is not a 
rejected face with respect to the chosen adaptation policy; 
its source and target vertices are valid (provided 
they exist, of course); the previous of its next halfedge is itself 
and the next of its previous halfedge is itself. 
*/ 
bool is_valid(); 

/// @}

}; /* end Halfedge */
} /* end namespace CGAL */
CONVINFO Missing include from ../doc_html/Voronoi_diagram_2/Voronoi_diagram_2_ref/Vertex.h

namespace CGAL {
CONVERROR Additional namespace Voronoi_diagram_2<DG,AT,AP>:: required
/*!
\ingroup PkgVoronoiDiagramAdaptor2

The class `Vertex` is the Voronoi vertex class provided by the class 
`Voronoi_diagram_2<DG,AT,AP>` class. Below we present its interface. 

\models `DefaultConstructible`
\models `CopyConstructible`
\models `Assignable`, 
\models `EqualityComparable`
\models `LessThanComparable`

\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>` 
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>::Halfedge` 
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>::Face` 
\sa `DelaunayGraph_2` 

*/

class Vertex {
public:

/// \name Types 
/// @{

/*! 
A type for the halfedges of the Voronoi diagram. 
*/ 
typedef Hidden_type Halfedge; 

/*! 
A type for the faces of the Voronoi diagram. 
*/ 
typedef Hidden_type Face; 

/*! 
Handle for the vertices of the Voronoi diagram. 
*/ 
typedef Hidden_type Vertex_handle; 

/*! 
Handle for the faces of the Voronoi diagram. 
*/ 
typedef Hidden_type Face_handle; 

/*! 
Handle for the halfedges of the Voronoi 
diagram. 
*/ 
typedef Hidden_type Halfedge_handle; 

/*! 
A type for the point represented by the 
vertex. 
*/ 
typedef Hidden_type Point_2; 

/*! 
A type for sizes. 
*/ 
typedef Hidden_type size_type; 

/*! 
A type for a bidirectional 
circulator that allows to traverse all incident halfedges, i.e., all 
halfedges that have the vertex as their target. The value 
type of the circulator is 
`CGAL::Voronoi_diagram_2<DG,AT,AP>::Halfedge` and is convertible to 
`Halfedge_handle`. 
*/ 
typedef Hidden_type Halfedge_around_vertex_circulator; 

/*! 
A type for the Delaunay graph. It is a 
model of the `DelaunayGraph_2` concept. 
*/ 
typedef Hidden_type Delaunay_graph; 

/*! 
A type for the handle of the dual face. 
*/ 
typedef Delaunay_graph::Face_handle Delaunay_face_handle; 

/*! 
A type for the vertex handles in the Delaunay graph. 
*/ 
typedef Delaunay_graph::Vertex_handle Delaunay_vertex_handle; 

/// @} 

/// \name Access Methods 
/// @{

/*! 
Returns an incident halfedge 
that has `v` as its target. 
*/ 
Halfedge_handle halfedge(); 

/*! 
Returns the in-degree of the vertex, 
i.e. the number of halfedges that have `v` as their target. 
*/ 
size_type degree(); 

/*! 
Returns the point represented by the vertex. 
*/ 
Point_2 point(); 

/*! 
Returns a handle to the corresponding dual face in the 
Delaunay graph. 
*/ 
Delaunay_face_handle dual(); 

/*! 
Returns a handle to the vertex in the Delaunay graph corresponding to 
the \f$ (i+1)\f$-th generating site of the Voronoi vertex. 
\pre `i` must be smaller or equal to 2. 
*/ 
Delaunay_vertex_handle site(unsigned int i); 

/*! 
Returns a bidirectional circulator that allows the traversal of the 
halfedges that have `v` as their target. Applying 
`operator++` (resp. `operator-`) to this circulator returns 
the next incident halfedge in the counterclockwise (resp. clockwise) sense. 
*/ 
Halfedge_around_vertex_circulator incident_halfedges(); 

/// @} 

/// \name Predicate Methods 
/// @{

/*! 
Returns `true` 
if the halfedge `e` is incident to `v`. 
*/ 
bool is_incident_edge(Halfedge_handle e); 

/*! 
Returns `true` 
if the face `f` is incident to `v`. 
*/ 
bool is_incident_face(Face_handle e); 

/*! 
Returns `true` if the following 
conditions are met: the dual face is not an infinite face; all 
incident halfedges have the vertex as their target. 
*/ 
bool is_valid(); 

/// @}

}; /* end Vertex */
} /* end namespace CGAL */
