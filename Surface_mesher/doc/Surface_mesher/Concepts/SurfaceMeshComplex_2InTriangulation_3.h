
/*!
\ingroup PkgSurfaceMesher3Concepts
\cgalConcept

The concept `SurfaceMeshComplex_2InTriangulation_3` describes a data structure
designed to represent a two dimensional pure complex
embedded in a three dimensional triangulation.

A <I>complex</I> is a set \f$ C\f$ of faces such that:

- any subface of a face in \f$ C\f$ is a face of \f$ C\f$

- two faces of \f$ C\f$ are disjoint or share a common subface

The complex is <I>two dimensional</I>, if its faces have dimension at most
two. It is <I>pure</I> if any face in the complex is a subface
of some face of maximal dimension.
Thus, a two dimensional pure complex is a set of facets
together with their edges and vertices.
A two dimensional pure complex embedded
in a three dimensional triangulation
is a subset of the facets
of this triangulation, together with their edges and vertices.

The concept `SurfaceMeshComplex_2InTriangulation_3` is particularly suited to handle
surface meshes obtained as the restriction to a surface of
a three dimensional Delaunay triangulation.
A model of this concept is a type to be plugged as first template
parameter in the
function template `CGAL::make_surface_mesh()`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr>}
\cgalHasModelsEnd

\sa `CGAL::make_surface_mesh()`

*/

class SurfaceMeshComplex_2InTriangulation_3 {
public:

/// \name Types
/// `SurfaceMeshComplex_2InTriangulation_3` provides the following types.
/// @{

/*!
The type of the
embedding 3D triangulation. Must be a model of
`SurfaceMeshTriangulation_3`.
*/
typedef unspecified_type Triangulation;

/*!
The type of
the embedding triangulation vertex handles.
*/
typedef Triangulation::Vertex_handle Vertex_handle;

/*!
The type of
the embedding triangulation cell handles.
*/
typedef Triangulation::Cell_handle Cell_handle;

/*!
The type of
the embedding triangulation facets.
*/
typedef Triangulation::Facet Facet;

/*!
The type of
the embedding triangulation edges.
*/
typedef Triangulation::Edge Edge;

/*!
Size type (an
unsigned integral type)
*/
typedef Triangulation::size_type size_type;

/*!
A type to describe the status of a face (facet, edge, or vertex) with respect to
the 2D pure complex. A `NOT_IN_COMPLEX` face does not belong to the 2D complex.
Facets can only be `NOT_IN_COMPLEX` or `REGULAR`
depending on whether they belong to the 2D complex on not.
Edges and vertices can be `NOT_IN_COMPLEX`, `BOUNDARY`,
`REGULAR` or `SINGULAR`.
An edge in the complex is
`BOUNDARY`,
`REGULAR`, or `SINGULAR`,
if it is incident to respectively 1, 2, or 3 or
more facets in the complex.
The status of a vertex is determined by
the adjacency graph of the facets of the 2D complex
incident to that vertex.
The vertex of the 2D complex is `BOUNDARY`, if this adjacency graph
is a simple path, it is `REGULAR`, if the adjacency graph is cyclic,
and `SINGULAR` in any other case.

*/
enum Face_status {NOT_IN_COMPLEX, BOUNDARY, REGULAR,
SINGULAR};

/*!
An iterator type to visit the facets
of the 2D complex.
*/
typedef unspecified_type Facet_iterator;

/*!
An iterator type to visit the
edges of the 2D complex.
*/
typedef unspecified_type Edge_iterator;

/*!
An iterator type to visit
vertices of the 2D complex.
*/
typedef unspecified_type Vertex_iterator;

/*!
An iterator type to visit the
boundary edges of the 2D complex.
*/
typedef unspecified_type Boundary_edges_iterator;

/// @}

/// \name Creation
/// @{

/*!
Builds an empty 2D complex embedded in the triangulation `t3`
*/
SurfaceMeshComplex_2InTriangulation_3(Triangulation& t3);

// /*!
// Builds a 2D complex embedded in the triangulation `t3`,
// including in the 2D complex the facets of `t3` for
// which the predicate `select` returns `true`.

// The type `FacetSelector` must be
// a function object with an operator to select facets:
// `bool operator()(Facet f);`.
// */
// template < class FacetSelector>
// SurfaceMeshComplex_2InTriangulation_3(Triangulation& t3,
// FacetSelector select);

/// @}

/// \name Member access
/// @{

/*!
Returns the reference to the triangulation.
*/
Triangulation& triangulation();

/// @}

/// \name Modifications
/// @{

/*!
Adds facet `f` to the 2D complex.
*/
void add_to_complex(Facet f);

/*!
Adds facet `(c,i)` to the 2D complex.
*/
void add_to_complex(Cell_handle c, int i);

/*!
Removes facet `f` from the 2D complex.
*/
void remove_from_complex(Facet f);

/*!
Removes facet `(c,i)` from the 2D complex.
*/
void remove_from_complex(Cell_handle c, int i);

/// @}

/// \name Queries
/// Queries on the status of individual face with respect to the 2D complex.
/// @{

/*!
Returns the number of facets that belong to the 2D complex.
*/
size_type number_of_facets() const;

/*!
Returns the status of the facet `f` with respect to the 2D complex.
*/
Face_status face_status(Facet f);

/*!
Returns the status of the facet `(c,i)` with respect to the 2D complex.
*/
Face_status face_status(Cell_handle c, int i);

/*!
Returns the status of edge `e` in the 2D complex.
*/
Face_status face_status(Edge e);

/*!
Returns the status of edge `(c,i,j)` in the 2D complex.
*/
Face_status face_status(Cell_handle c, int
i, int j);

/*!
Returns the status of vertex `v` in the 2D complex.
*/
Face_status face_status(Vertex_handle v);

/*!
Returns `true`, if the facet `f` belongs to the 2D complex.
*/
bool is_in_complex(Facet f);

/*!
Returns `true`, if the facet `(c,i)` belongs to the 2D complex.
*/
bool is_in_complex(Cell_handle c, int i);

/*!
Returns `true`, if the edge `e` belongs to the 2D complex.
*/
bool is_in_complex(Edge e);

/*!
Returns `true`, if the edge `(c,i,j)` belongs to the 2D complex.
*/
bool is_in_complex(Cell_handle c, int i, int j);

/*!
Returns `true`, if the vertex `v` belongs to the 2D complex.
*/
bool is_in_complex(Vertex_handle v);

/*!
Returns true if the status of vertex `v` is `REGULAR` or `BOUNDARY`.
\pre All the edges of the complex incident to `v` are `REGULAR` or `BOUNDARY`.
*/
bool is_regular_or_boundary_for_vertices (Vertex_handle v);

/// @}

/// \name Traversal of the complex
/// The data structure provides iterators to visit the facets, edges
/// and vertices of the complex. All those iterators are bidirectional
/// and non mutable.
/// @{

/*!
Returns an iterator with value type `Facet` to visit the facets
of the 2D complex.
*/
Facet_iterator facets_begin();

/*!
Returns the past the end iterator for the above iterator.
*/
Facet_iterator facets_end();

/*!
Returns an iterator with value type `Edge` to visit the
edges of the 2D complex which are not isolated.
*/
Edge_iterator edges_begin();

/*!
Returns the past the end iterator for the above iterator.
*/
Edge_iterator edges_end();

/*!
Returns an iterator with value type `Edge` to visit the
boundary edges of the complex.
*/
Boundary_edges_iterator boundary_edges_begin();

/*!
Returns the past the end iterator for the above iterator.
*/
Boundary_edges_iterator boundary_edges_end();

/*!
Returns an iterator with value type `Vertex_handle` to visit the
vertices of the 2D complex.
*/
Vertex_iterator vertices_begin();

/*!
Returns the past the end iterator for the above iterator.
*/
Vertex_iterator vertices_end();

/*!
Copies the facets of the complex incident to `v` to the output
iterator `facets`.
Returns the resulting output iterator.
\pre `c2t3.triangulation().dimension() == 3`, `v != Vertex_handle()`, `c2t3.triangulation().is_vertex(v)`.
*/
template <class OutputIterator>
OutputIterator
incident_facets(Vertex_handle v, OutputIterator facets);

/// @}

/// \name
/// The following function is the basic function to walk on the 2D complex
/// @{

/*!
Returns the facet of the complex which is the neighbor of
the facet `f` opposite to the vertex with index `j` of
`f`.
The vertices of the facet `f = (cell c, i)` are numbered
(0,1,2) (according to the `vertex_triple_index(i,j)` member function
of `Triangulation_3`)
in such a way that facet `f` is oriented by the
outward normal of tetraedra `c`.
If there is no such neighbor, or if the edge is singular the functions returns `Facet()`.
*/
Facet neighbor(Facet f, int j);

/*!
Returns the facet of the complex which is the neighbor of
the facet `f` opposite to the vertex with index `j` of `f`.
See above.
*/
Facet neighbor(Cell_handle c, int i, int j);

/// @}

}; /* end SurfaceMeshComplex_2InTriangulation_3 */

