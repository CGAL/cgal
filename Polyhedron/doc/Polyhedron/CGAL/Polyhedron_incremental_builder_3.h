
namespace CGAL {

/*!
\ingroup PkgPolyhedronRef

The auxiliary class `Polyhedron_incremental_builder_3` supports the incremental
construction of polyhedral surfaces, which is for example convenient
when constructing polyhedral surfaces from file formats, such as the
%Object File Format (OFF) \cgalCite{cgal:p-gmgv16-96},
OpenInventor \cgalCite{cgal:w-impoo-94} or
VRML \cgalCite{cgal:bpp-vrml-95}, \cgalCite{cgal:vrmls-97}.
`Polyhedron_incremental_builder_3` needs access to the internal halfedge data
structure of type `HDS` of the polyhedral surface. It is intended
to be used within a modifier, see `Modifier_base`.

The incremental builder might be of broader interest for other uses of
the halfedge data structures, but it is specifically bound to the
definition of polyhedral surfaces given here. During construction all
conditions of polyhedral surfaces are checked and in case of violation
an error status is set. A diagnostic message will be issued to
`std::cerr` if the `verbose` flag has been set at construction
time.

The incremental construction starts with a list of all point
coordinates and concludes with a list of all facet polygons. Edges are
not explicitly specified. They are derived from the vertex incidence
information provided from the facet polygons. The polygons are given as a
sequence of vertex indices. The halfedge data structure `HDS` must
support vertices (i.e., `Supports_halfedge_vertex` \f$ \equiv\f$
`Tag_true`). Vertices and facets can be added in arbitrary order
as long as a call to `add_vertex_to_facet()` refers only to a
vertex index that is already known. Some methods return already
handles to vertices, facets, and halfedges newly constructed. They can
be used to initialize additional fields, however, the incidences in
the halfedge-data structure are not stable and are not allowed to be
changed.

The incremental builder can work in two modes: `RELATIVE_INDEXING` (the
default), in which a polyhedral surface already contained in the
halfedge data structure is ignored and all indices are relative to the
newly added surface, or `ABSOLUTE_INDEXING`, in which all indices are
absolute indices including an already existing polyhedral surface. The
former mode allows to create easily independent connected components,
while the latter mode allows to to continue the construction of an
existing surface, the absolute indexing allows to address existing
vertices when creating new facets.

\sa `CGAL::Polyhedron_3<Traits>`
\sa `HalfedgeDS`
\sa `CGAL::Modifier_base`

\cgalHeading{Example}

A modifier class creates a new triangle in the halfedge data structure
using the incremental builder.

\cgalExample{Polyhedron/polyhedron_prog_incr_builder.cpp}

*/
template< typename HDS >
class Polyhedron_incremental_builder_3 {
public:

/// \name Types
/// @{

/*!
halfedge data structure `HDS`.
*/
typedef unspecified_type HalfedgeDS;

/*!
point type of the vertex.
*/
typedef unspecified_type Point_3;

/*!
size type.
*/
typedef unspecified_type size_type;

/*!

*/
typedef typename HalfedgeDS::Vertex_handle Vertex_handle;

/*!

*/
typedef typename HalfedgeDS::Halfedge_handle Halfedge_handle;

/*!

*/
typedef typename HalfedgeDS::Face_handle Facet_handle;

/// @}

/// \name Constants
/// @{

/*!
two different indexing modes.
*/
enum { RELATIVE_INDEXING, ABSOLUTE_INDEXING};

/// @}

/// \name Creation
/// @{

/*!
stores a reference to the halfedge data structure `hds` of a
polyhedral surface in its internal state. An existing polyhedral
surface in `hds` remains unchanged. The incremental builder
appends the new polyhedral surface. If `verbose` is `true`,
diagnostic messages will be printed to `cerr` in case of
malformed input data.
*/
Polyhedron_incremental_builder_3(HDS& hds,
bool verbose = false);

/// @}

/*!
\name Surface Creation

To build a polyhedral surface, the following regular expression gives
the correct and allowed order and nesting of method calls from this
section:

\code
begin_surface ( add_vertex  | ( begin_facet  add_vertex_to_facet  end_facet ) ) end_surface
\endcode
*/
/// @{

/*!
starts the construction. `v` is the number of new vertices
to expect, `f` the number of new facets, and `h` the number of
new halfedges. If `h` is unspecified (`== 0`) it is estimated using
Euler's equation (plus 5% for the so far unknown holes and genus of
the object). These values are used to reserve space in the
halfedge data structure `hds`. If the representation supports
insertion these values do not restrict the class of constructible
polyhedra. If the representation does not support insertion the
object must fit into the reserved sizes.

If `mode` is set to `ABSOLUTE_INDEXING` the incremental builder
uses absolute indexing and the vertices of the old polyhedral surface
can be used in new facets (needs preprocessing time linear in the
size of the old surface). Otherwise relative indexing is used
starting with new indices for the new construction.
*/
void begin_surface( size_type v, size_type f, size_type h = 0,
int mode = RELATIVE_INDEXING);

/*!

adds a new vertex for `p` and returns its handle.
*/
Vertex_handle add_vertex( const Point_3& p);

/*!

starts a new facet and returns its handle.
*/
Facet_handle begin_facet();

/*!
adds a vertex with
index  `i` to the current facet. The first point added with
`add_vertex()` has the index 0 if `mode` was set to
`RELATIVE_INDEXING`, otherwise the first vertex in the
referenced `hds` has the index 0.
*/
void add_vertex_to_facet( size_type i);

/*!
ends a newly constructed facet.
Returns the handle to the halfedge incident to the new facet that points
to the vertex added first. The halfedge can be safely used to traverse
the halfedge cycle around the new facet.
*/
Halfedge_handle end_facet();

/*!
ends the construction.
*/
void end_surface();

/// @}

/// \name Additional Operations
/// @{

/*!

is a synonym for `begin_facet()`, a call to `add_vertex_to_facet()` for each
value in the range `[first,beyond)`, and a call to `end_facet()`.
Returns the return value of `end_facet()`.
\pre The value type of `InputIterator` is `std::size_t`. All indices must refer to vertices already added.
*/
template <class InputIterator>
Halfedge_handle add_facet( InputIterator first, InputIterator beyond);

/*!

returns `true` if a facet described by the vertex indices in the range
`[first,beyond)` can be successfully inserted, e.g., with
`add_facet(first,beyond)`.
\pre The value type of `InputIterator` is `std::size_t`. All indices must refer to vertices already added.
*/
template <class InputIterator>
bool test_facet( InputIterator first, InputIterator beyond);

/*!

returns handle for the vertex of index `i`, or `Vertex_handle` if
there is no `i`-th vertex.
*/
Vertex_handle vertex( std::size_t i);

/*!
returns error status of the builder.
*/
bool error() const;

/*!
undoes all changes made to the halfedge
data structure since the last `begin_surface()` in relative
indexing, and deletes the whole surface in absolute indexing.
It needs a new call to `begin_surface()` to start inserting again.
*/
void rollback();

/*!
returns
`true` if unconnected vertices are detected. If `verbose` was set to
`true` (see the constructor above) debug information about the
unconnected vertices is printed.
*/
bool check_unconnected_vertices();

/*!
returns
`true` if all unconnected vertices could be removed successfully.
This happens either if no unconnected vertices had appeared or if the
halfedge data structure supports the removal of individual elements.
*/
bool remove_unconnected_vertices();

/// @}

}; /* end Polyhedron_incremental_builder_3 */
} /* end namespace CGAL */
