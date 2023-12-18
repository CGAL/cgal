namespace CGAL {

/*!
\ingroup PkgHalfedgeDS_Items

The class `HalfedgeDS_min_items` is a model of the `HalfedgeDSItems`
concept. It defines types for vertices, halfedges, and faces that
declare the minimal required incidences for a `HalfedgeDS`, which
are the `next()` and the `opposite()` member function for
halfedges.

\cgalModels{HalfedgeDSItems}

\sa `CGAL::HalfedgeDS_items_2`
\sa `CGAL::Polyhedron_items_3`
\sa `HalfedgeDS<Traits,Items,Alloc>`
\sa `PolyhedronItems_3`
\sa `CGAL::HalfedgeDS_vertex_min_base<Refs>`
\sa `CGAL::HalfedgeDS_halfedge_min_base<Refs>`
\sa `CGAL::HalfedgeDS_face_min_base<Refs>`

\cgalHeading{Example}

The following example shows the canonical implementation of the
`CGAL::HalfedgeDS_min_items` class. It uses the base classes for the
item types that are provided in the library.

\code{.cpp}

struct HalfedgeDS_min_items {
template < class Refs, class Traits>
struct Vertex_wrapper {
typedef CGAL::HalfedgeDS_vertex_min_base< Refs> Vertex;
};
template < class Refs, class Traits>
struct Halfedge_wrapper {
typedef CGAL::HalfedgeDS_halfedge_min_base< Refs> Halfedge;
};
template < class Refs, class Traits>
struct Face_wrapper {
typedef CGAL::HalfedgeDS_face_min_base< Refs> Face;
};
};

\endcode

*/

class HalfedgeDS_min_items {
public:

}; /* end HalfedgeDS_min_items */
} /* end namespace CGAL */
