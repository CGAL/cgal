namespace CGAL {

/*!
\ingroup PkgHalfedgeDS_Items

The class `HalfedgeDS_items_2` is a model of the `HalfedgeDSItems` concept.
It uses the default types for vertices, halfedges, and faces that
declare all incidences supported by a `HalfedgeDS`. The vertex
also contains a point of type `Traits::Point_2`, where `Traits`
is the template argument of the corresponding `HalfedgeDS`.

\cgalModels `HalfedgeDSItems`

\sa `CGAL::HalfedgeDS_min_items`
\sa `CGAL::Polyhedron_items_3`
\sa `HalfedgeDS<Traits,Items,Alloc>`
\sa `PolyhedronItems_3`
\sa `CGAL::HalfedgeDS_vertex_base<Refs>`
\sa `CGAL::HalfedgeDS_halfedge_base<Refs>`
\sa `CGAL::HalfedgeDS_face_base<Refs>`

\cgalHeading{Example}

The following example shows the canonical implementation of the
`HalfedgeDS_items_2` class. It uses the base classes for the item types that
are provided in the library.

\code{.cpp}

struct HalfedgeDS_items_2 {
template < class Refs, class Traits>
struct Vertex_wrapper {
typedef typename Traits::Point_2 Point;
typedef CGAL::HalfedgeDS_vertex_base< Refs, Tag_true, Point> Vertex;
};
template < class Refs, class Traits>
struct Halfedge_wrapper {
typedef CGAL::HalfedgeDS_halfedge_base< Refs> Halfedge;
};
template < class Refs, class Traits>
struct Face_wrapper {
typedef CGAL::HalfedgeDS_face_base< Refs> Face;
};
};

\endcode

The following example shows a class definition for a new items class
derived from the `HalfedgeDS_items_2` class. It replaces the `Face_wrapper`
with a new definition of a face that contains a member variable for
color. The new face makes use of the face base class provided in the
library.

\code{.cpp}

// A face type with a color member variable.
template <class Refs>
struct My_face : public CGAL::HalfedgeDS_face_base<Refs> {
CGAL::IO::Color color;
My_face() {}
My_face( CGAL::IO::Color c) : color(c) {}
};

// An items type using my face.
struct My_items : public CGAL::HalfedgeDS_items_2 {
template <class Refs, class Traits>
struct Face_wrapper {
typedef My_face<Refs> Face;
};
};

\endcode

*/

class HalfedgeDS_items_2 {
public:

}; /* end HalfedgeDS_items_2 */
} /* end namespace CGAL */
