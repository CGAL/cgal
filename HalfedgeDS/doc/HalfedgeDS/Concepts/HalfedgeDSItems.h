/*!
\ingroup PkgHDSConcepts
\cgalConcept

The concept `HalfedgeDSItems` wraps the three item types - vertex, 
halfedge, and face - for a halfedge data structure. A `HalfedgeDSItems` 
contains three member class templates named `Vertex_wrapper`, 
`Halfedge_wrapper`, and `Face_wrapper`, each with two template 
parameters, `Refs` and `Traits`. `Refs` requires an 
instantiated halfedge data structure `HalfedgeDS` as argument, 
`Traits` is a geometric traits class supplied by the class that 
uses the halfedge data structure as internal representation. `Traits` 
is not used by the halfedge data structure itself. These three member 
class templates provide a local type named `Vertex`, `Halfedge`, 
and `Face` respectively. The requirements on these 
types are described on the manual pages of the concepts `HalfedgeDSVertex`, 
`HalfedgeDSHalfedge`, and `HalfedgeDSFace` respectively. 

\cgalHasModel CGAL::HalfedgeDS_min_items 
\cgalHasModel CGAL::HalfedgeDS_items_2 
\cgalHasModel CGAL::Polyhedron_items_3 

\sa `HalfedgeDS<Traits,Items,Alloc>` 
\sa `HalfedgeDSVertex` 
\sa `HalfedgeDSHalfedge` 
\sa `HalfedgeDSFace` 
\sa `PolyhedronItems_3` 
\sa `CGAL::HalfedgeDS_vertex_base<Refs>` 
\sa `CGAL::HalfedgeDS_halfedge_base<Refs>` 
\sa `CGAL::HalfedgeDS_face_base<Refs>` 

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

See `CGAL::HalfedgeDS_items_2` for an example implementation.

*/

class HalfedgeDSItems {
public:

/// \name Types 
/// @{

/*!
model of `HalfedgeDSVertex`. 
*/ 
typedef Vertex_wrapper<Refs,Traits>::Vertex Vertex;

/*!
model of `HalfedgeDSHalfedge`. 
*/ 
typedef Halfedge_wrapper<Refs,Traits>::Halfedge Halfedge;

/*!
model of `HalfedgeDSFace`. 
*/ 
typedef Face_wrapper<Refs,Traits>::Face Face;

/// @}

}; /* end HalfedgeDSItems */
