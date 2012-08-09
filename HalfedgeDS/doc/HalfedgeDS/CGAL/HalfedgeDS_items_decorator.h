namespace CGAL {

/*!
\ingroup PkgHDS

The classes `CGAL::HalfedgeDS_items_decorator<HDS>`, 
`CGAL::HalfedgeDS_decorator<HDS>`, and 
`CGAL::HalfedgeDS_const_decorator<HDS>` provide additional functions 
to examine and to modify a halfedge data structure `HDS`. The class 
`CGAL::HalfedgeDS_items_decorator<HDS>` provides additional functions 
for vertices, halfedges, and faces of a halfedge data structure 
without knowing the containing halfedge data structure. The class 
`CGAL::HalfedgeDS_decorator<HDS>` stores a reference to the halfedge 
data structure and provides functions that modify the halfedge data 
structure, for example Euler-operators. The class 
`CGAL::HalfedgeDS_const_decorator<HDS>` stores a const reference to 
the halfedge data structure. It contains non-modifying functions, for 
example the test for validness of the data structure. 

All these additional functions take care of the different capabilities 
a halfedge data structure may have or may not have. The functions 
evaluate the type tags of the halfedge data structure to decide on the 
actions. If a particular feature is not supported nothing is done. 
Note that for example the creation of new halfedges is mandatory for 
all halfedge data structures and will not appear here again. 
\sa `CGAL::HalfedgeDS_decorator<HDS>` 
\sa `CGAL::HalfedgeDS_const_decorator<HDS>` 

Example 
-------------- 

The following program fragment illustrates how a refined halfedge 
class for a polyhedron can make use of the `find_prev()` member 
function to implement a `prev()` member function that works 
regardless of whether the halfedge data structure `HDS` provides a 
`prev()` member function for its halfedges or not. In the case that not, 
the implementation given here runs in time proportional to the size of the 
incident face. For const-correctness a second implementation with signature 
`Halfedge_const_handle prev() const;` is needed. 

Note also the use of the static member function `halfedge_handle()` 
of the halfedge data structure. It converts a pointer to the halfedge 
into a halfedge handle. This conversion encapsulates possible 
adjustments for hidden data members in the true halfedge type, such as 
linked-list pointers. 

\code{.cpp} 

struct Polyhedron_halfedge { 
// ... 
Halfedge_handle prev() { 
CGAL::HalfedgeDS_items_decorator<HDS> decorator; 
return decorator.find_prev( HDS::halfedge_handle(this)); 
} 
}; 

\endcode 

*/
template< typename HDS >
class HalfedgeDS_items_decorator {
public:

/// \name Types 
/// The respective `const_handle`'s and `const_iterator`'s are 
/// available as well.
/// @{

/*! 
halfedge data structure. 
*/ 
typedef Hidden_type HalfedgeDS; 

/*! 
traits class. 
*/ 
typedef Hidden_type Traits; 

/*! 
vertex type of `HalfedgeDS`. 
*/ 
typedef Hidden_type Vertex; 

/*! 
halfedge type of `HalfedgeDS`. 
*/ 
typedef Hidden_type Halfedge; 

/*! 
face type of `HalfedgeDS`. 
*/ 
typedef Hidden_type Face; 

/*! 

*/ 
typedef Hidden_type Vertex_handle; 

/*! 

*/ 
typedef Hidden_type Halfedge_handle; 

/*! 

*/ 
typedef Hidden_type Face_handle; 

/*! 

*/ 
typedef Hidden_type Vertex_iterator; 

/*! 

*/ 
typedef Hidden_type Halfedge_iterator; 

/*! 

*/ 
typedef Hidden_type Face_iterator; 

/*! 

*/ 
typedef Hidden_type size_type; 

/*! 

*/ 
typedef Hidden_type difference_type; 

/*! 

*/ 
typedef Hidden_type iterator_category; 

/*! 

*/ 
typedef Hidden_type Supports_vertex_halfedge; 

/*! 

*/ 
typedef Hidden_type Supports_halfedge_prev; 

/*! 

*/ 
typedef Hidden_type Supports_halfedge_vertex; 

/*! 

*/ 
typedef Hidden_type Supports_halfedge_face; 

/*! 

*/ 
typedef Hidden_type Supports_face_halfedge; 

/*! 

*/ 
typedef Hidden_type Supports_removal; 

/// @} 

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
HalfedgeDS_items_decorator(); 

/// @} 

/// \name Access Functions 
/// Corresponding member functions for `const_handle`'s are provided as well. 
/// @{

/*! 
returns the incident halfedge of \f$ v\f$ if supported, 
`Halfedge_handle()` otherwise. 
*/ 
Halfedge_handle get_vertex_halfedge( Vertex_handle v); 

/*! 
returns the incident vertex of \f$ h\f$ if supported, `Vertex_handle()` 
otherwise. 
*/ 
Vertex_handle get_vertex(Halfedge_handle h); 

/*! 
returns the previous halfedge of \f$ h\f$ if supported, 
`Halfedge_handle()` otherwise. 
*/ 
Halfedge_handle get_prev(Halfedge_handle h); 

/*! 
returns the previous halfedge of \f$ h\f$. Uses the `prev()` method 
if supported or performs a search around the face using `next()`. 
*/ 
Halfedge_handle find_prev(Halfedge_handle h); 

/*! 
returns the previous halfedge of \f$ h\f$. Uses the `prev()` method 
if supported or performs a search around the vertex using `next()`. 
*/ 
Halfedge_handle find_prev_around_vertex(Halfedge_handle h); 

/*! 
returns the incident face of \f$ h\f$ if supported, 
`Face_handle()` otherwise. 
*/ 
Face_handle get_face(Halfedge_handle h); 

/*! 
returns the incident halfedge of \f$ f\f$ if supported, 
`Halfedge_handle()` otherwise. 
*/ 
Halfedge_handle get_face_halfedge( Face_handle f); 

/// @} 

/// \name Modifying Functions (Composed) 
/// @{

/*! 
makes `h->opposite()` the successor of \f$ h\f$. 
*/ 
void close_tip( Halfedge_handle h) const; 

/*! 
makes `h->opposite()` the successor of \f$ h\f$ and sets the 
incident vertex of \f$ h\f$ to \f$ v\f$. 
*/ 
void close_tip( Halfedge_handle h, Vertex_handle v) const; 

/*! 
inserts the tip of the edge \f$ h\f$ into the halfedges around the vertex 
pointed to by \f$ v\f$. Halfedge `h->opposite()` is the new successor of 
\f$ v\f$ and `h->next()` will be set to `v->next()`. The vertex of \f$ h\f$ 
will be set to the vertex \f$ v\f$ refers to if vertices are supported. 
*/ 
void insert_tip( Halfedge_handle h, Halfedge_handle v) const; 

/*! 
removes the edge `h->next()->opposite()` from the halfedge 
circle around the vertex referred to by \f$ h\f$. The new successor 
halfedge of \f$ h\f$ will be `h->next()->opposite()->next()`. 
*/ 
void remove_tip( Halfedge_handle h) const; 

/*! 
inserts the halfedge \f$ h\f$ between \f$ f\f$ and `f->next()`. 
The face of \f$ h\f$ will be the one \f$ f\f$ refers to if faces 
are supported. 
*/ 
void insert_halfedge( Halfedge_handle h, Halfedge_handle f) const; 

/*! 
removes edge `h->next()` from the halfedge circle around 
the face referred to by \f$ h\f$. The new successor of \f$ h\f$ will be 
`h->next()->next()`. 
*/ 
void remove_halfedge( Halfedge_handle h) const; 

/*! 
loops around the vertex incident to \f$ h\f$ and sets all vertex 
pointers to \f$ v\f$. \pre `h != Halfedge_handle()`. 
*/ 
void set_vertex_in_vertex_loop( Halfedge_handle h, 
Vertex_handle v) const; 

/*! 
loops around the face incident to \f$ h\f$ and sets all face 
pointers to \f$ f\f$. \pre `h != Halfedge_handle()`. 
*/ 
void set_face_in_face_loop( Halfedge_handle h, Face_handle f) const; 

/*! 
performs an edge flip. It returns \f$ h\f$ after rotating the edge \f$ h\f$ one 
vertex in the direction of the face orientation. 
\pre `h != Halfedge_handle()` and both incident faces of \f$ h\f$ are triangles. 
*/ 
Halfedge_handle flip_edge( Halfedge_handle h) const; 

/// @} 

/// \name Modifying Functions (Primitives) 
/// @{

/*! 
sets the incident halfedge of \f$ v\f$ to \f$ g\f$. 
*/ 
void set_vertex_halfedge( Vertex_handle v, Halfedge_handle g) const; 

/*! 
sets the incident halfedge of the vertex incident to \f$ h\f$ to \f$ h\f$. 
*/ 
void set_vertex_halfedge( Halfedge_handle h) const; 

/*! 
sets the incident vertex of \f$ h\f$ to \f$ v\f$. 
*/ 
void set_vertex( Halfedge_handle h, Vertex_handle v) const; 

/*! 
sets the previous link of \f$ h\f$ to \f$ g\f$. 
*/ 
void set_prev( Halfedge_handle h, Halfedge_handle g) const; 

/*! 
sets the incident face of \f$ h\f$ to \f$ f\f$. 
*/ 
void set_face( Halfedge_handle h, Face_handle f) const; 

/*! 
sets the incident halfedge of \f$ f\f$ to \f$ g\f$. 
*/ 
void set_face_halfedge( Face_handle f, Halfedge_handle g) const; 

/*! 
sets the incident halfedge of the face incident to \f$ h\f$ to \f$ h\f$. 
*/ 
void set_face_halfedge( Halfedge_handle h) const; 

/// @}

}; /* end HalfedgeDS_items_decorator */
} /* end namespace CGAL */
