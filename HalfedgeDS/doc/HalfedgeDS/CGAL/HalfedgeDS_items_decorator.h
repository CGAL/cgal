namespace CGAL {

/*!
\ingroup PkgHDS_Decorators

The class 
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

\cgalHeading{Example}

The following program fragment illustrates how a refined halfedge 
class for a polyhedron can make use of the `find_prev()` member 
function to implement a `prev()` member function that works 
regardless of whether the halfedge data structure `HDS` provides a 
`prev()` member function for its halfedges or not. In the case that not, 
the implementation given here runs in time proportional to the size of the 
incident face. For const-correctness a second implementation with signature 
`Halfedge_const_handle prev() const;` is needed. 

Note also the use of the static member function `HalfedgeDS::halfedge_handle()` 
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
typedef unspecified_type HalfedgeDS; 

/*!
traits class. 
*/ 
typedef unspecified_type Traits; 

/*!
vertex type of `HalfedgeDS`. 
*/ 
typedef unspecified_type Vertex; 

/*!
halfedge type of `HalfedgeDS`. 
*/ 
typedef unspecified_type Halfedge; 

/*!
face type of `HalfedgeDS`. 
*/ 
typedef unspecified_type Face; 

/*!

*/ 
typedef unspecified_type Vertex_handle; 

/*!

*/ 
typedef unspecified_type Halfedge_handle; 

/*!

*/ 
typedef unspecified_type Face_handle; 

/*!

*/ 
typedef unspecified_type Vertex_iterator; 

/*!

*/ 
typedef unspecified_type Halfedge_iterator; 

/*!

*/ 
typedef unspecified_type Face_iterator; 

/*!

*/ 
typedef unspecified_type size_type; 

/*!

*/ 
typedef unspecified_type difference_type; 

/*!

*/ 
typedef unspecified_type iterator_category; 

/*!

*/ 
typedef unspecified_type Supports_vertex_halfedge; 

/*!

*/ 
typedef unspecified_type Supports_halfedge_prev; 

/*!

*/ 
typedef unspecified_type Supports_halfedge_vertex; 

/*!

*/ 
typedef unspecified_type Supports_halfedge_face; 

/*!

*/ 
typedef unspecified_type Supports_face_halfedge; 

/*!

*/ 
typedef unspecified_type Supports_removal; 

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
returns the incident halfedge of `v` if supported, 
`Halfedge_handle()` otherwise. 
*/ 
Halfedge_handle get_vertex_halfedge( Vertex_handle v); 

/*!
returns the incident vertex of `h` if supported, `Vertex_handle()` 
otherwise. 
*/ 
Vertex_handle get_vertex(Halfedge_handle h); 

/*!
returns the previous halfedge of `h` if supported, 
`Halfedge_handle()` otherwise. 
*/ 
Halfedge_handle get_prev(Halfedge_handle h); 

/*!
returns the previous halfedge of `h`. Uses the `prev()` method 
if supported or performs a search around the face using `next()`. 
*/ 
Halfedge_handle find_prev(Halfedge_handle h); 

/*!
returns the previous halfedge of `h`. Uses the `prev()` method 
if supported or performs a search around the vertex using `next()`. 
*/ 
Halfedge_handle find_prev_around_vertex(Halfedge_handle h); 

/*!
returns the incident face of `h` if supported, 
`Face_handle()` otherwise. 
*/ 
Face_handle get_face(Halfedge_handle h); 

/*!
returns the incident halfedge of `f` if supported, 
`Halfedge_handle()` otherwise. 
*/ 
Halfedge_handle get_face_halfedge( Face_handle f); 

/// @} 

/// \name Modifying Functions (Composed) 
/// @{

/*!
makes `h->%opposite()` the successor of `h`. 
*/ 
void close_tip( Halfedge_handle h) const; 

/*!
makes `h->%opposite()` the successor of `h` and sets the 
incident vertex of `h` to `v`. 
*/ 
void close_tip( Halfedge_handle h, Vertex_handle v) const; 

/*!
inserts the tip of the edge `h` into the halfedges around the vertex 
pointed to by `v`. Halfedge `h->%opposite()` is the new successor of 
`v` and `h->%next()` will be set to `v->next()`. The vertex of `h` 
will be set to the vertex `v` refers to if vertices are supported. 
*/ 
void insert_tip( Halfedge_handle h, Halfedge_handle v) const; 

/*!
removes the edge `h->%next()->%opposite()` from the halfedge 
circle around the vertex referred to by `h`. The new successor 
halfedge of `h` will be `h->%next()->%opposite()->%next()`. 
*/ 
void remove_tip( Halfedge_handle h) const; 

/*!
inserts the halfedge `h` between `f` and `f->%next()`. 
The face of `h` will be the one `f` refers to if faces 
are supported. 
*/ 
void insert_halfedge( Halfedge_handle h, Halfedge_handle f) const; 

/*!
removes edge `h->%next()` from the halfedge circle around 
the face referred to by `h`. The new successor of `h` will be 
`h->%next()->%next()`. 
*/ 
void remove_halfedge( Halfedge_handle h) const; 

/*!
loops around the vertex incident to `h` and sets all vertex 
pointers to `v`. \pre `h != Halfedge_handle()`. 
*/ 
void set_vertex_in_vertex_loop( Halfedge_handle h, 
Vertex_handle v) const; 

/*!
loops around the face incident to `h` and sets all face 
pointers to `f`. \pre `h != Halfedge_handle()`. 
*/ 
void set_face_in_face_loop( Halfedge_handle h, Face_handle f) const; 

/*!
performs an edge flip. It returns `h` after rotating the edge `h` one 
vertex in the direction of the face orientation. 
\pre `h != Halfedge_handle()` and both incident faces of `h` are triangles. 
*/ 
Halfedge_handle flip_edge( Halfedge_handle h) const; 

/// @} 

/// \name Modifying Functions (Primitives) 
/// @{

/*!
sets the incident halfedge of `v` to `g`. 
*/ 
void set_vertex_halfedge( Vertex_handle v, Halfedge_handle g) const; 

/*!
sets the incident halfedge of the vertex `v` to `h`. 
*/ 
void set_vertex_halfedge( Halfedge_handle h) const; 

/*!
sets the incident vertex of `h` to `v`. 
*/ 
void set_vertex( Halfedge_handle h, Vertex_handle v) const; 

/*!
sets the previous link of `h` to `g`. 
*/ 
void set_prev( Halfedge_handle h, Halfedge_handle g) const; 

/*!
sets the incident face of `h` to `f`. 
*/ 
void set_face( Halfedge_handle h, Face_handle f) const; 

/*!
sets the incident halfedge of `f` to `g`. 
*/ 
void set_face_halfedge( Face_handle f, Halfedge_handle g) const; 

/*!
sets the incident halfedge of the face `f` to `h`. 
*/ 
void set_face_halfedge( Halfedge_handle h) const; 

/// @}

}; /* end HalfedgeDS_items_decorator */
} /* end namespace CGAL */
