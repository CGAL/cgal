/*!
\ingroup PkgHDSConcepts
\cgalConcept

The concept of a halfedge data structure (abbreviated as `HalfedgeDS`, or 
`HDS` for template parameters) defines an edge-centered data structure 
capable of maintaining incidence information of vertices, edges, and 
faces, for example for planar maps or polyhedral surfaces. It is a 
combinatorial data structure, geometric interpretation is added by 
classes built on top of the halfedge data structure. 

The data structure defined here is known as the 
FE-structure \cgalCite{w-ebdss-85}, as 
halfedges \cgalCite{m-ism-88}, \cgalCite{cgal:bfh-mgedm-95} or as the doubly connected edge 
list (DCEL) \cgalCite{bkos-cgaa-97}, although the original reference for 
the DCEL \cgalCite{mp-fitcp-78} describes a different data structure. The 
halfedge data structure can also be seen as one of the variants of the 
quad-edge data structure \cgalCite{gs-pmgsc-85}. In general, the quad-edge 
data can represent non-orientable 2-manifolds, but the variant here is 
restricted to orientable 2-manifolds only. An overview and comparison 
of these different data structures together with a thorough 
description of the design implemented here can be found 
in \cgalCite{k-ugpdd-99}. 

Each edge is represented by two halfedges with opposite orientations. 
Each halfedge can store a reference to an incident face and an 
incident vertex. For each face and each vertex an incident halfedge 
is stored. Reduced variants of the halfedge data structure can omit 
some of these incidences, for example the reference to halfedges in 
vertices or the storage of vertices at all. See 
Figure \ref figureOptionalMethods 
for the incidences, the mandatory and optional member functions 
possible for vertices, halfedges, and faces. 

\anchor figureOptionalMethods 
\image html hds_optional_small.png "The three classes Vertex, Halfedge, and Face of the halfedge data structure. Member functions with shaded background are mandatory. The others are optionally supported."
\image latex hds_optional_small.png "The three classes Vertex, Halfedge, and Face of the halfedge data structure. Member functions with shaded background are mandatory. The others are optionally supported."

A `HalfedgeDS` organizes the internal storage of its items. Examples 
are a list-based or a vector-based storage. The `HalfedgeDS` exhibits 
most of the characteristics of the container class used internally, 
for example the iterator category. A vector resizes 
automatically when a new item exceeds the reserved space. Since 
resizing is an expensive operation for a `HalfedgeDS` in general and 
only possible in a well defined state of the data structure (no 
dangling handles), it must be called explicitly in advance for a 
`HalfedgeDS` before inserting new items beyond the current capacity. 
Classes built on top of a `HalfedgeDS` are advised to call the 
`reserve()` member function before creating new items. 

\cgalHeading{Parameters}

A `HalfedgeDS` is a class template and will be used as argument for 
other class templates, for example `CGAL::Polyhedron_3`. The 
template parameters to instantiate the `HalfedgeDS` will be provided by 
this other class template. Therefore, the three template parameters 
and their meaning are mandatory. We distinguish between the template 
`HalfedgeDS` and an instantiation of it. 

`Traits` is a traits class that will be passed to the 
item types in `Items`. It will not be used in `HalfedgeDS` itself. `Items` is a model of the `HalfedgeDSItems` concept. 
`Alloc` is an allocator that fulfills all requirements 
of allocators for \stl container classes. The `rebind` 
mechanism from `Alloc` will be used to create appropriate 
allocators internally. A default argument is mandatory for 
`Alloc`, for example, the macro `CGAL_ALLOCATOR(int)` 
from the `<CGAL/memory.h>` header file can be used as default 
allocator. 

\cgalHasModel CGAL::HalfedgeDS_default 
\cgalHasModel CGAL::HalfedgeDS_list 
\cgalHasModel CGAL::HalfedgeDS_vector 

\sa `HalfedgeDSItems` 
\sa `CGAL::Polyhedron_3<Traits>` 
\sa `CGAL::HalfedgeDS_vertex_base<Refs>` 
\sa `CGAL::HalfedgeDS_halfedge_base<Refs>` 
\sa `CGAL::HalfedgeDS_face_base<Refs>` 
\sa `CGAL::HalfedgeDS_items_decorator<HDS>` 
\sa `CGAL::HalfedgeDS_decorator<HDS>` 
\sa `CGAL::HalfedgeDS_const_decorator<HDS>` 

*/
template< typename Traits, typename Items, typename Alloc >
class HalfedgeDS {
public:

/// \name Types 
/// @{

/*!
traits class. 
*/ 
typedef unspecified_type Traits; 

/*!
model of `HalfedgeDSItems` concept. 
*/ 
typedef unspecified_type Items; 

/*!
size type. 
*/ 
typedef unspecified_type size_type; 

/*!
difference type. 
*/ 
typedef unspecified_type difference_type; 

/*!
iterator category for all iterators. 
*/ 
typedef unspecified_type iterator_category; 

/*!
allocator type `Alloc`. 
*/ 
typedef unspecified_type allocator_type; 

/*!
model of `HalfedgeDSVertex` concept. 
*/ 
typedef unspecified_type Vertex; 

/*!
model of `HalfedgeDSHalfedge` concept. 
*/ 
typedef unspecified_type Halfedge; 

/*!
model of `HalfedgeDSFace` concept. 
*/ 
typedef unspecified_type Face; 

/// @}

/*! \name Handles and Iterators
The following handles and iterators have appropriate non-mutable 
counterparts, i.e., `const_handle` and `const_iterator`. The 
mutable types are assignable to their non-mutable counterparts. The 
iterators are assignable to the respective handle types. Wherever the 
handles appear in function parameter lists, the corresponding 
iterators can be used as well. \note The handle types must have 
a default constructor that creates a unique and always the same handle 
value. It will be used in analogy to `NULL` for pointers. 
*/
/// @{

/*!
handle to vertex. 
*/ 
typedef unspecified_type Vertex_handle; 

/*!
handle to halfedge. 
*/ 
typedef unspecified_type Halfedge_handle; 

/*!
handle to face. 
*/ 
typedef unspecified_type Face_handle; 

/*!
iterator over all vertices. 
*/ 
typedef unspecified_type Vertex_iterator; 

/*!
iterator over all halfedges. 
*/ 
typedef unspecified_type Halfedge_iterator; 

/*!
iterator over all faces. 
*/ 
typedef unspecified_type Face_iterator; 

/// @} 

/* \name Types for Tagging Optional Features 
\cgalAdvancedBegin
The following types are equal to either `CGAL::Tag_true` or 
`CGAL::Tag_false`, depending on whether the named feature is 
supported or not. 

The following dependencies among these options must be regarded: 

Vertices are supported \f$ \Longleftrightarrow\f$ 
`Supports_halfedge_vertex` \f$ \equiv\f$ `CGAL::Tag_true`. 

Faces are supported \f$ \Longleftrightarrow\f$ 
`Supports_halfedge_face` \f$ \equiv\f$ `CGAL::Tag_true`. 

`Supports_vertex_halfedge` \f$ \equiv\f$ `CGAL::Tag_true` \f$ \Longrightarrow\f$ 
`Supports_halfedge_vertex` \f$ \equiv\f$ `CGAL::Tag_true`. 

`Supports_vertex_point` \f$ \equiv\f$ `CGAL::Tag_true` \f$ \Longrightarrow\f$ 
`Supports_halfedge_vertex` \f$ \equiv\f$ `CGAL::Tag_true`. 

`Supports_face_halfedge` \f$ \equiv\f$ `CGAL::Tag_true` \f$ \Longrightarrow\f$ 
`Supports_halfedge_face` \f$ \equiv\f$ `CGAL::Tag_true`. 
\cgalAdvancedEnd
*/
/// @{

/*!
`Vertex::halfedge()`. 
*/ 
typedef unspecified_type Supports_vertex_halfedge; 

/*!
`Halfedge::prev()`. 
*/ 
typedef unspecified_type Supports_halfedge_prev; 

/*!
`Halfedge::vertex()`. 
*/ 
typedef unspecified_type Supports_halfedge_vertex; 

/*!
`Halfedge::face()`. 
*/ 
typedef unspecified_type Supports_halfedge_face; 

/*!
`Face::halfedge()`. 
*/ 
typedef unspecified_type Supports_face_halfedge; 

/*!
removal of individual elements. 
*/ 
typedef unspecified_type Supports_removal; 

/// @}

/*! \name Static Member Functions

\cgalAdvancedBegin
When writing an items type, such as a user defined vertex, certain 
functions need to create a handle but knowing only a pointer, for 
example, the `this`-pointer. The following static member functions 
of `HalfedgeDS` create such a corresponding handle for an item type 
from a pointer. This conversion encapsulates possible adjustments for 
hidden data members in the true item type, such as linked-list 
pointers. Note that the user provides item types with the 
`Items` template argument, which may differ from the `Vertex`, 
`Halfedge`, and `Face` types defined in `HalfedgeDS`. If they 
differ, they are derived from the user provided item types. We denote the 
user item types with `Vertex_base`, `Halfedge_base`, and 
`Face_base` in the following. The fully qualified name for 
`Vertex_base` would be for example - assuming that the type `Self` 
refers to the instantiated `HalfedgeDS` - 

\code
typedef typename Items::template Vertex_wrapper<Self,Traits> Vertex_wrapper;
typedef typename Vertex_wrapper::Vertex Vertex_base;
\endcode

Implementing these functions relies on the fundamental assumption that 
an iterator (or handle) of the internally used container class can be 
constructed from a pointer of a contained item only. This is true and 
controlled by us for `CGAL::In_place_list`. It is true for the 
`std::vector` of major \stl distributions, but not necessarily 
guaranteed. We might switch to an internal implementation if need 
arises.
\cgalAdvancedEnd
*/
/// @{

/*!

*/ 
static Vertex_handle vertex_handle( Vertex_base* v); 

/*!

*/ 
static Vertex_const_handle vertex_handle( const Vertex_base* v); 

/*!

*/ 
static Halfedge_handle halfedge_handle( Halfedge_base* h); 

/*!

*/ 
static Halfedge_const_handle halfedge_handle( const Halfedge_base* h); 

/*!

*/ 
static Face_handle face_handle( Face_base* f); 

/*!

*/ 
static Face_const_handle face_handle( const Face_items* f); 

/// @} 

/// \name Creation 
/// @{

/*!
empty halfedge data structure. 
*/ 
HalfedgeDS(); 

/*!
storage reserved for `v` vertices, `h` halfedges, and `f` faces. 
*/ 
HalfedgeDS( size_type v, size_type h, size_type f); 

/*!
copy constructor. \pre `hds2` contains no dangling handles. 
*/ 
HalfedgeDS( const HalfedgeDS<Traits,Items,Alloc>& hds2); 

/*!
assignment operator. \pre `hds2` contains no dangling handles. 
*/ 
HalfedgeDS<Traits,Items,Alloc>& 
operator=( const HalfedgeDS<Traits,Items,Alloc>& hds2); 

/*!
reserves storage for `v` vertices, `h` halfedges, and `f` faces. 
If all capacities are already greater or equal than the requested sizes 
nothing happens. Otherwise, the halfedge data structure will be resized and all handles, 
iterators and circulators are invalid. 

\pre If resizing is necessary the halfedge data structure contains no dangling handles. 
*/ 
void reserve( size_type v, size_type h, size_type f); 

/// @} 

/// \name Access Member Functions 
/// @{

/*!
number of vertices. 
*/ 
Size size_of_vertices() const; 

/*!
number of halfedges. 
*/ 
Size size_of_halfedges() const; 

/*!
number of faces. 
*/ 
Size size_of_faces() const; 

/*!
space reserved for vertices. 
*/ 
Size capacity_of_vertices() const; 

/*!
space reserved for halfedges. 
*/ 
Size capacity_of_halfedges() const; 

/*!
space reserved for faces. 
*/ 
Size capacity_of_faces() const; 

/*!
bytes used for the halfedge data structure. 
*/ 
size_t bytes() const; 

/*!
bytes reserved for the halfedge data structure. 
*/ 
size_t bytes_reserved() const; 

/*!
allocator object. 
*/ 
allocator_type get_allocator() const; 

/*!
iterator over all vertices. 
*/ 
Vertex_iterator vertices_begin(); 

/*!

*/ 
Vertex_iterator vertices_end(); 

/*!
iterator over all halfedges 
*/ 
Halfedge_iterator halfedges_begin(); 

/*!

*/ 
Halfedge_iterator halfedges_end(); 

/*!
iterator over all faces. 
*/ 
Face_iterator faces_begin(); 

/*!

*/ 
Face_iterator faces_end(); 

/// @} 

/*! \name Insertion 
Note that the vertex-related and the face-related member functions may 
not be provided for a `HalfedgeDS` that does not support vertices or 
faces respectively. 
*/
/// @{

/*!
appends a copy of `v` to the halfedge data structure. Returns a handle of the new vertex. 
*/ 
Vertex_handle vertices_push_back( const Vertex& v); 

/*!
appends a copy of `h` and a copy of `g` to the halfedge data structure and makes them 
opposite to each other. Returns a handle of the copy of `h`. 
*/ 
Halfedge_handle edges_push_back( const Halfedge& h, 
const Halfedge& g); 

/*!
appends a copy of `h` and a copy of  \link HalfedgeDSHalfedge::opposite() `h->opposite()`\endlink  
to the halfedge data structure and makes them opposite to each other. Returns a handle of the copy of `h`. 
\pre  \link HalfedgeDSHalfedge::opposite() `h->opposite()`\endlink  denotes a halfedge. 
*/ 
Halfedge_handle edges_push_back( const Halfedge& h); 

/*!
appends a copy of `f` to the halfedge data structure. Returns a handle of the new face. 
*/ 
Face_handle faces_push_back( const Face& f); 

/// @} 

/*! \name Removal 
Erasing single elements is optional and indicated with the type tag 
`Supports_removal`. The three `pop_back()` and the `clear()` member 
functions are mandatory. If vertices or faces are not supported 
for a `HalfedgeDS` the three `pop_back()` and the `clear()` member 
functions must be provided as null operations. 
*/
/// @{

/*!
removes the first vertex if vertices are supported and 
`Supports_removal` \f$ \equiv\f$ `CGAL::Tag_true`. 
*/ 
void vertices_pop_front(); 

/*!
removes the last vertex. 
*/ 
void vertices_pop_back(); 

/*!
removes the vertex `v` if vertices are supported and 
`Supports_removal` \f$ \equiv\f$ `CGAL::Tag_true`. 
*/ 
void vertices_erase( Vertex_handle v); 

/*!
removes the range of vertices `[first,last)` if vertices 
are supported and `Supports_removal` \f$ \equiv\f$ `CGAL::Tag_true`. 
*/ 
void vertices_erase( Vertex_handle first, Vertex_handle last); 

/*!
removes the first two halfedges if 
`Supports_removal` \f$ \equiv\f$ `CGAL::Tag_true`. 
*/ 
void edges_pop_front(); 

/*!
removes the last two halfedges. 
*/ 
void edges_pop_back(); 

/*!
removes the pair of halfedges `h` and \link HalfedgeDSHalfedge::opposite() `h->opposite()`\endlink 
if `Supports_removal` \f$ \equiv\f$ `CGAL::Tag_true`. 
*/ 
void edges_erase( Halfedge_handle h); 

/*!
removes the range of edges `[first},last)` if 
`Supports_removal` \f$ \equiv\f$ `CGAL::Tag_true`. 
*/ 
void edges_erase( Halfedge_handle first, Halfedge_handle last); 

/*!

removes the first face if faces are supported and 
`Supports_removal` \f$ \equiv\f$ `CGAL::Tag_true`. 
*/ 
void faces_pop_front(); 

/*!

removes the last face. 
*/ 
void faces_pop_back(); 

/*!

removes the face `f` if faces are supported and 
`Supports_removal` \f$ \equiv\f$ `CGAL::Tag_true`. 
*/ 
void faces_erase( Face_handle f); 

/*!

removes the range of faces `[first,last)` if faces are 
supported and `Supports_removal` \f$ \equiv\f$ `CGAL::Tag_true`. 
*/ 
void faces_erase( Face_handle first, Face_handle last); 

/*!
removes all vertices. 
*/ 
void vertices_clear(); 

/*!
removes all halfedges. 
*/ 
void edges_clear(); 

/*!
removes all faces. 
*/ 
void faces_clear(); 

/*!
removes all elements. 
*/ 
void clear(); 

/// @} 

/*! \name Operations with Border Halfedges 
\cgalAdvancedBegin
The following notion of <I>border halfedges</I> is particular useful 
where the halfedge data structure is used to model surfaces with 
boundary, i.e., surfaces with missing faces or open regions. Halfedges 
incident to an open region are called <I>border halfedges</I>. A 
halfedge is a <I>border edge</I> if the halfedge itself or its 
opposite halfedge is a border halfedge. The only requirement to work 
with border halfedges is that the 
`Halfedge` class provides a member function `HalfedgeDSHalfedge::is_border()` 
returning a `bool`. Usually, the halfedge data structure 
supports faces and the value of the default constructor of the face 
handle will indicate a border halfedge, but this may not be the only 
possibility. The `HalfedgeDSHalfedge::is_border()` predicate divides the edges into 
two classes, the border edges and the non-border edges. The 
following normalization reorganizes the sequential storage of the 
edges such that the non-border edges precede the border edges, and 
that for each border edge the latter of the two halfedges is a 
border halfedge (the first one might be a border halfedge too). The 
normalization stores the number of border halfedges, as well as the 
halfedge iterator where the border edges start at, within the 
halfedge data structure. These values will be invalid after further 
halfedge insertions or removals and changes in the border status of 
a halfedge. There is no automatic update required. 
\cgalAdvancedEnd
*/
/// @{

/*!
sorts halfedges such that the non-border edges precede the 
border edges. For each border edge that is incident to a face, 
the halfedge iterator will reference the halfedge incident to the 
face right before the halfedge incident to the open region. 
*/ 
void normalize_border(); 

/*!
number of border halfedges. An edge with no incident face 
counts as two border halfedges. 
\pre `normalize_border()` has been called and no halfedge insertion or removal and no change in border status of the halfedges have occurred since then. 
*/ 
Size size_of_border_halfedges() const; 

/*!
number of border edges. If `size_of_border_edges()` is equal 
to `size_of_border_halfedges()` all border edges are incident to 
a face on one side and to an open region on the other side. 
\pre `normalize_border()` has been called and no halfedge insertion or removal and no change in border status of the halfedges have occurred since then. 
*/ 
Size size_of_border_edges() const; 

/*!
halfedge iterator starting with the border edges. The range 
[`halfedges_begin(), border_halfedges_begin()`) denotes 
all non-border edges. The range 
[`border_halfedges_begin(), halfedges_end()`) denotes all 
border edges. 
\pre `normalize_border()` has been called and no halfedge insertion or removal and no change in border status of the halfedges have occurred since then. 
*/ 
Halfedge_iterator border_halfedges_begin(); 

/// @}

}; /* end HalfedgeDS */
