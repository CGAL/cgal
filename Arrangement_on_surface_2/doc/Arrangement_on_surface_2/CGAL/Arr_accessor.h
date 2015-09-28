
namespace CGAL {

/*!
\ingroup PkgArrangement2

\anchor arr_refarr_access 

`Arr_accessor` provides an access to some of the private 
`Arrangement` functions. Users may use these functions to achieve 
more efficient programs when they have the exact topological information 
required by the specialized functions. 

It is however highly recommended to be very careful when using the 
accessor functions that modify the arrangement. As we have just mentioned, 
these functions have very specific requirement on their input on one hand, 
and perform no preconditions on the other hand, so providing incorrect 
topological input may invalidate the arrangement. 

*/
template< typename Arrangement >
class Arr_accessor {
public:

/// \name Types 
/// @{

/*!
the type of the associated arrangement. 
*/ 
typedef unspecified_type Arrangement_2; 

/*!
the point type. 
*/ 
typedef typename Arrangement_2::Point_2 Point_2; 

/*!
the \f$ x\f$-monotone curve type. 
*/ 
typedef typename Arrangement_2::X_monotone_curve_2 X_monotone_curve_2; 

/*!

*/ 
typedef typename Arrangement_2::Vertex_handle Vertex_handle; 

/*!

*/ 
typedef typename Arrangement_2::Halfedge_handle Halfedge_handle; 

/*!

*/ 
typedef typename Arrangement_2::Face_handle Face_handle; 

/*!
represents the boundary of a connected component (CCB). 
*/ 
typedef typename Arrangement_2::Ccb_halfedge_circulator Ccb_halfedge_circulator; 

/// @} 

/// \name Creation 
/// @{

/*!
constructs an accessor attached to the given arrangement 
`arr`. 
*/ 
Arr_accessor(Arrangement_2& arr); 

/// @} 

/// \name Accessing the Notification Functions 

/// @{

/*!
notifies the arrangement observer that a global change is going to take 
place (for the usage of the global functions that operate on 
arrangements). 
*/ 
void notify_before_global_change (); 

/*!
notifies the arrangement observer that a global change has taken 
place (for the usage of the global functions that operate on 
arrangements). 
*/ 
void notify_after_global_change (); 

/// @} 

/// \name Arrangement Predicates 
/// @{

/*!
locates a place for the curve `c` around the vertex `v` 
and returns a halfedge whose target is `v`, where c should be 
inserted between this halfedge and the next halfedge around `v` 
in a clockwise order. 
*/ 
Halfedge_handle 
locate_around_vertex (Vertex_handle v, 
const X_monotone_curve_2& c) const; 

/*!
counts the number of edges along the path from `e1` to `e2`. 
In case the two halfedges do not belong to the same connected component, 
the function returns (-1). 
*/ 
int halfedge_distance (Halfedge_const_handle e1, 
Halfedge_const_handle e2) const; 

/*!
determines whether a new halfedge we are about to create, which is to be associated 
with the curve `c` and directed from `pred1->target()` to 
`pred2->target()`, lies on the inner CCB of the new face that will be created, 
introducing this new edge. 
\pre `pred1->target()` and `pred2->target()` are associated with `c`'s endpoints. 
\pre `pred1` and `pred2` belong to the same connected component, such that a new face is created by connecting `pred1->target()` and `pred2->target()`. 
*/ 
bool is_inside_new_face (Halfedge_handle pred1, 
Halfedge_handle pred2, 
const X_monotone_curve_2& c) const; 

/*!
determines whether a given point lies within the region bounded by 
a boundary of the connected component that `he` belongs to. 
Note that if the function returns `true`, then `p` is contained within 
`he->face()` (but not on its boundary), or inside one of the holes inside this 
face, or it may coincide with an isolated vertex in this face. 
*/ 
bool point_is_in (const Point_2& p, 
Halfedge_const_handle he) const; 

/*!
determines whether `he` lies on the outer boundary of its incident 
face. 
*/ 
bool is_on_outer_boundary (Halfedge_const_handle he) const; 

/*!
determines whether `he` lies on the inner boundary of its incident 
face (that is, whether it lies on the boundary of one of the holes 
in this face). 
*/ 
bool is_on_inner_boundary (Halfedge_const_handle he) const; 

/// @} 

/// \name Arrangement Modifiers 
/// @{

/*!
creates a new vertex an associates it with the point `p`. 
\pre There is no existing vertex already associated with `p`. 
*/ 
Vertex_handle create_vertex (const Point_2& p); 

/*!
inserts the curve `c` as a new hole (inner component) of the face 
`f`, connecting the two isolated vertices `v1` and `v2`. 
`res` is the comparison result between these two end-vertices. 
The function returns a handle for one of the new halfedges 
corresponding to the inserted curve, directed from `v1` to `v2`. 
\pre `v1` and `v2` are associated with `c`'s endpoints, that they lie of `f`'s interior and that and that they have no incident edges. 
*/ 
Halfedge_handle insert_in_face_interior_ex (const X_monotone_curve_2& c, 
Face_handle f, 
Vertex_handle v1, 
Vertex_handle v2, 
Comparison_result res); 

/*!
inserts the curve `c` into the arrangement, such that one of its 
endpoints corresponds to an arrangement, which is the 
target vertex of the halfedge `pred`, such that `c` is inserted 
to the circular list of halfedges around `pred->target()` right 
between `pred` and its successor. The other end-vertex is given by 
an isolated vertex `v`, 
where `res` is the comparison result between `pred->target()` and `v`. 
The function returns a handle for one of the new halfedges directed from 
`pred->target()` to `v`. 
\pre `pred->target()` and `v` are associated with `c`'s endpoints and that and that `v` has no incident edges. 
*/ 
Halfedge_handle insert_from_vertex_ex (const X_monotone_curve_2& c, 
Halfedge_handle pred, 
Vertex_handle v, 
Comparison_result res); 

/*!
inserts the curve `c` into the arrangement, such that both `c`'s 
endpoints correspond to existing arrangement vertices, given by 
`pred1->target()` and `pred2->target()`. `res` is the comparison 
result between these two end-vertices. The function creates a 
new halfedge pair that connects the two vertices (with `pred1` and 
`pred2` indicate the exact place for these halfedges around the two 
target vertices) and returns a handle for the halfedge directed from 
`pred1->target()` to `pred2->target()`. 
The output flag `new_face` indicates whether a new face has been created 
following the insertion of the new curve. 
\pre `pred1->target()` and `pred2->target()` are associated with `c`'s endpoints and that if a new face is created, then `is_inside_new_face (pred1, pred2, c)` is `true`. 
*/ 
Halfedge_handle insert_at_vertices_ex (const X_monotone_curve_2& c, 
Halfedge_handle pred1, 
Halfedge_handle pred2, 
Comparison_result res, 
bool& new_face); 

/*!
inserts `v` as an isolated vertex inside `f`. 
\pre `v->point()` is contained in the interior of the given face. 
*/ 
void insert_isolated_vertex (Face_handle f, Vertex_handle v); 

/*!
moves the given hole from the interior of the face `f1` inside the 
face `f2`. 
\pre `hole` is currently contained in `f1` and should be moved to `f2`. 
*/ 
void move_hole (Face_handle f1, Face_handle f2, 
Ccb_halfedge_circulator hole); 

/*!
moves the given isolated vertex from the interior of the face `f1` 
inside the face `f2`. 
\pre `v` is indeed an isolated vertex currently contained in `f1` and should be moved to `f2`. 
*/ 
bool move_isolated_vertex (Face_handle f1, Face_handle f2, 
Vertex_handle v); 

/*!
relocates all holes and isolated vertices to their proper position 
immediately after a face has split due to the insertion of a new halfedge, 
namely after `insert_at_vertices_ex()` was invoked and indicated that 
a new face has been created. `he` is the halfegde returned by 
`insert_at_vertices_ex()`, 
such that `he->twin()->face` is the face that has just been split and 
`he->face()` is the newly created face. 
*/ 
void relocate_in_new_face (Halfedge_handle he); 

/*!
relocates all holes in a new face, as detailed above. 
*/ 
void relocate_holes_in_new_face (Halfedge_handle he); 

/*!
relocates all isolated vertices in a new face, as detailed above. 
*/ 
void relocate_isolated_vertices_in_new_face (Halfedge_handle he); 

/*!
modifies the point associated with the vertex `v` (the point may be 
geometrically different than the one currently associated with `v`). 
The function returns a handle to the modified vertex (same as `v`). 
\pre No other arrangement vertex is already associated with `p`. 
\pre The topology of the arrangement does not change after the vertex point is modified.
*/ 
Vertex_handle modify_vertex_ex (Vertex_handle v, 
const Point_2& p); 

/*!
modifies the \f$ x\f$-monotone curve associated with the edge `e` (the 
curve `c` may be geometrically different than the one currently 
associated with `e`). 
The function returns a handle to the modified edge (same as `e`). 
\pre The interior of `c` is disjoint from all existing arrangement vertices and edges. 
*/ 
Halfedge_handle modify_edge_ex (Halfedge_handle e, 
const X_monotone_curve_2& c); 

/*!
splits a given edge into two at the split point `p`, and associate the 
x-monotone curves `c1` and `c2` with the resulting edges, such that 
`c1` connects `he->source()` with `p` and `c2` connects 
`p` with `he->target()`. The function return a handle to the split 
halfedge directed from `he->source()` to the split point `p`. 
\pre The endpoints of `c1` and `c2` correspond to `p` and to `he`'s end-vertices, as indicated above. 
*/ 
Halfedge_handle split_edge_ex (Halfedge_handle he, 
const Point_2& p, 
const X_monotone_curve_2& c1, 
const X_monotone_curve_2& c2); 

/*!
splits a given edge into two at by the vertex `v`, and associate the 
x-monotone curves `c1` and `c2` with the resulting edges, such that 
`c1` connects `he->source()` with `v` and `c2` connects 
`v` with `he->target()`. The function return a handle to the split 
halfedge directed from `he->source()` to `v`. 
\pre The endpoints of `c1` and `c2` correspond to `v` and to `he`'s end-vertices, as indicated above. It is also assumed that `v` has no incident edges. 
*/ 
Halfedge_handle split_edge_ex (Halfedge_handle he, 
Vertex_handle v, 
const X_monotone_curve_2& c1, 
const X_monotone_curve_2& c2); 

/*!
removes the edge `he` from the arrangement, such that if the edge removal causes 
the creation of a new hole, `he->target()` lies on the boundary of this hole. 
The flags `remove_source` and `remove_target` indicate whether the end-vertices 
of `he` should be removed as well, in case they have no other incident edges. 
If the operation causes two faces to merge, the merged face is returned. 
Otherwise, the face to which the edge was incident is returned. 
*/ 
Face_handle remove_edge_ex (Halfedge_handle he, 
bool remove_source = true, 
bool remove_target = true); 

/// @}

}; /* end Arr_accessor */
} /* end namespace CGAL */
