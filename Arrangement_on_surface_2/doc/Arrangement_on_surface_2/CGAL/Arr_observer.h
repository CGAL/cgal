
namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2Ref

\anchor arr_refarr_obs

`Arr_observer` serves as an abstract base class for all observer
classes that are attached to an arrangement instance of type `Arrangement`
and receive notifications from the arrangement.
This base class handles the attachment of the
observer to a given arrangement instance or to the detachment of the
observer from this arrangement instance. It also gives a default empty
implementation to all notification functions that are invoked by the
arrangement to notify the observer on local or global changes it undergoes.
The notification functions are all virtual functions, so they can be
overridden by the concrete observer classes that inherit from
`Arr_observer`.

In order to implement a concrete arrangement observer-class, one simply
needs to derive from `Arr_observer` and override the relevant
notification functions. For example, if only face-split events are of
interest, it is sufficient to override just `before_split_face()`
(or just `after_split_face()`).

*/
template< typename Arrangement >
class Arr_observer {
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
In particular, holes are represented by a circulator for their outer CCB.
*/
typedef typename Arrangement_2::Ccb_halfedge_circulator Ccb_halfedge_circulator;

/// @}

/// \name Creation
/// @{

/*!
constructs an observer that is unattached to any arrangement instance.
*/
Arr_observer();

/*!
constructs an observer and attaches it to the given arrangement
`arr`.
*/
Arr_observer(Arrangement_2& arr);

/// @}

/// \name Modifiers
/// @{

/*!
attaches the observer to the given arrangement `arr`.
*/
void attach (Arrangement_2& arr);

/*!
detaches the observer from its arrangement.
*/
void detach ();

/// @}

/// \name Notifications on Global Arrangement Operations
/// @{

/*!
issued just before the attached arrangement is assigned with the contents
of another arrangement `arr`.
*/
virtual void before_assign (const Arrangement_2& arr);

/*!
issued immediately after the attached arrangement has been assigned with
the contents of another arrangement.
*/
virtual void after_assign ();

/*!
issued just before the attached arrangement is cleared.
*/
virtual void before_clear ();

/*!
issued immediately after the attached arrangement has been cleared, so it
now consists only of a the unbounded face `uf`.
*/
virtual void after_clear (Face_handle uf);

/*!
issued just before a global function starts to modify the attached
arrangement. It is guaranteed that no queries (especially no
point-location queries) are issued until the termination of the global
function is indicated by `after_global_change()`.
*/
virtual void before_global_change ();

/*!
issued immediately after a global function has stopped modifying the
attached arrangement.
*/
virtual void after_global_change ();

/// @}

/// \name Notifications on Attachment or Detachment
/// @{

/*!
issued just before the observer is attached to the arrangement instance
`arr`.
*/
virtual void before_attach (const Arrangement_2& arr);

/*!
issued immediately after the observer has been attached to an
arrangement instance.
*/
virtual void after_attach ();

/*!
issued just before the observer is detached from its arrangement instance.
*/
virtual void before_detach ();

/*!
issued immediately after the observer has been detached from its
arrangement instance.
*/
virtual void after_detach ();

/// @}

/// \name Notifications on Local Changes in the Arrangement
/// @{

/*!
issued just before a new vertex that corresponds to the point `p`
is created.
*/
virtual void before_create_vertex (const Point_2& p);

/*!
issued immediately after a new vertex `v` has been created.
Note that the vertex still has no incident edges and is not connected
to any other vertex.
*/
virtual void after_create_vertex (Vertex_handle v);

/*!
issued just before a new vertex at infinity is created, `cv` is the
curve incident to the surface boundary, `ind` is the relevant curve-end,
`ps_x` is the boundary condition of the vertex in \f$ x\f$ and `ps_y`
is the boundary condition of the vertex in \f$ y\f$.
*/
virtual void before_create_boundary_vertex (const X_monotone_curve_2& cv ,
Arr_curve_end ind ,
Arr_parameter_space ps_x ,
Arr_parameter_space ps_y );

/*!
issued immediately after a new vertex `v` has been created.
Note that the vertex still has no incident edges and is not connected
to any other vertex.
*/
virtual void after_create_boundary_vertex (Vertex_handle v);

/*!
issued just before a new edge that corresponds to the \f$ x\f$-monotone curve
`c` and connects the vertices `v1` and `v2` is created.
*/
virtual void before_create_edge (const X_monotone_curve_2& c,
Vertex_handle v1,
Vertex_handle v2);

/*!
issued immediately after a new edge `e` has been created.
The halfedge that is sent to this function is always directed from
`v1` to `v2` (see above).
*/
virtual void after_create_edge (Halfedge_handle e);

/*!
issued just before a vertex `v` is modified to be associated with
the point `p`.
*/
virtual void before_modify_vertex (Vertex_handle v,
const Point_2& p);

/*!
issued immediately after an existing vertex `v` has been modified.
*/
virtual void after_modify_vertex (Vertex_handle v);

/*!
issued just before an edge `e` is modified to be associated with
the \f$ x\f$-monotone curve `c`.
*/
virtual void before_modify_edge (Halfedge_handle e,
const X_monotone_curve_2& c);

/*!
issued immediately after an existing edge `e` has been modified.
*/
virtual void after_modify_edge (Halfedge_handle e);

/*!
issued just before an edge `e` is split into two edges that should
be associated with the \f$ x\f$-monotone curves `c1` and `c2`. The
vertex `v` corresponds to the split point, and will be used to
separate the two resulting edges.
*/
virtual void before_split_edge (Halfedge_handle e,
Vertex_handle v,
const X_monotone_curve_2& c1,
const X_monotone_curve_2& c2);

/*!
issued immediately after an existing edge has been split into the two
given edges `e1` and `e2`.
*/
virtual void after_split_edge (Halfedge_handle e1,
Halfedge_handle e2);

/*!
issued just before a fictitious edge `e` is split into two. The
vertex at infinity `v` corresponds to the split point, and will be
used to separate the two resulting edges.
*/
virtual void before_split_fictitious_edge (Halfedge_handle e,
Vertex_handle v);

/*!
issued immediately after an existing fictitious edge has been split into
the two given fictitious edges `e1` and `e2`.
*/
virtual void after_split_fictitious_edge (Halfedge_handle e1,
Halfedge_handle e2);

/*!
issued just before a face `f` is split into two, as a result of
the insertion of the edge `e` into the arrangement.
*/
virtual void before_split_face (Face_handle f,
Halfedge_handle e);

/*!
issued immediately after the existing face `f1` has been split
into the two faces `f1` and `f2`. The flag
`is_hole` designates whether `f2` forms a hole inside `f1`.
*/
virtual void after_split_face (Face_handle f1,
Face_handle f2,
bool is_hole);

/*!
issued just before outer ccb `h` inside a face `f` is split into
two, as a result of the removal of the edge `e` from the arrangement.
*/
virtual void before_split_outer_ccb (Face_handle f,
Ccb_halfedge_circulator h,
Halfedge_handle e);

/*!
issued immediately after outer ccb the face `f` has been split,
resulting in the two holes `h1` and `h2`.
*/
virtual void after_split_outer_ccb (Face_handle f,
Ccb_halfedge_circulator h1,
Ccb_halfedge_circulator h2);

/*!
issued just before inner ccb `h` inside a face `f` is split into
two, as a result of the removal of the edge `e` from the arrangement.
*/
virtual void before_split_inner_ccb (Face_handle f,
Ccb_halfedge_circulator h,
Halfedge_handle e);

/*!
issued immediately after inner ccb the face `f` has been split,
resulting in the two holes `h1` and `h2`.
*/
virtual void after_split_inner_ccb (Face_handle f,
Ccb_halfedge_circulator h1,
Ccb_halfedge_circulator h2);

/*!
issued just before the edge `e` is inserted as a new outer ccb inside
the face `f`.
*/
virtual void before_add_outer_ccb (Face_handle f,
Halfedge_handle e);

/*!
issued immediately after a new outer ccb `h` has been created. The
outer ccb always consists of a single pair of twin halfedges.
*/
virtual void after_add_outer_ccb (Ccb_halfedge_circulator h);

/*!
issued just before the edge `e` is inserted as a new inner ccb inside
the face `f`.
*/
virtual void before_add_inner_ccb (Face_handle f,
Halfedge_handle e);

/*!
issued immediately after a new inner ccb `h` has been created. The
inner ccb always consists of a single pair of twin halfedges.
*/
virtual void after_add_inner_ccb (Ccb_halfedge_circulator h);

/*!
issued just before the vertex `v` is inserted as an isolated
vertex inside the face `f`.
*/
virtual void before_add_isolated_vertex (Face_handle f,
Vertex_handle v);

/*!
issued immediately after the vertex `v` has been set as an
isolated vertex.
*/
virtual void after_add_isolated_vertex (Vertex_handle v);

/*!
issued just before the two edges `e1` and `e2` are merged to
form a single edge that will be associated with the \f$ x\f$-monotone curve
`c`.
*/
virtual void before_merge_edge (Halfedge_handle e1,
Halfedge_handle e2,
const X_monotone_curve_2& c);

/*!
issued immediately after two edges have been merged to form the edge
`e`.
*/
virtual void after_merge_edge (Halfedge_handle e);

/*!
issued just before the two fictitious edges `e1` and `e2` are
merged to form a single fictitious edge.
*/
virtual void before_merge_fictitious_edge (Halfedge_handle e1,
Halfedge_handle e2);

/*!
issued immediately after two fictitious edges have been merged to form
the fictitious edge `e`.
*/
virtual void after_merge_fictitious_edge (Halfedge_handle e);

/*!
issued just before the two edges `f1` and `f2` are merged to
form a single face, following the removal of the edge `e` from the
arrangement.
*/
virtual void before_merge_face (Face_handle f1,
Face_handle f2,
Halfedge_handle e);

/*!
issued immediately after two faces have been merged to form the face
`f`.
*/
virtual void after_merge_face (Face_handle f);

/*!
issued just before two outer ccbs `h1` and `h2` inside the face
`f` are merged to form a single connected component, following the
insertion of the edge `e` into the arrangement.
*/
virtual void before_merge_outer_ccb (Face_handle f,
Ccb_halfedge_circulator h1,
Ccb_halfedge_circulator h2,
Halfedge_handle e);

/*!
issued immediately after two outer ccbs have been merged to form a single
outer ccb `h` inside the face `f`.
*/
virtual void after_merge_outer_ccb (Face_handle f,
Ccb_halfedge_circulator h);

/*!
issued just before two inner ccbs `h1` and `h2` inside the face
`f` are merged to form a single connected component, following the
insertion of the edge `e` into the arrangement.
*/
virtual void before_merge_inner_ccb (Face_handle f,
Ccb_halfedge_circulator h1,
Ccb_halfedge_circulator h2,
Halfedge_handle e);

/*!
issued immediately after two inner ccbs have been merged to form a single
inner ccb `h` inside the face `f`.
*/
virtual void after_merge_inner_ccb (Face_handle f,
Ccb_halfedge_circulator h);

/*!
issued just before the outer ccb `h` is moved from one face to another.
This can happen if the face `to_f` containing the outer ccb has just
been split from `from_f`.
*/
virtual void before_move_outer_ccb (Face_handle from_f,
Face_handle to_f,
Ccb_halfedge_circulator h);

/*!
issued immediately after the outer ccb `h` has been moved to a new face.
*/
virtual void after_move_outer_ccb (Ccb_halfedge_circulator h);

/*!
issued just before the inner ccb `h` is moved from one face to another.
This can happen if the face `to_f` containing the inner ccb has just
been split from `from_f`.
*/
virtual void before_move_inner_ccb (Face_handle from_f,
Face_handle to_f,
Ccb_halfedge_circulator h);

/*!
issued immediately after the inner ccb `h` has been moved to a new face.
*/
virtual void after_move_inner_ccb (Ccb_halfedge_circulator h);

/*!
issued just before the isolated vertex `v` is moved from one face
to another.
This can happen if the face `to_f` containing the isolated vertex
has just been split from `from_f`.
*/
virtual void before_move_isolated_vertex (Face_handle from_f,
Face_handle to_f,
Vertex_handle v);

/*!
issued immediately after the isolated vertex `v` has been moved to a
new face.
*/
virtual void after_move_isolated_vertex (Vertex_handle v);

/*!
issued just before the vertex `v` is removed from the arrangement.
*/
virtual void before_remove_vertex (Vertex_handle v);

/*!
issued immediately after a vertex has been removed (and deleted)
from the arrangement.
*/
virtual void after_remove_vertex ();

/*!
issued just before the edge `e` is removed from the arrangement.
*/
virtual void before_remove_edge (Halfedge_handle e);

/*!
issued immediately after an edge has been removed (and deleted)
from the arrangement.
*/
virtual void after_remove_edge ();

/*!
issued just before the outer ccb `f` is removed from inside the
face `f`.
*/
virtual void before_remove_outer_ccb (Face_handle f,
Ccb_halfedge_circulator h);

/*!
issued immediately after a outer ccb has been removed (and deleted)
from inside the face `f`.
*/
virtual void after_remove_outer_ccb (Face_handle f);

/*!
issued just before the inner ccb `f` is removed from inside the
face `f`.
*/
virtual void before_remove_inner_ccb (Face_handle f,
Ccb_halfedge_circulator h);

/*!
issued immediately after a inner ccb has been removed (and deleted)
from inside the face `f`.
*/
virtual void after_remove_inner_ccb (Face_handle f);

/// @}

}; /* end Arr_observer */
} /* end namespace CGAL */
