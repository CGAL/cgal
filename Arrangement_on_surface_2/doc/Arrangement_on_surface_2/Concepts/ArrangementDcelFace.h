/*!
 * \ingroup PkgArrangementOnSurface2ConceptsDCEL
 * \cgalConcept
 *
 * A face record in a \dcel data structure. A face represents a region, which
 * may have outer and inner boundaries. A boundary conists of a chain of
 * incident halfedges, referred to as a Connected Component of the Boundary
 * (CCB). A face may be unbounded. Otherwise, it has one or more outer CCBs. A
 * face may also be bounded by inner CCBs, and it may contain isolated vertices
 * in its interior. A planar face may have only one outer CCBs and its inner
 * CCBs are referred to as holes.
 *
 * \sa `ArrangementDcel`
 * \sa `ArrangementDcelVertex`
 * \sa `ArrangementDcelHalfedge`
 */

class ArrangementDcelFace {
public:

/// \name Types
/// The non-mutable iterators `Outer_ccb_const_iterator`,
/// `Inner_ccb_const_iterator`, `Hole_const_iterator`, and
/// `Isolated_vertex_const_iterator` are also defined.
/// @{

/*! the corresponding \dcel vertex type. */
typedef unspecified_type Vertex;

/*! the corresponding \dcel halfedge type. */
typedef unspecified_type Halfedge;

/*! a bidirectional iterator over the outer CCBs of the face. Its value-type
 * is `Halfedge*`.
 */
typedef unspecified_type Outer_ccb_iterator;

/*! a bidirectional iterator over the inner CCBs of the face. Its value-type
 * is `Halfedge*`.
 */
typedef unspecified_type Inner_ccb_iterator;

/*! a bidirectional iterator over the holes (i.e., inner CCBs) of the face. Its
 * value-type is `Halfedge*`.
 */
typedef unspecified_type Hole_iterator;

/*! a bidirectional iterator over the isolated vertices in inside the face.
 * Its value-type is `Vertex*`.
 */
typedef unspecified_type Isolated_vertex_iterator;

/// @}

/// \name Creation
/// @{

/*! default constructor. */
Arr_dcel_face();

/*! assigns `f` with the contents of the `other` face. */
void assign(const Self& other);

/// @}

/// \name Access Functions
/// All functions below also have `const` counterparts, returning
/// non-mutable pointers or iterators:
/// @{

/*! determines whether the face is unbounded. */
bool is_unbounded() const;

/*! obtains an incident halfedge along the outer boundaries of the face.  If `f`
 * has no outer boundary, the function returns `nullptr`.
 */
Halfedge* halfedge();

/*! obtains the number of outer CCBs of `f`. In case of planar arrangement
 * this is either 0 or 1.
 */
size_t number_of_outer_ccbs() const;

/*! obtains a begin iterator for the outer CCBs of `f`. */
Outer_ccb_iterator outer_ccbs_begin();

/*! obtains a past-the-end iterator for the outer CCBs of `f`. */
Outer_ccb_iterator outer_ccbs_end();

/*! obtains the number of inner CCBs of `f`. */
size_t number_of_inner_ccbs() const;

/*! obtains a begin iterator for the inner CCBs of `f`. */
Inner_ccb_iterator inner_ccbs_begin();

/*! obtains a past-the-end iterator for the inner CCBs of `f`. */
Inner_ccb_iterator inner_ccbs_end();

/*! obtains the number of holes (i.e., inner CCBs) inside `f`. */
size_t number_of_holes() const;

/*! obtains a begin-iterator for the holes (i.e., inner CCBs) of `f`. */
Hole_iterator holes_begin();

/*! obtains a past-the-end iterator for the holes (i.e., inner CCBs) of `f`. */
Hole_iterator holes_end();

/*! obtains the number of isolated vertices inside `f`. */
size_t number_of_isolated_vertices() const;

/*! obtains a begin-iterator for the isolated vertices inside `f`. */
Isolated_vertex_iterator isolated_vertices_begin();

/*! obtains a past-the-end iterator for the isolated vertices inside `f`. */
Isolated_vertex_iterator isolated_vertices_end();

/// @}

/// \name Modifiers
/// @{

/*! sets the face as unbounded (if `flag` is `true`), or as a bounded face
 * (if it is `false`).
 */
void set_unbounded(bool flag);

/*! sets the incident halfedge. */
void set_halfedge(Halfedge* e);

/*! adds `e` as an outer CCB of `f`. */
void add_outer_ccb(Halfedge* e);

/*! removes the outer CCB that `it` points to from `f`. */
void erase_outer_ccb(Outer_ccb_iterator it);

/*! adds `e` as an inner CCB of `f`. */
void add_inner_ccb(Halfedge* e);

/*! removes the inner CCB that `it` points to from `f`. */
void erase_inner_ccb(Inner_ccb_iterator it);

/*! adds `e` as a hole (i.e., inner CCB) of `f`. */
void add_hole(Halfedge* e);

/*! removes the hole (i.e., inner CCB) that `it` points to from `f`. */
void erase_hole(Hole_iterator it);

/*! adds `v` as an isolated vertex inside `f`. */
void add_isolated_vertex(Vertex* v);

/*! removes the isolated vertex that `it` points to from inside `f`. */
void erase_isolated_vertex(Isolated_vertex_iterator it);

/// @}

}; /* end ArrangementDcelFace */
