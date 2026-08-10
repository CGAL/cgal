/*!
 * \ingroup PkgEnvelope2Concepts
 * \cgalConcept
 *
 * This concept defines the representation of an envelope diagram of a set of
 * planar curve. The <I>envelope diagram</I> is a subdivision of the \f$
 * x\f$-axis into 0-dimensional cells (<I>vertices</I>) and 1-dimensional cells
 * (<I>edges</I>), such that the identity of the curves that induce the lower
 * envelope (or the upper envelope) over each cell is fixed.
 *
 * \cgalRefines{ArrangementOnCurve_1}
 *
 * A vertex in an envelope diagram is therefore associated with a point on the
 * envelope, and corresponds to either a curve endpoint or to an intersection
 * point of two (or more) curves. Therefore each vertex is associated with a set
 * of \f$ x\f$-monotone curves that induce the envelope over this point. Each
 * vertex is incident to two edges, one lying to its left and the other to its
 * right.
 *
 * An edge in the envelope diagram represents a continuous portion of the
 * \f$x\f$-axis, and is associated with a (possibly empty) set of curves that
 * induce the envelope over this portion of the \f$x\f$-axis. An edge may be
 * bounded by two vertices, one to its left and the other to its right. However,
 * the diagram contains two unbounded edges, its <I>leftmost</I> edge,
 * representing the interval \f$(-\infty, x_l)\f$, and its <I>rightmost</I>
 * edge, representing the interval \f$(x_r, \infty)\f$, where \f$x_l\f$ and
 * \f$x_r\f$ are the \f$ x\f$-coordinates of the leftmost and the rightmost
 * vertices in the diagram, respectively.  Note that a diagram may contain no
 * vertices at all, in which case it comprises a single edge.
 *
 * Note that any model of the `EnvelopeDiagram_1` concept must define a
 * geometric traits class, which in turn defines the `Point_2` and
 * `X_monotone_curve_2` types defined with the diagram features.
 */

class EnvelopeDiagram_1 {
public:

/// \name Types
/// @{

/*! the geometric traits class.
 */
typedef unspecified_type Traits_2;

/*! the point type.
 */
typedef Traits_2::Point_2 Point_2;

/*! the \f$ x\f$-monotone curve type.
 */
typedef Traits_2::X_monotone_curve_2 X_monotone_curve_2;

/*! the size type (convertible to `size_t`).
 */
typedef unspecified_type Size;

/*! an iterator for the \f$ x\f$-monotone curves that induce a diagram feature, with value type `X_monotone_curve_2`.
 */
typedef unspecified_type Curve_const_iterator;

/*! a descriptor of a diagram vertex.
 */
typedef unspecified_type Vertex_descriptor;

/*! a non-mutable descriptor of a diagram vertex.
 */
typedef unspecified_type Vertex_const_descriptor;

/*! a descriptor of a diagram edge.
 */
typedef unspecified_type Edge_descriptor;

/*! a non-mutable descriptor of a diagram edge.
 */
typedef unspecified_type Edge_const_descriptor;

/// @}

/// \name Creation
/// @{

/*! constructs an empty diagram containing one unbounded edge,
 * which corresponds to the entire plane and has no \f$x\f$-monotone
 * curves that are associated with it.
 */
EnvelopeDiagram_1();

/*! copy constructor.
 */
Envelope_diagram_1(const Self& other);

/// @}

/// \name Accessors
/// @{

/*! obtains the leftmost edge of the diagram. (A non-const version is also available.)
 */
Edge_const_descriptor leftmost() const;

/*! obtains the rightmost edge of the diagram. (A non-const version is also available.)
 */
Edge_const_descriptor rightmost() const;

/// @}

/// \name Modifiers
/// @{

/*! sets the leftmost edge of the diagram to be `e`.
 */
void set_leftmost(Edge_const_descriptor e);

/*! sets the rightmost edge of the diagram to be `e`.
 */
void set_rightmost(Edge_const_descriptor e);

/*! creates a new diagram vertex associated with a given point.
 */
Vertex_descriptor new_vertex(const Point_2& p);

/*! creates a new diagram edge.
 */
Edge_descriptor new_edge();

/*! deletes the given vertex.
 */
void delete_vertex(Vertex_descriptor v);

/*! deletes the given edge.
 */
void delete_edge(Edge_descriptor e);

/*! clears the diagram; leavs only one unbounded edge.
 */
void clear();

/// @}

/// \name Curve Data Accessors
/// @{

/*! obtains the point associated with a given vertex.
 */
const Point_2& point(Vertex_const_descriptor v) const;

/*! obtains the diagram curves associated with a given vertex. (A non-const version is also available.)
 */
const Curve_container& vertex_curves(Vertex_descriptor v) const;

/*! obtains the diagram curves associated with a given edge. (A non-const version is also available.)
 */
const Curve_container& edge_curves(Edge_descriptor e) const;

/*! obtains the number of curves associated with a given vertex.
 */
Size number_of_vertex_curves(Vertex_const_descriptor v) const;

/*! obtains the number of curves associated with a given edge.
 */
Size number_of_edge_curves(Edge_const_descriptor e) const;

/*! determines whether there are no curves associated with a given vertex.
 */
bool empty_vertex_curves(Vertex_const_descriptor v) const;

/*! determines whether there are no curves associated with a given edge.
 */
bool empty_edge_curves(Edge_const_descriptor e) const;

/*! obtains the first curve assocuayed with a given vertex.
 * \pre `empty_vertex_curves(v)` is `false`.
 */
const X_monotone_curve_2& vertex_curve(Vertex_const_descriptor v) const;

/*! obtains the first curve associated with a given edge.
 * \pre `empty_vertex_curves(v)` is `false`.
 */
const X_monotone_curve_2& edge_curve(Edge_const_descriptor e) const;

/// @}

/// \name Curve Data Modifiers
/// @{

/*! sets the point associated with a fiven vertex.
 */
void set_point(Vertex_descriptor v, const Point_2& p);

/*! adds an \f$x\f$-monotone curve to the list of curves associated with a given vertex.
 */
void add_vertex_curve(Vertex_descriptor v, const X_monotone_curve_2& xcv);

/*! adds an \f$x\f$-monotone curve to the list of curves associated with a given edge.
 */
void add_edge_curve(Edge_descriptor e, const X_monotone_curve_2& xcv);

/*! adds a range of \f$x\f$-monotone curves to the list of curves associated with a given vertex.
 */
void add_edge_curves(Edge_descriptor e, Curve_const_iterator begin, Curve_const_iterator end);

/*! adds a range of \f$x\f$-monotone curves to the list of curves associated with a given edge.
 */
void add_vertex_curves(Vertex_descriptor v, Curve_const_iterator begin, Curve_const_iterator end);

/*! clears the curves associated with a given vertex.
 */
void clear_vertex_curves(Vertex_descriptor v);

/*! clears the curves associated with a given edge.
 */
void clear_edge_curves(Edge_descriptor e);

/// @}

}; /* end EnvelopeDiagram_1 */
