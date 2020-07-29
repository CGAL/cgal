/*!
\ingroup PkgEnvelope2Concepts
\cgalConcept

An edge record in an envelope diagram, which represents a continuous portion
of the \f$ x\f$-axis. It is associated with a (possibly empty) set of curves that
induce the envelope over this portion of the \f$ x\f$-axis. Note that all curves
in this set overlap over the interval represented by the edge.

\sa `EnvelopeDiagram_1`
\sa `EnvelopeDiagramVertex`

*/

class EnvelopeDiagramEdge {
public:

/// \name Types
/// @{

/*!
the size type (convertible to `size_t`).
*/
typedef unspecified_type Size;

/*!
the corresponding diagram-vertex type.
*/
typedef unspecified_type Vertex;

/*!
the \f$ x\f$-monotone curve type.
*/
typedef unspecified_type X_monotone_curve_2;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
EnvelopeDiagramEdge();

/*!
copy constructor.
*/
EnvelopeDiagramEdge (const EnvelopeDiagramEdge& other);

/// @}

/// \name Access Functions
/// @{

/*!
returns the number of \f$ x\f$-monotone curves associated with `e`.
*/
Size number_of_curves () const;

/*!
returns whether `e` represents an empty interval - namely, whether the set of \f$ x\f$-monotone curves associated with it is empty.
*/
bool is_empty () const;

/*!
returns a representative \f$ x\f$-monotone curve associated with `e`.
\pre `e` does not represent an empty interval.
*/
const X_monotone_curve_2& curve () const;

/*!
returns an iterator for the first \f$ x\f$-monotone curve associated with `e`.
*/
Curve_const_iterator curves_begin () const;

/*!
returns a past-the-end iterator for the \f$ x\f$-monotone curves associated with `e`.
*/
Curve_const_iterator curves_end () const;

/*!
returns the vertex lying to `e`'s left.
\pre `e` is not the leftmost edge in the diagram.
*/
Vertex_const_handle left() const;

/*!
returns the vertex lying to `e`'s right.
\pre `e` is not the rightmost edge in the diagram.
*/
Vertex_const_handle right() const;

/// @}

/// \name Modifiers
/// @{

/*!
clears the set of curves associated with `e`.
*/
void clear_curves();

/*!
adds the \f$ x\f$-monotone curve `cv` to the set of curves associated with `e`.
*/
void add_curve (const X_monotone_curve_2& cv);

/*!
adds the given range of \f$ x\f$-monotone curves to the set of curves associated with `e`.
*/
void add_curves (Curve_const_iterator begin,
Curve_const_iterator end);

/*!
sets the vertex lying to the left of `e` to be `v`.
*/
void set_left (Vertex_const_handle v);

/*!
sets the vertex lying to the right of `e` to be `v`.
*/
void set_right (Vertex_const_handle v);

/// @}

}; /* end EnvelopeDiagramEdge */
