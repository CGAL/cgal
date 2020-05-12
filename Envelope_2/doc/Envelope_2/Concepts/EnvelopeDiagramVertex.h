/*!
\ingroup PkgEnvelope2Concepts
\cgalConcept

A vertex record in an envelope diagram. It is always associated with a point
on the lower (upper) envelope of a non-empty set of curves. A vertex is also
associated with a set of \f$ x\f$-monotone curves that induce the envelope
over this point. It is incident to two edges, one lying to its
left and the other to its right.

\sa `EnvelopeDiagram_1`
\sa `EnvelopeDiagramEdge`

*/

class EnvelopeDiagramVertex {
public:

/// \name Types
/// @{

/*!
the size type (convertible to `size_t`).
*/
typedef unspecified_type Size;

/*!
the corresponding diagram-edge type.
*/
typedef unspecified_type Edge;

/*!
the point type associated with the vertex.
*/
typedef unspecified_type Point_2;

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
EnvelopeDiagramVertex();

/*!
copy constructor.
*/
EnvelopeDiagramVertex (const EnvelopeDiagramVertex& other);

/*!
constructs a vertex associated with the given point `p`.
*/
EnvelopeDiagramVertex (const Point_2& p);

/// @}

/// \name Access Functions
/// @{

/*!
returns the point associated with `v`.
*/
const Point_2& point () const;

/*!
returns the number of \f$ x\f$-monotone curves associated with `v`.
*/
Size number_of_curves () const;

/*!
returns an iterator for the first \f$ x\f$-monotone curve associated with `v`.
*/
Curve_const_iterator curves_begin () const;

/*!
returns a past-the-end iterator for the \f$ x\f$-monotone curves associated with `v`.
*/
Curve_const_iterator curves_end () const;

/*!
returns the edge lying to `v`'s left.
*/
Edge_const_handle left() const;

/*!
returns the edge lying to `v`'s right.
*/
Edge_const_handle right() const;

/// @}

/// \name Modifiers
/// @{

/*!
associates the point `p` with `v`.
*/
void set_point (const Point_2& p);

/*!
clears the set of curves associated with `v`.
*/
void clear_curves();

/*!
adds the \f$ x\f$-monotone curve `cv` to the set of curves associated with `v`.
*/
void add_curve (const X_monotone_curve_2& cv);

/*!
adds the given range of \f$ x\f$-monotone curves to the set of curves associated with `v`.
*/
void add_curves (Curve_const_iterator begin,
Curve_const_iterator end);

/*!
sets the edge lying to the left of `v` to be `e`.
*/
void set_left (Edge_const_handle e);

/*!
sets the edge lying to the right of `v` to be `e`.
*/
void set_right (Edge_const_handle e);

/// @}

}; /* end EnvelopeDiagramVertex */
