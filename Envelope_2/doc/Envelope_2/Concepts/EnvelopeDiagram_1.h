/*!
\ingroup PkgEnvelope2Concepts
\cgalConcept

This concept defines the representation of an envelope diagram of a set 
of planar curve. The <I>envelope diagram</I> is a subdivision of the \f$ x\f$-axis 
into 0-dimensional cells (<I>vertices</I>) and 1-dimensional cells 
(<I>edges</I>), such that the identity of the curves that induce the lower 
envelope (or the upper envelope) over each cell is fixed. 

A vertex in an envelope diagram is therefore associated with a point 
on the envelope, and corresponds to either a curve endpoint 
or to an intersection point of two (or more) curves. Therefore each vertex 
is associated with a set of \f$ x\f$-monotone curves that induce the envelope 
over this point. Each vertex is incident to two edges, one lying to its 
left and the other to its right. 

An edge in the envelope diagram represents a continuous portion of the 
\f$ x\f$-axis, and is associated with a (possibly empty) set of curves that 
induce the envelope over this portion of the \f$ x\f$-axis. An edge may be bounded 
by two vertices, one to its left and the other to its right. However, the 
diagram contains two unbounded edges, its <I>leftmost</I> edge, representing 
the interval \f$ (-\infty, x_l)\f$, and its <I>rightmost</I> edge, representing the 
interval \f$ (x_r, \infty)\f$, where \f$ x_l\f$ and \f$ x_r\f$ are the \f$ x\f$-coodinates of 
the leftmost and the rightmost vertices in the diagram, respectively. 
Note that a diagram may contain no vertices at all, in which case it 
comprises a single edge. 

Note that any model of the `EnvelopeDiagram_1` concept must define a geometric 
traits class, which in turn defines the `Point_2` and 
`X_monotone_curve_2` types defined with the diagram features. 

\sa `EnvelopeDiagramVertex` 
\sa `EnvelopeDiagramEdge` 

*/

class EnvelopeDiagram_1 {
public:

/// \name Types 
/// @{

/*!
the geometric traits class. 
*/ 
typedef unspecified_type Traits_2; 

/*!
the point type. 
*/ 
typedef Traits_2::Point_2 Point_2; 

/*!
the \f$ x\f$-monotone curve type. 
*/ 
typedef Traits_2::X_monotone_curve_2 X_monotone_curve_2; 

/*!
the size type (convertible to `size_t`). 
*/ 
typedef unspecified_type Size; 

/*!
an iterator for the \f$ x\f$-monotone curves that induce a diagram feature, with value type `X_monotone_curve_2`. 
*/ 
typedef unspecified_type Curve_const_iterator; 

/*!
the vertex type, a model of the concept `EnvelopeDiagramVertex`. 
*/ 
typedef unspecified_type Vertex; 

/*!
the edge type, a model of the concept `EnvelopeDiagramEdge`. 
*/ 
typedef unspecified_type Edge; 

/*!
a handle to a diagram vertex. 
*/ 
typedef unspecified_type Vertex_handle; 

/*!
a non-mutable handle to a diagram vertex. 
*/ 
typedef unspecified_type Vertex_const_handle; 

/*!
a handle to a diagram edge. 
*/ 
typedef unspecified_type Edge_handle; 

/*!
a non-mutable handle to a diagram edge. 
*/ 
typedef unspecified_type Edge_const_handle; 

/// @} 

/// \name Creation 
/// @{

/*!
constructs an empty diagram containing one unbounded edge, 
which corresponds to the entire plane and has no \f$ x\f$-monotone 
curves that are associated with it. 
*/ 
EnvelopeDiagram_1(); 

/*!
copy constructor. 
*/ 
Envelope_diagram_1 (const Self& other); 

/// @} 

/// \name Access Functions 
/// @{

/*!
returns the leftmost edge of the diagram (a non-const version is also available). 
*/ 
Edge_const_handle leftmost() const; 

/*!
returns the rightmost edge of the diagram (a non-const version is also available). 
*/ 
Edge_const_handle rightmost() const; 

/// @} 

/// \name Modifiers 
/// @{

/*!
sets the leftmost edge of the diagram to be `e`. 
*/ 
void set_leftmost (Edge_const_handle e); 

/*!
sets the rightmost edge of the diagram to be `e`. 
*/ 
void set_rightmost (Edge_const_handle e); 

/*!
creates a new diagram vertex, associated with the point `p`. 
*/ 
Vertex_handle new_vertex (const Point_2& p); 

/*!
creates a new diagram edge. 
*/ 
Edge_handle new_edge (); 

/*!
deletes the given vertex `v`. 
*/ 
void delete_vertex (Vertex_handle v); 

/*!
deletes the given edge `e`. 
*/ 
void delete_edge (Edge_handle e); 

/// @}

}; /* end EnvelopeDiagram_1 */
