
/*!
\ingroup PkgAlphaShapes2Concepts
\cgalConcept

The concept `AlphaShapeFace_2` describes the requirements for the base face of an alpha shape.

\cgalRefines `TriangulationFaceBase_2`, if the underlying triangulation of the alpha shape is a Delaunay triangulation.
\cgalRefines `RegularTriangulationFaceBase_2`, if the underlying triangulation of the alpha shape is a regular triangulation.
\cgalRefines `Periodic_2TriangulationFaceBase_2`, if the underlying triangulation of the alpha shape is a periodic triangulation.

\cgalHasModel `CGAL::Alpha_shape_face_base_2` (templated with the appropriate triangulation face base class).

*/
class AlphaShapeFace_2 {
public:

/// \name Types
/// @{

/*!
A container type to get (and put) the three special values
(\f$ \alpha_1, \alpha_2, \alpha_3\f$) associated with an alpha shape edge.
*/
typedef unspecified_type Interval_3;

/*!
A coordinate type.
The type must provide a copy constructor, assignment, comparison
operators, negation, multiplication, division and allow the
declaration and initialization with a small integer constant
(cf. requirements for number types). An obvious choice would be
coordinate type of the point class
*/
typedef unspecified_type FT;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
AlphaShapeFace_2();

/*!
constructor setting the incident vertices.
*/
AlphaShapeFace_2(const Vertex_handle& v0, const Vertex_handle& v1, const Vertex_handle& v2);

/*!
constructor setting the incident vertices and the neighboring faces.
*/
AlphaShapeFace_2(const Vertex_handle& v0, const Vertex_handle& v1, const Vertex_handle& v2, const Face_handle& n0, const Face_handle& n1, const Face_handle& n2);

/// @}

/// \name Access Functions
/// @{

/*!
returns the interval associated with the edge indexed with \f$ i\f$, which contains
three alpha values
\f$ \alpha_1 \leq\alpha_2 \leq\alpha_3\f$, such as for
\f$ \alpha\f$ between \f$ \alpha_1\f$ and \f$ \alpha_2\f$, the edge indexed with \f$ i\f$ is
attached but singular,
for \f$ \alpha\f$ between \f$ \alpha_2\f$ and \f$ \alpha_3\f$, the edge is regular, and for \f$ \alpha\f$
greater than \f$ \alpha_3\f$, the edge is interior.
*/
Interval_3 get_ranges(const int& i);

/*!
return the alpha value, under which the alpha shape contains the
face.
*/
FT get_alpha();

/// @}

/// \name Modifiers
/// @{

/*!
sets the interval associated with the edge indexed with \f$ i\f$, which contains three
alpha values
\f$ \alpha_1 \leq\alpha_2 \leq\alpha_3\f$, such as for
\f$ \alpha\f$ between \f$ \alpha_1\f$ and \f$ \alpha_2\f$, the edge indexed with \f$ i\f$ is
attached but singular,
for \f$ \alpha\f$ between \f$ \alpha_2\f$ and \f$ \alpha_3\f$, the edge is regular, and for \f$ \alpha\f$
greater than \f$ \alpha_3\f$, the edge is interior.
*/
void set_ranges(const int& i, const Interval_3& V);

/*!
sets the alpha value, under which the alpha shape contains the
face.
*/
void set_alpha(FT A);

/// @}

}; /* end AlphaShapeFace_2 */

