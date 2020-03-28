
/*!
\ingroup PkgTriangulation3Concepts
\cgalConcept

The base vertex of a regular triangulation must be a model of the concept
`RegularTriangulationVertexBase_3`, which refines the concept `TriangulationDSVertexBase_3`
by adding a geometric point member.

\cgalRefines `TriangulationDSVertexBase_3`

\cgalHasModel `CGAL::Regular_triangulation_vertex_base_3`

\sa `RegularTriangulationCellBase_3`

*/

class RegularTriangulationVertexBase_3 {
public:

/// \name Types
/// @{

/*!
Must be the same as the point type `RegularTriangulationTraits_3::Weighted_point_3`
defined by the geometric traits class of the triangulation.
*/
typedef unspecified_type Point;

/// @}

/// \name Creation
/// @{

/*!
Constructs a vertex whose geometric embedding is point `p`.
*/
RegularTriangulationVertexBase_3(const Point & p);

/*!
Constructs a vertex embedding the point `p` and pointing to cell `c`.
*/
RegularTriangulationVertexBase_3(const Point & p, Cell_handle c);

/// @}

/// \name Access Functions
/// @{

/*!
Returns the point.
*/
Point point() const;

/// @}

/// \name Setting
/// @{

/*!
Sets the point.
*/
void set_point(Point p);

/// @}

/// \name I/O
/// @{

/*!
Inputs the non-combinatorial information given by the vertex:
the point and other possible information.
*/
istream& operator>> (istream& is, RegularTriangulationVertexBase_3 & v);

/*!
Outputs the non-combinatorial information given by the vertex: the
point and other possible information.
*/
ostream& operator<< (ostream& os, const RegularTriangulationVertexBase_3 & v);

/// @}

}; /* end RegularTriangulationVertexBase_3 */

