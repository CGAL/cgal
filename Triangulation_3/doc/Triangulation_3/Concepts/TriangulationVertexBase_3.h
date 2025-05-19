
/*!
\ingroup PkgTriangulation3Concepts
\cgalConcept

The vertex base used by the geometric triangulation must store a point.
We list here the additional requirements compared to a vertex base usable
for the triangulation data structure.

\cgalRefines{TriangulationDSVertexBase_3}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Triangulation_vertex_base_3}
\cgalHasModels{CGAL::Triangulation_vertex_base_with_info_3}
\cgalHasModelsEnd

\sa `TriangulationCellBase_3`

*/

class TriangulationVertexBase_3 {
public:

/// \name Types
/// @{

/*!
Must be the same as the point type `TriangulationTraits_3::Point_3`
defined by the geometric traits class of the triangulation.
*/
typedef unspecified_type Point;

/// @}

/// \name Creation
/// @{

/*!
Constructs a vertex whose geometric embedding is point `p`.
*/
TriangulationVertexBase_3(const Point & p);

/*!
Constructs a vertex embedding the point `p` and pointing to cell `c`.
*/
TriangulationVertexBase_3(const Point & p, Cell_handle c);

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
istream& operator>> (istream& is, TriangulationVertexBase_3 & v);

/*!
Outputs the non-combinatorial information given by the vertex: the
point and other possible information.
*/
ostream& operator<< (ostream& os, const TriangulationVertexBase_3 & v);

/// @}

}; /* end TriangulationVertexBase_3 */

