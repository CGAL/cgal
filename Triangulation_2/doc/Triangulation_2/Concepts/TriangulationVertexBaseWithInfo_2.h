
/*!
\ingroup PkgTriangulation2Concepts
\cgalConcept

A type model of this concept can be used as vertex base by a triangulation
and provides an additional information storage.

\cgalRefines `TriangulationVertexBase_2`

\cgalHasModel `CGAL::Triangulation_vertex_base_with_info_2`

*/

class TriangulationVertexBaseWithInfo_2 {
public:

/// \name Types
/// @{

/*!
A type which is `DefaultConstructible` and `Assignable`.
*/
typedef unspecified_type Info;

/// @}

/// \name Access Functions
/// @{

/*!
returns a const reference to the object of type `Info` stored in the
vertex.
*/
const Info& info() const;

/*!
returns a reference to the object of type `Info` stored in the vertex.
*/
Info& info();

/// @}

}; /* end TriangulationVertexBaseWithInfo_2 */

