
/*!
\ingroup PkgTriangulation3Concepts
\cgalConcept

A type model of this concept can be used as cell base by a triangulation
and provides an additional information storage.

\cgalRefines{TriangulationCellBase_3}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Triangulation_cell_base_with_info_3}
\cgalHasModelsEnd

*/

class TriangulationCellBaseWithInfo_3 {
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
Returns a const reference to the object of type `Info` stored in the
vertex.
*/
const Info& info() const;

/*!
Returns a reference to the object of type `Info` stored in the vertex.
*/
Info& info();

/// @}

}; /* end TriangulationCellBaseWithInfo_3 */

