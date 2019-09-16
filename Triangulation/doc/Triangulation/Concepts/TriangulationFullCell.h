
/*!
\ingroup PkgTriangulationsConcepts
\cgalConcept

The concept `TriangulationFullCell` describes the requirements on the type used by the
class `CGAL::Triangulation<TriangulationTraits_, TriangulationDataStructure_>`, and its derived classes, to
represent a full cell.

\cgalRefines `TriangulationDSFullCell` We only list below the
additional specific requirements of `TriangulationFullCell`.

\cgalHasModel `CGAL::Triangulation_full_cell<TriangulationTraits_, TriangulationDSFullCell_>`

\sa `CGAL::Triangulation_full_cell<TriangulationTraits_, Data, TriangulationDSFullCell_>`
\sa `TriangulationVertex`
\sa `CGAL::Triangulation<TriangulationTraits_, TriangulationDataStructure_>`

*/

class TriangulationFullCell {
public:
/// \name Input/Output
/// These operators can be used directly and are called by the I/O
/// operator of class `Triangulation`.
/// @{

/*!
Inputs additional information stored in the full cell.
*/
std::istream & operator>>(std::istream & is, TriangulationFullCell & c);

/*!
Outputs additional information stored in the full cell.
*/
std::ostream & operator<<(std::ostream & os, const TriangulationFullCell & c);

/// @}

}; /* end TriangulationFullCell */
