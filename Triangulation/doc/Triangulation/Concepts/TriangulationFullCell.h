
/*!
\ingroup PkgTriangulationsConcepts
\cgalconcept

The concept `TriangulationFullCell` describes the requirements on the type used by the 
class `Triangulation<TriangulationTraits, TriangulationDataStructure>`, and its derived classes, to 
represent a full cell. 

\refines ::TriangulationDSFullCell 
We only list below the additional specific requirements of \refines ::TriangulationFullCell. 

\hasModel CGAL::Triangulation_full_cell<TriangulationTraits, TriangulationDSFullCell> 

Input/Output 
-------------- 

These operators can be used directly and are called by the I/O 
operator of class `Triangulation`. 

\sa `Triangulation_full_cell<TriangulationTraits, TriangulationDSFullCell>` 
\sa `TriangulationVertex` 
\sa `Triangulation<TriangulationTraits, TriangulationDataStructure>` 

*/

class TriangulationFullCell {
public:

CONVERROR: ccFunction inside class or concept, try to relate 
/*! 
Inputs additional information stored in the full cell. 
\relates TriangulationFullCell 
*/ 
istream & operator>>(istream & is, TriangulationFullCell & c); 

CONVERROR: ccFunction inside class or concept, try to relate 
/*! 
Outputs additional information stored in the full cell. 
\relates TriangulationFullCell 
*/ 
ostream & operator<<(ostream & os, const TriangulationFullCell & c); 

/// @}

}; /* end TriangulationFullCell */

