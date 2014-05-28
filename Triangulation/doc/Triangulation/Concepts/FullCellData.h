
/*!
\ingroup PkgTriangulationsConcepts
\cgalConcept

The concept `FullCellData` describes the requirements on the type which 
is used to mark some full cells, during modifications of the triangulation data
structure.

\sa `TriangulationDataStructure`
\sa `TriangulationDSFullCell`
*/

class FullCellData 
{
public:
/*!
Clear all data.
*/
void clear();
/*!
Mark the full cell as visited.
*/
void mark_visited();
/*!
Mark the full cell as not visited.
*/
void clear_visited();
/*!
Returns `true` if the full cell is <b>not</b> marked as visited, `false` otherwise.
*/
bool is_clear()   const;
/*!
Returns `true` if the full cell is marked as visited, `false` otherwise.
*/
bool is_visited() const;

}; /* end FullCellData */
