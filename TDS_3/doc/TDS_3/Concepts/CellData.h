
/*!
\ingroup PkgTDS3Concepts

\cgalConcept

The concept `CellData` describes the requirements on the type which
is used to mark some cells during modifications of the 3D triangulation data
structure.

\sa `TriangulationDataStructure_3`
\sa `TriangulationDSCellBase_3`
*/

class CellData
{
public:

  /*!
  Clear all data.
  */
  void clear();

  /*!
  Mark the cell as "in conflict".
  */
  void mark_in_conflict();

  /*!
  Mark the cell as "on boundary".
  */
  void mark_on_boundary();

  /*!
  Mark the cell as processed.
  */
  void mark_processed();

  /*!
  Returns `true` if the cell is not marked as "in conflict", "on boundary", or "processed", and `false` otherwise.
  */
  void is_clear();

  /*!
  Returns `true` if the cell is marked as "in conflict", `false` otherwise.
  */
  void is_in_conflict();

  /*!
  Returns `true` if the cell is marked as "on boundary", `false` otherwise.
  */
  void is_on_boundary();

  /*!
  Returns `true` if the cell is marked as processed, `false` otherwise.
  */
  bool processed() const;
};
