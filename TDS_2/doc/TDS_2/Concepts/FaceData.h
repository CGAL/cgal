
/*!
\ingroup PkgTDS2Concepts

\cgalConcept

The concept `FaceData` describes the requirements on the type which
is used to mark some faces during modifications of the 2D triangulation data
structure.

\sa `TriangulationDataStructure_2`
\sa `TriangulationDSFaceBase_2`
*/

class FaceData
{
public:

  /*!
  Clear all data.
  */
  void clear();

  /*!
  Mark the face as "in conflict".
  */
  void mark_in_conflict();

  /*!
  Mark the face as "on boundary".
  */
  void mark_on_boundary();

  /*!
  Mark the face as processed.
  */
  void mark_processed();

  /*!
  Returns `true` if the face is not marked as "in conflict", "on boundary", or "processed", and `false` otherwise.
  */
  void is_clear();

  /*!
  Returns `true` if the face is marked as "in conflict", `false` otherwise.
  */
  void is_in_conflict();

  /*!
  Returns `true` if the face is marked as "on boundary", `false` otherwise.
  */
  void is_on_boundary();

  /*!
  Returns `true` if the face is marked as processed, `false` otherwise.
  */
  bool processed() const;
};
