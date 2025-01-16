/*!
\ingroup PkgTetrahedralRemeshingConcepts
\cgalConcept

\cgalRefines{SimplicialMeshCellBase_3}

Cell base concept to be used in the triangulation type given to the function `CGAL::tetrahedral_isotropic_remeshing()`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3}
\cgalHasModelsEnd

*/
class RemeshingCellBase_3
{
  /*!
  These functions are used internally by tetrahedral remeshing.
  The class should provide storage, accessors and modificators for a cache value for sliverity.*/
  /// @{

  /*!
  */
  void set_sliver_value(double value);

  /*!
  */
  double sliver_value() const;

  /*!
  */
  bool is_cache_valid() const;

  /*!
  */
  void reset_cache_validity() const;

  /// @}
};
