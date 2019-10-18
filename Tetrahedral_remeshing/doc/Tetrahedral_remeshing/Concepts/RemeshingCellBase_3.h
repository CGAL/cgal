/// \ingroup PkgTetrahedralRemeshingConcepts
/// \cgalConcept
///
/// The concept `RemeshingCellBase_3` defines the requirements for the cell base
/// used in the triangulation given as input to the remeshing algorithm
///
/// \cgalRefines `TriangulationCellBaseWithInfo_3`, `CopyConstructible`
/// \cgalHasModel `CGAL::Tetrahedral_remeshing::Remeshing_cell_base`.


class RemeshingCellBase_3 {
public:
  /// Subdomain index
  typedef unspecified_type Subdomain_index;
  /// Surface patch index
  typedef unspecified_type Surface_patch_index;

  /// @name Operations
  /// @{
    /// Returns the index of the input subdomain of the triangulation
    /// that contains the cell.
    /// Cells with a non-zero `Subdomain_index` are considered as the "inside"
    /// of the domain to be remeshed
    const Subdomain_index& subdomain_index() const;

    /// Sets the subdomain index of the cell. 
    void set_subdomain_index(const Subdomain_index& si);

    /// returns `Surface_patch_index` of facet `i`. 
    const Surface_patch_index surface_patch_index(const int&) const;

    /// sets `Surface_patch_index` of facet `i` to `index`
    void set_surface_patch_index(const int i, const Surface_patch_index&)

    /// Returns `true` if the facet `i` lies on a surface patch
    bool is_facet_on_surface(const int& i) const;

  /// @}
};
