/// \ingroup PkgTetrahedralRemeshingConcepts
/// \cgalConcept
///
/// The concept `RemeshingVertexBase_3` defines the requirements for the vertex base
/// used in the triangulation given as input to the remeshing algorithm
///
/// \cgalRefines `TriangulationVertexBase_3`, `CopyConstructible`
/// \cgalHasModel `CGAL::Tetrahedral_remeshing::Remeshing_vertex_base`.


class RemeshingVertexBase_3 {
public:

  /// @name Operations
  /// @{

  /// Returns the dimension of the lowest dimensional face of the input 3D
  /// complex that contains the vertex
  int in_dimension() const;

  /// Sets the dimension of the lowest dimensional face of the input 3D complex
  /// that contains the vertex
  void set_dimension(const int dimension);

  /// Returns the number of incident facets
  std::size_t number_of_incident_facets() const;

  /// Returns the number of subdomains to which belong incident cells
  std::size_t number_of_incident_subdomains() const;

  /// Internal function that invalidates cache data stored for performance
  void invalidate_cache();
  /// Internal function that sets cache data stored for performance
  void set_cache(const std::size_t i, const std::size_t j);


  /// @}
};
