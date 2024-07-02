/*!
\ingroup PkgTetrahedralRemeshingConcepts
\cgalConcept

Sizing field functional, to be used as second parameter
of the function `CGAL::tetrahedral_isotropic_remeshing()`,
used to control the size of the mesh elements.

This concept is
equivalent to `MeshDomainField_3`, used in \cgal \ref PkgMesh3 package,
so they can be used interchangeably.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Uniform_sizing_field}
\cgalHasModels{CGAL::Adaptive_remeshing_sizing_field}
\cgalHasModelsEnd

*/
class RemeshingSizingField_3
{
public:
  /// \name Types
  /// @{
  /*!
  Numerical type.
  */
  typedef unspecified_type FT;

  /*!
  Point type.
  */
  typedef unspecified_type Point_3;

  /*!
  %Index type for points. Must match the type
  `RemeshingVertexBase_3::Index`.
  */
  typedef unspecified_type Index;
  /// @}

  /*! \name Operations
  The field value may depend on the query point location and/or
  on the input subcomplex including the query point.
  */
  /// @{

  /*!
  returns the value of the sizing field at the point `p`,
  assumed to be included in the input subcomplex with dimension `dim`
  and mesh subcomplex index `index`.
  */
  template<typename Index>
  FT operator()(const Point_3& p, const int dim, const Index& index) const;

  /// @}
}; /*end RemeshingSizingField_3*/
