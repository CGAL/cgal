/*!
\ingroup PkgHyperbolicSurfaceTriangulation2Concepts
\cgalConcept

\cgalRefines{HyperbolicDelaunayTriangulationTraits_2}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Hyperbolic_surface_traits_2}
\cgalHasModelsEnd
*/
class HyperbolicSurfacesTraits_2 {
public:
  /// \name Types
  /// @{
  /*!
          Field type.
  */
  typedef unspecified_type                          FT;
  /*!
          Represents a point in the Poincaré disk model of the hyperbolic plane.
  */
  typedef unspecified_type          Hyperbolic_point_2;
  /*!
          Represents a complex number over a field: must be a model of ComplexWithoutSqrt.
  */
  typedef unspecified_type                                   Complex;
  /// @}
};
