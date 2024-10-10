/*!
\ingroup PkgHyperbolicSurfaceTriangulation2Concepts
\cgalConcept

\cgalRefines{HyperbolicDelaunayTriangulationTraits_2}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Hyperbolic_surface_traits_2}
\cgalHasModelsEnd
*/
class HyperbolicSurfaceTraits_2 {
public:
  /// \name Types
  /// @{
  /*!
          represents a complex number over a field: must be a model of ComplexNumber.
  */
  typedef unspecified_type                                   Complex;
  /// @}
};
