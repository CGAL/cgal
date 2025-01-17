/*!
\ingroup PkgHyperbolicSurfaceTriangulation2Concepts
\cgalConcept

This traits class must have a type for complex numbers.

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
          represents a complex number, model of
          `ComplexNumber`, over the field `HyperbolicSurfaceTraits_2::FT` for its real and
          imaginary parts.
  */
  typedef unspecified_type  Complex;
  /// @}
};
