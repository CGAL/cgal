namespace CGAL{

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2TraitsClasses

\tparam HyperbolicTraitsClass must be a model of `HyperbolicDelaunayTriangulationTraits_2`.

\cgalModels{HyperbolicSurfaceTraits_2}
*/
template<class HyperbolicTraitsClass>
class Hyperbolic_surface_traits_2 : public HyperbolicTraitsClass {
public:
  /// \name Types
  /// @{
  /*!
          represents a complex number, model of
          `ComplexNumber`, over the field `HyperbolicSurfaceTraits_2::FT` for its real and
          imaginary parts.
  */
  typedef unspecified_type                                   Complex;
  /// @}
};

}; // namespace CGAL
