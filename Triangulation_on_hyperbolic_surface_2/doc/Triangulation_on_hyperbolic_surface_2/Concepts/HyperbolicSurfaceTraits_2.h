/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2Concepts
  \cgalConcept

  This traits class must have a type for complex numbers.

  \cgalRefines{HyperbolicDelaunayTriangulationTraits_2}

  \cgalHasModelsBegin
  \cgalHasModels{CGAL::Hyperbolic_surface_traits_2}
  \cgalHasModelsEnd
*/
class HyperbolicSurfaceTraits_2
{
public:
  /// \name Types
  /// @{

  /*!
    Field number type.
  */
  typedef unspecified_type  FT;

  /*!
    represents a complex number, model of `ComplexNumber`,
    over the field `HyperbolicSurfaceTraits_2::FT` for its real and imaginary parts.
  */
  typedef unspecified_type  Complex;

  /// @}

  /// \name Distance computation
  /// @{
  /*!
    \return the hyperbolic cosine of the hyperbolic distance between u and v in the Poincaré disk.
   */
   template <typename Number = FT, typename Point = Hyperbolic_point_2>
     static Number cosh_hd(Point const& u, Hyperbolic_point_2 const& v);


  /// @}
};
