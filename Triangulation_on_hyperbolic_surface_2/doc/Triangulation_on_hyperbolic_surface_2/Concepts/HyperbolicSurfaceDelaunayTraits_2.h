/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2Concepts
  \cgalConcept


  \cgalRefines{HyperbolicSurfaceTriangulationTraits_2}

  \cgalHasModelsBegin
  \cgalHasModels{CGAL::Hyperbolic_surface_Delaunay_traits_2}
  \cgalHasModelsEnd
*/
class HyperbolicSurfaceDelaunayTraits_2
{
public:
 /// \name Construction Types
  /// @{

   /*! A constructor object. Must provide the function operator
     `Hyperbolic_Voronoi_point_2 operator()(Hyperbolic_point_2 p,
     Hyperbolic_point_2 q, Hyperbolic_point_2 r),` which constructs an
     approximation of the hyperbolic circumcenter of the triangle with vertices
     `p, q`, and `r`. The approximation is tuned by the function
     `approximation_precision()`, if its value is `prec`, the coordinates are
     approximations at relative precision \f$ 2^{-53\times prec}\f$ of the exact
     circumcenter.
        */
  typedef unspecified_type     Construct_approximate_hyperbolic_circumcenter_2;
	/// @}

/// \name Operations
  /// @{
  Construct_approximate_hyperbolic_circumcenter_2
    construct_approximate_hyperbolic_circumcenter_2_object() const;

  /*!   sets the approximation precision to `p`.   */
  void approximation_precision(unsigned p);
  /*!    gets the approximation precision   */
   unsigned approximation_precision() const;

   /// @}

};
