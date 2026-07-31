namespace CGAL {

/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2TraitsClasses

  \tparam HyperbolicTraits must be a model of `HyperbolicDelaunayTriangulationTraits_2`.

  \cgalModels{HyperbolicSurfaceDelaunayTraits_2}
*/
template<class HyperbolicTraits>
class Hyperbolic_surface_Delaunay_traits_2 : public HyperbolicTraits {};

/// \name Approximation
  /// @{
  /*!
    This model of for the constrution of approximate hyperpolic circumcenter works as follows.
    If the number type of the point coordinates is `CGAL::Gmpq, the approximation precision is chosen via the function `approximation_precision()`. If another number type is used, the approximation uses the function `CGAL::to_double()` and there are no general guarantees whatsoever. When the number type is ``CGAL::Gmpq`, the integer `p` given to the function `set_circumcenter_approximation_precision()` means that computations are done with p times the double precision. More precisely each coordinate of an exact circumcenter is rounded to a fixed precision floating-point number of type CGAL::Gmpfr before being converted to a CGAL::Gmpq. The precision of the CGAL::Gmpfr numbers involved is p times 53 bits, where 53 bits is the precision of a double.
   */
   Construct_approximate_hyperbolic_circumcenter_2
     construct_approximate_hyperbolic_circumcenter_2_object() const;

 /*!   \return the approximation precision.  */
   unsigned approximation_precision() const {return approx_precision;}
   /*!
     If the number type of the point coordinates is `CGAL::Gmpq`, this function sets the precision for approximations to `p`. In particular, it is used in the function `CGAL::HyperbolicSurfaceDelaunayTraitsClass::Construct_approximate_hyperbolic_circumcenter_2` where the precision of computations is `p`x 53 bits, where 53 is the precision of a double.
  */
 void approximation_precision(unsigned p) {approx_precision = p;}

  /// @}
} // namespace CGAL
