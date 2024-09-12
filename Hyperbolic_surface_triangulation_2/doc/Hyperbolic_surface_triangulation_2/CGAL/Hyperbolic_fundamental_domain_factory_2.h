namespace CGAL{

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Factory class, whose only purpose is to construct random fundamental domains of
closed orientable hyperbolic surfaces.

The method `make_hyperbolic_fundamental_domain_g2` constructs such a domain for
a surface of genus 2.

\tparam Traits is the traits class and must be a model of `HyperbolicSurfaceTraits_2` (default model: `Hyperbolic_surface_traits_2`).
*/
template<class Traits>
class Hyperbolic_fundamental_domain_factory_2{
public:
  /// \name Creation
  /// @{
  /*!
    Constructor.
  */
  Hyperbolic_fundamental_domain_factory_2();
  /// @}

  /// \name Generation of a domain in genus 2.
  /// @{
  /*!
    randomly generates a convex domain of a closed orientable hyperbolic surface
   of genus two from a seed.
  */
  Hyperbolic_fundamental_domain_2<Traits> make_hyperbolic_fundamental_domain_g2(unsigned int seed);
  /// @}

};

}; // namespace CGAL
