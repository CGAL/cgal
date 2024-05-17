// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL{

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Factory class, whose purpose is to generate some convex domains of surfaces of genus two.

\tparam Traits is the traits class and must be a model of `HyperbolicSurfacesTraits_2` (default model : `Hyperbolic_surfaces_traits_2`).
*/
template<class Traits>
class Hyperbolic_fundamental_domain_factory_2{
public:
  /// \name Creation
  /// @{
  /*!
    Default constructor, the seed is used for generation.
  */
  Hyperbolic_fundamental_domain_factory_2(unsigned int seed);
  /// @}

  /// \name Generation of a domain
  /// @{
  /*!
    Randomly generates a convex domain of a closed orientable hyperbolic surface of genus two.
  */
  Hyperbolic_fundamental_domain_2<Traits> generate_domain_g2();
  /// @}

};

}; // namespace CGAL
