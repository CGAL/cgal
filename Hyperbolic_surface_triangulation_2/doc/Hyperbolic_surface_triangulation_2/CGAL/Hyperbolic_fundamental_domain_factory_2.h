// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Factory class, whose only purpose is to randomly generate some convex domains of closed orientable hyperbolic surfaces of genus two.

\tparam GeometricTraits_2 is the geometric traits class and must be a model of `HyperbolicSurfacesTraits_2`.
*/
template<class GeometricTraits_2>
class Hyperbolic_fundamental_domain_factory_2{
public:
  /// \name Types
  /// @{
  typedef GeometricTraits_2                                 Geometric_traits_2;
  /// @}

  /// \name Creation
  /// @{
  /*!
    Default constructor. The seed is used for random generation
  */
  Hyperbolic_fundamental_domain_factory_2(unsigned int seed);
  /// @}

  /// \name Generation of a domain
  /// @{
  /*!
    Randomly generates a convex domain of a closed orientable hyperbolic surface of genus two.
  */
  Hyperbolic_fundamental_domain_2<GeometricTraits_2> generate_domain_g2();
  /// @}

};

} //namespace CGAL
