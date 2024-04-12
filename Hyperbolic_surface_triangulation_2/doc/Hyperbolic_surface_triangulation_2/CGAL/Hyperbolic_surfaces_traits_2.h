// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2TraitsClasses

Traits class offered by \cgal as a default model for the traits concept `HyperbolicSurfacesTraits_2`.

\tparam HyperbolicTraitsClass must be a model of `HyperbolicDelaunayTriangulationTraits_2`.

The default value for `HyperbolicTraitsClass` is `CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<CGAL::Circular_kernel_2<CGAL::Cartesian<CGAL::Gmpq>,CGAL::Algebraic_kernel_for_circles_2_2<CGAL::Gmpq>>>`,
 which provides exact constructions and predicates over rational numbers.

\cgalModels{HyperbolicSurfacesTraits_2}
*/
template<class HyperbolicTraitsClass>
class Hyperbolic_surfaces_traits_2 : public HyperbolicTraitsClass {
public:
  /// \name Types
  /// @{
  typedef typename HyperbolicTraitsClass::FT                          FT;
  typedef typename HyperbolicTraitsClass::Hyperbolic_point_2          Hyperbolic_point_2;
  typedef Complex_without_sqrt<FT>                                    Complex;
  /// @}
};


} //namespace CGAL
