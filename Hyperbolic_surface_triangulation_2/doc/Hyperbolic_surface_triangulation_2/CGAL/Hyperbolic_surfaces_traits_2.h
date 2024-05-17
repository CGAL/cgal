// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL{

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2TraitsClasses

\tparam HyperbolicTraitsClass must be a model of `HyperbolicDelaunayTriangulationTraits_2`.

\cgalModels{HyperbolicSurfacesTraits_2}
*/
template<class HyperbolicTraitsClass>
class Hyperbolic_surfaces_traits_2 : public HyperbolicTraitsClass {
};

}; // namespace CGAL
