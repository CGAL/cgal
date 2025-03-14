namespace CGAL {

/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2TraitsClasses

  \tparam HyperbolicTraits must be a model of `HyperbolicDelaunayTriangulationTraits_2`.

  \cgalModels{HyperbolicSurfaceTraits_2}
*/
template<class HyperbolicTraits>
class Hyperbolic_surface_traits_2 : public HyperbolicTraits {};

} // namespace CGAL
