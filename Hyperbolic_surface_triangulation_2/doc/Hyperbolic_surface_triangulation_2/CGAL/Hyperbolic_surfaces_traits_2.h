namespace CGAL{

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2TraitsClasses

\tparam HyperbolicTraitsClass must be a model of `HyperbolicDelaunayTriangulationTraits_2`.

\cgalModels{HyperbolicSurfacesTraits_2}
*/
template<class HyperbolicTraitsClass>
class Hyperbolic_surface_traits_2 : public HyperbolicTraitsClass {
};

}; // namespace CGAL
