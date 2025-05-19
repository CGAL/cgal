
namespace CGAL {
namespace Surface_mesh_simplification {

/*
\ingroup PkgSurfaceMeshSimplificationRef

The class `Bounded_distance_placement` is a model for the concept `GetPlacement`.

This placement class is a wrapper around another (so-called <em>base</em>) placement class.
The position of a vertex resulting from the contraction of an edge is obtained by first querying
the base placement class, and checking whether this tentative position is not too far
(according to a user-provided distance bound) from the input mesh.
If it is too far, the position is rejected and no position is returned; otherwise,
the position is returned.

\tparam BasePlacement must be a model of `GetPlacement`.
\tparam GeomTraits must be a model of `Kernel` and be identical to the traits specified
                  in the named parameters of the function `edge_collapse()` (if specified).

The distance check is performed using an AABB tree and this class thus depends on the package \ref PkgAABBTree.

\cgalModels{GetPlacement}

*/
template<class BasePlacement, class GeomTraits>
class Bounded_distance_placement
  : public BasePlacement
{
public:
  //
  typedef typename GeomTraits::FT FT;

  // \name Creation
  //
  // @{

  // The distance bound `d` is used to control that during simplification,
  // no vertex has a distance to the input that would be greater than `d`.
  Bounded_distance_placement(const FT d, const BasePlacement& base_placement = BasePlacement());

  // @}

  // \name Operations
  // @{

  // Returns the placement computed by `base_placement`, provided the distance between the input
  // and this placement is smaller than `d`. Otherwise, nothing is returned.
  std::optional<typename Edge_profile::Point> operator()(const Edge_profile& profile) const;

  // @}
};

} // namespace Surface_Mesh_Simplification
} // namespace CGAL
