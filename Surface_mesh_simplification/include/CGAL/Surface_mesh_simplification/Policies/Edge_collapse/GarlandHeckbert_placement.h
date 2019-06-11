#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_PLACEMENT_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_PLACEMENT_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/Lindstrom_Turk_core.h>


#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_core.h>

namespace CGAL {
namespace Surface_mesh_simplification {

template<class TM_>
class GarlandHeckbert_placement
{
public:
  typedef TM_                                                            TM;


  typedef typename internal::GarlandHeckbertCore<TM> GHC;
  typedef typename GHC::garland_heckbert_map_type    garland_heckbert_map_type;
  typedef typename GHC::Matrix4x4                    Matrix4x4;
  typedef typename GHC::Row4                         Row4;
  typedef typename GHC::Col4                         Col4;
  typedef typename GHC::FT                           FT;

  typedef typename boost::optional<FT>               Optional_FT;

  GarlandHeckbert_placement(const garland_heckbert_map_type& aCostMatrices)
    : mCostMatrices(aCostMatrices)
  {}

  template <typename Profile>
  boost::optional<typename Profile::Point> operator()(const Profile& aProfile) const
  {
    const Matrix4x4 combinedMatrix = GHC::combine_matrices(
                  mCostMatrices.at(aProfile.v0()),
                  mCostMatrices.at(aProfile.v1())
                );

    const Col4 p0 = std::move(GHC::point_to_homogenous_column(aProfile.p0()));
    const Col4 p1 = std::move(GHC::point_to_homogenous_column(aProfile.p1()));

    const Col4 opt = GHC::optimal_point(combinedMatrix, p0, p1);

    boost::optional<typename Profile::Point> pt;
    pt = typename Profile::Point(opt(0) / opt(3), opt(1) / opt(3), opt(2) / opt(3));

    return pt;
  }

private:
  const garland_heckbert_map_type& mCostMatrices;
};
} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_PLACEMENT_H
