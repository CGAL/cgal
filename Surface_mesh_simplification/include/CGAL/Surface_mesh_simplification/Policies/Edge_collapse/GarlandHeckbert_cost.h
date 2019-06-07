#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_COST_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_COST_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_core.h>

namespace CGAL {
namespace Surface_mesh_simplification {

template<class TM_>
class GarlandHeckbert_cost
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

  GarlandHeckbert_cost(const garland_heckbert_map_type& aCostMatrices)
    : mCostMatrices(aCostMatrices) {
  }

  template <typename Profile>
  boost::optional<typename Profile::FT>
  operator()(const Profile& aProfile,
             const boost::optional<typename Profile::Point>& aPlacement) const
  {
    Matrix4x4 combinedMatrix = GHC::combine_matrices(
                  mCostMatrices.at(aProfile.v0()),
                  mCostMatrices.at(aProfile.v1())
                );

    Col4 pt;
    pt << (*aPlacement).x(), (*aPlacement).y(), (*aPlacement).z(), 1;


    Optional_FT cost = (pt.transpose() * combinedMatrix * pt)(0,0);

    return cost;
  }

private:
  const garland_heckbert_map_type& mCostMatrices;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_COST_H
