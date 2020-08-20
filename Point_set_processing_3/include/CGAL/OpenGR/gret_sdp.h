#ifndef CGAL_OPENGR_GRET_SDP_H
#define CGAL_OPENGR_GRET_SDP_H

#include <CGAL/license/Point_set_processing_3.h>

#if defined(CGAL_LINKED_WITH_OPENGR) || defined(DOXYGEN_RUNNING)

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/range/iterator_range.hpp>
#include "boost/tuple/tuple.hpp"

#include <gr/algorithms/GRET_SDP.h>
#include <gr/accelerators/gret_sdp/wrappers.h>
#include <gr/shared.h>
#include <memory>

#include <memory>

namespace CGAL {

namespace OpenGR {

template <typename Kernel, class PointRange, class CorrespondencesRange, class NamedParameters, class TransformRange>
void computeRegistrationTransformations(const std::vector<PointRange>& point_clouds, const CorrespondencesRange& correspondences, 
                                        const NamedParameters& np, TransformRange& transformations)
{
    using Scalar = typename Kernel::FT;
    using GR_PointType = gr::Point3D<Scalar> ;
    using GR_IndexedPointType = std::pair<GR_PointType, std::size_t>;
    using GR_MatcherType = gr::GRET_SDP<GR_PointType, gr::DummyTransformVisitor, gr::GRET_SDP_Options>;
    using GR_Options = typename GR_MatcherType::OptionsType;
    using GR_PatchType = std::vector<GR_IndexedPointType>;
    using GR_VectorType = typename GR_MatcherType::VectorType;
    
    using PointMap = typename CGAL::GetPointMap<PointRange, NamedParameters>::type;

    using parameters::choose_parameter;
    using parameters::get_parameter;
    
    const int num_point_clouds = point_clouds.size();
    const int num_global_coordinates = correspondences.size();

    PointMap point_map = choose_parameter(get_parameter(np, internal_np::point_map), PointMap());

    // compute patches
    std::vector<GR_PatchType> patches(num_point_clouds);
    for (size_t i = 0; i < num_global_coordinates; i++){
        for(const auto& correspondence : correspondences[i]){
            const auto& p = get (point_map, point_clouds[correspondence.first][correspondence.second]);
            GR_PointType out(p.x(), p.y(), p.z());
            patches[correspondence.first].emplace_back(out ,i);
        }
    }
    
    // initialize matcher
    gr::Utils::Logger logger(gr::Utils::Verbose);
    GR_Options options;
    GR_MatcherType matcher(options, logger);

    // register patches
    gr::DummyTransformVisitor tr_visitor;
    matcher.template RegisterPatches<gr::MOSEK_WRAPPER<Scalar>>(patches, num_global_coordinates, tr_visitor);
    
    // get transformations
    using GR_TrafoType = typename GR_MatcherType::MatrixType;
    std::vector<GR_TrafoType> gr_transformations;
    matcher.getTransformations(gr_transformations);

    // convert to cgal affin transformations
    transformations.reserve(num_point_clouds);
    for(const GR_TrafoType& gr_trafo : gr_transformations)
        transformations.emplace_back(
            gr_trafo.coeff(0,0), gr_trafo.coeff(0,1), gr_trafo.coeff(0,2), gr_trafo.coeff(0,3),
            gr_trafo.coeff(1,0), gr_trafo.coeff(1,1), gr_trafo.coeff(1,2), gr_trafo.coeff(1,3),
            gr_trafo.coeff(2,0), gr_trafo.coeff(2,1), gr_trafo.coeff(2,2), gr_trafo.coeff(2,3)
        );
}

template <typename Kernel, class PointRange, class CorrespondencesRange, class NamedParameters>
void registerPointClouds(   const std::vector<PointRange>& point_clouds, const CorrespondencesRange& correspondences, 
                            const NamedParameters& np, PointRange& registered_points)
{
    // compute registration transformations
    std::vector<typename Kernel::Aff_transformation_3> transformations;
    computeRegistrationTransformations<Kernel>(point_clouds, correspondences, np, transformations);

    // add transformed point clouds to registered_points
    for (size_t i = 0; i < point_clouds.size(); i++){
        for (size_t j = 0; j < point_clouds[i].size(); j++){
            registered_points.push_back(point_clouds[i][j].transform(transformations[i]));
        }
    }
}

} } // end of namespace CGAL::OpenGR

#endif // CGAL_LINKED_WITH_OPENGR

#endif // CGAL_OPENGR_GRET_SDP_H