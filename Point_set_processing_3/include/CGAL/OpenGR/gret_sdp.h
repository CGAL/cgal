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



/**
   \ingroup PkgPointSetProcessing3Algorithms

   Computes the registration of the point clouds in `point_clouds` using the correspondences provided 
   in `correspondences` and stores the corresponding affine transformations in `transformations`.

   Registration is computed using the GRET-SDP algorithm that is described in [this paper](https://arxiv.org/abs/1306.5226)
   
   \note This function requires the \ref thirdpartyOpenGR library.

   \warning 

   \tparam PointRange is a model of `Range`. The value type of its iterator is
   the key type of the named parameter `point_map` in `NamedParameters`.
   \tparam CorrespondencesRange is a model of `Range`. The value type of its iterator is
   is another range of correspondences. The value type of the iterator of a correspondences range 
   is `std::pair<size_t,size_t>`.
   \tparam TransformRange is a model of `Range`. The value type of its iterator is `Kernel::Aff_transformation_3`
   where `Kernel` is the same kernel used to construct the value type of PointRange `Kernel_traits<typename PointRange::value_type>::Kernel`. ???
   The value type of its iterator could be any class with a constructor that accepts 12 doubles, ints,... ???

   \param point_clouds vector of input point ranges of the point clouds to be registered.
   \param correspondences input range of correspondences. Each entry in correspondences is a range itself that stores the indexes of the points that are 
   considered to correspond to the same global coordinate. A points indexes are stores using a `std::pair<size_t,size_t>` where `first` equals the point cloud index 
   and `second` equals the index within this point cloud that it belongs to (indexes with respect to `point_clouds`).  
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.
   \param transformations ouput range of the affine transformations where the i'th transformation should be applied to the i'th point cloud in order to register them.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` whose key type
     is the value type of the iterator of `PointRange` and whose value type is
     `geom_traits::Point_3`.  If this parameter is omitted,
     `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
   \cgalNamedParamsEnd

*/
template <class PointRange, class CorrespondencesRange, class NamedParameters, class TransformRange>
void compute_registration_transformations(const std::vector<PointRange>& point_clouds, const CorrespondencesRange& correspondences, 
                                        const NamedParameters& np, TransformRange& transformations)
{
    using PointMap = typename CGAL::GetPointMap<PointRange, NamedParameters>::type;
    using Kernel = typename CGAL::Kernel_traits<typename PointMap::value_type>::Kernel;
    using Scalar = typename Kernel::FT;
    using GR_PointType = gr::Point3D<Scalar> ;
    using GR_IndexedPointType = std::pair<GR_PointType, std::size_t>;
    using GR_MatcherType = gr::GRET_SDP<GR_PointType, gr::DummyTransformVisitor, gr::GRET_SDP_Options>;
    using GR_Options = typename GR_MatcherType::OptionsType;
    using GR_PatchType = std::vector<GR_IndexedPointType>;
    using GR_VectorType = typename GR_MatcherType::VectorType;
    

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


/**
   \ingroup PkgPointSetProcessing3Algorithms

   Computes the registration of the point clouds in `point_clouds` using the correspondences provided 
   in `correspondences` and stores the registered point cloud in `registered_points`.

   Registration is computed using the GRET-SDP algorithm that is described in [this paper](https://arxiv.org/abs/1306.5226)
   
   \note This function requires the \ref thirdpartyOpenGR library.

   \warning 

   \tparam PointRange is a model of `Range`. The value type of its iterator is
   the key type of the named parameter `point_map` in `NamedParameters`.
   \tparam CorrespondencesRange is a model of `Range`. The value type of its iterator is
   is another range of correspondences. The value type of the iterator of a correspondences range 
   is `std::pair<size_t,size_t>`.

   \param point_clouds vector of input point ranges of the point clouds to be registered.
   \param correspondences input range of correspondences. Each entry in correspondences is a range itself that stores the indexes of the points that are 
   considered to correspond to the same global coordinate. A points indexes are stores using a `std::pair<size_t,size_t>` where `first` equals the point cloud index 
   and `second` equals the index within this point cloud that it belongs to (indexes with respect to `point_clouds`).  
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.
   \param registered_points ouput point range of the registered point clouds.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` whose key type
     is the value type of the iterator of `PointRange` and whose value type is
     `geom_traits::Point_3`.  If this parameter is omitted,
     `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
   \cgalNamedParamsEnd

*/
template <class PointRange, class CorrespondencesRange, class NamedParameters>
void register_point_clouds( const std::vector<PointRange>& point_clouds, const CorrespondencesRange& correspondences, 
                            const NamedParameters& np, PointRange& registered_points)
{
    namespace PSP = CGAL::Point_set_processing_3;
    // property map types
    typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::type PointMap;
    typedef typename PSP::GetNormalMap<PointRange, NamedParameters>::type NormalMap;

    using Kernel = typename CGAL::Kernel_traits<typename PointMap::value_type>::Kernel;
    using Pwn = typename PointRange::value_type;
    using parameters::choose_parameter;
    using parameters::get_parameter;

    // compute registration transformations
    std::vector<typename Kernel::Aff_transformation_3> transformations;
    compute_registration_transformations(point_clouds, correspondences, np, transformations);




    PointMap point_map = choose_parameter(get_parameter(np, internal_np::point_map), PointMap());
    NormalMap normal_map = choose_parameter(get_parameter(np, internal_np::normal_map), NormalMap());

    // add transformed point clouds to registered_points
    for (size_t i = 0; i < point_clouds.size(); i++){
        for (size_t j = 0; j < point_clouds[i].size(); j++){
            auto point = get(point_map, point_clouds[i][j]).transform(transformations[i]);
            auto normal = get(normal_map, point_clouds[i][j]);
            Pwn pwn;
            put(point_map, pwn, point);
            put(normal_map, pwn, normal);
            registered_points.push_back(pwn);
        }
    }
}

} } // end of namespace CGAL::OpenGR

#endif // CGAL_LINKED_WITH_OPENGR

#endif // CGAL_OPENGR_GRET_SDP_H