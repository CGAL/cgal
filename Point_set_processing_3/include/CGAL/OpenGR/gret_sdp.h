#ifndef CGAL_OPENGR_GRET_SDP_H
#define CGAL_OPENGR_GRET_SDP_H

#include <CGAL/license/Point_set_processing_3.h>

#if defined(CGAL_LINKED_WITH_OPENGR) || defined(DOXYGEN_RUNNING)

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

// currently used for CGAL_range_and_pmaps_to_opengr_point3d_range
#include <CGAL/OpenGR/compute_registration_transformation.h>

#include <boost/range/iterator_range.hpp>

#include <gr/algorithms/GRET_SDP.h>
#include <gr/accelerators/MOSEKWrapper.h>
#include <gr/shared.h>
#include <memory>

namespace CGAL {

namespace OpenGR {

    namespace internal {


        template <typename Scalar, typename InputRange, typename PointMap, typename IndexMap, typename VectorMap>
        struct CGAL_range_and_pmaps_to_opengr_indexed_point3d_range
        {
        typedef typename InputRange::const_iterator::value_type argument_type;
        typedef typename gr::Point3D<Scalar> point_type;
        typedef std::pair<point_type, int> result_type;
        typedef typename point_type::VectorType vector_type;

        PointMap point_map;
        IndexMap index_map;
        VectorMap normal_map;

        CGAL_range_and_pmaps_to_opengr_indexed_point3d_range (PointMap point_map, IndexMap index_map, VectorMap normal_map)
            : point_map (point_map), index_map(index_map), normal_map (normal_map)
        { }

        result_type operator() (const argument_type& arg) const
        {
            const auto& p = get (point_map, arg);
            const auto& n = get (normal_map, arg);
            const auto& i = get (index_map, arg);

            point_type out (p.x(), p.y(), p.z());
            out.set_normal ( vector_type(n.x(), n.y(), n.z()) );

            return std::make_pair(out, i);
        }
        };


        template <typename InputRange, typename Scalar, typename PointMap, typename IndexMap, typename VectorMap>
        struct CGAL_range_and_pmaps_range_to_opengr_point3d_PatchRange {
            typedef typename InputRange::value_type argument_type;
            typedef typename argument_type::value_type value_type;


            // typedef typename boost::iterator_range<<boost::iterators::transform_iterator<
            // CGAL_range_and_pmaps_to_opengr_point3d_range<Scalar, argument_type, PointMap, VectorMap>,
            // argument_type::iterator>> return_type;

            CGAL_range_and_pmaps_to_opengr_indexed_point3d_range<Scalar, argument_type, PointMap, IndexMap, VectorMap>
            unary_function;

            CGAL_range_and_pmaps_range_to_opengr_point3d_PatchRange(PointMap point_map, IndexMap index_map, VectorMap vector_map) 
                :  unary_function (point_map, index_map, vector_map) 
                { }

            auto operator() (const argument_type& range) const {
                return boost::make_iterator_range(
                    boost::make_transform_iterator (range.begin(), unary_function),
                    boost::make_transform_iterator (range.end(),   unary_function));
                }
        };
    }

    template <typename Scalar>
    using GRET_SDPOptions = typename gr::GRET_SDP<gr::Point3D<Scalar>, gr::DummyTransformVisitor, gr::GRET_SDP_Options>::OptionsType;

    template <class Scalar>
    class GRET_SDP {
        using GrMatcherType = gr::GRET_SDP<gr::Point3D<Scalar>, gr::DummyTransformVisitor, gr::GRET_SDP_Options>;
        using GrOptionsType = typename GrMatcherType::OptionsType;
        

        public:
        GRET_SDP();

        template <class PatchRange, class NamedParameters>
        void registerPatches(const PatchRange& patches, const int n, const NamedParameters& np);

        // returns transformations
        template<typename TrRange>
        void getTransformations(TrRange& transformations);

        // returns registered points
        template<typename PointRange>
        void getRegisteredPatches(PointRange& registered_points);

        private:
        std::unique_ptr<GrMatcherType> gr_matcher;

        template <class PatchRange, class PointMap, class IndexMap, class VectorMap>
        void registerPatches(const PatchRange& patches, const int n, PointMap point_map, IndexMap index_map, VectorMap vector_map);

        template <class NamedParameters>
        GrOptionsType GetOptionsFromNamedParameters(const NamedParameters& np);
    };

    template <class Scalar>
    GRET_SDP<Scalar>::GRET_SDP() {};

    template <class Scalar>
    template <class NamedParameters>
    typename GRET_SDP<Scalar>::GrOptionsType 
    GRET_SDP<Scalar>::GetOptionsFromNamedParameters(const NamedParameters& np) {
        GrOptionsType options;
        return options;
    };

    template <class Scalar>
    template <class PatchRange, class NamedParameters>
    void GRET_SDP<Scalar>::registerPatches(const PatchRange& patches, const int n, const NamedParameters& np)
    {
        namespace PSP = CGAL::Point_set_processing_3;
        typedef typename PatchRange::value_type PointRange;
        using parameters::choose_parameter;
        using parameters::get_parameter;

        // property map types
        typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::type PointMap;
        typedef typename PSP::GetNormalMap<PointRange, NamedParameters>::type NormalMap;

        //typedef typename PSP::GetK<PointRange, NamedParameters>::Kernel Kernel;
        PointMap point_map = choose_parameter(get_parameter(np, internal_np::point_map), PointMap());
        NormalMap normal_map = choose_parameter(get_parameter(np, internal_np::normal_map), NormalMap());
        auto index_map = get_parameter(np, internal_np::vertex_index);

        // using GrMatcherType = gr::GRET_SDP<gr::Point3D<typename Kernel::FT>, gr::DummyTransformVisitor, gr::GRET_SDP_Options>;
        GrOptionsType options = GetOptionsFromNamedParameters(np);
        gr::Utils::Logger logger(gr::Utils::NoLog);

        gr_matcher = std::unique_ptr<GrMatcherType>(new GrMatcherType(options, logger));

        registerPatches(patches, n, point_map, index_map, normal_map);
    }
    
    // returns transformations
    template <class Scalar>
    template<typename TrRange>
    void GRET_SDP<Scalar>::getTransformations(TrRange& transformations){

    }

    // returns registered points
    template <class Scalar>
    template<typename PointRange>
    void GRET_SDP<Scalar>::getRegisteredPatches(PointRange& registered_points){
        
    }
    template <class Scalar>
    template <class PatchRange, class PointMap, class IndexMap, class VectorMap>
    void GRET_SDP<Scalar>::registerPatches(const PatchRange& patches, const int n, PointMap point_map, IndexMap index_map, VectorMap vector_map){
        // unary function that converty CGAL PatchRange to OpenGR PatchRange
        internal::CGAL_range_and_pmaps_range_to_opengr_point3d_PatchRange<PatchRange,Scalar, PointMap, IndexMap, VectorMap>
        unary_function(point_map, index_map, vector_map);

        // construct gr patches
        auto gr_patches = boost::make_iterator_range(
        boost::make_transform_iterator (patches.begin(), unary_function),
        boost::make_transform_iterator (patches.end(),   unary_function));

        gr::DummyTransformVisitor tr_visitor;

        GrMatcherType matcher(GrOptionsType(), gr::Utils::Logger(gr::Utils::NoLog));
        matcher.template RegisterPatches<gr::MOSEK_WRAPPER<Scalar>>(gr_patches, n, tr_visitor);

        gr_matcher->template RegisterPatches<gr::MOSEK_WRAPPER<Scalar>>(gr_patches, n, tr_visitor);

    }


} } // end of namespace CGAL::OpenGR

#endif // CGAL_LINKED_WITH_OPENGR

#endif // CGAL_OPENGR_GRET_SDP_H