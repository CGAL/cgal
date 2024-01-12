// Author(s)     : Roberto M. Dyke
#ifndef CGAL_COVERAGE_FEATURE_H
#define CGAL_COVERAGE_FEATURE_H
#include <CGAL/license/Classification.h>

#include <CGAL/Classification/Feature_base.h>

#include <CGAL/Surface_mesh.h>

#include <limits>
struct Vertex { float x, y, z; };
struct Face { int a, b, c; };

#include <CGAL/Surface_mesh.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Cartesian_converter.h>

#ifndef CGAL_COVERAGE_FEATURE_DISPLAY_PROGRESS
#define CGAL_COVERAGE_FEATURE_DISPLAY_PROGRESS false
#endif
#if CGAL_COVERAGE_FEATURE_DISPLAY_PROGRESS
#include <boost/timer/progress_display.hpp>
#endif

#ifdef CGAL_LINKED_WITH_TBB
#include <thread>
#endif


//#undef CGAL_LINKED_WITH_EMBREE
#ifdef CGAL_LINKED_WITH_EMBREE
#include <embree3/rtcore.h>
#else
// CGAL AABB tree
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#endif

#include <CGAL/Real_timer.h>

namespace CGAL {
namespace   Classification {
namespace     Feature {
    /*!
    \ingroup PkgClassificationFeatures

    %Feature based on local coverage w.r.t. the wrapped point cloud.
    The coverage of a point quantifies the number of angles that a
    given point can "see out" from.

    Its default name is "coverage".

    \tparam GeomTraits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator` and its value type is the key type of
    `PointMap`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `GeomTraits::Point_3`.

  */
    namespace AW3 = CGAL::Alpha_wraps_3;
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel_Alpha;
    typedef Kernel_Alpha::Point_3 Point_Alpha;
    using Point_container = std::vector<Point_Alpha>;
    using Mesh = CGAL::Surface_mesh<Point_Alpha>;
#ifndef CGAL_LINKED_WITH_EMBREE
    //typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
    //typedef CGAL::AABB_traits<Kernel_Alpha, Primitive> AABB_triangle_traits;
    //typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
    //typedef Kernel_Alpha::Ray_3 Ray;
    //typedef std::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;
#endif

    template <typename GeomTraits, typename PointRange, typename PointMap, typename Feature_generator, typename ConcurrencyTag>
    void ComputeCoverage(const PointRange& input, Feature_set& features, Feature_generator& generator, const double alpha, const double offset) {
        /*!
        \brief constructs the feature and appends it to a feature set (makes incorporation of feature a little tidier).

        \param input point range.
        \param feature set to append Coverage feature to.
        \param feature generator that determines the feature scale(s).
        \param alpha the wrap fidelity.
        \param offset the distance between the wrap and the point cloud.
        */
        Mesh wrap;

        #if CGAL_COVERAGE_FEATURE_DISPLAY_PROGRESS
        std::cerr << "Converting point cloud data" << std::endl;
        #endif
        // convert to a kernel that is more stable for Alpha Wrap
        Point_container points;
        for (auto& point : input.points()) {
            Point_Alpha pt(point.x(), point.y(), point.z());
            //Point_Alpha pt(point);
            points.push_back(pt);
        }

        // construct the wrap
        #if CGAL_COVERAGE_FEATURE_DISPLAY_PROGRESS
        std::cerr << "Wrapping points..." << std::endl;
        #endif
        CGAL::Real_timer t;
        t.start();
        CGAL::alpha_wrap_3(points, alpha, offset, wrap);
        t.stop();
        std::cerr << "Result: " << num_vertices(wrap) << " vertices, " << num_faces(wrap) << " faces " << std::endl;

        std::cerr << "Took " << t.time() << " s" << std::endl;

        using Coverage_ = CGAL::Classification::Feature::Coverage<GeomTraits, PointRange, PointMap, Mesh, ConcurrencyTag>;
        for (std::size_t i = 0; i < generator.number_of_scales(); ++i) {
            std::cerr << "  scale: " << i << ", radius: " << generator.radius_neighbors(i) << std::endl;
            features.add_with_scale_id<Coverage_>(i, input, input.point_map(), wrap, generator.radius_neighbors(i));
        }
    }

    template <typename GeomTraits, typename PointRange, typename PointMap, typename Mesh, typename ConcurrencyTag>
    class Coverage : public CGAL::Classification::Feature_base
    {
        using MeshKernel = typename Mesh::Point::R::Kernel;
        using MeshPoint = typename Mesh::Point;
        using Point_to_Mesh = CGAL::Cartesian_converter<GeomTraits, MeshKernel>;
        Point_to_Mesh to_mesh;

#ifndef CGAL_LINKED_WITH_EMBREE
        typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
        typedef CGAL::AABB_traits<MeshKernel, Primitive> AABB_triangle_traits;
        typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
        typedef typename MeshKernel::Ray_3 Ray;
        typedef std::optional<typename Tree::template Intersection_and_primitive_id<Ray>::Type> Ray_intersection;
#endif

        using FloatMap = typename PointRange::template Property_map<float>;
        PointMap point_map;
        std::vector<float> coverage;

    public:
    /*!
    \brief constructs the feature.

    \note By default the AABB tree package is used. Alternatively, if Embree is included in the project, then Embree is used as it executes faster.

    \tparam GeomTraits must be a model of `Kernel`
    \tparam PointRange is a model of `Range`. The value type of
    its iterator is the key type of the named parameter `point_map`.
    \tparam PointMap is a model of `ReadablePropertyMap` whose value type is `GeomTraits::Point_3`
    \tparam ConcurrencyTag enables sequential versus parallel algorithm. Possible values are `Sequential_tag`,
    `Parallel_tag`, and `Parallel_if_available_tag`.
    \tparam Mesh is a surface mesh.

    \param input point range.
    \param point_map property map to access the input points.
    \param wrap a mesh that is offset from the `input` point set.
    \param tnear distance from from a point (along a ray) that the ray is actually shot from. This parameter
    controls the scale of the feature.
    \param tfar maximum distance a ray is shot.
    \param n_rays the number of rays shot per point.
  */

#ifdef CGAL_LINKED_WITH_EMBREE
        Coverage(const PointRange& input, PointMap point_map, const Mesh& wrap, float tnear = 0.f, int n_rays = 25) : point_map(point_map)
        {
            this->set_name("coverage_" + std::to_string(tnear));

            tnear = tnear * .5f;

            unsigned int n = wrap.number_of_vertices();
            unsigned int m = wrap.number_of_faces();

            RTCDevice device = rtcNewDevice(NULL);
            RTCScene scene = rtcNewScene(device);

            // Convert surface_mesh data into Embree's format
            RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
            Vertex* vertices = (Vertex*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * sizeof(float), n);
            {
                int i = 0;
                for (const auto& point : wrap.points()) {
                    vertices[i].x = point.x();
                    vertices[i].y = point.y();
                    vertices[i].z = point.z();
                    ++i;
                }
            }
            Face* faces = (Face*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(int), m);
            {
                int i = 0;
                for (const auto& fd : wrap.faces()) {
                    auto& hd = wrap.halfedge(fd);
                    faces[i].a = wrap.source(hd);
                    hd = wrap.next(hd);
                    faces[i].b = wrap.source(hd);
                    hd = wrap.next(hd);
                    faces[i].c = wrap.source(hd);
                    hd = wrap.next(hd);
                    ++i;
                }
            }

            rtcCommitGeometry(geom);
            rtcAttachGeometry(scene, geom);
            rtcReleaseGeometry(geom);
            rtcCommitScene(scene);

            RTCIntersectContext context;
            rtcInitIntersectContext(&context);

            // construct rays
            int n_thetas = floor(sqrt(n_rays / 2.f) * 2.f);
            int n_phis = ceil(sqrt(n_rays / 2.f) * 2.f);
            int N_rays = n_thetas * n_phis;

#if CGAL_COVERAGE_FEATURE_DISPLAY_PROGRESS
            std::cerr << "Shooting rays (Embree)..." << std::endl;
#ifdef BOOST_TIMER_PROGRESS_DISPLAY_HPP_INCLUDED
            boost::timer::progress_display progress(input.number_of_points());
#else
            CGAL::Real_timer t;
            t.reset();
            t.start();
#endif
#endif
            coverage.reserve(input.number_of_points());
            CGAL::for_each<ConcurrencyTag>
                (CGAL::make_counting_range<std::size_t>(0, input.size()),
                    [&](const std::size_t& s) -> bool
                    {
                        RTCRay* rays = new RTCRay[N_rays];
                        {
                            int k = 0;
                            for (int i = 0; i < n_thetas; ++i) {
                                float theta = (i + 1.f) / (n_thetas + 1.f) * CGAL_PI;
                                for (int j = 0; j < n_phis; ++j) {
                                    //rays[k].org_x = 0.f; rays[k].org_y = 0.f; rays[k].org_z = 0.f;
                                    float phi = (j + 1.f) / (n_thetas + 1.f) * 2.f * CGAL_PI;
                                    rays[k].dir_x = (float)(sin(theta) * cos(phi));
                                    rays[k].dir_y = (float)(sin(theta) * sin(phi));
                                    rays[k].dir_z = (float)(cos(theta));
                                    rays[k].tnear = tnear;
                                    rays[k].tfar = 1e+9f;
                                    ++k;
                                }
                            }
                        }


                        const auto& point = input.point(s);
                        // prepare rays
                        for (int i = 0; i < N_rays; ++i) {
                            rays[i].org_x = point.x(); rays[i].org_y = point.y(); rays[i].org_z = point.z();
                            rays[i].tfar = 1e+9f;
                        }

                        rtcOccluded1M(scene, &context, rays, N_rays, sizeof(RTCRay));

                        int n_occluded = 0;
                        for (int i = 0; i < N_rays; ++i)
                            if (rays[i].tfar < 0)
                                n_occluded++;

                        coverage[s] = n_occluded / float(N_rays);

#if CGAL_COVERAGE_FEATURE_DISPLAY_PROGRESS
#ifdef BOOST_TIMER_PROGRESS_DISPLAY_HPP_INCLUDED
                        ++progress;
#endif
#endif

                        delete[] rays;
                        return true;
                    });
#if CGAL_COVERAGE_FEATURE_DISPLAY_PROGRESS
#ifndef BOOST_TIMER_PROGRESS_DISPLAY_HPP_INCLUDED
            t.stop();
            std::cout << "Took " << t.time() << " s." << std::endl;
#endif
#endif

        }
#else
        Coverage(const PointRange& input, PointMap point_map, const Mesh& wrap, float tnear = 0.f, int n_rays = 25) : point_map(point_map)
        {
            

            this->set_name("coverage_" + std::to_string(tnear));

            tnear = tnear * .5f;

            unsigned int n = input.number_of_points();
            unsigned int m = wrap.number_of_faces();

            int n_thetas = floor(sqrt(n_rays / 2.f) * 2.f);
            int n_phis = ceil(sqrt(n_rays / 2.f) * 2.f);
            int N_rays = n_thetas * n_phis;

            // construct AABB tree
            Tree tree(wrap.faces_begin(), wrap.faces_end(), wrap);
            tree.accelerate_distance_queries();

#ifdef CGAL_COVERAGE_FEATURE_DISPLAY_PROGRESS
            std::cout << "Shooting rays (AABB)...";
            CGAL::Real_timer t;
            t.reset();
            t.start();
#endif

            std::vector<MeshKernel::Vector_3> dirs;
            dirs.resize(N_rays);
            int k = 0;
            for (int i = 0; i < n_thetas; ++i) {
                float theta = (i + 1.f) / (n_thetas + 1.f) * CGAL_PI;
                for (int j = 0; j < n_phis; ++j) {
                    //rays[k].org_x = 0.f; rays[k].org_y = 0.f; rays[k].org_z = 0.f;
                    float phi = (j + 1.f) / (n_thetas + 1.f) * 2.f * CGAL_PI;
                    dirs[k] = MeshKernel::Vector_3(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
                    ++k;
                }
            }
            tree.accelerate_distance_queries();
            coverage.resize(n);


            CGAL::for_each<ConcurrencyTag>
                (CGAL::make_counting_range<std::size_t>(0, input.size()),
                    [&](const std::size_t& s) -> bool
                    {
                        MeshPoint point = to_mesh(input.point(s));
                        int n_occluded = 0;
                        for (int k = 0; k < N_rays; ++k) {
                            const MeshPoint org = point + dirs[k] * tnear;
                            Ray ray(org, dirs[k]);
                            auto intersection = tree.any_intersection(ray);
                            if (intersection && std::get_if<MeshPoint>(&(intersection->first))) {
                                const MeshPoint* p = std::get_if<MeshPoint>(&(intersection->first));
                                auto dist = CGAL::sqrt(CGAL::squared_distance(*p, org));
                                n_occluded++;
                            }
                        }
                        coverage[s] = n_occluded / float(N_rays);
                        return true;
                    });

#ifdef CGAL_COVERAGE_FEATURE_DISPLAY_PROGRESS
            t.stop();
            std::cout << "done." << std::endl;
            std::cout << "Took " << t.time() << " s." << std::endl;
#endif
        }
#endif // CGAL_LINKED_WITH_EMBREE
        /// \cond SKIP_IN_MANUAL
        float value(std::size_t pt_index) {
            return coverage[pt_index];
        }
        /// \endcond
    };


} // namespace Feature
} // namespace Classification
} // namespace CGAL

#endif // CGAL_COVERAGE_FEATURE_H