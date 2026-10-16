#if defined (_MSC_VER) && !defined (_WIN64)
#pragma warning(disable:4244) // boost::number_distance::distance()
// converts 64 to 32 bits integers
#endif

#define COMPUTE_COVERAGE false
#define COMPUTE_SKELETON false
#define COMPUTE_COMPACTNESS false
#define COMPUTE_BETTI true
#define COMPUTE_FRACTAL_DIMENSIONALITY false

// Features
#include <CGAL/Classification/Feature/Fractal_dimensionality.h>
#include <CGAL/Classification/Feature/Coverage.h>
#include <CGAL/Classification/Feature/Skeletonise.h>
#include <CGAL/Classification/Feature/Compactness.h>
#include <CGAL/Classification/Feature/Betti_numbers.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#ifdef BOOST_TIMER_PROGRESS_DISPLAY_HPP_INCLUDED
#define CGAL_CLASSIFICATION_VERBOSE true
#include <boost/timer/progress_display.hpp>
#endif

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Real_timer.h>
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
typedef Point_set::Point_map Pmap;
namespace Classification = CGAL::Classification;
typedef Classification::Feature_set                                             Feature_set;
typedef Classification::Point_set_feature_generator<Kernel, Point_set, Pmap>    Feature_generator;

using Neighborhood = Classification::Point_set_neighborhood<Kernel, Point_set, Pmap>;
using Neighbor_query = typename Neighborhood::Sphere_neighbor_query;

#include <CGAL/Surface_mesh.h>
#include <CGAL/alpha_wrap_3.h>
namespace AW3 = CGAL::Alpha_wraps_3;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel_Alpha;
typedef Kernel_Alpha::Point_3 Point_Alpha;
using Point_container = std::vector<Point_Alpha>;
using Mesh = CGAL::Surface_mesh<Point_Alpha>;

int main(int argc, char** argv)
{
    if (argc == 1) {
        std::cout << "Usage: " << argv[0] << " filename scale alpha offset" << std::endl;
    }
    std::string filename = (argc > 1) ? argv[1] : "C:/Users/rdyke/Documents/Datasets/ROAD-AI/Chambon/Chambon_Scan_Riegl_20210712.ply";
    filename = "C:/Users/rdyke/Documents/Datasets/ROAD-AI/Chambon/Chambon_Scan_Riegl_20210712_small.ply";
    //filename = "C:/Users/rdyke/Documents/Research/ROAD-AI/Betti\ number/wrap_intersect_smaller_2_comp.ply";
    //filename = "C:/Users/rdyke/Documents/Datasets/ROAD-AI/Talus_CereMap3D/Run1-Falaise-Sud-training_data.ply";
    filename = CGAL::data_file_path("meshes/b9.ply");

    std::ifstream in(filename.c_str(), std::ios::binary);
    Point_set pts;
    std::cerr << "Reading input" << std::endl;
    in >> pts;
    if (pts.number_of_points() == 0) {
        std::cerr << "Error: no vertices found." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "No. of vertices: " << pts.number_of_points() << std::endl;

    std::cout << "Properties found:" << std::endl;
    for (const auto& prop : pts.properties_and_types()) {
        std::cout << "  " << prop.first << std::endl;
    }

    const float radius_size = (argc > 2) ? atof(argv[2]) : 0.5f; //0.1; specify radius of neighborhoods (default: 10cm)
    const float voxel_size = radius_size / 3.f; // re-scale for CGAL's feature generator

    // Wrap point cloud
    // Compute the alpha and offset values
    const double relative_alpha = (argc > 3) ? std::stod(argv[3]) : 0.1f;
    const double relative_offset = (argc > 4) ? std::stod(argv[4]) : 2.f;
    std::cout << "relative alpha = " << relative_alpha << " relative offset = " << relative_offset << std::endl;
    double alpha = radius_size / relative_alpha; 
    double offset = radius_size / relative_offset;
    std::cout << "absolute alpha = " << alpha << " absolute offset = " << offset << std::endl;

    if (radius_size < offset) {
        std::cerr << "Warning: The smallest radius_size (" << radius_size << "cm) is smaller than the offset (" << offset << "cm). This is too small to be useful." << std::endl;
    }

    CGAL::Real_timer t;
    t.start();

    // convert to a kernel that is more stable for Alpha Wrap
    Point_container points;
    for (const auto& point : pts.points()) {
        Point_Alpha pt(point.x(), point.y(), point.z());
        points.push_back(pt);
    }

    // construct the wrap
    Mesh wrap;
    if (COMPUTE_BETTI || COMPUTE_COMPACTNESS || COMPUTE_COVERAGE || COMPUTE_SKELETON) {
        CGAL::alpha_wrap_3(points, alpha, offset, wrap);
        std::cout << "Result: " << num_vertices(wrap) << " vertices, " << num_faces(wrap) << " faces, " << std::endl;

        CGAL::IO::write_polygon_mesh("wrap.ply", wrap, CGAL::parameters::stream_precision(17));
        std::cout << "Wrap saved" << std::endl;
    }
    t.stop();
    std::cout << "Took " << t.time() << " s" << std::endl;


    std::cout << "Computing advanced features..." << std::endl;

    Feature_set features;
    std::cerr << "Initialising feature generator...";

    t.reset(); t.start();
    Feature_generator generator(pts, pts.point_map(), 5, voxel_size);  // using 5 scales
    t.stop();
    std::cout << "done in " << t.time() << " second(s)" << std::endl;

    std::cout << "Neighbourhood radii: " << generator.radius_neighbors(0);
    for (std::size_t i = 1; i < generator.number_of_scales(); ++i)
        std::cout << ", " << generator.radius_neighbors(i);
    std::cout << std::endl;

    if (COMPUTE_COVERAGE) {
        std::cout << "Computing coverage feature..." << std::endl;
        using MyCoverage = CGAL::Classification::Feature::Coverage<Kernel, Point_set, Pmap, Mesh, CGAL::Parallel_if_available_tag>;// <typename GeomTraits, typename PointRange, typename PointMap, typename Feature_generator, typename ConcurrencyTag>
        t.reset(); t.start();
        for (std::size_t i = 0; i < generator.number_of_scales(); ++i) {
            std::cout << "  scale: " << i << std::endl;
            features.add_with_scale_id<MyCoverage>(i, pts, pts.point_map(), wrap, generator.radius_neighbors(i));
        }
        t.stop();
        std::cout << "done in " << t.time() << " second(s)" << std::endl;
    }

    // construct the wrap
    if (COMPUTE_SKELETON) {
        std::cout << "Computing skeletonize feature..." << std::endl;
        t.reset(); t.start();
        using MySkeleton = CGAL::Classification::Feature::Skeleton<Kernel, Point_set, Pmap, Mesh>;
        features.add_multidimensional_feature<MySkeleton>(2, pts, pts.point_map(), wrap);
        std::cout << "done in " << t.time() << " second(s)" << std::endl;
    }

    if (COMPUTE_COMPACTNESS) {
        std::cout << "Computing compactness feature..." << std::endl;
        using MyCompactness = CGAL::Classification::Feature::Compactness<Kernel, Point_set, Point_set::Point_map, Mesh, CGAL::Parallel_if_available_tag>;
        t.reset(); t.start();
        auto bbox = CGAL::bounding_box(CGAL::make_transform_iterator_from_property_map(pts.begin(), pts.point_map()),
            CGAL::make_transform_iterator_from_property_map(pts.end(), pts.point_map()));
        using Planimetric_grid = Classification::Planimetric_grid<Kernel, Point_set, Pmap>;
        for (std::size_t i = 0; i < generator.number_of_scales(); ++i) {
            std::cout << "  scale: " << i << std::endl;
            Planimetric_grid grid(pts, pts.point_map(), bbox, generator.grid_resolution(i));
            features.add_with_scale_id<MyCompactness>(i, pts, pts.point_map(), grid, generator.radius_neighbors(i), wrap);
        }
        std::cout << "done in " << t.time() << " second(s)" << std::endl;
    }

    if (COMPUTE_BETTI) {
        std::cout << "Computing betti feature..." << std::endl;
        using MyBettiNumbers = CGAL::Classification::Feature::Betti_numbers<Kernel, Point_set, Point_set::Point_map, Mesh, CGAL::Parallel_if_available_tag>;
        t.reset(); t.start();
        auto bbox = CGAL::bounding_box(CGAL::make_transform_iterator_from_property_map(pts.begin(), pts.point_map()),
            CGAL::make_transform_iterator_from_property_map(pts.end(), pts.point_map()));
        using Planimetric_grid = Classification::Planimetric_grid<Kernel, Point_set, Pmap>;
        for (std::size_t i = 0; i < generator.number_of_scales(); ++i) {
            std::cout << "  scale: " << i << std::endl;
            Planimetric_grid grid(pts, pts.point_map(), bbox, generator.grid_resolution(i));
            features.add_multidimensional_feature_with_scale_id<MyBettiNumbers>(i, 3, pts, pts.point_map(), grid, generator.radius_neighbors(i), wrap);
        }
        std::cout << "done in " << t.time() << " second(s)" << std::endl;
    }

    if (COMPUTE_FRACTAL_DIMENSIONALITY) {
        std::cout << "Computing fractal feature..." << std::endl;
        using MyFractalDimensionality = CGAL::Classification::Feature::Fractal_dimensionality<Kernel, Point_set, Point_set::Point_map>;
        t.reset(); t.start();
        for (std::size_t i = 0; i < generator.number_of_scales(); ++i) {
            features.add_with_scale_id<MyFractalDimensionality>(i, pts, pts.point_map(), generator.radius_neighbors(i));
        }
        std::cout << "done in " << t.time() << " second(s)" << std::endl;
    }


    // Export features (adds features to output file)
    std::cout << "exporting features...";
    for (auto feature : features) {
        Point_set::Property_map<float> prop = pts.add_property_map<float>("scalar_" + feature->name(), 0).first;
        for (size_t i = 0; i < pts.size(); ++i) {
            prop[i] = feature->value(i);
        }
    }
    std::cout << "done" << std::endl;

    //CGAL::IO::write_PLY("computed_features.ply", pts, CGAL::parameters::use_binary_mode(true));
#ifdef CGAL_LINKED_WITH_EMBREE
    CGAL::IO::write_PLY("computed_features_EMBREE.ply", pts, CGAL::parameters::use_binary_mode(false));
#else
    CGAL::IO::write_PLY("computed_features_AABB.ply", pts, CGAL::parameters::use_binary_mode(false));
#endif
    std::cout << "Completed successfully" << std::endl;
    return EXIT_SUCCESS;
}