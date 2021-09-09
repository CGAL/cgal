#ifndef POLYFIT_TEST_FRAMEWORK_H_
#define POLYFIT_TEST_FRAMEWORK_H_

#include <CGAL/IO/read_points.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Shape_detection/Efficient_RANSAC.h>
#include <CGAL/Polygonal_surface_reconstruction.h>
#include <CGAL/Timer.h>

#include <string>

template <typename Reconstruction, typename SMesh, typename MIP_Solver>
int reconstruct (Reconstruction& algo, SMesh& model, const MIP_Solver&)
{
        CGAL::Timer t;
        t.start();
  if (!algo.template reconstruct<MIP_Solver>(model)) {
                std::cerr << " Failed: " << algo.error_message() << std::endl;
                return EXIT_FAILURE;
        }
  std::cout << " Done. Time: " << t.time() << " sec. ";

        if (model.is_valid()) {
    std::cout << "\tReconstructed model has " << model.number_of_faces() << " faces" << std::endl;
    return EXIT_SUCCESS;
        }
        else {
    std::cout << "\tReconstructed model is not valid. Reconstruction maybe failed?" << std::endl;
    return EXIT_FAILURE;
        }
}

template <typename Reconstruction, typename SMesh>
int reconstruct (Reconstruction&, SMesh&, const int&)
{
  std::cout << " No solver. Nothing happened." << std::endl;
  return EXIT_SUCCESS;
}


// This function enables to run the Polygonal Surface Reconstruction algorithm with different
//    - kernels(Simple_cartesian, EPICK)
//    - solvers(GLPK, SCIP)
//    - use/ignore provided planar segmentation
//    - input file formats (pwn, ply). For ply format, a property named "segment_index"
//                stores the plane index for each point(-1 if the point is not assigned to a plane).

template <typename Kernel, typename MIP_Solver>
int reconstruct(const std::string& input_file, bool force_extract_planes)
{
    typedef typename Kernel::Point_3                                Point;
    typedef typename Kernel::Vector_3                                Vector;
    typedef        CGAL::Polygonal_surface_reconstruction<Kernel>                Polygonal_surface_reconstruction;
    typedef CGAL::Surface_mesh<Point>                                Surface_mesh;

    // Point with normal, and plane index
    typedef boost::tuple<Point, Vector, int>                                                        PNI;
    typedef std::vector<PNI>                                                                                        Point_vector;
    typedef CGAL::Nth_of_tuple_property_map<0, PNI>                                                Point_map;
    typedef CGAL::Nth_of_tuple_property_map<1, PNI>                                                Normal_map;
    typedef CGAL::Nth_of_tuple_property_map<2, PNI>                                                Plane_index_map;

    typedef CGAL::Shape_detection::Efficient_RANSAC_traits<Kernel, Point_vector, Point_map, Normal_map>     Traits;

    typedef CGAL::Shape_detection::Efficient_RANSAC<Traits>             Efficient_ransac;
    typedef CGAL::Shape_detection::Plane<Traits>                                                Plane;
    typedef CGAL::Shape_detection::Point_to_shape_index_map<Traits>     Point_to_shape_index_map;

        Point_vector points;

        // Loads point set from a file.
        std::ifstream input_stream(input_file.c_str());
    if (!input_stream) {
        std::cerr << " Error: cannot read file " << input_file << std::endl;
        return EXIT_FAILURE;
    }
        std::cout << "\t\t\tLoading point cloud: " << input_file << "...";

        std::string extension = input_file.substr(input_file.find_last_of('.'));

        bool input_has_planes = false;

        CGAL::Timer t;
        t.start();
        if (extension == ".pwn") {
          if (!CGAL::IO::read_XYZ(input_stream,
                                  std::back_inserter(points),
                                  CGAL::parameters::point_map(Point_map()).normal_map(Normal_map())))
          {
            std::cerr << " Error: cannot read file " << input_file << std::endl;
            return EXIT_FAILURE;
          }
          else
            std::cout << " Done. " << points.size() << " points. Time: " << t.time() << " sec." << std::endl;
        }
        else if (extension == ".ply") {
        if (!CGAL::IO::read_PLY_with_properties(
                    input_stream,
                    std::back_inserter(points),
                    CGAL::make_ply_point_reader(Point_map()),
                    CGAL::make_ply_normal_reader(Normal_map()),
                    std::make_pair(Plane_index_map(), CGAL::PLY_property<int>("segment_index"))))
        {
          std::cerr << " Error: cannot read file " << input_file << std::endl;
          return EXIT_FAILURE;
        }
        else
          std::cout << " Done. " << points.size() << " points. Time: " << t.time() << " sec." << std::endl;

        int max_plane_index = 0;
        for (std::size_t i = 0; i < points.size(); ++i) {
          int plane_index = points[i].template get<2>();
          if (plane_index > max_plane_index)
            max_plane_index = plane_index;
        }

                int num_planes = max_plane_index + 1;
                // if fewer than 4 planes provided, we consider it as not provided.
                if (num_planes > 3) {
                        input_has_planes = true;
                        std::cout << "\t\t\tInput has " << num_planes << " planes" << std::endl;
                }
        }

        if (!input_has_planes || force_extract_planes) {
                if (input_has_planes && force_extract_planes)
                        std::cout << "\t\t\tExtracting planes (forced)...";
                else if (!input_has_planes)
                        std::cout << "\t\t\tExtracting planes...";

                // Shape detection
                Efficient_ransac ransac;
                ransac.set_input(points);
        ransac.template add_shape_factory<Plane>();

                t.reset();
                ransac.detect();

        typename Efficient_ransac::Plane_range planes = ransac.planes();
                std::size_t num_planes = planes.size();

                std::cout << " Done. " << num_planes << " planes extracted. Time: " << t.time() << " sec." << std::endl;

                // store the plane index of each point as the third element of the tuple.
                Point_to_shape_index_map shape_index_map(points, planes);
                for (std::size_t i = 0; i < points.size(); ++i) {
                        // Use the get function from the property map that accesses the 3rd element of the tuple.
                        int plane_index = get(shape_index_map, i);
            points[i].template get<2>() = plane_index;
                }
        }

        //////////////////////////////////////////////////////////////////////////

        std::cout << "\t\t\tGenerating candidate faces...";
        t.reset();

        Polygonal_surface_reconstruction algo(
                points,
                Point_map(),
                Normal_map(),
                Plane_index_map()
        );

        std::cout << " Done. Time: " << t.time() << " sec." << std::endl;

        //////////////////////////////////////////////////////////////////////////

        Surface_mesh model;

        std::cout << "\t\t\tReconstructing...";

  return reconstruct (algo, model, MIP_Solver());
}


#endif        // POLYFIT_TEST_FRAMEWORK_H_
