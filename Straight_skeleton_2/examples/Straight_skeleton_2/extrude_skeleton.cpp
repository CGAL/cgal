#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/draw_straight_skeleton_2.h>
#include <CGAL/draw_surface_mesh.h>

#include "CGAL/input_helpers.h" // polygon reading, random polygon with weights generation
#include <CGAL/extrude_skeleton.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Real_timer.h>

#include <boost/shared_ptr.hpp>

#include <iostream>
#include <unordered_map>
#include <vector>

namespace SS = CGAL::CGAL_SS_i;
namespace PMP = CGAL::Polygon_mesh_processing;

// Kernel choice:
// EPICK: Robust and fast
// EPECK_with_sqrt: Exact and slow
// EPECK: More robust, and less slow than EPECK_with_sqrt

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
// using K = CGAL::Exact_predicates_exact_constructions_kernel;
// using K = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;

using FT = K::FT;
using Point_2 = K::Point_2;
using Segment_2 = K::Segment_2;
using Line_2 = K::Line_2;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;

using Polygon_2 = CGAL::Polygon_2<K>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

using Straight_skeleton_2 = CGAL::Straight_skeleton_2<K>;
using Straight_skeleton_2_ptr = boost::shared_ptr<Straight_skeleton_2>;

using Mesh = CGAL::Surface_mesh<Point_3>;

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  const int argc_check = argc - 1;

  char* poly_filename = nullptr;
  char* speeds_filename = nullptr;

  FT height = FT{(std::numeric_limits<double>::max)()};
  bool use_angles = false; // whether the input is SLS edge weights, or taper angles
  bool flip_weights = false; // takes the opposite for weights, and the complement for angles

  // below is only used for random weight generation
  double min_weight = 1., max_weight = 10.;
  int seed = std::time(nullptr);

  for(int i = 1; i < argc; ++i)
  {
    if(!strcmp("-h", argv[i]) || !strcmp("--help", argv[i]) || !strcmp("-?", argv[i]))
    {
      std::cout << "Usage: " << argv[0] << "[options].\n"
        "Options:\n"
        "   -i <input_filename>: input polygon filename.\n"
        "   -t <value>: height. Must be non-zero.\n"
        "   -a <angles_filename>: angles. Format: one angle per line, a space to separate borders.\n"
        "   -w <weights_filename>: weights. Format: one weight per line, a space to separate borders.\n"
        " Note: -w and -a are exclusive.\n"
                << std::endl;

      return EXIT_FAILURE;
    } else if(!strcmp("-i", argv[i]) && i < argc_check) {
      poly_filename = argv[++i];
    } else if(!strcmp("-w", argv[i])) {
      if(speeds_filename != nullptr)
      {
        std::cerr << "Error: -w and -a are exclusive." << std::endl;
        return EXIT_FAILURE;
      }
      speeds_filename = argv[++i];
      use_angles = false;
    } else if(!strcmp("-a", argv[i])) {
      if(speeds_filename != nullptr)
      {
        std::cerr << "Error: -w and -a are exclusive." << std::endl;
        return EXIT_FAILURE;
      }
      speeds_filename = argv[++i];
      use_angles = true;
    } else if(!strcmp("-t", argv[i]) && i < argc_check) {
      height = std::stod(argv[++i]);
    } else if(!strcmp("-mw", argv[i]) && i < argc_check) {
      min_weight = std::stod(argv[++i]);
    } else if(!strcmp("-Mw", argv[i]) && i < argc_check) {
      max_weight = std::stod(argv[++i]);
    } else if(!strcmp("-f", argv[i]) && i < argc) {
      flip_weights = true;
    } else if(!strcmp("-s", argv[i]) && i < argc_check) {
      seed = std::stoi(argv[++i]);
    }
  }

  if(CGAL::is_zero(height))
  {
    std::cerr << "Error: height must be non-zero" << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::Real_timer timer;
  timer.start();

  Polygon_with_holes_2 pwh;
  if(poly_filename == nullptr)
  {
    pwh = generate_random_polygon<Polygon_with_holes_2>(seed);
  }
  else if(!read_input_polygon(poly_filename, pwh) || pwh.outer_boundary().is_empty())
  {
    std::cerr << "Error: failure during polygon read" << std::endl;
    return EXIT_FAILURE;
  }

#ifdef CGAL_SLS_OUTPUT_FILES
  std::ofstream out_poly("input.dat");
  out_poly.precision(17);
  out_poly << pwh;
  out_poly.close();
#endif

  // read segment speeds (angles or weights)
  std::vector<std::vector<FT> > speeds;
  if(speeds_filename == nullptr)
    generate_random_weights(pwh, min_weight, max_weight, seed, speeds);
  else
    read_segment_speeds(speeds_filename, speeds);

  if(flip_weights)
  {
    if(use_angles)
    {
      for(auto& contour_speeds : speeds)
        for(FT& a : contour_speeds)
          a = 180 - a;
    }
    else
    {
      for(auto& contour_speeds : speeds)
        for(FT& w : contour_speeds)
          w = -w;
    }
  }

  timer.stop();
  std::cout << "Reading input(s) took " << timer.time() << " s." << std::endl;

  // End of I/O, do some slope preprocessing and check the validity of the input(s)
  // -----------------------------------------------------------------------------------------------

  timer.reset();
  timer.start();

  Mesh sm;
  if(use_angles)
    CGAL::extrude_skeleton(pwh, sm, CGAL::parameters::angles(speeds).maximum_height(height));
  else
    CGAL::extrude_skeleton(pwh, sm, CGAL::parameters::weights(speeds).maximum_height(height));

  timer.stop();
  std::cout << "Extrusion computation took " << timer.time() << " s." << std::endl;

  CGAL::draw(sm);

  CGAL::IO::write_polygon_mesh("extruded_skeleton.off", sm, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}
