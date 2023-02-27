#define CGAL_SLS_DEBUG_DRAW

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/extrude_skeleton.h>
#include "print.h"
#include <CGAL/draw_surface_mesh.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Real_timer.h>

#include <boost/shared_ptr.hpp>

#include <fstream>
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

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

bool read_dat_polygon(const char* filename,
                      Polygon_with_holes_2& p)
{
  std::ifstream in(filename);
  if(!in)
  {
    std::cerr << "Error: Could not read " << filename << std::endl;
    return false;
  }

  bool is_number_of_CC_in_input = false;
  if(CGAL::IO::internal::get_file_extension(filename) == "poly")
  {
    is_number_of_CC_in_input = true;
  }

  std::vector<Polygon_2> polys;

  auto read_polygon = [&in, &polys](int i) -> void
  {
    std::vector<Point_2> poly;

    int v_count = 0;
    in >> v_count;
    for(int j=0; j<v_count && in; ++j)
    {
      double x = 0., y = 0.;
      in >> x >> y;
      poly.push_back(Point_2(x, y));
    }

    if(poly.size() >= 3)
    {
      bool is_simple = CGAL::is_simple_2(poly.begin(), poly.end(), K());
      if(!is_simple)
        std::cerr << "Warning: input polygon not simple (hopefully it is strictly simple...)" << std::endl;

      CGAL::Orientation expected = (i == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE);

      const double area = CGAL::to_double(CGAL::polygon_area_2(poly.begin(), poly.end(), K()));
      CGAL::Orientation orientation = area > 0 ? CGAL::COUNTERCLOCKWISE : area < 0 ? CGAL::CLOCKWISE : CGAL::COLLINEAR;

      if(orientation == expected)
        polys.emplace_back(poly.begin(), poly.end());
      else
        polys.emplace_back(poly.rbegin(), poly.rend());
    }
  };

  if(is_number_of_CC_in_input)
  {
    int ccb_count = 0;
    in >> ccb_count;
    for(int i=0; i<ccb_count && in; ++i)
      read_polygon(i);
  }
  else
  {
    int i = 0;
    while(in)
      read_polygon(i++);
  }

  if(polys.empty())
  {
    std::cerr << "Error: empty input?" << std::endl;
    return false;
  }

  std::cout <<"Polygon with border of size: " << polys[0].size() << std::endl;
  if(polys.size() > 1)
    std::cout << polys.size() - 1 << " hole(s)" << std::endl;

  p = Polygon_with_holes_2(polys[0]);
  for(std::size_t i=0; i<polys.size()-1; ++i)
    p.add_hole(polys[i+1]);

  return true;
}

bool read_input_polygon(const char* filename,
                        Polygon_with_holes_2& p)
{
  std::string ext = CGAL::IO::internal::get_file_extension(filename);
  if(ext == "dat")
  {
    return read_dat_polygon(filename, p);
  }
  else
  {
    std::cerr << "Error: unknown file extension: " << ext << std::endl;
    return false;
  }
}

bool read_segment_speeds(const char* filename,
                        std::vector<std::vector<FT> >& weights)
{
  std::ifstream in(filename);
  if(!in)
  {
    std::cerr << "Error: Could not read " << filename << std::endl;
    return false;
  }

  std::vector<FT> border_weights;

  std::string line;
  while(getline(in, line))
  {
    if(line.empty())
    {
      weights.push_back(border_weights);
      border_weights.clear();
    }

    std::istringstream iss(line);
    double w;
    if(iss >> w)
      border_weights.push_back(w);
  }

  // in case the last line is not empty
  if(!border_weights.empty())
    weights.push_back(border_weights);

  return true;
}

// create a random 10-gon in a square
Polygon_with_holes_2 generate_random_polygon(const std::size_t seed)
{
  typedef CGAL::Random_points_in_square_2<Point_2> Point_generator;

  CGAL::Random rnd(seed);

  Polygon_2 poly;
  CGAL::random_polygon_2(10, std::back_inserter(poly), Point_generator(0.5, rnd));
  return Polygon_with_holes_2{poly};
}

void generate_random_weights(const Polygon_with_holes_2& p,
                             const double min_weight,
                             const double max_weight,
                             const unsigned int seed,
                             std::vector<std::vector<FT> >& speeds)
{
  CGAL::Random rnd(seed);
  std::cout << "Seed is " << rnd.get_seed() << std::endl;

  CGAL_assertion(max_weight > 1);

  auto prev = [](const auto& it, const auto& container)
  {
    return it == container.begin() ? std::prev(container.end()) : std::prev(it);
  };

  auto next = [](const auto& it, const auto& container)
  {
    return (it == std::prev(container.end())) ? container.begin() : std::next(it);
  };

  auto generate_range_weights = [prev, next, min_weight, max_weight, &rnd](const auto& c)
  {
    using Container = typename std::remove_reference<decltype(c)>::type;
    using Iterator = typename Container::const_iterator;

    std::map<Iterator, std::size_t /*rnd weight*/> weight;

    // start somewhere not collinear
    Iterator start_it;
    for(Iterator it=c.begin(); it<c.end(); ++it)
    {
      // the edge is [prev_1 ; it], check for collinearity with the previous edge [prev_2; prev_1]
      auto prev_2 = prev(prev(it, c), c);
      auto prev_1 = prev(it, c);
      Segment_2 s0 {*prev_2, *prev_1}, s1 {*prev_1, *it};
      if(!SS::are_edges_orderly_collinear(s0, s1))
        start_it = it;
    }

    CGAL_assertion(start_it != Iterator()); // all collinear is impossible

    Iterator it=start_it, end=start_it;
    do
    {
      auto prev_2 = prev(prev(it, c), c);
      auto prev_1 = prev(it, c);
      Segment_2 s0 {*prev_2, *prev_1}, s1 {*prev_1, *it};
      if(SS::are_edges_orderly_collinear(s0, s1))
      {
        CGAL_assertion(weight.count(prev_1) != 0);
        weight[it] = weight[prev_1];
      }
      else
      {
        CGAL_assertion(weight.count(it) == 0);
        weight[it] = rnd.get_double(min_weight, max_weight);
      }

      it = next(it, c);
    }
    while(it != end);

    std::vector<FT> weights;
    for(auto it=c.begin(); it<c.end(); ++it)
      weights.push_back(weight[it]);

    return weights;
  };

  speeds.push_back(generate_range_weights(p.outer_boundary()));
  for(auto hit=p.holes_begin(); hit!=p.holes_end(); ++hit)
    speeds.push_back(generate_range_weights(*hit));
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  const int argc_check = argc - 1;

  char* poly_filename = nullptr;
  char* speeds_filename = nullptr;

  FT height = FT{std::numeric_limits<double>::max()};
  bool use_angles = false; // whether the input is SLS edge weights, or taper angles
  bool flip_weights = false; // takes the opposite for weights, and the complement for angles

  // below is only used for random weight generation
  double min_weight = 1., max_weight = 10.;
  std::size_t seed = std::time(nullptr);

  for(int i = 1; i < argc; ++i)
  {
    if(!strcmp("-h", argv[i]) || !strcmp("--help", argv[i]) || !strcmp("-?", argv[i]))
    {
      std::cout << "Usage: " << argv[0] << "[options].\n"
        "Options:\n"
        "   -i <input_filename>: input polygon filename.\n"
        "   -t <value>: time (== height). Must be strictly positive.\n"
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

  if(height <= FT(0))
  {
    std::cerr << "Error: height/offset/time must be strictly positive; it is " << height << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::Real_timer timer;
  timer.start();

  Polygon_with_holes_2 pwh;
  if(poly_filename == nullptr)
  {
    pwh = generate_random_polygon(seed);
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
    extrude_skeleton(pwh, height, sm, CGAL::parameters::angles(speeds));
  else
    extrude_skeleton(pwh, height, sm, CGAL::parameters::weights(speeds));

  timer.stop();
  std::cout << "Extrusion computation took " << timer.time() << " s." << std::endl;

  CGAL::draw(sm);

  CGAL::IO::write_polygon_mesh("extruded_skeleton.off", sm, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}