// Copyright (c) 2023 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).

#ifndef CGAL_STRAIGHT_SKELETON_EXTRUSION_INPUT_HELPERS_H
# define CGAL_STRAIGHT_SKELETON_EXTRUSION_INPUT_HELPERS_H

#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>

#include <CGAL/Kernel_traits.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/IO/helpers.h>
#include <CGAL/enum.h>

#include <fstream>
#include <iostream>
#include <vector>

template <typename PolygonWithHoles>
bool read_dat_polygon(const char* filename,
                      PolygonWithHoles& p)
{
  using Polygon_2 = typename PolygonWithHoles::Polygon_2;
  using Point_2 = typename Polygon_2::Point_2;
  using K = typename CGAL::Kernel_traits<Point_2>::type;

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

  p = PolygonWithHoles(polys[0]);
  for(std::size_t i=0; i<polys.size()-1; ++i)
    p.add_hole(polys[i+1]);

  return true;
}

template <typename PolygonWithHoles>
bool read_input_polygon(const char* filename,
                        PolygonWithHoles& p)
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

template <typename FT>
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
template <typename PolygonWithHoles>
PolygonWithHoles generate_random_polygon(const int seed)
{
  using Polygon_2 = typename PolygonWithHoles::Polygon_2;
  using Point_2 = typename Polygon_2::Point_2;

  typedef CGAL::Random_points_in_square_2<Point_2> Point_generator;

  CGAL::Random rnd(seed);

  Polygon_2 poly;
  CGAL::random_polygon_2(10, std::back_inserter(poly), Point_generator(0.5, rnd));
  return PolygonWithHoles{poly};
}

template <typename PolygonWithHoles, typename FT>
void generate_random_weights(const PolygonWithHoles& p,
                             const double min_weight,
                             const double max_weight,
                             const int seed,
                             std::vector<std::vector<FT> >& speeds)
{
  using Polygon_2 = typename PolygonWithHoles::Polygon_2;
  using Point_2 = typename Polygon_2::Point_2;
  using K = typename CGAL::Kernel_traits<Point_2>::type;
  using Segment_2 = typename K::Segment_2;

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

    std::map<Iterator, FT> weight;

    // start somewhere not collinear
    Iterator start_it;
    CGAL_assertion_code(bool found = false);
    for(Iterator it=c.begin(); it<c.end(); ++it)
    {
      // the edge is [prev_1 ; it], check for collinearity with the previous edge [prev_2; prev_1]
      auto prev_2 = prev(prev(it, c), c);
      auto prev_1 = prev(it, c);
      Segment_2 s0 {*prev_2, *prev_1}, s1 {*prev_1, *it};
      if(!CGAL::CGAL_SS_i::are_edges_orderly_collinear(s0, s1)){
        start_it = it;
        CGAL_assertion_code(found = true);
        break;
      }
    }

    CGAL_assertion(found); // all collinear is impossible

    Iterator it=start_it, end=start_it;
    do
    {
      auto prev_2 = prev(prev(it, c), c);
      auto prev_1 = prev(it, c);
      Segment_2 s0 {*prev_2, *prev_1}, s1 {*prev_1, *it};
      if(CGAL::CGAL_SS_i::are_edges_orderly_collinear(s0, s1))
      {
        CGAL_assertion(weight.count(prev_1) != 0);
        weight[it] = weight[prev_1];
      }
      else
      {
        CGAL_assertion(weight.count(it) == 0);
        weight[it] = FT(rnd.get_double(min_weight, max_weight));
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

#endif // CGAL_STRAIGHT_SKELETON_EXTRUSION_INPUT_HELPERS_H
