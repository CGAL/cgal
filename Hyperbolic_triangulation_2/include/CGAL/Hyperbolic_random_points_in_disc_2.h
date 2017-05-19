// Copyright (c) 2010-2016  INRIA Sophia Antipolis, INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: 
// $Id: 
// 
//
// Author(s)     : Mikhail Bogdanov
//                 Monique Teillaud <Monique.Teillaud@inria.fr>

#ifndef HYPERBOLIC_RANDOM_POINTS_IN_DISC_2
#define HYPERBOLIC_RANDOM_POINTS_IN_DISC_2

#include <boost/math/special_functions/acosh.hpp>

#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>

// Euclidean radius vector to hyperbolic raidus vector
template<class FT>
FT r_h(FT r_e)
{
  return boost::math::acosh((1 + r_e*r_e)/(1 - r_e*r_e));
}

// hyperbolic raidus vector to Euclidean radius vector
template<class FT>
FT r_e(FT r_h)
{
  FT dist = std::tanh(r_h/FT(2));
  
  if(dist > 0) {
    return dist;
  }
  return -dist;
}

// if seed = -1, then the seed will get a random value.
template<class Gt>
void Hyperbolic_random_points_in_disc_2(std::vector<typename Gt::Point_2>& output, int nb, int seed = 1, typename Gt::FT e = 0.0001)
{
  typedef typename Gt::FT FT;
  typedef typename Gt::Point_2 Point_2;
  typedef typename Gt::Vector_2 Vector_2;
  
  FT re = FT(1) - e;
  FT rh = r_h(re);
  
  typedef CGAL::Creator_uniform_2<FT, Point_2> Creator;
  
  CGAL::Random rand;
  if (seed != -1) {
    rand = CGAL::Random(seed);
  }
  CGAL::Random_points_in_disc_2<Point_2, Creator> in_Euclidean_disk(rh, rand);
  
  std::vector<Point_2> pts;
  pts.reserve(nb);
  for(int i = 0; i < nb ; i++) {
    pts.push_back(*in_Euclidean_disk);
    in_Euclidean_disk++;
  }
  
  for(int i = 0; i < nb ; i++) {
    Vector_2 v = Vector_2(Point_2(0, 0), pts[i]);
    
    FT sq_dist = v.squared_length();
    FT dist = CGAL::sqrt(sq_dist);
    
    FT dist_in_disc = r_e(dist);
    
    output.push_back(Point_2((pts[i].x()*dist_in_disc)/dist, (pts[i].y()*dist_in_disc)/dist));
  }
}

template<class Gt>
void Random_points_in_disc_2(std::vector<typename Gt::Point_2>& output, int nb, int seed = 1, typename Gt::FT e = 0.0001)
{
  typedef typename Gt::FT FT;
  typedef typename Gt::Point_2 Point_2;
  
  FT re = FT(1) - e;
  
  typedef CGAL::Creator_uniform_2<FT, Point_2> Creator;
  CGAL::Random rand(seed);
  CGAL::Random_points_in_disc_2<Point_2, Creator> in_Euclidean_disk(re, rand);
  
  for(int i = 0; i < nb ; i++) {
    output.push_back(*in_Euclidean_disk);
    in_Euclidean_disk++;
  }
}

#endif // HYPERBOLIC_RANDOM_POINTS_IN_DISC_2
