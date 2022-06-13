// Copyright (c) 2022 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot
//
//******************************************************************************
//
//******************************************************************************


#ifndef CGAL_MESH_3_TRIPLE_LINES_H
#define CGAL_MESH_3_TRIPLE_LINES_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/number_utils.h>

#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>

#include <vector>

namespace CGAL
{
namespace Mesh_3
{

template<typename P, typename Functor>
std::vector<P> create_polyline(const double start,
    const double end,
    P starting_point,
    P ending_point,
    Functor f,
    const int prec)
{
    using Polyline = std::vector<P>;
    Polyline poly;
    poly.reserve(prec + 1);
    poly.push_back(std::move(starting_point));
#ifdef CGAL_DEBUG_TRIPLE_LINES
    std::cerr << "new poly\n";
    std::cerr << "  poly.push_back(" << poly.back() << "\n";
#endif // CGAL_DEBUG_TRIPLE_LINES
    const double step = (end - start) / prec;
    const double stop = end - step / 2;
    const bool step_is_positive = (step > 0);
    for (double t = start + step;
        (step_is_positive ? t < stop : t > stop);
        t += step)
    {
        poly.push_back(f(t));
#ifdef CGAL_DEBUG_TRIPLE_LINES
        std::cerr << "  poly.push_back(" << poly.back() << ")\n";
#endif // CGAL_DEBUG_TRIPLE_LINES
    }
    poly.push_back(std::move(ending_point));
#ifdef CGAL_DEBUG_TRIPLE_LINES
    std::cerr << "  poly.push_back(" << poly.back() << ")\n";
#endif // CGAL_DEBUG_TRIPLE_LINES
    return poly;
}

template <typename P, typename Functor>
std::vector<P> create_polyline(const double start,
    const double end,
    P starting_point,
    Functor f,
    const int prec)
{
    return create_polyline(start, end, starting_point, f(end), f, prec);
}

template <typename P, typename Functor>
std::vector<P> create_polyline(const double start,
    const double end,
    Functor f,
    const int prec)
{
    return create_polyline<P>(start, end, f(start), f(end), f, prec);
}


/*
      c6 ------- c7
     /|         /|
    / |        / |
   /  |       /  |
  c4 ------- c5  |
  |   c2 ----|-- c3
  |  /       |  /
  | /        | /
  |/         |/
  c0 ------- c1

  ci = color at the corner

 The coordinates of the unit cube are given by:

     for(int z = 0; z <=1 ; ++z)
       for(int y = 0; y <=1 ; ++y)
         for(int x = 0; x <=1 ; ++x)
           cube.emplace_back(x, y, z);
*/


//
// One internal corner
//

  // 00001221
  //   corner is (1/2,1/2,2/3)
  //   curve_1  : x=1/2, y=[0,1/2], z=2/3
  //   curve_1' : x=1/2, y=[1/2,1], z=2/3
  //   curve_2  : x=[0,1/2], y=1/2, z=2/3
  //   curve_2' : x=[1/2,1], y=1/2, z=2/3
  //   curve_3  : x=1/2, y=1/2, z=[2/3,1]
template<typename P>
std::vector<std::vector<P>> poly00001221(const int /* no sampling for segments */)
{
    P corner{ 1. / 2, 1. / 2, 2. / 3 };
    P      a{ 1. / 2,    0, 2. / 3 };
    P      b{ 1. / 2,    1, 2. / 3 };
    P      c{ 0,    1. / 2, 2. / 3 };
    P      d{ 1,    1. / 2, 2. / 3 };
    P      e{ 1. / 2,    1. / 2, 1 };

    return
    {
      { a, corner}, // segment curve_1
      { corner, b}, // segment curve_1'
      { corner, c}, // segment curve_2
      { corner, d}, // segment curve_2'
      { corner, e}, // segment curve_3
    };
}
// 00111202
//   corner is (1/2,1/2,2/3)
//   curve_1  : x=1/2,y=[0,1/2],z=2/3
//   curve_1' : x=1/2,y=[1/2,1],z=2/3
//   curve_2  : x=1/(3*z),y=1/2,z=[2/3,1]
//   ADDED   curve_2' : x=1/(3*z),y=1/2,z=[1/3,2/3]
//   curve_3  : x=(2*z−1)/z,y=1/2, z=[1/2,2/3]
//   REMOVED curve_3' : x=(2*z−1)/z,y=1/2, z=[2/3,1]
template<typename P>
std::vector<std::vector<P>> poly00111202(const int prec = 10)
{
    P      a{ 1. / 2, 0.  , 2. / 3 };
    P corner{ 1. / 2, 1. / 2, 2. / 3 };
    P      b{ 1. / 2, 1.  , 2. / 3 };

    auto f_curve_2 = [](double z) { return P(1. / (3 * z), 1. / 2, z); };
    auto f_curve_3 = [](double z) { return P((2. * z - 1) / z, 1. / 2, z); };

    return
    {
      { a, corner}, // segment curve_1
      { corner, b}, // segment curve_1'
      create_polyline(2. / 3, 1   , corner, f_curve_2, prec), // curve_2
      create_polyline(2. / 3, 1. / 3, corner, f_curve_2, prec), // curve_2' ADDED
      create_polyline(2. / 3, 1. / 2, corner, f_curve_3, prec), // curve_3
      // create_polyline(2./3, 1.  , corner, f_curve_3, prec), // curve_3' REMOVED
    };
}

// 01101001
//   corner is (1/2,1/2,1/2)
//   curve_1  : x=1/2, y=1/2, z=[0,1/2]
//   curve_1' : x=1/2, y=1/2, z=[1/2,1]
//   curve_2  : x=1/2, y=[0,1/2], z=1/2
//   curve_2' : x=1/2, y=[1/2,1], z=1/2
//   curve_3  : x=[0,1/2], y=1/2, z=1/2
//   curve_3' : x=[1/2,1], y=1/2, z=1/2
template<typename P>
std::vector<std::vector<P>> poly01101001(const int /* no sampling for segments */)
{
    P corner{ 1. / 2, 1. / 2, 1. / 2 };
    P      a{ 1. / 2, 1. / 2,    0 };
    P      b{ 1. / 2, 1. / 2,    1 };
    P      c{ 1. / 2,    0, 1. / 2 };
    P      d{ 1. / 2,    1, 1. / 2 };
    P      e{ 0,    1. / 2, 1. / 2 };
    P      f{ 1,    1. / 2, 1. / 2 };

    return
    {
      { a, corner}, // segment curve_1
      { corner, b}, // segment curve_1'
      { corner, c}, // segment curve_2
      { corner, d}, // segment curve_2'
      { corner, e}, // segment curve_3
      { corner, f}, // segment curve_3'
    };
}

//
// Two curves
//
template<typename P>
std::vector<std::vector<P>> poly00011022(const int prec = 10)
{
    // x = (3*z^2-2*z)/(3*z^2-1), y = 1/(3*z), z = [1/3,1/2]
    // x = (3*z^2-2*z)/(3*z^2-1), y = 1/(3*z), z = [2/3,1]
    auto f = [](double z) {return P(z * (3 * z - 2) / (3 * z * z - 1),
        1 / (3 * z),
        z); };
    return {
      create_polyline<P>(1. / 3, 1. / 2, f, prec),
      create_polyline<P>(2. / 3, 1.  , f, prec),
    };
}

// 00011221
//   curve_1 : x = ]1/2,1], y = (3 * x * x - sqrt(9 * x * x * x * x - 30 * x * x * x + 45 * x * x - 24 * x + 4) + 3 * x - 2)/(6 * x * (2 * x - 1)), z = 1./(3*x+3*y-6*x*y)
//   point limit x = 1/2, y = 0, z = 2/3
//   curve_2 : x = ]0,1/3], y = (3 * x * x + sqrt(9 * x * x * x * x - 30 * x * x * x + 45 * x * x - 24 * x + 4) + 3 * x - 2)/(6 * x * (2 * x - 1)), z = 1./(3*x+3*y-6*x*y)
//   point limit x = 0, y = 1/2, z = 2/3
template<typename P>
std::vector<std::vector<P>> poly00011221(const int prec = 10)
{
    auto sq_exp = [](double x) {
        return sqrt(9 * x * x * x * x - 30 * x * x * x + 45 * x * x - 24 * x + 4);
    };
    auto y1 = [sq_exp](double x) { return (3 * x * x - sq_exp(x) + 3 * x - 2) / (6 * x * (2 * x - 1)); };
    P corner1 = { 1. / 2, 0, 2. / 3 };
    auto y2 = [sq_exp](double x) { return (3 * x * x + sq_exp(x) + 3 * x - 2) / (6 * x * (2 * x - 1)); };
    P corner2 = { 0, 1. / 2, 2. / 3 };
    auto z = [](double x, double y) { return 1. / (3 * x + 3 * y - 6 * x * y); };
    return {
      create_polyline(1. / 2, 1, corner1,
                      [y1, z](double x) { return P(x, y1(x), z(x, y1(x))); },
                      prec),
      create_polyline(0, 1. / 3, corner2,
                      [y2, z](double x) { return P(x, y2(x), z(x, y2(x))); },
                      prec)
    };
}

// 00011222
//   curve_1 : x = ]1/2,1[, y = (3 * x * x - sqrt(9 * x * x * x * x - 18 * x * x * x + 25 *x * x - 16 * x + 4) + x - 2)/(6 * (x - 1) * x), z = 1./(3*x+3*y-3*x*y)
//   point limit x = 1/2, y = 1, z=1/3
//   point limit x = 1, y = 1/2, z=1/3
//   curve_2 : x = ]0,1/2], y = (3 * x * x + sqrt(9 * x * x * x * x - 18 * x * x * x + 25 *x * x - 16 * x + 4) + x - 2)/(6 * (x - 1) * x), z = 1./(3*x+3*y-3*x*y)
//   point limit x = 0, y=1/2, z=2/3
template<typename P>
std::vector<std::vector<P>> poly00011222(const int prec = 10)
{
    auto sq_exp = [](double x) {
        return sqrt(9 * x * x * x * x - 18 * x * x * x + 25 * x * x - 16 * x + 4);
    };
    auto y1 = [sq_exp](double x) { return (3 * x * x - sq_exp(x) + x - 2) / (6 * (x - 1) * x); };
    P corner1 = { 1. / 2, 1, 1. / 3 };
    P corner1_bis = { 1, 1. / 2, 1. / 3 };
    auto y2 = [sq_exp](double x) { return (3 * x * x + sq_exp(x) + x - 2) / (6 * (x - 1) * x); };
    P corner2 = { 0, 1. / 2, 2. / 3 };
    auto z = [](double x, double y) { return 1. / (3 * x + 3 * y - 3 * x * y); };
    return {
      create_polyline(1. / 2, 1, corner1, corner1_bis,
                      [y1, z](double x) { return P(x, y1(x), z(x, y1(x))); },
                      prec),
      create_polyline(0, 1. / 2, corner2,
                      [y2, z](double x) { return P(x, y2(x), z(x, y2(x))); },
                      prec)
    };
}

// 00121200
//   curve_1  : x = 1/2, y = (3*z-2)/(6*z-3), z = [0,1/3]
//   curve_1' : x = 1/2, y = (3*z-2)/(6*z-3), z = [2/3,1]
template<typename P>
std::vector<std::vector<P>> poly00121200(const int prec = 10)
{
    auto y = [](double z) { return (3 * z - 2) / (6 * z - 3); };
    return {
      create_polyline<P>(0, 1. / 3,
                      [y](double z) { return P(1. / 2, y(z), z); },
                      prec),
      create_polyline<P>(2. / 3, 1,
                      [y](double z) { return P(1. / 2, y(z), z); },
                      prec)
    };
}

// 00121221
//   curve_1 : x = 1/2, y = (3*z-2)/(3*z-3), z = [0,2/3]
//   curve_2 : x = 1/2, y = z/(3*z-1), z = [1/2,1]
template<typename P>
std::vector<std::vector<P>> poly00121221(const int prec = 10)
{
    auto y1 = [](double z) { return (3 * z - 2) / (3 * z - 3); };
    auto y2 = [](double z) { return z / (3 * z - 1); };
    return {
      create_polyline<P>(0, 2. / 3,
                      [y1](double z) { return P(1. / 2, y1(z), z); },
                      prec),
      create_polyline<P>(1. / 2, 1,
                      [y2](double z) { return P(1. / 2, y2(z), z); },
                      prec)
    };
}

// 00122100
//   curve_1  : x = 1/2, y = (3*z-2)/(6*z-3), z = [0,1/3]
//   curve_1' : x = 1/2, y = (3*z-2)/(6*z-3), z = [2/3,1]
template<typename P>
std::vector<std::vector<P>> poly00122100(const int prec = 10)
{
    auto y = [](double z) { return (3 * z - 2) / (6 * z - 3); };
    return {
      create_polyline<P>(0, 1. / 3,
                      [y](double z) { return P(1. / 2, y(z), z); },
                      prec),
      create_polyline<P>(2. / 3, 1,
                      [y](double z) { return P(1. / 2, y(z), z); },
                      prec)
    };
}

// 00122101
//   curve_1 : x = [1/3,1/2] , y = (-sqrt(24 x^3 - 35 x^2 + 18 x - 3) - 5 x + 3)/(6 (x - 1)^2), x^2 - 3 x + 1!=0,
//                             z = (-sqrt(24 x^3 - 35 x^2 + 18 x - 3) - 7 x + 3)/(6 (x^2 - 3 x + 1))
//   curve_2 : x = [1./2,1[ ,  y = ( sqrt(24 x^3 - 35 x^2 + 18 x - 3) - 5 x + 3)/(6 (x - 1)^2),
//                             z = ( sqrt(24 x^3 - 35 x^2 + 18 x - 3) - 7 x + 3)/(6 (x^2 - 3 x + 1))
//   point limit of curve_2 when x -> 1 x = 1, y = 1/2, z = 1/3
template<typename P>
std::vector<std::vector<P>> poly00122101(const int prec = 10)
{
    auto sq_exp = [](double x) {
        return sqrt(24 * x * x * x - 35 * x * x + 18 * x - 3);
    };
    auto y1 = [sq_exp](double x)
      { return (-sq_exp(x) - 5 * x + 3) / (6 * (x - 1) * (x - 1)); };
    auto z1 = [sq_exp](double x)
      { return (-sq_exp(x) - 7 * x + 3) / (6 * (x * x - 3 * x + 1)); };
    auto y2 = [sq_exp](double x)
      { return (sq_exp(x) - 5 * x + 3) / (6 * (x - 1) * (x - 1)); };
    auto z2 = [sq_exp](double x)
      { return (sq_exp(x) - 7 * x + 3) / (6 * (x * x - 3 * x + 1)); };
    P corner{ 1., .5, 1. / 3 };
    return {
      create_polyline<P>(1. / 3, 1. / 2,
                      [y1, z1](double x) { return P(x, y1(x), z1(x)); },
                      prec),
      create_polyline<P>(1., 1. / 2, corner,
                      [y2, z2](double x) { return P(x, y2(x), z2(x)); },
                      prec),
    };
}

//
// One curve
//
template<typename P>
std::vector<std::vector<P>> poly00000012(const int prec = 10)
{
    // curve : x = 1/2, y = 2/(3*z), z = [2/3,1]
    return { create_polyline<P>(2. / 3, 1.,
                             [](double z) { return P(0.5, 2. / (3. * z), z); },
                             prec) };
}

// 00000112
//   x =[1/2,1], y = x/(3 * x - 1), z = (3 * x - 1)/(3 * x * x)
template<typename P>
std::vector<std::vector<P>> poly00000112(const int prec = 10)
{
    return { create_polyline<P>(1. / 2, 1,
                             [](double x) { return P(x, x / (3 * x - 1), (3 * x - 1) / (3 * x * x)); },
                             prec) };
}

// 00000121
//   curve : x = 1/(3*z), y = 1/(3*z-1), z = [2/3,1]
template<typename P>
std::vector<std::vector<P>> poly00000121(const int prec = 10)
{
    return { create_polyline<P>(2. / 3, 1,
                             [](double z) { return P(1 / (3 * z), 1 / (3 * z - 1), z); },
                             prec) };
}

// 00001112
//   curve : x = 1/(2*y), y = [1/2, 1],z = 2/3
template<typename P>
std::vector<std::vector<P>> poly00001112(const int prec = 10)
{
    return { create_polyline<P>(1. / 2, 1,
                             [](double y) { return P(1 / (2 * y), y, 2. / 3); },
                             prec) };
}

// 00001122
//   curve : x = [0,1], y = 1/2, z = 2/3
template<typename P>
std::vector<std::vector<P>> poly00001122(const int prec = 10)
{
    return { create_polyline<P>(0, 1,
                             [](double x) { return P(x, 1. / 2, 2. / 3); },
                             prec) };
}

// 00010121
//   curve : x =y * z / (z+y), y = ((3 * z * z - 1) - sqrt(CGAL::square(1 - 3 * z * z) - 12 * (z - 1) * z * z))/(6 * (z - 1) * z), z=[1,1/2[
//   point limit =  (1/3, 1/2, 1)
template<typename P>
std::vector<std::vector<P>> poly00010121(const int prec = 10)
{
    auto y = [](double z) { return ((3 * z * z - 1) - sqrt(CGAL::square(1 - 3 * z * z) - 12 * (z - 1) * z * z)) / (6 * (z - 1) * z); };
    auto x = [](double y, double z) { return y * z / (z + y); };
    P corner(1. / 3, 1. / 2, 1);
    return { create_polyline<P>(1, 1. / 2, corner,
                             [x, y](double z) { return P(x(y(z), z), y(z), z); },
                             prec) };
}

// 00010122
//   curve : x = z/(3*z^2-2*z+1), y = 1/(3*z), z = [1/3,1]
template<typename P>
std::vector<std::vector<P>> poly00010122(const int prec = 10)
{
    return { create_polyline<P>(1. / 3, 1,
                             [](double z) { return P(z / (3 * z * z - 2 * z + 1), 1. / (3 * z), z); },
                             prec) };
}
// 00011002
//   curve : x = (y+1)/(4*y-1), y = [2/3,1], z = 1/2
template<typename P>
std::vector<std::vector<P>> poly00011002(const int prec = 10)
{
    return { create_polyline<P>(2. / 3, 1,
                             [](double y) { return P((y + 1) / (4 * y - 1), y, 1. / 2); },
                             prec) };
}
// 00011012
//   curve : x = (3*z^2-2*z+1)/(3*z^2), y = z/(3*z^2-2*z+1), z = [1/2, 1]
template<typename P>
std::vector<std::vector<P>> poly00011012(const int prec = 10)
{
    return { create_polyline<P>(1. / 2, 1,
                             [](double z) { return P((3 * z * z - 2 * z + 1) / (3 * z * z), z / (3 * z * z - 2 * z + 1), z); },
                             prec) };
}
// 00011110
//   curve : x = 1/(2*y), y = [1/2,1], z = 1/2
template<typename P>
std::vector<std::vector<P>> poly00011110(const int prec = 10)
{
    return { create_polyline<P>(1. / 2, 1,
                             [](double y) { return P(1. / (2 * y), y, 1. / 2); },
                             prec) };
}
// 00011120
//   curve : x = (3*z^2-2*z)/(3*z^2-1), y = (3*z^2-1)/(6*z^2-3*z), z = [2/3,1]
template<typename P>
std::vector<std::vector<P>> poly00011120(const int prec = 10)
{
    return { create_polyline<P>(2. / 3, 1,
                             [](double z) { return P((3 * z * z - 2 * z) / (3 * z * z - 1), (3 * z * z - 1) / (6 * z * z - 3 * z), z); },
                             prec) };
}
// 00011121
//   curve : x = (3*z^2-2*z)/(3*z^2-z-1), y = (3*z^2-z-1)/(3*z^2-3*z), z =[1/2,2./3]
template<typename P>
std::vector<std::vector<P>> poly00011121(const int prec = 10)
{
    return { create_polyline<P>(1. / 2, 2. / 3,
                             [](double z) { return P((3 * z * z - 2 * z) / (3 * z * z - z - 1), (3 * z * z - z - 1) / (3 * z * z - 3 * z), z); },
                             prec) };
}
// 00011122
//   curve : x = (3*z*z-2*z)/(z-1),y = 1/(3*z),z = [1/3, 2/3]
//
template<typename P>
std::vector<std::vector<P>> poly00011122(const int prec = 10)
{
    return { create_polyline<P>(1. / 3,2. / 3,
                             [](double z) { return P((3 * z * z - 2 * z) / (z - 1), 1 / (3 * z), z); },
                             prec) };
}

// 00011220
//   curve : x=[0,1/2], y = (2 * x - 1)/(3 * x - 2), z = (2 * (x^2 - 2 * x + 1))/(5 * x^2 - 7 * x + 3)
template<typename P>
std::vector<std::vector<P>> poly00011220(const int prec = 10)
{
    return { create_polyline<P>(0, 1. / 2,
                             [](double x) { return P(x, (2 * x - 1) / (3 * x - 2), (2 * (x * x - 2 * x + 1)) / (5 * x * x - 7 * x + 3)); },
                             prec) };
}
// 00012002
//   curve_1 : x = [2/3,1], y = (3 * x*x + sqrt(9 * x*x*x*x - 24 * x*x*x + 30 * x*x - 12 * x + 1) - 1)/(6 * x * (2 * x - 1)), z = 1 - 1./(3*x*y)
template<typename P>
std::vector<std::vector<P>> poly00012002(const int prec = 10)
{
    auto y = [](double x) { return (3 * x * x + sqrt(9 * x * x * x * x - 24 * x * x * x + 30 * x * x - 12 * x + 1) - 1) / (6 * x * (2 * x - 1)); };
    return { create_polyline<P>(2. / 3, 1,
                             [y](double x) { return P(x, y(x), 1 - 1. / (3 * x * y(x))); },
                             prec) };
}
// 00012012
//   curve_1 : x = [0, 1/2[ , y = 1./(-6*x*z+3*x+3*z), z = (3 *x*x + sqrt(9 *x*x*x*x - 18 *x*x*x + 25 *x*x - 16 * x + 4) - 7 * x + 2)/(6 *(x - 1) * (2 * x - 1))
//   curve_1 : x = ]1/2,1[, y = 1./(-6*x*z+3*x+3*z), z = (3 *x*x + sqrt(9 *x*x*x*x - 18 *x*x*x + 25 *x*x - 16 * x + 4) - 7 * x + 2)/(6 *(x - 1) * (2 * x - 1))
//  point limit x = 1, y = 2/3, z = 1/2
//  point limit x = 1/2, y= 2./3, z = 2./3);
//
template<typename P>
std::vector<std::vector<P>> poly00012012(const int prec = 10)
{
    auto y = [](double x, double z) { return 1. / (-6 * x * z + 3 * x + 3 * z); };
    auto z = [](double x) { return (3 * x * x + sqrt(9 * x * x * x * x - 18 * x * x * x + 25 * x * x - 16 * x + 4) - 7 * x + 2) / (6 * (x - 1) * (2 * x - 1)); };
    P corner1(1. / 2, 2. / 3, 2. / 3);
    P corner2(1, 2. / 3, 1. / 2);
    return { create_polyline<P>(1. / 2, 0, corner1,
                             [y, z](double x) { return P(x, y(x, z(x)), z(x)); },
                             prec),
             create_polyline<P>(1. / 2, 1, corner1, corner2,
                                      [y, z](double x) { return P(x, y(x, z(x)), z(x)); },
                                      prec)
    };
}
// 00012021
//   curve : x = (3*z-1)/(3*z), y = z/(3*z-1), z = [1/2,1]
template<typename P>
std::vector<std::vector<P>> poly00012021(const int prec = 10)
{
    return { create_polyline<P>(1. / 2, 1,
                             [](double z) { return P((3 * z - 1) / (3 * z), z / (3 * z - 1), z); },
                             prec) };
}
// 00012110
//   curve : x=]0,1/2], y = (3 * x * x + sqrt(9 * x * x * x * x - 18 * x * x * x + 25 * x * x - 16 * x + 4) + x - 2)/(6 * (x - 1) * x), z = (3*x*y-1)/(9*x*y-3*x-3*y)
//   point limit : x = 0, y = 1/2, z = 2/3
//
template<typename P>
std::vector<std::vector<P>> poly00012110(const int prec = 10)
{
    auto y = [](double x) { return (3 * x * x + sqrt(9 * x * x * x * x - 18 * x * x * x + 25 * x * x - 16 * x + 4) + x - 2) / (6 * (x - 1) * x); };
    auto z = [](double x, double y) { return (3 * x * y - 1) / (9 * x * y - 3 * x - 3 * y); };
    P corner(0, 1. / 2, 2. / 3);
    return { create_polyline<P>(0, 1. / 2, corner,
                             [y, z](double x) { return P(x, y(x), z(x, y(x))); },
                             prec) };
}
// 00012112
//   x = ]0,1/2[, y = (3 * x*x + sqrt(9 * x * x * x * x - 36 * x * x * x + 40 * x * x - 20 * x + 4) + 2 * x - 2)/(6 * x * (2 * x - 1)), z = (3*x*y-1)/(9*x*y-3*x-3*y)
//   point limit x = 0, y = 1/2, z = 2/3
//   point limit x = 1/2, y = 0, z = 2/3
//
template<typename P>
std::vector<std::vector<P>> poly00012112(const int prec = 10)
{
    auto y = [](double x) { return (3 * x * x + sqrt(9 * x * x * x * x - 36 * x * x * x + 40 * x * x - 20 * x + 4) + 2 * x - 2) / (6 * x * (2 * x - 1)); };
    auto z = [](double x, double y) { return (3 * x * y - 1) / (9 * x * y - 3 * x - 3 * y); };
    P corner1(0, 1. / 2, 2. / 3);
    P corner2(1. / 2, 0, 2. / 3);
    return { create_polyline<P>(0, 1. / 2, corner1, corner2,
                             [y, z](double x) { return P(x, y(x), z(x, y(x))); },
                             prec) };
}

// 00012120
//   curve : x = (3*z-1)/(3*z), y = (3*z^2-2*z)/(6*z^2-5*z+1), z = [2/3,1]
template<typename P>
std::vector<std::vector<P>> poly00012120(const int prec = 10)
{
    return { create_polyline<P>(2. / 3, 1,
                             [](double z) { return P((3 * z - 1) / (3 * z), (3 * z * z - 2 * z) / (6 * z * z - 5 * z + 1), z); },
                             prec) };
}
//
// 00012121
//  curve : x = (3*z-1)/(3*z), y = (3*z^2-2*z)/(3*z^2-4*z+1), z = [1/2,2/3]
template<typename P>
std::vector<std::vector<P>> poly00012121(const int prec = 10)
{
    return { create_polyline<P>(1. / 2, 2. / 3,
                             [](double z) { return P((3 * z - 1) / (3 * z), (3 * z * z - 2 * z) / (3 * z * z - 4 * z + 1), z); },
                             prec) };
}
// 00012122
//  curve : x = (6*z^2-6*z+1)/(3*z^2-3*z), y = (3*z^2-2*z)/(6*z^2-6*z+1), z = [1/3,2/3]
template<typename P>
std::vector<std::vector<P>> poly00012122(const int prec = 10)
{
    auto x = [](double z) { return (6 * z * z - 6 * z + 1) / (3 * z * z - 3 * z); };
    auto y = [](double z) { return (3 * z * z - 2 * z) / (6 * z * z - 6 * z + 1); };
    return { create_polyline<P>(1. / 3, 2. / 3,
                             [x,y](double z) { return P(x(z), y(z), z); },
                             prec) };
}
// 00012221
//   curve : x = 1/(3*y),y = [1/3,1],z = 1/2
template<typename P>
std::vector<std::vector<P>> poly00012221(const int prec = 10)
{
    return { create_polyline<P>(1. / 3, 1,
                             [](double y) { return P(1. / (3 * y),y , 1. / 2); },
                             prec) };
}

// 00111100
//   curve : x = [0,1], y = 1/2, z = 1/2
template<typename P>
std::vector<std::vector<P>> poly00111100(const int prec = 10)
{
    return { create_polyline<P>(0, 1,
                             [](double x) { return P(x , 1. / 2, 1. / 2); },
                             prec) };
}

// 00111102
//   curve : x = (2*z-1)/(3*z^2-z), y = (3*z-1)/(6*z-3), z = [2/3,1]
template<typename P>
std::vector<std::vector<P>> poly00111102(const int prec = 10)
{
    return { create_polyline<P>(2. / 3, 1,
                             [](double z) { return P((2 * z - 1) / (3 * z * z - z), (3 * z - 1) / (6 * z - 3), z); },
                             prec) };
}

// 00111220
//   segment 1/2 0 2/3 1/2 1 2/3
template<typename P>
std::vector<std::vector<P>> poly00111220(const int /*not needed for a segment*/)
{
    return { { P(1. / 2, 0, 2. / 3), P(1. / 2, 1, 2. / 3) } };
}

// 00121201
//   curve_1 : x =[1/2, 2/3], y = (-sqrt(-24 * x^3 + 37 * x^2 - 20 * x + 4) + 5 * x - 2)/(6 * x^2), z = ( sqrt(-24 * x^3 + 37 * x^2 - 20 * x + 4) + 5 * x - 2)/(6 * x^2)
//   curve_2 : x =[1/2, 2/3], y = ( sqrt(-24 * x^3 + 37 * x^2 - 20 * x + 4) + 5 * x - 2)/(6 * x^2), z = (-sqrt(-24 * x^3 + 37 * x^2 - 20 * x + 4) + 5 * x - 2)/(6 * x^2)
//   point 0 1/2 1/2
template<typename P>
std::vector<std::vector<P>> poly00121201(const int prec = 10)
{
    auto sq_exp = [](double x) {
        return sqrt(-24 * x * x * x + 37 * x * x - 20 * x + 4);
    };
    auto y = [sq_exp](double x) { return (-sq_exp(x) + 5 * x - 2) / (6 * x * x); };
    auto z = [sq_exp](double x) { return (sq_exp(x) + 5 * x - 2) / (6 * x * x); };
    P corner{ 2. / 3, y(2. / 3), z(2. / 3) };
    return {
      create_polyline<P>(2. / 3, 1. / 2, corner,
                      [y, z](double x) { return P(x, y(x), z(x)); },
                      prec),
      create_polyline<P>(2. / 3, 1. / 2, corner,
                      [y, z](double x) { return P(x, z(x), y(x)); },
                      prec),
                      // { P(0., .5, .5), P(0., .5, .5) }
    };
}

// 00000123
// curve_1 : x = 3./4. - sqrt(z*(9.*z - 8))/(4*z), y = 3./4. - sqrt(z*(9.*z - 8))/(4*z), z = [8/9, 1]
// curve 2 : x = 3./4. + sqrt(z*(9.*z - 8))/(4*z), y = 3./4. + sqrt(z*(9.*z - 8))/(4*z), z = [8/9, 1]
// curve 3 : x =1/2, y = -(z-2.)/(2.*z), z = [2/3, 1]
// curve 4 : x = -(z/2.-1)/z, y = 1./2., z =[2/3, 1]
template<typename P>
std::vector<std::vector<P>> poly00000123(const int prec = 10)
{
  return {
  create_polyline(8./9, 1., P(3./4, 3./4, 8./9),
                  [](double z) { return P(  3./4. - sqrt(z*(9.*z - 8))/(4*z), 3./4. - sqrt(z*(9.*z - 8))/(4*z), z ); },
                  prec),
  create_polyline(8./9, 1., P(3./4, 3./4, 8./9),
                  [](double z) { return P(  3./4. + sqrt(z*(9.*z - 8))/(4*z), 3./4. + sqrt(z*(9.*z - 8))/(4*z), z ); },
                  prec),

  create_polyline(2./3, 1., P(1./2, 1, 2./3),
                  [](double z) { return P(  1./2, -(z-2.)/(2.*z), z); },
                  prec),
  create_polyline(2./3, 1., P(1, 1./2, 2./3),
                  [](double z) { return P(-(z-2.)/(2.*z), 1./2, z); },
                  prec),
};
}


//00001123
//
//
template<typename P>
std::vector<std::vector<P>> poly00001123(const int /*prec*/ = 10)
{
return {


};
}

//00001223
//
//
template<typename P>
std::vector<std::vector<P>> poly00001223(const int /*prec*/ = 10)
{
return {


};
}

//00010123
//
//
template<typename P>
std::vector<std::vector<P>> poly00010123(const int /*prec*/ = 10)
{
return {


};
}

//00010230
//
//
template<typename P>
std::vector<std::vector<P>> poly00010230(const int /*prec*/ = 10)
{
return {


};
}

//00010231
//
//
template<typename P>
std::vector<std::vector<P>> poly00010231(const int /*prec*/ = 10)
{
return {


};
}

//00011023
//
//
template<typename P>
std::vector<std::vector<P>> poly00011023(const int /*prec*/ = 10)
{
return {


};
}

//00011123
//
//
template<typename P>
std::vector<std::vector<P>> poly00011123(const int /*prec*/ = 10)
{
return {


};
}

//00011223
//
//
template<typename P>
std::vector<std::vector<P>> poly00011223(const int /*prec*/ = 10)
{
return {


};
}

//00011230
//
//
template<typename P>
std::vector<std::vector<P>> poly00011230(const int /*prec*/ = 10)
{
return {


};
}

//00011231
//
//
template<typename P>
std::vector<std::vector<P>> poly00011231(const int /*prec*/ = 10)
{
return {


};
}

//00011232
//
//
template<typename P>
std::vector<std::vector<P>> poly00011232(const int /*prec*/ = 10)
{
return {


};
}

//00012003
//
//
template<typename P>
std::vector<std::vector<P>> poly00012003(const int /*prec*/ = 10)
{
return {


};
}

//00012013
//
//
template<typename P>
std::vector<std::vector<P>> poly00012013(const int /*prec*/ = 10)
{
return {


};
}

//00012023
//
//
template<typename P>
std::vector<std::vector<P>> poly00012023(const int /*prec*/ = 10)
{
return {


};
}

//00012033
//
//
template<typename P>
std::vector<std::vector<P>> poly00012033(const int /*prec*/ = 10)
{
return {


};
}

//00012113
//
//
template<typename P>
std::vector<std::vector<P>> poly00012113(const int /*prec*/ = 10)
{
return {


};
}

//00012123
//
//
template<typename P>
std::vector<std::vector<P>> poly00012123(const int /*prec*/ = 10)
{
return {


};
}

//00012130
//
//
template<typename P>
std::vector<std::vector<P>> poly00012130(const int /*prec*/ = 10)
{
return {


};
}

//00012131
//
//
template<typename P>
std::vector<std::vector<P>> poly00012131(const int /*prec*/ = 10)
{
return {


};
}

//00012132
//
//
template<typename P>
std::vector<std::vector<P>> poly00012132(const int /*prec*/ = 10)
{
return {


};
}

//00012133
//
//
template<typename P>
std::vector<std::vector<P>> poly00012133(const int /*prec*/ = 10)
{
return {


};
}

//00012223
//
//
template<typename P>
std::vector<std::vector<P>> poly00012223(const int /*prec*/ = 10)
{
return {


};
}

//00012230
//
//
template<typename P>
std::vector<std::vector<P>> poly00012230(const int /*prec*/ = 10)
{
return {


};
}

//00012231
//
//
template<typename P>
std::vector<std::vector<P>> poly00012231(const int /*prec*/ = 10)
{
return {


};
}

//00012232
//
//
template<typename P>
std::vector<std::vector<P>> poly00012232(const int /*prec*/ = 10)
{
return {


};
}

//00012233
//
//
template<typename P>
std::vector<std::vector<P>> poly00012233(const int /*prec*/ = 10)
{
return {


};
}

//00012330
//
//
template<typename P>
std::vector<std::vector<P>> poly00012330(const int /*prec*/ = 10)
{
return {


};
}

//00012331
//
//
template<typename P>
std::vector<std::vector<P>> poly00012331(const int /*prec*/ = 10)
{
return {


};
}

//00012332
//
//
template<typename P>
std::vector<std::vector<P>> poly00012332(const int /*prec*/ = 10)
{
return {


};
}

//00012333
//
//
template<typename P>
std::vector<std::vector<P>> poly00012333(const int /*prec*/ = 10)
{
return {


};
}

//00111123
//
//
template<typename P>
std::vector<std::vector<P>> poly00111123(const int /*prec*/ = 10)
{
return {


};
}

//00111203
//
//
template<typename P>
std::vector<std::vector<P>> poly00111203(const int /*prec*/ = 10)
{
return {


};
}

//00111223
//
//
template<typename P>
std::vector<std::vector<P>> poly00111223(const int /*prec*/ = 10)
{
return {


};
}

//00111230
//
//
template<typename P>
std::vector<std::vector<P>> poly00111230(const int /*prec*/ = 10)
{
return {


};
}

//00111232
//
//
template<typename P>
std::vector<std::vector<P>> poly00111232(const int /*prec*/ = 10)
{
return {


};
}

//00111233
//
//
template<typename P>
std::vector<std::vector<P>> poly00111233(const int /*prec*/ = 10)
{
return {


};
}

//00112233
//
//
template<typename P>
std::vector<std::vector<P>> poly00112233(const int /*prec*/ = 10)
{
return {


};
}

//00112323
//
//
template<typename P>
std::vector<std::vector<P>> poly00112323(const int /*prec*/ = 10)
{
return {


};
}

//00112332
//
//
template<typename P>
std::vector<std::vector<P>> poly00112332(const int /*prec*/ = 10)
{
return {


};
}

//00121203
//
//
template<typename P>
std::vector<std::vector<P>> poly00121203(const int /*prec*/ = 10)
{
return {


};
}

//00121223
//
//
template<typename P>
std::vector<std::vector<P>> poly00121223(const int /*prec*/ = 10)
{
return {


};
}

//00121233
//
//
template<typename P>
std::vector<std::vector<P>> poly00121233(const int /*prec*/ = 10)
{
return {


};
}

//00121300
//
//
template<typename P>
std::vector<std::vector<P>> poly00121300(const int /*prec*/ = 10)
{
return {


};
}

//00121301
//
//
template<typename P>
std::vector<std::vector<P>> poly00121301(const int /*prec*/ = 10)
{
return {


};
}

//00121302
//
//
template<typename P>
std::vector<std::vector<P>> poly00121302(const int /*prec*/ = 10)
{
return {


};
}

//00121320
//
//
template<typename P>
std::vector<std::vector<P>> poly00121320(const int /*prec*/ = 10)
{
return {


};
}

//00121321
//
//
template<typename P>
std::vector<std::vector<P>> poly00121321(const int /*prec*/ = 10)
{
return {


};
}

//00121323
//
//
template<typename P>
std::vector<std::vector<P>> poly00121323(const int /*prec*/ = 10)
{
return {


};
}

//00122103
//
//
template<typename P>
std::vector<std::vector<P>> poly00122103(const int /*prec*/ = 10)
{
return {


};
}

//00122113
//
//
template<typename P>
std::vector<std::vector<P>> poly00122113(const int /*prec*/ = 10)
{
return {


};
}

//00122133
//
//
template<typename P>
std::vector<std::vector<P>> poly00122133(const int /*prec*/ = 10)
{
return {


};
}

//00122300
//
//
template<typename P>
std::vector<std::vector<P>> poly00122300(const int /*prec*/ = 10)
{
return {


};
}

//00122301
//
//
template<typename P>
std::vector<std::vector<P>> poly00122301(const int /*prec*/ = 10)
{
return {


};
}

//00122302
//
//
template<typename P>
std::vector<std::vector<P>> poly00122302(const int /*prec*/ = 10)
{
return {


};
}

//00122313
//
//
template<typename P>
std::vector<std::vector<P>> poly00122313(const int /*prec*/ = 10)
{
return {


};
}

//00122331
//
//
template<typename P>
std::vector<std::vector<P>> poly00122331(const int /*prec*/ = 10)
{
return {


};
}

//01101023
//
//
template<typename P>
std::vector<std::vector<P>> poly01101023(const int /*prec*/ = 10)
{
return {


};
}

//01101223
//
//
template<typename P>
std::vector<std::vector<P>> poly01101223(const int /*prec*/ = 10)
{
return {


};
}

//01101231
//
//
template<typename P>
std::vector<std::vector<P>> poly01101231(const int /*prec*/ = 10)
{
return {


};
}

//01102332
//
//
template<typename P>
std::vector<std::vector<P>> poly01102332(const int /*prec*/ = 10)
{
return {


};
}

//01121223
//
//
template<typename P>
std::vector<std::vector<P>> poly01121223(const int /*prec*/ = 10)
{
return {


};
}

//01121230
//
//
template<typename P>
std::vector<std::vector<P>> poly01121230(const int /*prec*/ = 10)
{
return {


};
}

//01122330
//
//
template<typename P>
std::vector<std::vector<P>> poly01122330(const int /*prec*/ = 10)
{
return {


};
}

//01123023
//
//
template<typename P>
std::vector<std::vector<P>> poly01123023(const int /*prec*/ = 10)
{
return {


};
}

//01233210
//
//
template<typename P>
std::vector<std::vector<P>> poly01233210(const int /*prec*/ = 10)
{
return {


};
}


//00001234
//
//
template<typename P>
std::vector<std::vector<P>> poly00001234(const int /*prec*/ = 10)
{
return {


};
}

//00010234
//
//
template<typename P>
std::vector<std::vector<P>> poly00010234(const int /*prec*/ = 10)
{
return {


};
}

//00011234
//
//
template<typename P>
std::vector<std::vector<P>> poly00011234(const int /*prec*/ = 10)
{
return {


};
}

//00012034
//
//
template<typename P>
std::vector<std::vector<P>> poly00012034(const int /*prec*/ = 10)
{
return {


};
}

//00012134
//
//
template<typename P>
std::vector<std::vector<P>> poly00012134(const int /*prec*/ = 10)
{
return {


};
}

//00012234
//
//
template<typename P>
std::vector<std::vector<P>> poly00012234(const int /*prec*/ = 10)
{
return {


};
}

//00012334
//
//
template<typename P>
std::vector<std::vector<P>> poly00012334(const int /*prec*/ = 10)
{
return {


};
}

//00012340
//
//
template<typename P>
std::vector<std::vector<P>> poly00012340(const int /*prec*/ = 10)
{
return {


};
}

//00012341
//
//
template<typename P>
std::vector<std::vector<P>> poly00012341(const int /*prec*/ = 10)
{
return {


};
}

//00012342
//
//
template<typename P>
std::vector<std::vector<P>> poly00012342(const int /*prec*/ = 10)
{
return {


};
}

//00012343
//
//
template<typename P>
std::vector<std::vector<P>> poly00012343(const int /*prec*/ = 10)
{
return {


};
}

//00111234
//
//
template<typename P>
std::vector<std::vector<P>> poly00111234(const int /*prec*/ = 10)
{
return {


};
}

//00112234
//
//
template<typename P>
std::vector<std::vector<P>> poly00112234(const int /*prec*/ = 10)
{
return {


};
}

//00112324
//
//
template<typename P>
std::vector<std::vector<P>> poly00112324(const int /*prec*/ = 10)
{
return {


};
}

//00112334
//
//
template<typename P>
std::vector<std::vector<P>> poly00112334(const int /*prec*/ = 10)
{
return {


};
}

//00121234
//
//
template<typename P>
std::vector<std::vector<P>> poly00121234(const int /*prec*/ = 10)
{
return {


};
}

//00121304
//
//
template<typename P>
std::vector<std::vector<P>> poly00121304(const int /*prec*/ = 10)
{
return {


};
}

//00121324
//
//
template<typename P>
std::vector<std::vector<P>> poly00121324(const int /*prec*/ = 10)
{
return {


};
}

//00121340
//
//
template<typename P>
std::vector<std::vector<P>> poly00121340(const int /*prec*/ = 10)
{
return {


};
}

//00121341
//
//
template<typename P>
std::vector<std::vector<P>> poly00121341(const int /*prec*/ = 10)
{
return {


};
}

//00121342
//
//
template<typename P>
std::vector<std::vector<P>> poly00121342(const int /*prec*/ = 10)
{
return {


};
}

//00121344
//
//
template<typename P>
std::vector<std::vector<P>> poly00121344(const int /*prec*/ = 10)
{
return {


};
}

//00122134
//
//
template<typename P>
std::vector<std::vector<P>> poly00122134(const int /*prec*/ = 10)
{
return {


};
}

//00122304
//
//
template<typename P>
std::vector<std::vector<P>> poly00122304(const int /*prec*/ = 10)
{
return {


};
}

//00122314
//
//
template<typename P>
std::vector<std::vector<P>> poly00122314(const int /*prec*/ = 10)
{
return {


};
}

//00122324
//
//
template<typename P>
std::vector<std::vector<P>> poly00122324(const int /*prec*/ = 10)
{
return {


};
}

//00122334
//
//
template<typename P>
std::vector<std::vector<P>> poly00122334(const int /*prec*/ = 10)
{
return {


};
}

//00122344
//
//
template<typename P>
std::vector<std::vector<P>> poly00122344(const int /*prec*/ = 10)
{
return {


};
}

//00123400
//
//
template<typename P>
std::vector<std::vector<P>> poly00123400(const int /*prec*/ = 10)
{
return {


};
}

//00123401
//
//
template<typename P>
std::vector<std::vector<P>> poly00123401(const int /*prec*/ = 10)
{
return {


};
}

//00123414
//
//
template<typename P>
std::vector<std::vector<P>> poly00123414(const int /*prec*/ = 10)
{
return {


};
}

//00123421
//
//
template<typename P>
std::vector<std::vector<P>> poly00123421(const int /*prec*/ = 10)
{
return {


};
}

//00123423
//
//
template<typename P>
std::vector<std::vector<P>> poly00123423(const int /*prec*/ = 10)
{
return {


};
}

//01101234
//
//
template<typename P>
std::vector<std::vector<P>> poly01101234(const int /*prec*/ = 10)
{
return {


};
}

//01102334
//
//
template<typename P>
std::vector<std::vector<P>> poly01102334(const int /*prec*/ = 10)
{
return {


};
}

//01121234
//
//
template<typename P>
std::vector<std::vector<P>> poly01121234(const int /*prec*/ = 10)
{
return {


};
}

//01121340
//
//
template<typename P>
std::vector<std::vector<P>> poly01121340(const int /*prec*/ = 10)
{
return {


};
}

//01121341
//
//
template<typename P>
std::vector<std::vector<P>> poly01121341(const int /*prec*/ = 10)
{
return {


};
}

//01122034
//
//
template<typename P>
std::vector<std::vector<P>> poly01122034(const int /*prec*/ = 10)
{
return {


};
}

//01122334
//
//
template<typename P>
std::vector<std::vector<P>> poly01122334(const int /*prec*/ = 10)
{
return {


};
}

//01122340
//
//
template<typename P>
std::vector<std::vector<P>> poly01122340(const int /*prec*/ = 10)
{
return {


};
}

//01123024
//
//
template<typename P>
std::vector<std::vector<P>> poly01123024(const int /*prec*/ = 10)
{
return {


};
}

//01233214
//
//
template<typename P>
std::vector<std::vector<P>> poly01233214(const int /*prec*/ = 10)
{
return {


};
}

//00012345
//
//
template<typename P>
std::vector<std::vector<P>> poly00012345(const int /*prec*/ = 10)
{
return {


};
}

//00112345
//
//
template<typename P>
std::vector<std::vector<P>> poly00112345(const int /*prec*/ = 10)
{
return {


};
}

//00121345
//
//
template<typename P>
std::vector<std::vector<P>> poly00121345(const int /*prec*/ = 10)
{
return {


};
}

//00122345
//
//
template<typename P>
std::vector<std::vector<P>> poly00122345(const int /*prec*/ = 10)
{
return {


};
}

//00123405
//
//
template<typename P>
std::vector<std::vector<P>> poly00123405(const int /*prec*/ = 10)
{
return {


};
}

//00123415
//
//
template<typename P>
std::vector<std::vector<P>> poly00123415(const int /*prec*/ = 10)
{
return {


};
}

//00123425
//
//
template<typename P>
std::vector<std::vector<P>> poly00123425(const int /*prec*/ = 10)
{
return {


};
}

//00123455
//
//
template<typename P>
std::vector<std::vector<P>> poly00123455(const int /*prec*/ = 10)
{
return {


};
}

//01102345
//
//
template<typename P>
std::vector<std::vector<P>> poly01102345(const int /*prec*/ = 10)
{
return {


};
}

//01121345
//
//
template<typename P>
std::vector<std::vector<P>> poly01121345(const int /*prec*/ = 10)
{
return {


};
}

//01122345
//
//
template<typename P>
std::vector<std::vector<P>> poly01122345(const int /*prec*/ = 10)
{
return {


};
}

//01123045
//
//
template<typename P>
std::vector<std::vector<P>> poly01123045(const int /*prec*/ = 10)
{
return {


};
}

//01123445
//
//
template<typename P>
std::vector<std::vector<P>> poly01123445(const int /*prec*/ = 10)
{
return {


};
}

//01123453
//
//
template<typename P>
std::vector<std::vector<P>> poly01123453(const int /*prec*/ = 10)
{
return {


};
}

//01233245
//
//
template<typename P>
std::vector<std::vector<P>> poly01233245(const int /*prec*/ = 10)
{
return {


};
}

//00123456
//
//
template<typename P>
std::vector<std::vector<P>> poly00123456(const int /*prec*/ = 10)
{
return {


};
}

//01123456
//
//
template<typename P>
std::vector<std::vector<P>> poly01123456(const int /*prec*/ = 10)
{
return {


};
}

//01233456
//
//
template<typename P>
std::vector<std::vector<P>> poly01233456(const int /*prec*/ = 10)
{
return {


};
}

//01234567
//
//
template<typename P>
std::vector<std::vector<P>> poly01234567(const int /*prec*/ = 10)
{
return {


};
}

// Cube (begin definition)
using Cube = std::array<std::uint8_t, 8>;

inline constexpr Cube convert_to_cube(unsigned int n) {
  assert(n < (1 << 24));
  return {
    (std::uint8_t)((n & 070000000) >> 21), (std::uint8_t)((n & 007000000) >> 18),
    (std::uint8_t)((n & 000700000) >> 15), (std::uint8_t)((n & 000070000) >> 12),
    (std::uint8_t)((n & 000007000) >> 9), (std::uint8_t)((n & 000000700) >> 6),
    (std::uint8_t)((n & 000000070) >> 3), (std::uint8_t)((n & 000000007) >> 0),
  };
}

// User-defined literal operator.  Given an integer in octal notation, like
// 01234567, gives the cube with the same colors. For example, `01234567_c`
// is `Cube{0, 1, 2, 3, 4, 5, 6, 7}`.
inline constexpr Cube operator"" _c(unsigned long long n)
{
  assert(n < (1 << 24));
  return convert_to_cube(unsigned(n));
}

inline std::string config_name(const Cube cube) {
  std::stringstream filename_prefix_ss;
  for (int j = 0; j < 8; ++j) {
    filename_prefix_ss << int(cube[j]);
  }
  return filename_prefix_ss.str();
}
// Cube (end)


template<typename Point>
class Triple_line_extractor
{
    using P = Point;
    using Polyline = std::vector<P>;
    using Polylines = std::vector<Polyline>;

    typedef Polylines(*create_polylines_fct)(const int /* prec */);


public:
    boost::unordered_map<Cube, create_polylines_fct> create_polylines_fcts {
          // One internal corner
            { 00001221_c, CGAL::Mesh_3::poly00001221 } ,
            { 00111202_c, CGAL::Mesh_3::poly00111202 } ,
            { 01101001_c, CGAL::Mesh_3::poly01101001 } ,
          // Two curves
            { 00011022_c, CGAL::Mesh_3::poly00011022 } ,
            { 00011221_c, CGAL::Mesh_3::poly00011221 } ,
            { 00011222_c, CGAL::Mesh_3::poly00011222 } ,
            { 00121200_c, CGAL::Mesh_3::poly00121200 } ,
            { 00121221_c, CGAL::Mesh_3::poly00121221 } ,
            { 00122100_c, CGAL::Mesh_3::poly00122100 } ,
            { 00122101_c, CGAL::Mesh_3::poly00122101 } ,
          // One curve
            { 00000012_c, CGAL::Mesh_3::poly00000012 } ,
            { 00000112_c, CGAL::Mesh_3::poly00000112 } ,
            { 00000121_c, CGAL::Mesh_3::poly00000121 } ,
            { 00001112_c, CGAL::Mesh_3::poly00001112 } ,
            { 00001122_c, CGAL::Mesh_3::poly00001122 } ,
            { 00010121_c, CGAL::Mesh_3::poly00010121 } ,
            { 00010122_c, CGAL::Mesh_3::poly00010122 } ,
            { 00011002_c, CGAL::Mesh_3::poly00011002 } ,
            { 00011012_c, CGAL::Mesh_3::poly00011012 } ,
            { 00011110_c, CGAL::Mesh_3::poly00011110 } ,
            { 00011120_c, CGAL::Mesh_3::poly00011120 } ,
            { 00011121_c, CGAL::Mesh_3::poly00011121 } ,
            { 00011122_c, CGAL::Mesh_3::poly00011122 } ,
            { 00011220_c, CGAL::Mesh_3::poly00011220 } ,
            { 00012002_c, CGAL::Mesh_3::poly00012002 } ,
            { 00012012_c, CGAL::Mesh_3::poly00012012 } ,
            { 00012021_c, CGAL::Mesh_3::poly00012021 } ,
            { 00012110_c, CGAL::Mesh_3::poly00012110 } ,
            { 00012112_c, CGAL::Mesh_3::poly00012112 } ,
            { 00012120_c, CGAL::Mesh_3::poly00012120 } ,
            { 00012121_c, CGAL::Mesh_3::poly00012121 } ,
            { 00012122_c, CGAL::Mesh_3::poly00012122 } ,
            { 00012221_c, CGAL::Mesh_3::poly00012221 } ,
            { 00111100_c, CGAL::Mesh_3::poly00111100 } ,
            { 00111102_c, CGAL::Mesh_3::poly00111102 } ,
            { 00111220_c, CGAL::Mesh_3::poly00111220 } ,
            { 00121201_c, CGAL::Mesh_3::poly00121201 } ,
          // 4 color cases
            { 00000123_c, CGAL::Mesh_3::poly00000123 } ,
            { 00001123_c, CGAL::Mesh_3::poly00001123 } ,
            { 00001223_c, CGAL::Mesh_3::poly00001223 } ,
            { 00010123_c, CGAL::Mesh_3::poly00010123 } ,
            { 00010230_c, CGAL::Mesh_3::poly00010230 } ,
            { 00010231_c, CGAL::Mesh_3::poly00010231 } ,
            { 00011023_c, CGAL::Mesh_3::poly00011023 } ,
            { 00011123_c, CGAL::Mesh_3::poly00011123 } ,
            { 00011223_c, CGAL::Mesh_3::poly00011223 } ,
            { 00011230_c, CGAL::Mesh_3::poly00011230 } ,
            { 00011231_c, CGAL::Mesh_3::poly00011231 } ,
            { 00011232_c, CGAL::Mesh_3::poly00011232 } ,
            { 00012003_c, CGAL::Mesh_3::poly00012003 } ,
            { 00012013_c, CGAL::Mesh_3::poly00012013 } ,
            { 00012023_c, CGAL::Mesh_3::poly00012023 } ,
            { 00012033_c, CGAL::Mesh_3::poly00012033 } ,
            { 00012113_c, CGAL::Mesh_3::poly00012113 } ,
            { 00012123_c, CGAL::Mesh_3::poly00012123 } ,
            { 00012130_c, CGAL::Mesh_3::poly00012130 } ,
            { 00012131_c, CGAL::Mesh_3::poly00012131 } ,
            { 00012132_c, CGAL::Mesh_3::poly00012132 } ,
            { 00012133_c, CGAL::Mesh_3::poly00012133 } ,
            { 00012223_c, CGAL::Mesh_3::poly00012223 } ,
            { 00012230_c, CGAL::Mesh_3::poly00012230 } ,
            { 00012231_c, CGAL::Mesh_3::poly00012231 } ,
            { 00012232_c, CGAL::Mesh_3::poly00012232 } ,
            { 00012233_c, CGAL::Mesh_3::poly00012233 } ,
            { 00012330_c, CGAL::Mesh_3::poly00012330 } ,
            { 00012331_c, CGAL::Mesh_3::poly00012331 } ,
            { 00012332_c, CGAL::Mesh_3::poly00012332 } ,
            { 00012333_c, CGAL::Mesh_3::poly00012333 } ,
            { 00111123_c, CGAL::Mesh_3::poly00111123 } ,
            { 00111203_c, CGAL::Mesh_3::poly00111203 } ,
            { 00111223_c, CGAL::Mesh_3::poly00111223 } ,
            { 00111230_c, CGAL::Mesh_3::poly00111230 } ,
            { 00111232_c, CGAL::Mesh_3::poly00111232 } ,
            { 00111233_c, CGAL::Mesh_3::poly00111233 } ,
            { 00112233_c, CGAL::Mesh_3::poly00112233 } ,
            { 00112323_c, CGAL::Mesh_3::poly00112323 } ,
            { 00112332_c, CGAL::Mesh_3::poly00112332 } ,
            { 00121203_c, CGAL::Mesh_3::poly00121203 } ,
            { 00121223_c, CGAL::Mesh_3::poly00121223 } ,
            { 00121233_c, CGAL::Mesh_3::poly00121233 } ,
            { 00121300_c, CGAL::Mesh_3::poly00121300 } ,
            { 00121301_c, CGAL::Mesh_3::poly00121301 } ,
            { 00121302_c, CGAL::Mesh_3::poly00121302 } ,
            { 00121320_c, CGAL::Mesh_3::poly00121320 } ,
            { 00121321_c, CGAL::Mesh_3::poly00121321 } ,
            { 00121323_c, CGAL::Mesh_3::poly00121323 } ,
            { 00122103_c, CGAL::Mesh_3::poly00122103 } ,
            { 00122113_c, CGAL::Mesh_3::poly00122113 } ,
            { 00122133_c, CGAL::Mesh_3::poly00122133 } ,
            { 00122300_c, CGAL::Mesh_3::poly00122300 } ,
            { 00122301_c, CGAL::Mesh_3::poly00122301 } ,
            { 00122302_c, CGAL::Mesh_3::poly00122302 } ,
            { 00122313_c, CGAL::Mesh_3::poly00122313 } ,
            { 00122331_c, CGAL::Mesh_3::poly00122331 } ,
            { 01101023_c, CGAL::Mesh_3::poly01101023 } ,
            { 01101223_c, CGAL::Mesh_3::poly01101223 } ,
            { 01101231_c, CGAL::Mesh_3::poly01101231 } ,
            { 01102332_c, CGAL::Mesh_3::poly01102332 } ,
            { 01121223_c, CGAL::Mesh_3::poly01121223 } ,
            { 01121230_c, CGAL::Mesh_3::poly01121230 } ,
            { 01122330_c, CGAL::Mesh_3::poly01122330 } ,
            { 01123023_c, CGAL::Mesh_3::poly01123023 } ,
            { 01233210_c, CGAL::Mesh_3::poly01233210 } ,
          // 5 colors
            { 00001234_c, CGAL::Mesh_3::poly00001234 } ,
            { 00010234_c, CGAL::Mesh_3::poly00010234 } ,
            { 00011234_c, CGAL::Mesh_3::poly00011234 } ,
            { 00012034_c, CGAL::Mesh_3::poly00012034 } ,
            { 00012134_c, CGAL::Mesh_3::poly00012134 } ,
            { 00012234_c, CGAL::Mesh_3::poly00012234 } ,
            { 00012334_c, CGAL::Mesh_3::poly00012334 } ,
            { 00012340_c, CGAL::Mesh_3::poly00012340 } ,
            { 00012341_c, CGAL::Mesh_3::poly00012341 } ,
            { 00012342_c, CGAL::Mesh_3::poly00012342 } ,
            { 00012343_c, CGAL::Mesh_3::poly00012343 } ,
            { 00111234_c, CGAL::Mesh_3::poly00111234 } ,
            { 00112234_c, CGAL::Mesh_3::poly00112234 } ,
            { 00112324_c, CGAL::Mesh_3::poly00112324 } ,
            { 00112334_c, CGAL::Mesh_3::poly00112334 } ,
            { 00121234_c, CGAL::Mesh_3::poly00121234 } ,
            { 00121304_c, CGAL::Mesh_3::poly00121304 } ,
            { 00121324_c, CGAL::Mesh_3::poly00121324 } ,
            { 00121340_c, CGAL::Mesh_3::poly00121340 } ,
            { 00121341_c, CGAL::Mesh_3::poly00121341 } ,
            { 00121342_c, CGAL::Mesh_3::poly00121342 } ,
            { 00121344_c, CGAL::Mesh_3::poly00121344 } ,
            { 00122134_c, CGAL::Mesh_3::poly00122134 } ,
            { 00122304_c, CGAL::Mesh_3::poly00122304 } ,
            { 00122314_c, CGAL::Mesh_3::poly00122314 } ,
            { 00122324_c, CGAL::Mesh_3::poly00122324 } ,
            { 00122334_c, CGAL::Mesh_3::poly00122334 } ,
            { 00122344_c, CGAL::Mesh_3::poly00122344 } ,
            { 00123400_c, CGAL::Mesh_3::poly00123400 } ,
            { 00123401_c, CGAL::Mesh_3::poly00123401 } ,
            { 00123414_c, CGAL::Mesh_3::poly00123414 } ,
            { 00123421_c, CGAL::Mesh_3::poly00123421 } ,
            { 00123423_c, CGAL::Mesh_3::poly00123423 } ,
            { 01101234_c, CGAL::Mesh_3::poly01101234 } ,
            { 01102334_c, CGAL::Mesh_3::poly01102334 } ,
            { 01121234_c, CGAL::Mesh_3::poly01121234 } ,
            { 01121340_c, CGAL::Mesh_3::poly01121340 } ,
            { 01121341_c, CGAL::Mesh_3::poly01121341 } ,
            { 01122034_c, CGAL::Mesh_3::poly01122034 } ,
            { 01122334_c, CGAL::Mesh_3::poly01122334 } ,
            { 01122340_c, CGAL::Mesh_3::poly01122340 } ,
            { 01123024_c, CGAL::Mesh_3::poly01123024 } ,
            { 01233214_c, CGAL::Mesh_3::poly01233214 } ,
          // 6 colors
            { 00012345_c, CGAL::Mesh_3::poly00012345 } ,
            { 00112345_c, CGAL::Mesh_3::poly00112345 } ,
            { 00121345_c, CGAL::Mesh_3::poly00121345 } ,
            { 00122345_c, CGAL::Mesh_3::poly00122345 } ,
            { 00123405_c, CGAL::Mesh_3::poly00123405 } ,
            { 00123415_c, CGAL::Mesh_3::poly00123415 } ,
            { 00123425_c, CGAL::Mesh_3::poly00123425 } ,
            { 00123455_c, CGAL::Mesh_3::poly00123455 } ,
            { 01102345_c, CGAL::Mesh_3::poly01102345 } ,
            { 01121345_c, CGAL::Mesh_3::poly01121345 } ,
            { 01122345_c, CGAL::Mesh_3::poly01122345 } ,
            { 01123045_c, CGAL::Mesh_3::poly01123045 } ,
            { 01123445_c, CGAL::Mesh_3::poly01123445 } ,
            { 01123453_c, CGAL::Mesh_3::poly01123453 } ,
            { 01233245_c, CGAL::Mesh_3::poly01233245 } ,
          // 7 colors
            { 00123456_c, CGAL::Mesh_3::poly00123456 } ,
            { 01123456_c, CGAL::Mesh_3::poly01123456 } ,
            { 01233456_c, CGAL::Mesh_3::poly01233456 } ,
          // 8 colors
            { 01234567_c, CGAL::Mesh_3::poly01234567 }
    };

};//class Triple_line_extractor
}//namespace Mesh_3
}//namespace CGAL


#endif // CGAL_MESH_3_TRIPLE_LINES_H
