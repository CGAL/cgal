#ifndef TRIPLE_LIGNES_H
#define TRIPLE_LIGNES_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef Point_3 P;
typedef std::vector<Point_3> Polyline;
typedef std::vector<Polyline>   Polylines;

typedef Polylines (*create_polylines_fct)(const int /* prec */);


#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>

#include <CGAL/Mesh_3/triple_lines_extraction/cube.h>

//extern boost::unordered_map<Cube, create_polylines_fct> create_polylines_fcts;

template <typename Functor>
Polyline create_polyline(const double start,
    const double end,
    P starting_point,
    P ending_point,
    Functor f,
    const int prec)
{
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

template <typename Functor>
Polyline create_polyline(const double start,
    const double end,
    P starting_point,
    Functor f,
    const int prec)
{
    return create_polyline(start, end, starting_point, f(end), f, prec);
}

template <typename Functor>
Polyline create_polyline(const double start,
    const double end,
    Functor f,
    const int prec)
{
    return create_polyline(start, end, f(start), f(end), f, prec);
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
Polylines poly00001221(const int /* no sampling for segments */)
{
    P corner{ 1. / 2, 1. / 2, 2. / 3 };
    P      a{ 1. / 2,    0, 2. / 3 };
    P      b{ 1. / 2,    1, 2. / 3 };
    P      c{ 0,    1. / 2, 2. / 3 };
    P      d{ 1,    1. / 2, 2. / 3 };
    P      e{ 1,    1. / 2, 1 };

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
Polylines poly00111202(const int prec = 10)
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
Polylines poly01101001(const int /* no sampling for segments */)
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
Polylines poly00011022(const int prec = 10)
{
    // x = (3*z^2-2*z)/(3*z^2-1), y = 1/(3*z), z = [1/3,1/2]
    // x = (3*z^2-2*z)/(3*z^2-1), y = 1/(3*z), z = [2/3,1]
    auto f = [](double z) {return P(z * (3 * z - 2) / (3 * z * z - 1),
        1 / (3 * z),
        z); };
    return {
      create_polyline(1. / 3, 1. / 2, f, prec),
      create_polyline(2. / 3, 1.  , f, prec),
    };
}

// 00011221
//   curve_1 : x = ]1/2,1], y = (3 * x * x - sqrt(9 * x * x * x * x - 30 * x * x * x + 45 * x * x - 24 * x + 4) + 3 * x - 2)/(6 * x * (2 * x - 1)), z = 1./(3*x+3*y-6*x*y)
//   point limit x = 1/2, y = 0, z = 2/3
//   curve_2 : x = ]0,1/3], y = (3 * x * x + sqrt(9 * x * x * x * x - 30 * x * x * x + 45 * x * x - 24 * x + 4) + 3 * x - 2)/(6 * x * (2 * x - 1)), z = 1./(3*x+3*y-6*x*y)
//   point limit x = 0, y = 1/2, z = 2/3
Polylines poly00011221(const int prec = 10)
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
Polylines poly00011222(const int prec = 10)
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
Polylines poly00121200(const int prec = 10)
{
    auto y = [](double z) { return (3 * z - 2) / (6 * z - 3); };
    return {
      create_polyline(0, 1. / 3,
                      [y](double z) { return P(1. / 2, y(z), z); },
                      prec),
      create_polyline(2. / 3, 1,
                      [y](double z) { return P(1. / 2, y(z), z); },
                      prec)
    };
}

// 00121221
//   curve_1 : x = 1/2, y = (3*z-2)/(3*z-3), z = [0,2/3]
//   curve_2 : x = 1/2, y = z/(3*z-1), z = [1/2,1]
Polylines poly00121221(const int prec = 10)
{
    auto y1 = [](double z) { return (3 * z - 2) / (3 * z - 3); };
    auto y2 = [](double z) { return z / (3 * z - 1); };
    return {
      create_polyline(0, 2. / 3,
                      [y1](double z) { return P(1. / 2, y1(z), z); },
                      prec),
      create_polyline(1. / 2, 1,
                      [y2](double z) { return P(1. / 2, y2(z), z); },
                      prec)
    };
}

// 00122100
//   curve_1  : x = 1/2, y = (3*z-2)/(6*z-3), z = [0,1/3]
//   curve_1' : x = 1/2, y = (3*z-2)/(6*z-3), z = [2/3,1]
Polylines poly00122100(const int prec = 10)
{
    auto y = [](double z) { return (3 * z - 2) / (6 * z - 3); };
    return {
      create_polyline(0, 1. / 3,
                      [y](double z) { return P(1. / 2, y(z), z); },
                      prec),
      create_polyline(2. / 3, 1,
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
Polylines poly00122101(const int prec = 10)
{
    auto sq_exp = [](double x) {
        return sqrt(24 * x * x * x - 35 * x * x + 18 * x - 3);
    };
    auto y1 = [sq_exp](double x) { return (-sq_exp(x) - 5 * x + 3) / (6 * CGAL::square(x - 1)); };
    auto z1 = [sq_exp](double x) { return (-sq_exp(x) - 7 * x + 3) / (6 * (x * x - 3 * x + 1)); };
    auto y2 = [sq_exp](double x) { return (sq_exp(x) - 5 * x + 3) / (6 * CGAL::square(x - 1)); };
    auto z2 = [sq_exp](double x) { return (sq_exp(x) - 7 * x + 3) / (6 * (x * x - 3 * x + 1)); };
    P corner{ 1., .5, 1. / 3 };
    return {
      create_polyline(1. / 3, 1. / 2,
                      [y1, z1](double x) { return P(x, y1(x), z1(x)); },
                      prec),
      create_polyline(1., 1. / 2, corner,
                      [y2, z2](double x) { return P(x, y2(x), z2(x)); },
                      prec),
    };
}

//
// One curve
//
Polylines poly00000012(const int prec = 10)
{
    // curve : x = 1/2, y = 2/(3*z), z = [2/3,1]
    return { create_polyline(2. / 3, 1.,
                             [](double z) { return P(0.5, 2. / (3. * z), z); },
                             prec) };
}

// 00000112
//   x =[1/2,1], y = x/(3 * x - 1), z = (3 * x - 1)/(3 * x * x)
Polylines poly00000112(const int prec = 10)
{
    return { create_polyline(1. / 2, 1,
                             [](double x) { return P(x, x / (3 * x - 1), (3 * x - 1) / (3 * x * x)); },
                             prec) };
}

// 00000121
//   curve : x = 1/(3*z), y = 1/(3*z-1), z = [2/3,1]
Polylines poly00000121(const int prec = 10)
{
    return { create_polyline(2. / 3, 1,
                             [](double z) { return P(1 / (3 * z), 1 / (3 * z - 1), z); },
                             prec) };
}

// 00001112
//   curve : x = 1/(2*y), y = [1/2, 1],z = 2/3
Polylines poly00001112(const int prec = 10)
{
    return { create_polyline(1. / 2, 1,
                             [](double y) { return P(1 / (2 * y), y, 2. / 3); },
                             prec) };
}

// 00001122
//   curve : x = [0,1], y = 1/2, z = 2/3
Polylines poly00001122(const int prec = 10)
{
    return { create_polyline(0, 1,
                             [](double x) { return P(x, 1. / 2, 2. / 3); },
                             prec) };
}

// 00010121
//   curve : x =y * z / (z+y), y = ((3 * z * z - 1) - sqrt(CGAL::square(1 - 3 * z * z) - 12 * (z - 1) * z * z))/(6 * (z - 1) * z), z=[1,1/2[
//   point limit =  (1/3, 1/2, 1)
Polylines poly00010121(const int prec = 10)
{
    auto y = [](double z) { return ((3 * z * z - 1) - sqrt(CGAL::square(1 - 3 * z * z) - 12 * (z - 1) * z * z)) / (6 * (z - 1) * z); };
    auto x = [](double y, double z) { return y * z / (z + y); };
    P corner(1. / 3, 1. / 2, 1);
    return { create_polyline(1, 1. / 2, corner,
                             [x, y](double z) { return P(x(y(z), z), y(z), z); },
                             prec) };
}

// 00010122
//   curve : x = z/(3*z^2-2*z+1), y = 1/(3*z), z = [1/3,1]
Polylines poly00010122(const int prec = 10)
{
    return { create_polyline(1. / 3, 1,
                             [](double z) { return P(z / (3 * z * z - 2 * z + 1), 1. / (3 * z), z); },
                             prec) };
}
// 00011002
//   curve : x = (y+1)/(4*y-1), y = [2/3,1], z = 1/2
Polylines poly00011002(const int prec = 10)
{
    return { create_polyline(2. / 3, 1,
                             [](double y) { return P((y + 1) / (4 * y - 1), y, 1. / 2); },
                             prec) };
}
// 00011012
//   curve : x = (3*z^2-2*z+1)/(3*z^2), y = z/(3*z^2-2*z+1), z = [1/2, 1]
Polylines poly00011012(const int prec = 10)
{
    return { create_polyline(1. / 2, 1,
                             [](double z) { return P((3 * z * z - 2 * z + 1) / (3 * z * z), z / (3 * z * z - 2 * z + 1), z); },
                             prec) };
}
// 00011110
//   curve : x = 1/(2*y), y = [1/2,1], z = 1/2
Polylines poly00011110(const int prec = 10)
{
    return { create_polyline(1. / 2, 1,
                             [](double y) { return P(1. / (2 * y), y, 1. / 2); },
                             prec) };
}
// 00011120
//   curve : x = (3*z^2-2*z)/(3*z^2-1), y = (3*z^2-1)/(6*z^2-3*z), z = [2/3,1]
Polylines poly00011120(const int prec = 10)
{
    return { create_polyline(2. / 3, 1,
                             [](double z) { return P((3 * z * z - 2 * z) / (3 * z * z - 1), (3 * z * z - 1) / (6 * z * z - 3 * z), z); },
                             prec) };
}
// 00011121
//   curve : x = (3*z^2-2*z)/(3*z^2-z-1), y = (3*z^2-z-1)/(3*z^2-3*z), z =[1/2,2./3]
Polylines poly00011121(const int prec = 10)
{
    return { create_polyline(1. / 2, 2. / 3,
                             [](double z) { return P((3 * z * z - 2 * z) / (3 * z * z - z - 1), (3 * z * z - z - 1) / (3 * z * z - 3 * z), z); },
                             prec) };
}
// 00011122
//   curve : x = (3*z*z-2*z)/(z-1),y = 1/(3*z),z = [1/3, 2/3]
//
Polylines poly00011122(const int prec = 10)
{
    return { create_polyline(1. / 3,2. / 3,
                             [](double z) { return P((3 * z * z - 2 * z) / (z - 1), 1 / (3 * z), z); },
                             prec) };
}

// 00011220
//   curve : x=[0,1/2], y = (2 * x - 1)/(3 * x - 2), z = (2 * (x^2 - 2 * x + 1))/(5 * x^2 - 7 * x + 3)
Polylines poly00011220(const int prec = 10)
{
    return { create_polyline(0, 1. / 2,
                             [](double x) { return P(x, (2 * x - 1) / (3 * x - 2), (2 * (x * x - 2 * x + 1)) / (5 * x * x - 7 * x + 3)); },
                             prec) };
}
// 00012002
//   curve_1 : x = [2/3,1], y = (3 * x*x + sqrt(9 * x*x*x*x - 24 * x*x*x + 30 * x*x - 12 * x + 1) - 1)/(6 * x * (2 * x - 1)), z = 1 - 1./(3*x*y)
Polylines poly00012002(const int prec = 10)
{
    auto y = [](double x) { return (3 * x * x + sqrt(9 * x * x * x * x - 24 * x * x * x + 30 * x * x - 12 * x + 1) - 1) / (6 * x * (2 * x - 1)); };
    return { create_polyline(2. / 3, 1,
                             [y](double x) { return P(x, y(x), 1 - 1. / (3 * x * y(x))); },
                             prec) };
}
// 00012012
//   curve_1 : x = [0, 1/2[ , y = 1./(-6*x*z+3*x+3*z), z = (3 *x*x + sqrt(9 *x*x*x*x - 18 *x*x*x + 25 *x*x - 16 * x + 4) - 7 * x + 2)/(6 *(x - 1) * (2 * x - 1))
//   curve_1 : x = ]1/2,1[, y = 1./(-6*x*z+3*x+3*z), z = (3 *x*x + sqrt(9 *x*x*x*x - 18 *x*x*x + 25 *x*x - 16 * x + 4) - 7 * x + 2)/(6 *(x - 1) * (2 * x - 1))
//  point limit x = 1, y = 2/3, z = 1/2
//  point limit x = 1/2, y= 2./3, z = 2./3);
//
Polylines poly00012012(const int prec = 10)
{
    auto y = [](double x, double z) { return 1. / (-6 * x * z + 3 * x + 3 * z); };
    auto z = [](double x) { return (3 * x * x + sqrt(9 * x * x * x * x - 18 * x * x * x + 25 * x * x - 16 * x + 4) - 7 * x + 2) / (6 * (x - 1) * (2 * x - 1)); };
    P corner1(1. / 2, 2. / 3, 2. / 3);
    P corner2(1, 2. / 3, 1. / 2);
    return { create_polyline(1. / 2, 0, corner1,
                             [y, z](double x) { return P(x, y(x, z(x)), z(x)); },
                             prec),
             create_polyline(1. / 2, 1, corner1, corner2,
                                      [y, z](double x) { return P(x, y(x, z(x)), z(x)); },
                                      prec)
    };
}
// 00012021
//   curve : x = (3*z-1)/(3*z), y = z/(3*z-1), z = [1/2,1]
Polylines poly00012021(const int prec = 10)
{
    return { create_polyline(1. / 2, 1,
                             [](double z) { return P((3 * z - 1) / (3 * z), z / (3 * z - 1), z); },
                             prec) };
}
// 00012110
//   curve : x=]0,1/2], y = (3 * x * x + sqrt(9 * x * x * x * x - 18 * x * x * x + 25 * x * x - 16 * x + 4) + x - 2)/(6 * (x - 1) * x), z = (3*x*y-1)/(9*x*y-3*x-3*y)
//   point limit : x = 0, y = 1/2, z = 2/3
//
Polylines poly00012110(const int prec = 10)
{
    auto y = [](double x) { return (3 * x * x + sqrt(9 * x * x * x * x - 18 * x * x * x + 25 * x * x - 16 * x + 4) + x - 2) / (6 * (x - 1) * x); };
    auto z = [](double x, double y) { return (3 * x * y - 1) / (9 * x * y - 3 * x - 3 * y); };
    P corner(0, 1. / 2, 2. / 3);
    return { create_polyline(0, 1. / 2, corner,
                             [y, z](double x) { return P(x, y(x), z(x, y(x))); },
                             prec) };
}
// 00012112
//   x = ]0,1/2[, y = (3 * x*x + sqrt(9 * x * x * x * x - 36 * x * x * x + 40 * x * x - 20 * x + 4) + 2 * x - 2)/(6 * x * (2 * x - 1)), z = (3*x*y-1)/(9*x*y-3*x-3*y)
//   point limit x = 0, y = 1/2, z = 2/3
//   point limit x = 1/2, y = 0, z = 2/3
//
Polylines poly00012112(const int prec = 10)
{
    auto y = [](double x) { return (3 * x * x + sqrt(9 * x * x * x * x - 36 * x * x * x + 40 * x * x - 20 * x + 4) + 2 * x - 2) / (6 * x * (2 * x - 1)); };
    auto z = [](double x, double y) { return (3 * x * y - 1) / (9 * x * y - 3 * x - 3 * y); };
    P corner1(0, 1. / 2, 2. / 3);
    P corner2(1. / 2, 0, 2. / 3);
    return { create_polyline(0, 1. / 2, corner1, corner2,
                             [y, z](double x) { return P(x, y(x), z(x, y(x))); },
                             prec) };
}

// 00012120
//   curve : x = (3*z-1)/(3*z), y = (3*z^2-2*z)/(6*z^2-5*z+1), z = [2/3,1]
Polylines poly00012120(const int prec = 10)
{
    return { create_polyline(2. / 3, 1,
                             [](double z) { return P((3 * z - 1) / (3 * z), (3 * z * z - 2 * z) / (6 * z * z - 5 * z + 1), z); },
                             prec) };
}
//
// 00012121
//  curve : x = (3*z-1)/(3*z), y = (3*z^2-2*z)/(3*z^2-4*z+1), z = [1/2,2/3]
Polylines poly00012121(const int prec = 10)
{
    return { create_polyline(1. / 2, 2. / 3,
                             [](double z) { return P((3 * z - 1) / (3 * z), (3 * z * z - 2 * z) / (3 * z * z - 4 * z + 1), z); },
                             prec) };
}
// 00012122
//  curve : x = (6*z^2-6*z+1)/(3*z^2-3*z), y = (3*z^2-2*z)/(6*z^2-6*z+1), z = [1/3,2/3]
Polylines poly00012122(const int prec = 10)
{
    auto x = [](double z) { return (6 * z * z - 6 * z + 1) / (3 * z * z - 3 * z); };
    auto y = [](double z) { return (3 * z * z - 2 * z) / (6 * z * z - 6 * z + 1); };
    return { create_polyline(1. / 3, 2. / 3,
                             [x,y](double z) { return P(x(z), y(z), z); },
                             prec) };
}
// 00012221
//   curve : x = 1/(3*y),y = [1/3,1],z = 1/2
Polylines poly00012221(const int prec = 10)
{
    return { create_polyline(1. / 3, 1,
                             [](double y) { return P(1. / (3 * y),y , 1. / 2); },
                             prec) };
}

// 00111100
//   curve : x = [0,1], y = 1/2, z = 1/2
Polylines poly00111100(const int prec = 10)
{
    return { create_polyline(0, 1,
                             [](double x) { return P(x , 1. / 2, 1. / 2); },
                             prec) };
}

// 00111102
//   curve : x = (2*z-1)/(3*z^2-z), y = (3*z-1)/(6*z-3), z = [2/3,1]
Polylines poly00111102(const int prec = 10)
{
    return { create_polyline(2. / 3, 1,
                             [](double z) { return P((2 * z - 1) / (3 * z * z - z), (3 * z - 1) / (6 * z - 3), z); },
                             prec) };
}

// 00111220
//   segment 1/2 0 2/3 1/2 1 2/3
Polylines poly00111220(const int /*not needed for a segment*/)
{
    return { { P(1. / 2, 0, 2. / 3), P(1. / 2, 1, 2. / 3) } };
}

// 00121201
//   curve_1 : x =[1/2, 2/3], y = (-sqrt(-24 * x^3 + 37 * x^2 - 20 * x + 4) + 5 * x - 2)/(6 * x^2), z = ( sqrt(-24 * x^3 + 37 * x^2 - 20 * x + 4) + 5 * x - 2)/(6 * x^2)
//   curve_2 : x =[1/2, 2/3], y = ( sqrt(-24 * x^3 + 37 * x^2 - 20 * x + 4) + 5 * x - 2)/(6 * x^2), z = (-sqrt(-24 * x^3 + 37 * x^2 - 20 * x + 4) + 5 * x - 2)/(6 * x^2)
//   point 0 1/2 1/2
Polylines poly00121201(const int prec = 10)
{
    auto sq_exp = [](double x) {
        return sqrt(-24 * x * x * x + 37 * x * x - 20 * x + 4);
    };
    auto y = [sq_exp](double x) { return (-sq_exp(x) + 5 * x - 2) / (6 * x * x); };
    auto z = [sq_exp](double x) { return (sq_exp(x) + 5 * x - 2) / (6 * x * x); };
    P corner{ 2. / 3, y(2. / 3), z(2. / 3) };
    return {
      create_polyline(2. / 3, 1. / 2, corner,
                      [y, z](double x) { return P(x, y(x), z(x)); },
                      prec),
      create_polyline(2. / 3, 1. / 2, corner,
                      [y, z](double x) { return P(x, z(x), y(x)); },
                      prec),
                      // { P(0., .5, .5), P(0., .5, .5) }
    };
}


#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>

boost::unordered_map<Cube, create_polylines_fct> create_polylines_fcts{
    // One internal corner
      { 00001221_c, poly00001221 } ,
      { 00111202_c, poly00111202 } ,
      { 01101001_c, poly01101001 } ,
      // Two curves
      { 00011022_c, poly00011022 } ,
      { 00011221_c, poly00011221 } ,
      { 00011222_c, poly00011222 } ,
      { 00121200_c, poly00121200 } ,
      { 00121221_c, poly00121221 } ,
      { 00122100_c, poly00122100 } ,
      { 00122101_c, poly00122101 } ,
      // One curve
      { 00000012_c, poly00000012 } ,
      { 00000112_c, poly00000112 } ,
      { 00000121_c, poly00000121 } ,
      { 00001112_c, poly00001112 } ,
      { 00001122_c, poly00001122 } ,
      { 00010121_c, poly00010121 } ,
      { 00010122_c, poly00010122 } ,
      { 00011002_c, poly00011002 } ,
      { 00011012_c, poly00011012 } ,
      { 00011110_c, poly00011110 } ,
      { 00011120_c, poly00011120 } ,
      { 00011121_c, poly00011121 } ,
      { 00011122_c, poly00011122 } ,
      { 00011220_c, poly00011220 } ,
      { 00012002_c, poly00012002 } ,
      { 00012012_c, poly00012012 } ,
      { 00012021_c, poly00012021 } ,
      { 00012110_c, poly00012110 } ,
      { 00012112_c, poly00012112 } ,
      { 00012120_c, poly00012120 } ,
      { 00012121_c, poly00012121 } ,
      { 00012122_c, poly00012122 } ,
      { 00012221_c, poly00012221 } ,
      { 00111100_c, poly00111100 } ,
      { 00111102_c, poly00111102 } ,
      { 00111220_c, poly00111220 } ,
      { 00121201_c, poly00121201 }
};



#endif // TRIPLE_LIGNES_H
