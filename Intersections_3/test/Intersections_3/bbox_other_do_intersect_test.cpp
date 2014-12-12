// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stephane Tayeb
//

#include <string>

#include <CGAL/config.h>

#if defined(BOOST_MSVC)
#  pragma warning(disable:4244) // int to float conversion warning
#endif  
#include <CGAL/Exact_rational.h>

// leda_rational, or Gmpq, or Quotient<MP_float>
typedef CGAL::Exact_rational         Rational;
#include <CGAL/intersections.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/Kernel_traits.h>

#include <iomanip>

#include <boost/math/special_functions/next.hpp> // for nextafter


double random_in(const double a,
                 const double b)
{
  double r = rand() / (double)RAND_MAX;
  return a + (b - a) * r;
}

template <class K>
typename K::Point_3 random_point_in(const CGAL::Bbox_3& bbox)
{
  typedef typename K::FT FT;
  FT x = (FT)random_in(bbox.xmin(),bbox.xmax());
  FT y = (FT)random_in(bbox.ymin(),bbox.ymax());
  FT z = (FT)random_in(bbox.zmin(),bbox.zmax());
  return typename K::Point_3(x,y,z);
}


template <class T>
bool test_aux(const T& t,
              const std::string& name,
              const CGAL::Bbox_3& bbox,
              bool expected, bool /*exact_predicates*/ = false)
{
  bool b = CGAL::do_intersect(t,bbox);
  
  if ( b != expected )
    std::cerr << "ERROR: do_intersect(" << name
              << ") did not answer the expected result !" << std::endl;
  
  return (b == expected);
}

template <class K, 
          class FT>
bool test_case(const FT& px, const FT& py, const FT& pz,
               const FT& qx, const FT& qy, const FT& qz,
               FT bxmin, FT bymin, FT bzmin,
               FT bxmax, FT bymax, FT bzmax,
               const bool expected, 
               const bool exact_k = false,
               const bool exactness_issue = false,
               const bool change_signs = true,
               const bool swap_coords  = true,
               const bool opposite_seg = true,
               const bool translate    = true,
               const bool scale        = true)
{
  bool b = true;
  if(change_signs) {
    b &= test_case<K>( px,     py,     pz,
                       qx,     qy,     qz,
                       bxmin,  bymin,  bzmin,
                       bxmax,  bymax,  bzmax, expected, exact_k, exactness_issue, false);
    b &= test_case<K>(-px,     py,     pz,
                      -qx,     qy,     qz,
                      -bxmin,  bymin,  bzmin,
                      -bxmax,  bymax,  bzmax, expected, exact_k, exactness_issue, false);
    b &= test_case<K>( px,    -py,     pz,
                       qx,    -qy,     qz,
                       bxmin, -bymin,  bzmin,
                       bxmax, -bymax,  bzmax, expected, exact_k, exactness_issue, false);
    b &= test_case<K>( px,     py,    -pz,
                       qx,     qy,    -qz,
                       bxmin,  bymin, -bzmin,
                       bxmax,  bymax, -bzmax, expected, exact_k, exactness_issue, false);
    b &= test_case<K>(-px,    -py,     pz,
                      -qx,    -qy,     qz,
                      -bxmin, -bymin,  bzmin,
                      -bxmax, -bymax,  bzmax, expected, exact_k, exactness_issue, false);
    b &= test_case<K>( px,    -py,    -pz,
                       qx,    -qy,    -qz,
                       bxmin, -bymin, -bzmin,
                       bxmax, -bymax, -bzmax, expected, exact_k, exactness_issue, false);
    b &= test_case<K>(-px,     py,    -pz,
                      -qx,     qy,    -qz,
                      -bxmin,  bymin, -bzmin,
                      -bxmax,  bymax, -bzmax, expected, exact_k, exactness_issue, false);
    b &= test_case<K>(-px,    -py,    -pz,
                      -qx,    -qy,    -qz,
                      -bxmin, -bymin, -bzmin,
                      -bxmax, -bymax, -bzmax, expected, exact_k, exactness_issue, false);
  } else if(swap_coords) {
    // xyz
    b &= test_case<K>( px,     py,     pz,
                       qx,     qy,     qz,
                       bxmin,  bymin,  bzmin,
                       bxmax,  bymax,  bzmax, expected, exact_k, exactness_issue,
                       false, false);
    // xzy
    b &= test_case<K>( px,     pz,     py,
                       qx,     qz,     qy,
                       bxmin,  bzmin,  bymin,
                       bxmax,  bzmax,  bymax, expected, exact_k, exactness_issue,
                       false, false);
    // yxz 
    b &= test_case<K>( py,     px,     pz,
                       qy,     qx,     qz,
                       bymin,  bxmin,  bzmin,
                       bymax,  bxmax,  bzmax, expected, exact_k, exactness_issue,
                       false, false);
    // zxy
    b &= test_case<K>( pz,     px,     py,
                       qz,     qx,     qy,
                       bzmin,  bxmin,  bymin,
                       bzmax,  bxmax,  bymax, expected, exact_k, exactness_issue,
                       false, false);

    // yzx 
    b &= test_case<K>( py,     pz,     px,
                       qy,     qz,     qx,
                       bymin,  bzmin,  bxmin,
                       bymax,  bzmax,  bxmax, expected, exact_k, exactness_issue,
                       false, false);
    // zyx
    b &= test_case<K>( pz,     py,     px,
                       qz,     qy,     qx,
                       bzmin,  bymin,  bxmin,
                       bzmax,  bymax,  bxmax, expected, exact_k, exactness_issue,
                       false, false);
  } else if(opposite_seg) {
    b &= test_case<K>(px,     py,     pz,
                      qx,     qy,     qz,
                      bxmin,  bymin,  bzmin,
                      bxmax,  bymax,  bzmax, expected, exact_k, exactness_issue,
                      false, false, false);
    b &= test_case<K>(qx,     qy,     qz,
                      px,     py,     pz,
                      bxmin,  bymin,  bzmin,
                      bxmax,  bymax,  bzmax, expected, exact_k, exactness_issue,
                      false, false, false);
  } else if(translate) {
    b &= test_case<K>(px,     py,     pz,
                      qx,     qy,     qz,
                      bxmin,  bymin,  bzmin,
                      bxmax,  bymax,  bzmax, expected, exact_k, exactness_issue,
                      false, false, false, false);
    FT dx = 10, dy = 20, dz = 30;
    b &= test_case<K>(dx + px,     dy + py,     dz + pz,
                      dx + qx,     dy + qy,     dz + qz,
                      dx + bxmin,  dy + bymin,  dz + bzmin,
                      dx + bxmax,  dy + bymax,  dz + bzmax, expected, exact_k, exactness_issue,
                      false, false, false, false);    
    dx = (1 >> 10), dy = dx, dz = dx;
    b &= test_case<K>(dx + px,     dy + py,     dz + pz,
                      dx + qx,     dy + qy,     dz + qz,
                      dx + bxmin,  dy + bymin,  dz + bzmin,
                      dx + bxmax,  dy + bymax,  dz + bzmax, expected, exact_k, exactness_issue,
                      false, false, false, false);    
    dx = -(1 >> 10), dy = dx, dz = dx;
    b &= test_case<K>(dx + px,     dy + py,     dz + pz,
                      dx + qx,     dy + qy,     dz + qz,
                      dx + bxmin,  dy + bymin,  dz + bzmin,
                      dx + bxmax,  dy + bymax,  dz + bzmax, expected, exact_k, exactness_issue,
                      false, false, false, false);    
  } else if(scale) {
    b &= test_case<K>(px,     py,     pz,
                      qx,     qy,     qz,
                      bxmin,  bymin,  bzmin,
                      bxmax,  bymax,  bzmax, expected, exact_k, exactness_issue,
                      false, false, false, false, false);
    FT delta = 9;
    b &= test_case<K>(delta * qx,     delta * qy,     delta * qz,
                      delta * px,     delta * py,     delta * pz,
                      delta * bxmin,  delta * bymin,  delta * bzmin,
                      delta * bxmax,  delta * bymax,  delta * bzmax, expected, exact_k, exactness_issue,
                      false, false, false, false, false);
    delta = (1 << 10);
    b &= test_case<K>(delta * qx,     delta * qy,     delta * qz,
                      delta * px,     delta * py,     delta * pz,
                      delta * bxmin,  delta * bymin,  delta * bzmin,
                      delta * bxmax,  delta * bymax,  delta * bzmax, expected, exact_k, exactness_issue,
                      false, false, false, false, false);
    delta = (1 << 10);
    delta = 1/delta;
    b &= test_case<K>(delta * qx,     delta * qy,     delta * qz,
                      delta * px,     delta * py,     delta * pz,
                      delta * bxmin,  delta * bymin,  delta * bzmin,
                      delta * bxmax,  delta * bymax,  delta * bzmax, expected, exact_k, exactness_issue,
                      false, false, false, false, false);
    delta = 7;
    delta /= 2;
    b &= test_case<K>(delta * qx,     delta * qy,     delta * qz,
                      delta * px,     delta * py,     delta * pz,
                      delta * bxmin,  delta * bymin,  delta * bzmin,
                      delta * bxmax,  delta * bymax,  delta * bzmax, expected, exact_k, exactness_issue,
                      false, false, false, false, false);
    delta = 1;
    delta /= 8;
    b &= test_case<K>(delta * qx,     delta * qy,     delta * qz,
                      delta * px,     delta * py,     delta * pz,
                      delta * bxmin,  delta * bymin,  delta * bzmin,
                      delta * bxmax,  delta * bymax,  delta * bzmax, expected, exact_k, exactness_issue,
                      false, false, false, false, false);
  } else {
    using CGAL::do_intersect;
    using CGAL::Bbox_3;
    typedef typename K::Point_3 Point_3;
    typedef typename K::Segment_3 Segment_3;
    if(bxmin > bxmax) std::swap(bxmin, bxmax);
    if(bymin > bymax) std::swap(bymin, bymax);
    if(bzmin > bzmax) std::swap(bzmin, bzmax);
    if(do_intersect(Bbox_3(bxmin, bymin, bzmin,
                           bxmax, bymax, bzmax),
                    Segment_3(Point_3(px, py, pz),
                              Point_3(qx, qy, qz))) != expected)
    {
      if(!exactness_issue || exact_k) {
        b = false;
        CGAL::set_pretty_mode(std::cerr);
        std::cerr.precision(17);
        std::cerr << "Wrong result for do_intersect(" 
                  << Bbox_3(bxmin, bymin, bzmin,
                            bxmax, bymax, bzmax)
                  << ",\n"
                  << "                              "
                  << Segment_3(Point_3(px, py, pz),
                               Point_3(qx, qy, qz))
                  << ")\n"
                  << "  it should have been " << std::boolalpha << expected 
                  << std::endl;
      }
    }
  }
  return b;
}

template <class T>
void speed(const std::string& name)
{
  // types
  typedef typename CGAL::Kernel_traits<T>::Kernel K;
  typedef typename K::Point_3 Point;
  
  // speed
  double d1 = 1.44258699;
  CGAL::Bbox_3 bbox_small(-d1, -d1, -d1, d1, d1, d1);
  double d2 = 10;
  CGAL::Bbox_3 bbox_big(-d2,-d2,-d2,d2,d2,d2);
  
  std::vector<T> segment_vector;
  for ( int i = 0 ; i < 1e4 ; ++i )
  {
    Point source = random_point_in<K>(bbox_big);
    Point target = random_point_in<K>(bbox_big);
    
    segment_vector.push_back(T(source, target));
  }
  
  int nb_loops = 0;
  
  CGAL::Timer timer;
  timer.start();
  std::size_t success = 0;
  while ( timer.time() < 5. )
  {
    for ( typename std::vector<T>::iterator it = segment_vector.begin();
         it != segment_vector.end() ; ++it )
    {
      success += do_intersect(bbox_small, *it);
    }
    ++nb_loops;
  }
  timer.stop();

  std::cout << std::fixed << std::setprecision(1);
  std::cout << "\tDo_intersect(bbox, " << name << "): " 
            << (nb_loops*segment_vector.size()) / (timer.time()*1000) 
            << " computations / ms  " 
            << (success / ((0.+ nb_loops*segment_vector.size()) / 100))
            << "% of intersection" << std::endl;
}

template <class K>
void test_speed()
{
  typedef typename K::Segment_3 Segment;
  typedef typename K::Ray_3 Ray;
  typedef typename K::Line_3 Line;
  
  speed<Segment>("segment");
  speed<Ray>("ray");
  speed<Line>("line");
}

template <class K>
bool intensive_test(bool exact_predicates = true)
{
  bool b = true;

  // Test vertical segments
  for(double x = 4.; x <= 7.; x+=1.)
    for(double ymin = 0.; ymin <= 7.; ymin+=1.)
      for(double ymax = ymin; ymax <= 7.; ymax+=1.)
      {
        const bool expected = 
          x >= -5. && x<= 5. &&
          ymin <= 5. && ymax >= -5;
        b &= test_case<K>(x, ymin, 0.,
                          x, ymax, 0.,
                          -5., -5., -5., 5., 5., 5., expected, exact_predicates);
      }
  // Test slanted segments
  for(double x = -7.; x <= 6.; x+=1.)
    for(double y = -1.; y <= 6.; y+=1.)
    {
      const bool expected = 
        x >= -6. && x <= 5. &&
        y >= -6. && y <= 5. &&
        y <= x + 10. && y >= x - 10.;
      b &= test_case<K>(x, y, 0.,
                        x + 1., y + 1., 0.,
                        -5., -5., -5., 5., 5., 5., expected, exact_predicates);
    }
  for(double x = -9.; x <= 6.; x+=1.)
    for(double y = -3.; y <= 6.; y+=1.)
    {
      const bool expected = 
        x >= -8. && x <= 5. &&
        y >= -7. && y <= 5. &&
        3 * y <= 2 * x + 25. && 3 * y >= 2 * x - 25.;
      b &= test_case<K>(x, y, 0.,
                        x + 3., y + 2., 0.,
                        -5., -5., -5., 5., 5., 5., expected, exact_predicates);
    }
  return b;
}

template <class K>
bool test(bool exact_kernel = false)
{
	// types
  typedef typename K::FT FT;
  typedef typename K::Line_3 Line;
  typedef typename K::Point_3 Point;
  typedef typename K::Vector_3 Vector;
  typedef typename K::Segment_3 Segment;
  typedef typename K::Ray_3 Ray;
  typedef typename K::Line_3 Line;
  typedef typename K::Sphere_3 Sphere;
  typedef typename K::Plane_3 Plane;
  typedef typename K::Triangle_3 Triangle;
  
  CGAL::Bbox_3 bbox(1.0,1.0,1.0,10.0,50.0,100.0);
  
  Point p1(FT(0.), FT(0.), FT(0.));
  Point p2(FT(0.), FT(100.), FT(100.));
  Point p3(FT(100.), FT(0.), FT(100.));
  Point p4(FT(100.), FT(100.), FT(0.));
  Point p5(FT(0.), FT(0.), FT(100.));
  Point p6(FT(100.), FT(100.), FT(100.));
  Point p7(FT(-1.), FT(-1.), FT(-1.));
  
  Point pA(FT(5.), FT(0.), FT(0.));
  Point pB(FT(5.), FT(10.), FT(100.));
  Point pC(FT(5.), FT(10.), FT(-1.));
  Point pE(FT(5.), FT(10.), FT(0.));
  
  Segment s12(p1,p2);
  Segment s13(p1,p3);
  Segment s14(p1,p4);
  Segment s15(p1,p5);
  Segment s16(p1,p6);
  Segment s17(p1,p7);
  Segment s71(p7,p1);
  
  Segment sAB(pA,pB);
  Segment sBA(pB,pA);
  Segment sBC(pB,pC);
  Segment sCE(pC,pE);
  Segment sEC(pE,pC);
  
  bool b = test_aux(s12,"s12",bbox,false);
  b &= test_aux(s13,"s13",bbox,false);
  b &= test_aux(s14,"s14",bbox,false);
  b &= test_aux(s15,"s15",bbox,false);
  b &= test_aux(s16,"s16",bbox,true);
  b &= test_aux(s17,"s17",bbox,false);
  b &= test_aux(s71,"s71",bbox,false);
  
  b &= test_aux(sAB,"sAB",bbox,true);
  b &= test_aux(sBA,"sBA",bbox,true);
  b &= test_aux(sBC,"sBC",bbox,true);
  b &= test_aux(sCE,"sCE",bbox,false);
  b &= test_aux(sEC,"sEC",bbox,false);

  
  CGAL::Bbox_3 bbox_elem_1834(2, 1.81818, 0.166666,
                              2.18182, 2.18182, 0.333333);
  
  Point source_1834(2, 2, 1);
  Point target_1834(2, 2, 0.75);
  Segment segment_query_1834( source_1834, target_1834 );
  Point source_1834a(2, 2, 0.5);
  Point target_1834a(2, 2, 0.75);
  Segment segment_query_1834a( source_1834a, target_1834a );

  b &= test_aux(segment_query_1834,
                "segment_query_1834", bbox_elem_1834, false);
  b &= test_aux(segment_query_1834.opposite(),
                "segment_query_1834.opposite()", bbox_elem_1834, false);
  b &= test_aux(segment_query_1834a,
                "segment_query_1834a", bbox_elem_1834, false);
  b &= test_aux(segment_query_1834a.opposite(),
                "segment_query_1834a.opposite()", bbox_elem_1834, false);
  b &= test_case<K>(2., 2., 1.,
                    2., 2., 0.75,
                    1.81818, 2., 0.,
                    2., 2.18182, 0.333333, false); 
  b &= test_case<K>(2., 2., 0.5,
                    2., 2., 0.75,
                    1.81818, 2., 0.,
                    2., 2.18182, 0.333333, false); 

  CGAL::Bbox_3 bbox_elem_1834b(1.81818, 2, 0,
                               2, 2.18182, 0.333333);
  Segment segment_query_1834b(Point(2, 2, 1),
                              Point(2, 2, 0.75));
  b &= test_aux(segment_query_1834b,
                "segment_query_1834b", bbox_elem_1834b, false);
  b &= test_aux(segment_query_1834b.opposite(),
                "segment_query_1834b.opposite()", bbox_elem_1834b, false);
  b &= test_case<K>(2., 2., 1.,
                    2., 2., 0.75,
                    1.81818, 2., 0.,
                    2., 2.18182, 0.333333, false); 

  Ray r12(p1,p2);
  Ray r13(p1,p3);
  Ray r14(p1,p4);
  Ray r15(p1,p5);
  Ray r16(p1,p6);
  Ray r17(p1,p7);
  Ray r71(p7,p1);
  
  Ray rAB(pA,pB);
  Ray rBA(pB,pA);
  Ray rBC(pB,pC);
  Ray rCE(pC,pE);
  Ray rEC(pE,pC);
  
  b &= test_aux(r12,"r12",bbox,false);
  b &= test_aux(r13,"r13",bbox,false);
  b &= test_aux(r14,"r14",bbox,false);
  b &= test_aux(r15,"r15",bbox,false);
  b &= test_aux(r16,"r16",bbox,true);
  b &= test_aux(r17,"r17",bbox,false);
  b &= test_aux(r71,"r71",bbox,true);
  
  b &= test_aux(rAB,"rAB",bbox,true);
  b &= test_aux(rBA,"rBA",bbox,true);
  b &= test_aux(rBC,"rBC",bbox,true);
  b &= test_aux(rCE,"rCE",bbox,true);
  b &= test_aux(rEC,"rEC",bbox,false);
    
  Line l12(p1,p2);
  Line l13(p1,p3);
  Line l14(p1,p4);
  Line l15(p1,p5);
  Line l16(p1,p6);
  Line l17(p1,p7);
  Line l71(p7,p1);
  
  Line lAB(pA,pB);
  Line lBA(pB,pA);
  Line lBC(pB,pC);
  Line lCE(pC,pE);
  Line lEC(pE,pC);
  
  b &= test_aux(l12,"l12",bbox,false);
  b &= test_aux(l13,"l13",bbox,false);
  b &= test_aux(l14,"l14",bbox,false);
  b &= test_aux(l15,"l15",bbox,false);
  b &= test_aux(l16,"l16",bbox,true);
  b &= test_aux(l17,"l17",bbox,true);
  b &= test_aux(l71,"l71",bbox,true);
  
  b &= test_aux(lAB,"lAB",bbox,true);
  b &= test_aux(lBA,"lBA",bbox,true);
  b &= test_aux(lBC,"lBC",bbox,true);
  b &= test_aux(lCE,"lCE",bbox,true);
  b &= test_aux(lEC,"lEC",bbox,true);


  Sphere sphA_1(pA,1);
  Sphere sphB_1(pB,1);
  Sphere sphC_1(pC,1);
  Sphere sphE_1(pE,1);
  Sphere sph1_3(p1,3);
  
  b &= test_aux(sphA_1,"sphA_1",bbox,false);
  b &= test_aux(sphB_1,"sphB_1",bbox,true);
  b &= test_aux(sphC_1,"sphC_1",bbox,false);
  b &= test_aux(sphE_1,"sphE_1",bbox,true);
  b &= test_aux(sph1_3,"sph1_3",bbox,true);
  
  Plane Pl1(Point(1,1,1),Vector(0,0,1));
  Plane Pl2(Point(1,1,1),Vector(1,1,1));
  Plane Pl3(Point(5,5,5),Vector(1,2,-5));
  Plane Pl4(Point(11,50,100),Vector(1,1,1));
  
  b &= test_aux(Pl1,"Pl1",bbox,true);
  b &= test_aux(Pl2,"Pl2",bbox,true);
  b &= test_aux(Pl3,"Pl3",bbox,true);
  b &= test_aux(Pl4,"Pl4",bbox,false);
  
  Triangle t123(p1,p2,p3);
  Triangle t124(p1,p2,p4);
  Triangle t126(p1,p2,p6);
  Triangle t136(p1,p3,p6);
  Triangle tABC(pA,pB,pC);
  Triangle t1(Point(1,1,1),Point(0,0,0),Point(0,0,1));
  Triangle t2(Point(4,1,7),Point(8,1,99),Point(7,1,11));
  Triangle t3(Point(0,1,1),Point(0,0,0),Point(0,0,1));

  b &= test_aux(t123,"t123",bbox,true);
  b &= test_aux(t124,"t124",bbox,true);
  b &= test_aux(t126,"t126",bbox,true);
  b &= test_aux(t136,"t136",bbox,true);
  b &= test_aux(tABC,"tABC",bbox,true);
  b &= test_aux(t1,"t1",bbox,true);
  b &= test_aux(t2,"t2",bbox,true);
  b &= test_aux(t3,"t3",bbox,false);
  
  
  // Test more bboxes
  CGAL::Bbox_3 bbox2(-0.248143,-0.49325,0.0747943,-0.107021,-0.406955,0.151042);
  Segment seg2(Point(FT(0.10114),FT(0.23963),FT(0.0854394)),
               Point(FT(1.35831),FT(1.52921),FT(0.524127)) );
  
  CGAL::Bbox_3 bbox3(-0.489613, -0.333874, -0.154123, -0.100856, 0.0477374, 0.155364);
  Segment seg3(Point(FT(-0.123786),FT(0.0689497),FT(0.274589)),
               Point(FT(-0.203405),FT(-0.119905),FT(0.0661125)) );

  CGAL::Bbox_3 bbox4(-0.161409, -0.462146, -0.104619, -0.111159, -0.408844, -0.0594407);
  Segment seg4(Point(FT(0.0586705),FT(-0.276077),FT(0.0726862)),
               Point(FT(-1.12212),FT(-1.19194),FT(-1.25109)) );

  Ray ray2(seg2.source(), seg2.target());
  Ray ray3(seg3.source(), seg3.target());
  Ray ray4(seg4.source(), seg4.target());
  
  Line line2(seg2);
  Line line3(seg3);
  Line line4(seg4);
  
  b &= test_aux(seg2, "seg2", bbox2, false);
  b &= test_aux(seg3, "seg3", bbox3, true);
  b &= test_aux(seg4, "seg4", bbox4, false);
  
  b &= test_aux(ray2, "ray2", bbox2, false);
  b &= test_aux(ray3, "ray3", bbox3, true);
  b &= test_aux(ray4, "ray4", bbox4, false);
  
  b &= test_aux(line2, "line2", bbox2, false);
  b &= test_aux(line3, "line3", bbox3, true);
  b &= test_aux(line4, "line4", bbox4, false);
  
  // Use do_intersect(bbox,bbox)
  CGAL::do_intersect(bbox2,bbox4);

  b &= test_case<K>(1., 1., 0.,
                    1., 1., 1.,
                    0., 0., 0.,
                    1., 1., 1., true);
  b &= test_case<K>(0.5, 0.5, -0.5,
                    0.5, 0.5, 0.5,
                    -0.5, -0.5, -0.5,
                    0.5, 0.5, 0.5, true);
  float f = 0.5f;
  double d = boost::math::nextafter(f, f+1);
  double d2 = boost::math::nextafter(f, f-1);
  b &= test_case<K>(d, 0.5, -0.5,
                    d, 0.5, 0.5,
                    -0.5, -0.5, -0.5,
                    0.5, 0.5, 0.5, false, exact_kernel, true,
                    false, false, false, false, false);
  b &= test_case<K>(d2, 0.5, -0.5,
                    d, 0.5, 0.5,
                    -0.5, -0.5, -0.5,
                    0.5, 0.5, 0.5, true, exact_kernel, true,
                    false, false, false, false, false);

  b &= test_case<K>(1., 1., 0.,
                    2., 2., 2.,
                    0., 0., 0.,
                    1., 1., 1., true);
  b &= test_case<K>(1., 1., 1.,
                    1., 1., 1.,
                    0., 0., 0.,
                    1., 1., 1., true);
  b &= test_case<K>(0.9, 0.9, 0.9,
                    0.9, 0.9, 0.9,
                    0., 0., 0.,
                    1., 1., 1., true);
  b &= test_case<K>(0., 0., 0.,
                    0., 0., 0.,
                    0., 0., 0.,
                    1., 1., 1., true);
  b &= test_case<K>(0.1, 0., 0.1,
                    0.1, 0., 0.1,
                    0., 0., 0.,
                    1., 1., 1., true);
  b &= test_case<K>(0.1, -0.1, 0.1,
                    0.1, -0.1, 0.1,
                    0., 0., 0.,
                    1., 1., 1., false);
  b &= test_case<K>(0.1, 0.1, 0.1,
                    0.1, 0.1, 0.1,
                    0., 0., 0.,
                    1., 1., 1., true);
  b &= test_case<K>(1., 1., 1.1,
                    1., 1., 1.1,
                    0., 0., 0.,
                    1., 1., 1., false);
  return b;
}

template <typename K>
bool test_kernel(bool exact_predicates = true, K /*k*/ = K())
{
  bool b = test<K>(exact_predicates) &&
    intensive_test<K>(exact_predicates);
  test_speed<K>();
  return b;
}

int main()
{
  srand(0);
  std::cout << std::setprecision(5);
  
  bool b;
  std::cout << "Testing with Simple_cartesian<float>..." << std::endl ;
  b = test_kernel<CGAL::Simple_cartesian<float> >(false);
  
  std::cout << std::endl << "Testing with Simple_cartesian<double>..." << std::endl ;
  b &= test_kernel<CGAL::Simple_cartesian<double> >(true);

  std::cout << std::endl << "Testing with Simple_cartesian<Rational>..." << std::endl ;
  b &= test_kernel<CGAL::Simple_cartesian<Rational> >(true);

  std::cout << std::endl << "Testing with Cartesian<float>..." << std::endl ;
  b &= test_kernel<CGAL::Cartesian<float> >(false);
  
  std::cout << std::endl << "Testing with Cartesian<double>..." << std::endl ;
  b &= test_kernel<CGAL::Cartesian<double> >(true);
  
  std::cout << std::endl << "Testing with Filtered_kernel<Simple_cartesian<double> > without static filters..." << std::endl ;
  typedef CGAL::Filtered_kernel<CGAL::Simple_cartesian<double>, false> Fk_no_static;
  b &= test_kernel<Fk_no_static>();

  std::cout << std::endl << "Testing with Exact_predicates_inexact_constructions_kernel..." << std::endl ;
  b &= test_kernel<CGAL::Exact_predicates_inexact_constructions_kernel>();
  
  std::cout << std::endl << "Testing with Exact_predicates_exact_constructions_kernel..." << std::endl ;
  b &= test_kernel<CGAL::Exact_predicates_exact_constructions_kernel>();
  
  if ( b )
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}
