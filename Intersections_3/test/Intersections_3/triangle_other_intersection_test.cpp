// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb

#include <string>

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>




// -----------------------------------
// Kernels
// -----------------------------------
struct Sc_f : public CGAL::Simple_cartesian<float> {};
struct Sc_d : public CGAL::Simple_cartesian<double> {};
struct C_f : public CGAL::Cartesian<float> {};
struct C_d : public CGAL::Cartesian<double> {};
struct Epic : public CGAL::Exact_predicates_inexact_constructions_kernel {};
struct Epec : public CGAL::Exact_predicates_exact_constructions_kernel {};



// -----------------------------------
// Random intersection tests
// -----------------------------------

// Checkers (partial specialization for epic & epec)
template <class K>
struct Checker
{
  template <typename Query>
  void operator()(const Query& q, const typename K::Triangle_3& t) const
  {
    CGAL::Object result = CGAL::intersection(q, t);

    if ( ! result.empty() )
    {
      assert(   nullptr != CGAL::object_cast<typename K::Point_3>(&result)
             || nullptr != CGAL::object_cast<typename K::Segment_3>(&result));
    }
  }
};

template <>
struct Checker<Epic>
{
  typedef Epic K;

  template <typename Query>
  void operator()(const Query& q, const typename K::Triangle_3& t) const
  {
    CGAL::Object result = CGAL::intersection(q, t);

    if ( ! result.empty() )
    {
      assert( CGAL::do_intersect(q, t) );
      assert(   nullptr != CGAL::object_cast<typename K::Point_3>(&result)
             || nullptr != CGAL::object_cast<typename K::Segment_3>(&result));
    }
  }

  void operator()(const K::Line_3& l, const K::Triangle_3& t) const
  {
    CGAL::Object result = CGAL::intersection(l, t);

    if ( ! result.empty() )
    {
      // Here we can't check do_intersect, because there are constructions when
      // building points on line
      assert(   nullptr != CGAL::object_cast<K::Point_3>(&result)
             || nullptr != CGAL::object_cast<K::Segment_3>(&result));
    }
  }
};

template <>
struct Checker<Epec>
{
  typedef Epec K;

  template <typename Query>
  void operator()(const Query& q, const typename K::Triangle_3& t) const
  {
    typedef typename K::Point_3 Point_3;
    typedef typename K::Segment_3 Segment_3;

    CGAL::Object result = CGAL::intersection(q, t);

    if ( ! result.empty() )
    {
      assert( CGAL::do_intersect(q, t) );

      // Verify answer is correct
      const Point_3* p = CGAL::object_cast<Point_3>(&result);
      const Segment_3* s = CGAL::object_cast<Segment_3>(&result);

      assert(   (nullptr!=p && t.has_on(*p) && q.has_on(*p))
             || (nullptr!=s && t.has_on(s->source()) && t.has_on(s->target())
                         && q.has_on(s->source()) && q.has_on(s->target())) );
    }
    else
    {
      assert ( !CGAL::do_intersect(q, t) );
    }
  }
};



// random number generation
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

// random_test()
template <class K>
void random_test()
{
  typedef typename K::Point_3 Point;
  typedef typename K::Segment_3 Segment;
  typedef typename K::Ray_3 Ray;
  typedef typename K::Line_3 Line;
  typedef typename K::Triangle_3 Triangle;
  typedef typename K::Plane_3 Plane;

  typename K::Is_degenerate_3 is_degenerate = K().is_degenerate_3_object();

  Checker<K> check;

  double box_size = 1e12;
  CGAL::Bbox_3 bbox(-box_size,-box_size,-box_size,box_size,box_size,box_size);

  // Use 10 triangles, 100 queries for each triangle
  for ( int i=0 ; i<10 ; ++i )
  {
    Triangle t(random_point_in<K>(bbox),
               random_point_in<K>(bbox),
               random_point_in<K>(bbox));

    if ( is_degenerate(t) )
      continue;

    Plane p = t.supporting_plane();

    for ( int j=0 ; j<100 ; ++j )
    {
      Point a = random_point_in<K>(bbox);
      Point b = random_point_in<K>(bbox);

      Segment s (a,b);
      Ray r(a,b);
      Line l (a,b);

      if ( ! is_degenerate(s) )
        check(s,t);
      if ( ! is_degenerate(r) )
        check(r,t);
      if ( ! is_degenerate(l) )
        check(l,t);

      // Project points on triangle plane to have degenerate queries
      Point c = p.projection(a);
      Point d = p.projection(b);

      Segment s2 (c,d);
      Ray r2 (c,d);
      Line l2 (c,d);

      if ( ! is_degenerate(s2) )
        check(s2,t);
      if ( ! is_degenerate(r2) )
        check(r2,t);
      if ( ! is_degenerate(l2) )
        check(l2,t);
    }
  }
}



// -----------------------------------
// Precomputed results test
// -----------------------------------
template <class Triangle, class Query, class Result>
bool test_aux(bool is_kernel_exact,
              const Triangle t,
              const Query& q,
              const std::string& name,
              const Result& expected,
              double sq_espilon = 1e-20)
{
  CGAL::Object object = CGAL::intersection(t,q);
  const Result* pr = CGAL::object_cast<Result>(&object);

  if ( (nullptr != pr) &&
       (is_kernel_exact ?
        (expected == *pr) :
        CGAL::to_double(CGAL::squared_distance(expected, *pr)) < sq_espilon ) )
  {
    return true;
  }
  else
  {
    std::cout << "ERROR: intersection(" << name
    << ") did not answer the expected result !";

    if ( nullptr != pr )
      std::cout << " (answer: ["<< *pr << "])";

    std::cout << std::endl;
  }

  return false;
}

template <class K>
bool test(bool is_kernel_exact = true)
{
        // types
  typedef typename K::FT FT;
  typedef typename K::Line_3 Line;
  typedef typename K::Point_3 Point;
  typedef typename K::Segment_3 Segment;
  typedef typename K::Ray_3 Ray;
  typedef typename K::Line_3 Line;
  typedef typename K::Triangle_3 Triangle;

  /* -------------------------------------
  // Test data is something like that (in t supporting plane)
  // Triangle is (p1,p2,p3)
  //
  //       +E          +1
  //                 /   \
  //        +C     6+  +8  +4      +B
  //              /   9++7  \
  //            3+-------+5--+2
  //
  //         +F        +A
  ------------------------------------- */

  Point p1(FT(1.), FT(0.), FT(0.));
  Point p2(FT(0.), FT(1.), FT(0.));
  Point p3(FT(0.), FT(0.), FT(1.));

  Triangle t(p1,p2,p3);

  // Edges of t
  Segment s12(p1,p2);
  Segment s21(p2,p1);
  Segment s13(p1,p3);
  Segment s23(p2,p3);
  Segment s32(p3,p2);
  Segment s31(p3,p1);

  bool b = test_aux(is_kernel_exact,t,s12,"t-s12",s12);
  b &= test_aux(is_kernel_exact,t,s21,"t-s21",s21);
  b &= test_aux(is_kernel_exact,t,s13,"t-s13",s13);
  b &= test_aux(is_kernel_exact,t,s23,"t-s23",s23);

  // Inside points
  Point p4(FT(0.5), FT(0.5), FT(0.));
  Point p5(FT(0.), FT(0.75), FT(0.25));
  Point p6(FT(0.5), FT(0.), FT(0.5));
  Point p7(FT(0.25), FT(0.625), FT(0.125));
  Point p8(FT(0.5), FT(0.25), FT(0.25));

  Segment s14(p1,p4);
  Segment s41(p4,p1);
  Segment s24(p2,p4);
  Segment s42(p4,p2);
  Segment s15(p1,p5);
  Segment s25(p2,p5);
  Segment s34(p3,p4);
  Segment s35(p3,p5);
  Segment s36(p3,p6);
  Segment s45(p4,p5);
  Segment s16(p1,p6);
  Segment s26(p2,p6);
  Segment s62(p6,p2);
  Segment s46(p4,p6);
  Segment s48(p4,p8);
  Segment s56(p5,p6);
  Segment s65(p6,p5);
  Segment s64(p6,p4);
  Segment s17(p1,p7);
  Segment s67(p6,p7);
  Segment s68(p6,p8);
  Segment s86(p8,p6);
  Segment s78(p7,p8);
  Segment s87(p8,p7);

  b &= test_aux(is_kernel_exact,t,s14,"t-s14",s14);
  b &= test_aux(is_kernel_exact,t,s41,"t-s41",s41);
  b &= test_aux(is_kernel_exact,t,s24,"t-s24",s24);
  b &= test_aux(is_kernel_exact,t,s42,"t-s42",s42);
  b &= test_aux(is_kernel_exact,t,s15,"t-s15",s15);
  b &= test_aux(is_kernel_exact,t,s25,"t-s25",s25);
  b &= test_aux(is_kernel_exact,t,s34,"t-s34",s34);
  b &= test_aux(is_kernel_exact,t,s35,"t-s35",s35);
  b &= test_aux(is_kernel_exact,t,s36,"t-s36",s36);
  b &= test_aux(is_kernel_exact,t,s45,"t-s45",s45);
  b &= test_aux(is_kernel_exact,t,s16,"t-s16",s16);
  b &= test_aux(is_kernel_exact,t,s26,"t-s26",s26);
  b &= test_aux(is_kernel_exact,t,s62,"t-s62",s62);
  b &= test_aux(is_kernel_exact,t,s46,"t-s46",s46);
  b &= test_aux(is_kernel_exact,t,s65,"t-s65",s65);
  b &= test_aux(is_kernel_exact,t,s64,"t-s64",s64);
  b &= test_aux(is_kernel_exact,t,s48,"t-s48",s48);
  b &= test_aux(is_kernel_exact,t,s56,"t-s56",s56);
  b &= test_aux(is_kernel_exact,t,s17,"t-t17",s17);
  b &= test_aux(is_kernel_exact,t,s67,"t-t67",s67);
  b &= test_aux(is_kernel_exact,t,s68,"t-s68",s68);
  b &= test_aux(is_kernel_exact,t,s86,"t-s86",s86);
  b &= test_aux(is_kernel_exact,t,s78,"t-t78",s78);
  b &= test_aux(is_kernel_exact,t,s87,"t-t87",s87);

  // Outside points (in triangle plane)
  Point pA(FT(-0.5), FT(1.), FT(0.5));
  Point pB(FT(0.5), FT(1.), FT(-0.5));
  Point pC(FT(0.5), FT(-0.5), FT(1.));
  Point pE(FT(1.), FT(-1.), FT(1.));
  Point pF(FT(-1.), FT(0.), FT(2.));

  Segment sAB(pA,pB);
  Segment sBC(pB,pC);
  Segment s2E(p2,pE);
  Segment sE2(pE,p2);
  Segment s2A(p2,pA);
  Segment s6E(p6,pE);
  Segment sB8(pB,p8);
  Segment sC8(pC,p8);
  Segment s8C(p8,pC);
  Segment s1F(p1,pF);
  Segment sF6(pF,p6);

  b &= test_aux(is_kernel_exact,t,sAB,"t-sAB",p2);
  b &= test_aux(is_kernel_exact,t,sBC,"t-sBC",s46);
  b &= test_aux(is_kernel_exact,t,s2E,"t-s2E",s26);
  b &= test_aux(is_kernel_exact,t,sE2,"t-sE2",s62);
  b &= test_aux(is_kernel_exact,t,s2A,"t-s2A",p2);
  b &= test_aux(is_kernel_exact,t,s6E,"t-s6E",p6);
  b &= test_aux(is_kernel_exact,t,sB8,"t-sB8",s48);
  b &= test_aux(is_kernel_exact,t,sC8,"t-sC8",s68);
  b &= test_aux(is_kernel_exact,t,s8C,"t-s8C",s86);
  b &= test_aux(is_kernel_exact,t,s1F,"t-s1F",s13);
  b &= test_aux(is_kernel_exact,t,sF6,"t-sF6",s36);

  // Outside triangle plane
  Point pa(FT(0.), FT(0.), FT(0.));
  Point pb(FT(2.), FT(0.), FT(0.));
  Point pc(FT(1.), FT(0.), FT(1.));
  Point pe(FT(1.), FT(0.5), FT(0.5));

  Segment sab(pa,pb);
  Segment sac(pa,pc);
  Segment sae(pa,pe);
  Segment sa8(pa,p8);
  Segment sb2(pb,p2);

  b &= test_aux(is_kernel_exact,t,sab,"t-sab",p1);
  b &= test_aux(is_kernel_exact,t,sac,"t-sac",p6);
  b &= test_aux(is_kernel_exact,t,sae,"t-sae",p8);
  b &= test_aux(is_kernel_exact,t,sa8,"t-sa8",p8);
  b &= test_aux(is_kernel_exact,t,sb2,"t-sb2",p2);

  // -----------------------------------
  // ray queries
  // -----------------------------------
  // Edges of t
  Ray r12(p1,p2);
  Ray r21(p2,p1);
  Ray r13(p1,p3);
  Ray r23(p2,p3);

  b &= test_aux(is_kernel_exact,t,r12,"t-r12",s12);
  b &= test_aux(is_kernel_exact,t,r21,"t-r21",s21);
  b &= test_aux(is_kernel_exact,t,r13,"t-r13",s13);
  b &= test_aux(is_kernel_exact,t,r23,"t-r23",s23);

  // In triangle
  Point p9_(FT(0.), FT(0.5), FT(0.5));
  Point p9(FT(0.25), FT(0.375), FT(0.375));

  Ray r14(p1,p4);
  Ray r41(p4,p1);
  Ray r24(p2,p4);
  Ray r42(p4,p2);
  Ray r15(p1,p5);
  Ray r25(p2,p5);
  Ray r34(p3,p4);
  Ray r35(p3,p5);
  Ray r36(p3,p6);
  Ray r45(p4,p5);
  Ray r16(p1,p6);
  Ray r26(p2,p6);
  Ray r62(p6,p2);
  Ray r46(p4,p6);
  Ray r48(p4,p8);
  Ray r56(p5,p6);
  Ray r47(p4,p7);
  Ray r89(p8,p9);
  Ray r86(p8,p6);
  Ray r68(p6,p8);
  Segment r89_res(p8,p9_);

  b &= test_aux(is_kernel_exact,t,r14,"t-r14",s12);
  b &= test_aux(is_kernel_exact,t,r41,"t-r41",s41);
  b &= test_aux(is_kernel_exact,t,r24,"t-r24",s21);
  b &= test_aux(is_kernel_exact,t,r42,"t-r42",s42);
  b &= test_aux(is_kernel_exact,t,r15,"t-r15",s15);
  b &= test_aux(is_kernel_exact,t,r25,"t-r25",s23);
  b &= test_aux(is_kernel_exact,t,r34,"t-r34",s34);
  b &= test_aux(is_kernel_exact,t,r35,"t-r35",s32);
  b &= test_aux(is_kernel_exact,t,r36,"t-r36",s31);
  b &= test_aux(is_kernel_exact,t,r45,"t-r45",s45);
  b &= test_aux(is_kernel_exact,t,r16,"t-r16",s13);
  b &= test_aux(is_kernel_exact,t,r26,"t-r26",s26);
  b &= test_aux(is_kernel_exact,t,r62,"t-r62",s62);
  b &= test_aux(is_kernel_exact,t,r46,"t-r46",s46);
  b &= test_aux(is_kernel_exact,t,r48,"t-r48",s46);
  b &= test_aux(is_kernel_exact,t,r56,"t-r56",s56);
  b &= test_aux(is_kernel_exact,t,r47,"t-r47",s45);
  b &= test_aux(is_kernel_exact,t,r89,"t-t89",r89_res);
  b &= test_aux(is_kernel_exact,t,r68,"t-r68",s64);
  b &= test_aux(is_kernel_exact,t,r86,"t-r86",s86);


  // Outside points (in triangre prane)
  Ray rAB(pA,pB);
  Ray rBC(pB,pC);
  Ray r2E(p2,pE);
  Ray rE2(pE,p2);
  Ray r2A(p2,pA);
  Ray r6E(p6,pE);
  Ray rB8(pB,p8);
  Ray rC8(pC,p8);
  Ray r8C(p8,pC);
  Ray r1F(p1,pF);
  Ray rF6(pF,p6);

  b &= test_aux(is_kernel_exact,t,rAB,"t-rAB",p2);
  b &= test_aux(is_kernel_exact,t,rBC,"t-rBC",s46);
  b &= test_aux(is_kernel_exact,t,r2E,"t-r2E",s26);
  b &= test_aux(is_kernel_exact,t,rE2,"t-rE2",s62);
  b &= test_aux(is_kernel_exact,t,r2A,"t-r2A",p2);
  b &= test_aux(is_kernel_exact,t,r6E,"t-r6E",p6);
  b &= test_aux(is_kernel_exact,t,rB8,"t-rB8",s46);
  b &= test_aux(is_kernel_exact,t,rC8,"t-rC8",s64);
  b &= test_aux(is_kernel_exact,t,r8C,"t-r8C",s86);
  b &= test_aux(is_kernel_exact,t,r1F,"t-r1F",s13);
  b &= test_aux(is_kernel_exact,t,rF6,"t-rF6",s31);

  // Outside triangle plane
  Ray rab(pa,pb);
  Ray rac(pa,pc);
  Ray rae(pa,pe);
  Ray ra8(pa,p8);
  Ray rb2(pb,p2);

  b &= test_aux(is_kernel_exact,t,rab,"t-rab",p1);
  b &= test_aux(is_kernel_exact,t,rac,"t-rac",p6);
  b &= test_aux(is_kernel_exact,t,rae,"t-rae",p8);
  b &= test_aux(is_kernel_exact,t,ra8,"t-ra8",p8);
  b &= test_aux(is_kernel_exact,t,rb2,"t-rb2",p2);

  // -----------------------------------
  // Line queries
  // -----------------------------------
  // Edges of t
  Line l12(p1,p2);
  Line l21(p2,p1);
  Line l13(p1,p3);
  Line l23(p2,p3);

  b &= test_aux(is_kernel_exact,t,l12,"t-l12",s12);
  b &= test_aux(is_kernel_exact,t,l21,"t-l21",s21);
  b &= test_aux(is_kernel_exact,t,l13,"t-l13",s13);
  b &= test_aux(is_kernel_exact,t,l23,"t-l23",s23);

  // In triangle
  Line l14(p1,p4);
  Line l41(p4,p1);
  Line l24(p2,p4);
  Line l42(p4,p2);
  Line l15(p1,p5);
  Line l25(p2,p5);
  Line l34(p3,p4);
  Line l35(p3,p5);
  Line l36(p3,p6);
  Line l45(p4,p5);
  Line l16(p1,p6);
  Line l26(p2,p6);
  Line l62(p6,p2);
  Line l46(p4,p6);
  Line l48(p4,p8);
  Line l56(p5,p6);
  Line l47(p4,p7);
  Line l89(p8,p9);
  Line l86(p8,p6);
  Line l68(p6,p8);
  Segment l89_res(p1,p9_);


  b &= test_aux(is_kernel_exact,t,l14,"t-l14",s12);
  b &= test_aux(is_kernel_exact,t,l41,"t-l41",s21);
  b &= test_aux(is_kernel_exact,t,l24,"t-l24",s21);
  b &= test_aux(is_kernel_exact,t,l42,"t-l42",s12);
  b &= test_aux(is_kernel_exact,t,l15,"t-l15",s15);
  b &= test_aux(is_kernel_exact,t,l25,"t-l25",s23);
  b &= test_aux(is_kernel_exact,t,l34,"t-l34",s34);
  b &= test_aux(is_kernel_exact,t,l35,"t-l35",s32);
  b &= test_aux(is_kernel_exact,t,l36,"t-l36",s31);
  b &= test_aux(is_kernel_exact,t,l45,"t-l45",s45);
  b &= test_aux(is_kernel_exact,t,l16,"t-l16",s13);
  b &= test_aux(is_kernel_exact,t,l26,"t-l26",s26);
  b &= test_aux(is_kernel_exact,t,l62,"t-l62",s62);
  b &= test_aux(is_kernel_exact,t,l46,"t-l46",s46);
  b &= test_aux(is_kernel_exact,t,l48,"t-l48",s46);
  b &= test_aux(is_kernel_exact,t,l56,"t-l56",s56);
  b &= test_aux(is_kernel_exact,t,l47,"t-l47",s45);
  b &= test_aux(is_kernel_exact,t,l89,"t-t89",l89_res);
  b &= test_aux(is_kernel_exact,t,l68,"t-l68",s64);
  b &= test_aux(is_kernel_exact,t,l86,"t-l86",s46);


  // Outside points (in triangle plane)
  Line lAB(pA,pB);
  Line lBC(pB,pC);
  Line l2E(p2,pE);
  Line lE2(pE,p2);
  Line l2A(p2,pA);
  Line l6E(p6,pE);
  Line lB8(pB,p8);
  Line lC8(pC,p8);
  Line l8C(p8,pC);
  Line l1F(p1,pF);
  Line lF6(pF,p6);

  b &= test_aux(is_kernel_exact,t,lAB,"t-lAB",p2);
  b &= test_aux(is_kernel_exact,t,lBC,"t-lBC",s46);
  b &= test_aux(is_kernel_exact,t,l2E,"t-l2E",s26);
  b &= test_aux(is_kernel_exact,t,lE2,"t-lE2",s62);
  b &= test_aux(is_kernel_exact,t,l2A,"t-l2A",p2);
  b &= test_aux(is_kernel_exact,t,l6E,"t-l6E",s26);
  b &= test_aux(is_kernel_exact,t,lB8,"t-lB8",s46);
  b &= test_aux(is_kernel_exact,t,lC8,"t-lC8",s64);
  b &= test_aux(is_kernel_exact,t,l8C,"t-l8C",s46);
  b &= test_aux(is_kernel_exact,t,l1F,"t-l1F",s13);
  b &= test_aux(is_kernel_exact,t,lF6,"t-lF6",s31);

  // Outside triangle plane
  Line lab(pa,pb);
  Line lac(pa,pc);
  Line lae(pa,pe);
  Line la8(pa,p8);
  Line lb2(pb,p2);

  b &= test_aux(is_kernel_exact,t,lab,"t-lab",p1);
  b &= test_aux(is_kernel_exact,t,lac,"t-lac",p6);
  b &= test_aux(is_kernel_exact,t,lae,"t-lae",p8);
  b &= test_aux(is_kernel_exact,t,la8,"t-la8",p8);
  b &= test_aux(is_kernel_exact,t,lb2,"t-lb2",p2);


        return b;
}



// -----------------------------------
// Main
// -----------------------------------
int main()
{
  // -----------------------------------
  // Test intersection results
  // -----------------------------------
  std::cout << "Test precomputed intersection results" << std::endl;
  std::cout << "\tTesting with Simple_cartesian<float>..." << std::endl ;
  bool b = test<Sc_f>(false);

  std::cout << "\tTesting with Simple_cartesian<double>..." << std::endl ;
        b &= test<Sc_d>(false);

  std::cout << "\tTesting with Cartesian<float>..." << std::endl ;
        b &= test<C_f>(false);

  std::cout << "\tTesting with Cartesian<double>..." << std::endl ;
        b &= test<C_d>(false);

  std::cout << "\tTesting with Exact_predicates_inexact_constructions_kernel..." << std::endl ;
  b &= test<Epic>(false);

  std::cout << "\tTesting with Exact_predicates_exact_constructions_kernel..." << std::endl ;
  b &= test<Epec>(true);
  //test with a coplanar segment
  b &= !bool(CGAL::intersection(
    Epec::Segment_3(Epec::Point_3(0.125, 0, -0.125),Epec::Point_3(0.25, 0, -0.125) ),
    Epec::Triangle_3( Epec::Point_3(0.2500001, 0, -0.125),
                      Epec::Point_3(1.0278171, 0, -0.125) /* vertex 10*/,
                      Epec::Point_3(1.0278171, 0, -0.250001) /* vertex 9*/ ) ));
  // -----------------------------------
  // Test random intersection
  // -----------------------------------
  srand( static_cast<unsigned int>(time(nullptr)) );
  std::cout << std::endl << "Test random intersections" << std::endl;
  std::cout << "\tTesting with Simple_cartesian<float>..." << std::endl ;
  random_test<Sc_f>();

  std::cout << "\tTesting with Simple_cartesian<double>..." << std::endl ;
        random_test<Sc_d>();

  std::cout << "\tTesting with Cartesian<float>..." << std::endl ;
        random_test<C_f>();

  std::cout << "\tTesting with Cartesian<double>..." << std::endl ;
        random_test<C_d>();

  std::cout << "\tTesting with Exact_predicates_inexact_constructions_kernel..." << std::endl ;
  random_test<Epic>();

  std::cout << "\tTesting with Exact_predicates_exact_constructions_kernel..." << std::endl ;
  random_test<Epec>();

  if ( b ) {
    return EXIT_SUCCESS;
  } else {
    return EXIT_FAILURE;
  }
}
