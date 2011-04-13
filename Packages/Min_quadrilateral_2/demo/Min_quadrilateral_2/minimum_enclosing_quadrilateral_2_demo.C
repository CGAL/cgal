// ============================================================================
//
// Copyright (c) 1999, 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : minimum_enclosing_quadrilateral_2_demo.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Min_quadrilaterals $
// source        : oops.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch> and
//                 Emo Welzl <emo@inf.ethz.ch>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// Demo Program: Computing minimum enclosing quadrilaterals
// ============================================================================

#include <CGAL/basic.h>
#ifdef CGAL_USE_LEDA
#include <LEDA/basic.h>
#endif

#if defined(CGAL_USE_LEDA) && (__LEDA__ >= 400)

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/min_quadrilateral_2.h>
#include <vector>
#include <iostream>


//using CGAL::to_double;
#include <CGAL/Segment_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/geowin_support.h>
#include <CGAL/Timer.h>
#include <CGAL/leda_real.h>
#include <list>

#ifdef _MSC_VER
// the general templates from geowin support seem not to work here
leda_point
convert_to_leda(const CGAL::Point_2<CGAL::Cartesian<leda_real> >& obj)
{
  double x = CGAL::to_double(obj.x());
  double y = CGAL::to_double(obj.y());
  leda_point p(x,y);
  return p;
}
#endif // _MSC_VER


using CGAL::Polygon_traits_2;
using CGAL::Creator_uniform_2;
using CGAL::Random_points_in_square_2;
using CGAL::random_convex_set_2;
using CGAL::min_rectangle_2;
using CGAL::min_parallelogram_2;
using CGAL::min_strip_2;
using std::back_inserter;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::flush;
using std::list;
using CGAL::Timer;
using CGAL::convex_hull_points_2;
using CGAL::squared_distance;

#include <LEDA/list.h>

//typedef CGAL::Cartesian< double >                      R;
typedef CGAL::Cartesian< leda_real >                   R;
typedef R::Point_2                                     Point_2;
typedef R::Line_2                                      Line_2;
typedef Polygon_traits_2< R >                          P_traits;
typedef vector< Point_2 >                              Container;
typedef CGAL::Polygon_2< P_traits, Container >         Polygon_2;
typedef Creator_uniform_2< double, Point_2 >           Creator;
typedef Random_points_in_square_2< Point_2, Creator >  Point_generator;

struct Minimum_rectangle_2 :
public geowin_update< CGALPointlist, list< leda_polygon > > {
  void
  update(const CGALPointlist& ipts, list<leda_polygon>& poly)
  {
    poly.clear();

    // (possibly) convert to our number type
    Polygon_2 pts;
#ifdef _MSC_VER
    {
#endif
      typedef CGALPointlist::const_iterator Listconstiter;
      for (Listconstiter i = ipts.begin(); i != ipts.end(); ++i)
        pts.push_back(Point_2(i->x(), i->y()));
#ifdef _MSC_VER
    }
#endif

    // building convex polygon ...
    Polygon_2 p;
    convex_hull_points_2(
      pts.vertices_begin(), pts.vertices_end(), back_inserter(p));
    cout << "convex hull area " << CGAL::to_double(p.area()) << endl;
    CGAL_assertion(p.is_simple());
    CGAL_assertion(p.is_convex());

    Polygon_2 kg;
    Timer t;
    t.start();
    min_rectangle_2(
      p.vertices_begin(), p.vertices_end(), back_inserter(kg));
    t.stop();
    cout << "min_rectangle area "
         << CGAL::to_double(kg.area())
         << " [time: " << t.time() << " sec.]"
         << endl;
    {
      CGAL_assertion(kg.is_simple());
      CGAL_assertion(kg.is_convex());
      for (Polygon_2::iterator ii = p.vertices_begin();
           ii != p.vertices_end();
           ++ii)
        CGAL_assertion(!kg.has_on_unbounded_side(*ii));
      CGAL_assertion(kg.area() >= p.area());
      cout << "convex hull area 2 " << CGAL::to_double(p.area()) << endl;
    }
    leda_list< leda_point > HLP;
    for (int i = 0; i < kg.size(); ++i)
      HLP.append(convert_to_leda(kg.container()[i]));
    leda_polygon back(HLP);
    poly.push_back(back);
  } // update(pts, poly)
};
struct Minimum_parallelogram_2 :
public geowin_update< CGALPointlist, list< leda_polygon > > {
  void
  update(const CGALPointlist& ipts, list<leda_polygon>& poly)
  {
    poly.clear();

    // (possibly) convert to our number type
    Polygon_2 pts;
#ifdef _MSC_VER
    {
#endif
      typedef CGALPointlist::const_iterator Listconstiter;
      for (Listconstiter i = ipts.begin(); i != ipts.end(); ++i)
        pts.push_back(Point_2(i->x(), i->y()));
#ifdef _MSC_VER
    }
#endif

    // building convex polygon ...
    Polygon_2 p;
    convex_hull_points_2(
      pts.vertices_begin(), pts.vertices_end(), back_inserter(p));
    cout << "convex hull area " << CGAL::to_double(p.area()) << endl;
    CGAL_assertion(p.is_simple());
    CGAL_assertion(p.is_convex());

    Polygon_2 kg;
    Timer t;
    t.start();
    min_parallelogram_2(
      p.vertices_begin(), p.vertices_end(), back_inserter(kg));
    t.stop();
    cout << "min_parallelogram area "
         << CGAL::to_double(kg.area())
         << " [time: " << t.time() << " sec.]"
         << endl;
    {
      CGAL_assertion(kg.is_simple());
      CGAL_assertion(kg.is_convex());
      for (Polygon_2::iterator ii = p.vertices_begin();
           ii != p.vertices_end();
           ++ii)
        CGAL_assertion(!kg.has_on_unbounded_side(*ii));
      CGAL_assertion(kg.area() >= p.area());
      cout << "convex hull area 2 " << CGAL::to_double(p.area()) << endl;
    }
    leda_list< leda_point > HLP;
    for (int i = 0; i < kg.size(); ++i)
      HLP.append(convert_to_leda(kg.container()[i]));
    leda_polygon back(HLP);
    poly.push_back(back);
  } // update(pts, poly)
};
struct Minimum_strip_2 :
public geowin_update< CGALPointlist, list< leda_line > > {
  void
  update(const CGALPointlist& ipts, list<leda_line>& lines)
  {
    lines.clear();

    // (possibly) convert to our number type
    Polygon_2 pts;
#ifdef _MSC_VER
    {
#endif
      typedef CGALPointlist::const_iterator Listconstiter;
      for (Listconstiter i = ipts.begin(); i != ipts.end(); ++i)
        pts.push_back(Point_2(i->x(), i->y()));
#ifdef _MSC_VER
    }
#endif

    // building convex polygon ...
    Polygon_2 p;
    convex_hull_points_2(
      pts.vertices_begin(), pts.vertices_end(), back_inserter(p));

    typedef std::vector< Line_2 > Linelist;
    Linelist ll;
    Timer t;
    t.start();
    min_strip_2(
      p.vertices_begin(), p.vertices_end(), back_inserter(ll));
    t.stop();
    cout << "min_strip width^2 ";
    if (ll.size() > 1)
      cout << CGAL::to_double(squared_distance(ll[0], ll[1]));
    else
      cout << "undef";
    cout << " [time: " << t.time() << " sec.]" << endl;
    for (Linelist::iterator i = ll.begin(); i != ll.end(); ++i)
      lines.push_back(
        leda_line(
          leda_point(CGAL::to_double(i->point(0).x()),
                     CGAL::to_double(i->point(0).y())),
          leda_point(CGAL::to_double(i->point(1).x()),
                     CGAL::to_double(i->point(1).y()))));
  } // update(pts, lines)
};

int main()
{
  geowin_init_default_type((CGALPointlist*)0,
    leda_string("CGALPointList"));

  CGALPointlist L;

  GeoWin GW("CGAL - Minimum Enclosing Quadrilaterals");

  // build a new edit scene
  geo_scene my_scene = GW.new_scene(L);
  GW.set_point_style(my_scene, leda_disc_point);

  // add min strip
  Minimum_strip_2 minstrip;
  geo_scene strip_scene =
    GW.new_scene(minstrip,
                 my_scene,
                 leda_string("Minimum Width Strip"));
  GW.set_color(strip_scene, leda_green);
  GW.set_line_width(strip_scene, 3);
  GW.set_fill_color(strip_scene, leda_green2);
  GW.set_visible(strip_scene, true);

  // add min rectangle
  Minimum_rectangle_2 minrect;
  geo_scene rect_scene =
    GW.new_scene(minrect,
                 my_scene,
                 leda_string("Minimum Area Rectangle"));
  GW.set_color(rect_scene, leda_red);
  GW.set_line_width(rect_scene, 3);
  GW.set_fill_color(rect_scene, leda_orange);
  GW.set_visible(rect_scene, true);

  // add min parallelogram
  Minimum_parallelogram_2 minpara;
  geo_scene para_scene =
    GW.new_scene(minpara,
                 my_scene,
                 leda_string("Minimum Area Parallelogram"));
  GW.set_line_width(para_scene, 3);
  GW.set_color(para_scene, leda_blue);
  GW.set_fill_color(para_scene, leda_cyan);
  GW.set_visible(para_scene, true);

  GW.edit(my_scene);
  return 0;
}

#else // ! (CGAL_USE_LEDA && __LEDA__ >= 400)

#include <iostream>

int main()
{
  std::cerr << "LEDA >= 4.0 is required to run this demo."
            << std::endl;
  return 0;
}

#endif // (CGAL_USE_LEDA && __LEDA__ >= 400)
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

