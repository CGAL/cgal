// ======================================================================
//
// Copyright (c) The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-24 $
// release_date  : $CGAL_Date: 2000/12/29 $
//
// file          : include/CGAL/Snap_rounding_2.h
// package       : arr (1.73)
// maintainer    : Eli Packer <elip@post.tau.ac.il>
// author(s)     : Eli Packer
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_SR_2_H
#define CGAL_SR_2_H

#include <CGAL/leda_rational.h> 

#include <iostream>

#ifndef CGAL_ENUM_H
#include <CGAL/enum.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Random.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/intersection_2.h>

#endif

/*#ifndef SWEEP_TO_PRODUCE_PLANAR_MAP_SUBCURVES_H
#include <CGAL/sweep_to_produce_planar_map_subcurves.h>
#endif*/

//#include <CGAL/Sweep_line_tight_2.h>
#include <CGAL/Sweep_line_2.h>

//#ifndef CGAL_ARR_SEGMENT_TRAITS_2_H
#include <CGAL/Arr_segment_traits_2.h>
//#endif

#ifndef CGAL_ARR_POLYLINE_TRAITS_H
#include <CGAL/Arr_polyline_traits.h>
#endif

#ifdef ISR_DEBUG
#include <CGAL/IO/leda_window.h>
#endif

#include <list>
#include <set>

#include <CGAL/leda_real.h>

#include "../../include/CGAL/Snap_rounding_kd_2.h"

CGAL_BEGIN_NAMESPACE

#if defined(ISR_DEBUG)
typedef CGAL::Window_stream Window_stream;
#endif

template<class Rep_>
class Segment_data {

typedef Rep_                                Rep;
typedef typename Rep::FT                    NT;
typedef CGAL::Segment_2<Rep>                Segment_2;
typedef CGAL::Point_2<Rep>                  Point_2;

private:
  NT x1;
  NT y1;
  NT x2;
  NT y2;

public:
  Segment_data();
  Segment_data(NT inp_x1,NT inp_y1,NT inp_x2,NT inp_y2);
  void debug() const;
  NT get_x1();
  NT get_y1();
  NT get_x2();
  NT get_y2();
  inline void set_data(NT inp_x1,NT inp_y1,NT inp_x2,NT inp_y2);
  void determine_direction();
  bool equal(Segment_2 s);
  Segment_data(const Segment_data &other);
  // !!!! need operator =
};

template<class Rep_>
class Hot_Pixel {

typedef Rep_ Rep;
typedef typename Rep::FT NT;
typedef CGAL::Segment_2<Rep> Segment_2;
typedef CGAL::Point_2<Rep> Point_2;

private:
  // (x,y) is the center of the hot pixel
  NT x;
  NT y;
  NT pixel_size;
  Segment_2 *right_seg;
  Segment_2 *left_seg;
  Segment_2 *top_seg;
  Segment_2 *bot_seg;

public:
  void debug() const;
  template<class Out>
  void draw(Out &o) const;
  Hot_Pixel(NT inp_x,NT inp_y,NT inp_pixel_size);
  ~Hot_Pixel();
  inline NT get_x() const;
  inline NT get_y() const;
  bool intersect_left(Segment_2 &seg) const;
  bool intersect_right(Segment_2 &seg) const;
  bool intersect_bot(Segment_2 &seg) const;
  bool intersect_top(Segment_2 &seg) const;
  bool intersect(Segment_data<Rep> &seg) const;
};

// a function for compare two hot pixels for the set of hot pixels
template<class Rep_>
struct hot_pixel_auclidian_cmp
{
  bool operator()(const Hot_Pixel<Rep_> *h1,const Hot_Pixel<Rep_> *h2) const;
};

// a function for compare two hot pixels for the set of hot pixels a
// certain segment intersect
template<class Rep_>
struct hot_pixel_dir_cmp
{
  bool operator ()(const Hot_Pixel<Rep_> *h1,const Hot_Pixel<Rep_> *h2);
};


template<class Rep_>
class Snap_rounding_2 {

typedef CGAL::Arr_segment_traits_2<Rep_ > Traits;
typedef Rep_ Rep;
typedef typename Rep::FT NT;
typedef typename Traits::X_curve X_curve;
typedef typename Traits::Curve Curve;
typedef std::list<X_curve>              CurveContainer;
typedef typename CurveContainer::iterator            CurveContainerIter;

public:
  typedef CGAL::Segment_2<Rep> Segment_2;
  typedef CGAL::Point_2<Rep> Point_2;
  typedef std::list<Point_2> PointList;
  typedef typename std::list<std::list<CGAL::Point_2<Rep> > >
             Polylines_container;
  typedef typename Polylines_container::const_iterator Polyline_const_iterator;
  typedef typename Polylines_container::iterator Polyline_iterator;
  typedef typename std::list<CGAL::Point_2<Rep> > Points_container;
  typedef typename Points_container::const_iterator Point_const_iterator;
  typedef typename std::list<Segment_2> Segments_container;
  typedef typename Segments_container::const_iterator Segment_const_iterator;
  typedef typename Segments_container::iterator Segment_iterator;

  enum Direction {UP_RIGHT,UP_LEFT,DOWN_RIGHT,DOWN_LEFT,UP,DOWN,LEFT,
                  RIGHT,POINT_SEG};

  static Direction seg_dir;
  static bool erase_hp;
  static inline Direction get_direction() {return(seg_dir);}
  static inline void set_direction(Direction dir) {seg_dir = dir;}
  static inline bool get_erase_hp() {return(erase_hp);}

  // ctor
  Snap_rounding_2(Segment_const_iterator begin,
                  Segment_const_iterator end,
                  NT inp_pixel_size,bool inp_do_isr = true,
                  int inp_number_of_kd_trees = default_number_of_kd_trees);
  // ctor
  Snap_rounding_2(NT inp_pixel_size,bool inp_do_isr = true,
                  int inp_number_of_kd_trees = default_number_of_kd_trees);
  // cctor
  Snap_rounding_2(const Snap_rounding_2 &other);

#ifdef ISR_DEBUG
  template<class Out>
  void output_distances(Out &o);
#endif
  // !!!! change names to output and input
  const Polyline_const_iterator polylines_begin();
  const Polyline_const_iterator polylines_end();

  inline Segment_const_iterator segments_begin() const {
             return(seg_2_list.begin());}
  inline Segment_const_iterator segments_end() const {
    return(seg_2_list.end());}
  inline Segment_iterator segments_begin() {return(seg_2_list.begin());}
  inline Segment_iterator segments_end() {return(seg_2_list.end());}

  bool insert(Segment_2 seg);
  bool push_back(Segment_2 seg);
  template < class InputIterator >
    int insert(InputIterator first, InputIterator last);
  bool remove(Segment_2 seg);
  void clear();
  bool change_number_of_kd_trees(int inp_number_of_kd_trees);
  bool change_pixel_size(NT inp_pixel_size);
  void do_isr(bool inp_do_isr);

  template<class Out>
  void output(Out &o);

  /*#ifdef ISR_DEBUG
  void window_output(Window_stream &w,bool wait_for_click);
  #endif*/

private:
  // the next variable is for lazy evaluation:
  // it determines whether an isr/sr work has
  // to be done (at the beginning, after insertion, etc) 
  bool need_sr;

  static const int default_number_of_kd_trees = 5;

  std::set<Hot_Pixel<Rep> *,hot_pixel_auclidian_cmp<Rep> > hp_set;
  NT pixel_size,min_x,max_x,min_y,max_y;
  Segments_container seg_2_list;
  std::list<Segment_data<Rep> > seg_list;
  std::list<std::list<Point_2> > segments_output_list;
  int number_of_segments,number_of_kd_trees;
  Multiple_kd_tree<NT,Hot_Pixel<Rep> *> *mul_kd_tree;
  bool wheteher_to_do_isr;

  void find_hot_pixels_and_create_kd_trees();
  void find_intersected_hot_pixels(Segment_data<Rep> &seg,
                         std::set<Hot_Pixel<Rep> *,
                         hot_pixel_dir_cmp<Rep> > &hot_pixels_intersected_set,
                         int &number_of_intersections);

  void debug(std::set<Hot_Pixel<Rep>,hot_pixel_dir_cmp<Rep> > &s);
  void reroute_sr(std::set<Hot_Pixel<Rep> *,hot_pixel_dir_cmp<Rep> >
                  &inp_hot_pixels_intersected_set,std::list<Point_2>
                  &seg_output);
  void reroute_isr(std::set<Hot_Pixel<Rep> *,hot_pixel_dir_cmp<Rep> >
                   &inp_hot_pixels_intersected_set,std::list<Point_2>
                   &seg_output,int number_of_intersections,bool first_time);
  void iterate();
};

#if defined(ISR_DEBUG) || defined(TEST)
#include <CGAL/squared_distance_2.h>

#endif

#ifdef ISR_DEBUG
int max_rec = 1,cur_rec = -1,cur_max,needed_hp = 0,unneeded_hp = 0;
#elif defined XXXX
int needed_hp = 0,unneeded_hp = 0;
#endif

// ctor
template<class Rep_>
Segment_data<Rep_>::Segment_data() {}
template<class Rep_>
Segment_data<Rep_>::Segment_data(NT inp_x1,NT inp_y1,NT inp_x2,NT inp_y2) :
                    x1(inp_x1),y1(inp_y1),x2(inp_x2),y2(inp_y2) {}

// cctor
template<class Rep_>
Segment_data<Rep_>::Segment_data(const Segment_data &other)
{
  x1 = other.x1;
  y1 = other.y1;
  x2 = other.x2;
  y2 = other.y2;
}

/*template<class Rep_>
void Segment_data<Rep_>::debug() const {std::cerr << "Segment (" << x1 << ","
<< y1 << "):(" << x2 << ":" << y2 << ")\n";}*/

template<class Rep_>
typename Rep_::FT Segment_data<Rep_>::get_x1() {return(x1);}

template<class Rep_>
typename Rep_::FT Segment_data<Rep_>::get_y1() {return(y1);}

template<class Rep_>
typename Rep_::FT Segment_data<Rep_>::get_x2() {return(x2);}

template<class Rep_>
typename Rep_::FT Segment_data<Rep_>::get_y2() {return(y2);}

template<class Rep_>
inline void Segment_data<Rep_>::set_data(NT inp_x1,NT inp_y1,NT inp_x2,
            NT inp_y2)
{
   x1 = inp_x1;
   y1 = inp_y1;
   x2 = inp_x2;
   y2 = inp_y2;
}

template<class Rep_>
bool Segment_data<Rep_>::equal(Segment_2 s)
{
  return(
	 s.source().x() == x1 &
	 s.source().y() == y1 &
	 s.target().x() == x2 &
	 s.target().y() == y2);
}

template<class Rep_>
void Segment_data<Rep_>::determine_direction()
{
  if(x1 < x2) {
    if(y1 < y2)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::UP_RIGHT);
    else if(y1 == y2)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::RIGHT);
    else
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::DOWN_RIGHT);
  } else if(x1 == x2) {
    if(y1 < y2)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::UP);
    else if(y1 == y2)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::POINT_SEG);
    else
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::DOWN);
  } else {
    if(y1 < y2)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::UP_LEFT);
    else if(y1 == y2)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::LEFT);
    else
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::DOWN_LEFT);
  }
}


template<class Rep_>
void Hot_Pixel<Rep_>::debug() const {std::cerr << "Hot Pixel (" << x << ":"
                      << y << ")\n";}

template<class Rep_>
template<class Out> void Hot_Pixel<Rep_>::draw(Out &o) const
  {
    o << *right_seg;
    o << *left_seg;
    o << *top_seg;
    o << *bot_seg;
  }

  // intersection pixel
template<class Rep_>
Hot_Pixel<Rep_>::Hot_Pixel(NT inp_x,NT inp_y,NT inp_pixel_size) :
                           pixel_size(inp_pixel_size)
  {
    x = NT(floor((inp_x / pixel_size).to_double())) * pixel_size +
        pixel_size / 2.0;

    y = NT(floor((inp_y / pixel_size).to_double())) * pixel_size +
        pixel_size / 2.0;

    right_seg = new Segment_2(Point_2(x + pixel_size / 2.0,y -
                              pixel_size / 2.0),
                              Point_2(x + pixel_size / 2.0,y +
                              pixel_size / 2.0));
    left_seg = new Segment_2(Point_2(x - pixel_size / 2.0,y -
                             pixel_size / 2.0),
                             Point_2(x - pixel_size / 2.0,y +
                             pixel_size / 2.0));
    top_seg = new Segment_2(Point_2(x - pixel_size / 2.0,y +
                            pixel_size / 2.0),
                            Point_2(x + pixel_size / 2.0,y +
                            pixel_size / 2.0));
    bot_seg = new Segment_2(Point_2(x - pixel_size / 2.0,y -
                            pixel_size / 2.0),
                            Point_2(x + pixel_size / 2.0,y -
                            pixel_size / 2.0));
  }

template<class Rep_>
Hot_Pixel<Rep_>::~Hot_Pixel()
  {
    if(Snap_rounding_2<Rep_>::get_erase_hp()) {
      delete(right_seg);
      delete(left_seg);
      delete(top_seg);
      delete(bot_seg);
    }
  }

template<class Rep_>
inline typename Rep_::FT Hot_Pixel<Rep_>::get_x() const {return(x);}

template<class Rep_>
inline typename Rep_::FT Hot_Pixel<Rep_>::get_y() const {return(y);}

template<class Rep_>
bool Hot_Pixel<Rep_>::intersect_left(Segment_2 &seg) const
  {
    CGAL::Object result;    
    Point_2 p;
    Segment_2 s;

    result = CGAL::intersection(seg,*left_seg);

    if(CGAL::assign(p,result)) {
      NT tmp = y + pixel_size / 2.0;
      return(p.y() != tmp || Snap_rounding_2<Rep_>::get_direction() ==
             Snap_rounding_2<Rep_>::UP_LEFT && seg.source().y() != tmp ||
             Snap_rounding_2<Rep_>::get_direction() ==
             Snap_rounding_2<Rep_>::DOWN_RIGHT && seg.target().y() != tmp);
    } else if(CGAL::assign(s,result))
      return(true);
    else
      return(false);
  }


template<class Rep_>
bool Hot_Pixel<Rep_>::intersect_right(Segment_2 &seg) const
  {
    CGAL::Object result;    
    Point_2 p;
    Segment_2 s;
    
    result = CGAL::intersection(seg,*right_seg);

    if(CGAL::assign(p,result)) {
      // bottom right point was checked in intersect_bot

      NT tmp = y + pixel_size / 2.0;
      if(p.y() == tmp)
        return(Snap_rounding_2<Rep_>::get_direction() ==
               Snap_rounding_2<Rep_>::UP_RIGHT && seg.source().y() != tmp ||
               Snap_rounding_2<Rep_>::get_direction() ==
               Snap_rounding_2<Rep_>::DOWN_LEFT && seg.target().y() != tmp);
      else if(p.y() == y - pixel_size / 2.0)
        return(false);// was checked
      else
        return((Snap_rounding_2<Rep_>::get_direction() ==
                Snap_rounding_2<Rep_>::LEFT ||
                Snap_rounding_2<Rep_>::get_direction() ==
                Snap_rounding_2<Rep_>::DOWN_LEFT ||
                Snap_rounding_2<Rep_>::get_direction() ==
                Snap_rounding_2<Rep_>::UP_LEFT) && seg.target().x() !=
                x + pixel_size / 2.0 ||
                (Snap_rounding_2<Rep_>::get_direction() ==
                Snap_rounding_2<Rep_>::RIGHT ||
                Snap_rounding_2<Rep_>::get_direction() ==
                Snap_rounding_2<Rep_>::DOWN_RIGHT ||
                Snap_rounding_2<Rep_>::get_direction() ==
                Snap_rounding_2<Rep_>::UP_RIGHT) && seg.source().x() !=
                x + pixel_size / 2.0);
    } else
      return(false);
  }

template<class Rep_>
bool Hot_Pixel<Rep_>::intersect_bot(Segment_2 &seg) const
  {
    CGAL::Object result;
    Point_2 p;
    Segment_2 s;

    result = CGAL::intersection(seg,*bot_seg);

    if(CGAL::assign(p,result)) {
      NT tmp = x + pixel_size / 2.0;
      return(p.x() != tmp || Snap_rounding_2<Rep_>::get_direction() ==
             Snap_rounding_2<Rep_>::UP_LEFT && seg.target().x() != tmp ||
             Snap_rounding_2<Rep_>::get_direction() ==
             Snap_rounding_2<Rep_>::DOWN_RIGHT && seg.source().x() != tmp);
    } else if(CGAL::assign(s,result))
      return(true);
    else
      return(false);
  }

template<class Rep_>
bool Hot_Pixel<Rep_>::intersect_top(Segment_2 &seg) const
  {
    CGAL::Object result;
    Point_2 p;
    Segment_2 s;
    
    result = CGAL::intersection(seg,*top_seg);

    if(CGAL::assign(p,result)) {
      // corner points was checked in intersect_bot
      NT tar_y = seg.target().y(),sou_y = seg.source().y();

      if(p.x() == x - pixel_size / 2.0 || p.x() == x + pixel_size / 2.0)
        return(false);// were checked
      else
        return((Snap_rounding_2<Rep_>::get_direction() ==
               Snap_rounding_2<Rep_>::DOWN ||
               Snap_rounding_2<Rep_>::get_direction() ==
               Snap_rounding_2<Rep_>::DOWN_LEFT ||
               Snap_rounding_2<Rep_>::get_direction() ==
               Snap_rounding_2<Rep_>::DOWN_RIGHT) && tar_y != y +
               pixel_size / 2.0 ||
               (Snap_rounding_2<Rep_>::get_direction() ==
               Snap_rounding_2<Rep_>::UP ||
               Snap_rounding_2<Rep_>::get_direction() ==
               Snap_rounding_2<Rep_>::UP_LEFT ||
               Snap_rounding_2<Rep_>::get_direction() ==
               Snap_rounding_2<Rep_>::UP_RIGHT) && sou_y != y +
               pixel_size / 2.0);
    } else
    return(false);
  }

template<class Rep_>
bool Hot_Pixel<Rep_>::intersect(Segment_data<Rep_> &seg) const
  {
    Segment_2 s(Point_2(seg.get_x1(),seg.get_y1()),Point_2(seg.get_x2(),
                seg.get_y2()));

    return(intersect_bot(s) || intersect_left(s) || intersect_right(s) ||
           intersect_top(s));
  }

// a function for compare two hot pixels for the set of hot pixels
template<class Rep_>
bool hot_pixel_auclidian_cmp<Rep_>::operator()(const Hot_Pixel<Rep_> *h1,
     const Hot_Pixel<Rep_> *h2) const
  {
    return(h1->get_x() < h2->get_x() ||
         h1->get_x() == h2->get_x() && h1->get_y() < h2->get_y());
  }

// a function for compare two hot pixels for the set of hot pixels a certain
// segment intersect
template<class Rep_>
bool hot_pixel_dir_cmp<Rep_>::operator ()(const Hot_Pixel<Rep_> *h1,\
     const Hot_Pixel<Rep_> *h2) 
{
  return(
     // Point segment intersects only one pixel, thus ignored
    Snap_rounding_2<Rep_>::get_direction() ==
    Snap_rounding_2<Rep_>::UP_RIGHT &&
    (h1->get_x() < h2->get_x() || 
     h1->get_x() == h2->get_x() && h1->get_y() < h2->get_y()) ||
    Snap_rounding_2<Rep_>::get_direction() ==
    Snap_rounding_2<Rep_>::UP_LEFT &&
    (h1->get_x() > h2->get_x() || 
     h1->get_x() == h2->get_x() && h1->get_y() < h2->get_y()) ||
    Snap_rounding_2<Rep_>::get_direction() ==
    Snap_rounding_2<Rep_>::DOWN_RIGHT &&
    (h1->get_x() < h2->get_x() || 
     h1->get_x() == h2->get_x() && h1->get_y() > h2->get_y()) ||
    Snap_rounding_2<Rep_>::get_direction() ==
    Snap_rounding_2<Rep_>::DOWN_LEFT &&
    (h1->get_x() > h2->get_x() || 
     h1->get_x() == h2->get_x() && h1->get_y() > h2->get_y()) ||
    Snap_rounding_2<Rep_>::get_direction() == Snap_rounding_2<Rep_>::UP &&
    h1->get_y() < h2->get_y() ||
    Snap_rounding_2<Rep_>::get_direction() == Snap_rounding_2<Rep_>::DOWN &&
    h1->get_y() > h2->get_y() ||
    Snap_rounding_2<Rep_>::get_direction() == Snap_rounding_2<Rep_>::LEFT &&
    h1->get_x() > h2->get_x() ||
    Snap_rounding_2<Rep_>::get_direction() ==
    Snap_rounding_2<Rep_>::RIGHT &&
    h1->get_x() < h2->get_x());
}



template<class Rep_>
void Snap_rounding_2<Rep_>::find_hot_pixels_and_create_kd_trees()
  {
    Hot_Pixel<Rep_> *hp;
    typename std::list<Segment_data<Rep_> >::iterator iter1;
    CGAL::Object result;
    Point_2 p;
    std::list<std::pair<std::pair<NT,NT>,Hot_Pixel<Rep_> *> > hot_pixels_list;

    list<X_curve> segments;
    for(iter1 = seg_list.begin();iter1 != seg_list.end();++iter1)
      segments.push_back(X_curve(Point_2(iter1->get_x1(),iter1->get_y1()),
                                 Point_2(iter1->get_x2(),iter1->get_y2())));

    //    PM pm(new CGAL::Pm_naive_point_location<PM>);
    // sweep_to_construct_planar_map(segments.begin(), segments.end(), pm);
    std::list<X_curve>  subcurves;

    /*    sweep_to_produce_planar_map_subcurves(segments.begin(), 
					  segments.end(),  
					  traits, 
					  subcurves);*/

    /*// get subcurves with overlapping
//    CGAL::Sweep_line_tight_2<CurveContainerIter, Traits, Event, SubCurve> sl;
    CGAL::Sweep_line_2<CurveContainerIter, Traits> sl;
    sl.get_subcurves(segments.begin(), segments.end(),
    std::back_inserter(subcurves));*/

    // get intersection points (with endpoints)
    PointList mypointlist;
    //    CGAL::Sweep_line_tight_2<CurveContainerIter, Traits,
    //                             Event, SubCurve> sl;
    CGAL::Sweep_line_2<CurveContainerIter, Traits> sl;
    sl.get_intersection_points(segments.begin(), segments.end(),
                             std::back_inserter(mypointlist));

    for(typename std::list<Point_2>::const_iterator
            v_iter = mypointlist.begin();
	v_iter != mypointlist.end();++v_iter) {
      hp = new Hot_Pixel<Rep_>(v_iter->x(),v_iter->y(),pixel_size);
      hot_pixels_list.push_back(std::pair<std::pair<NT,NT>,Hot_Pixel<Rep_> *>(
            std::pair<NT,NT>(hp->get_x(),hp->get_y()),hp));
    }

    /*    for(list<X_curve>::iterator v_iter = subcurves.begin();
        v_iter != subcurves.end();
        ++v_iter) {
      hp = new Hot_Pixel<Rep_>(v_iter->source().x(),
                               v_iter->source().y(),
                               pixel_size);
      if(hp_set.insert(hp).second)
        hot_pixels_list.push_back(pair<pair<NT,NT>,Hot_Pixel<Rep_> *>(
            pair<NT,NT>(hp->get_x(),hp->get_y()),hp));
      hp = new Hot_Pixel<Rep_>(v_iter->target().x(),
                               v_iter->target().y(),
                               pixel_size);
      if(hp_set.insert(hp).second)
        hot_pixels_list.push_back(pair<pair<NT,NT>,Hot_Pixel<Rep_> *>(
            pair<NT,NT>(hp->get_x(),hp->get_y()),hp));
	    }*/

    // create kd multiple tree
    // create simple_list from seg_list
    std::list<std::pair<std::pair<NT,NT>,std::pair<NT,NT> > > simple_seg_list;
    for(typename std::list<Segment_data<Rep_> >::iterator iter =
        seg_list.begin();iter != seg_list.end();++iter) {
      std::pair<NT,NT> first(iter->get_x1(),iter->get_y1()),
   	          second(iter->get_x2(),iter->get_y2());
      simple_seg_list.push_back(std::pair<std::pair<NT,NT>,std::pair<NT,NT> >(
                                first,second));
    }

    mul_kd_tree = new Multiple_kd_tree<NT,Hot_Pixel<Rep_> *>(hot_pixels_list,
                  number_of_kd_trees,simple_seg_list);
  }

template<class Rep_>
void Snap_rounding_2<Rep_>::find_intersected_hot_pixels(Segment_data<Rep_>
                    &seg,
                    std::set<Hot_Pixel<Rep_> *,
                    hot_pixel_dir_cmp<Rep_> > &hot_pixels_intersected_set,
                    int &number_of_intersections)
  {
    typename std::list<Hot_Pixel<Rep_> *>::iterator iter;

    hot_pixels_intersected_set.clear();
    seg.determine_direction();
    number_of_intersections = 0;

    std::list<Hot_Pixel<Rep_> *> hot_pixels_list;
    mul_kd_tree->get_intersecting_points(hot_pixels_list,
	   Segment_2(Point_2(seg.get_x1(),seg.get_y1()),
           Point_2(seg.get_x2(),seg.get_y2())),pixel_size);

    for(iter = hot_pixels_list.begin();iter != hot_pixels_list.end();++iter) {
      if((*iter)->intersect(seg)) {

#if defined ISR_DEBUG
        ++needed_hp;
#endif
        hot_pixels_intersected_set.insert(*iter);
      }
#if defined ISR_DEBUG
        else
          ++unneeded_hp;
#endif
    }

    number_of_intersections = hot_pixels_intersected_set.size();
  }

template<class Rep_>
void debug(std::set<Hot_Pixel<Rep_>,hot_pixel_dir_cmp<Rep_> > &s)
  {
    std::cerr << "  Debugging inp_hot_pixels_intersected_set\n";
    for(typename std::set<Hot_Pixel<Rep_>,hot_pixel_dir_cmp<Rep_> >::
        iterator iter = s.begin();iter != s.end();++iter) {
      std::cerr << "    ";
      iter->debug();
    }
    std::cerr << "  Finish Debugging inp_hot_pixels_intersected_set\n";
  }

template<class Rep_>
void Snap_rounding_2<Rep_>::reroute_sr(std::set<Hot_Pixel<Rep_> *,
     hot_pixel_dir_cmp<Rep_> > &inp_hot_pixels_intersected_set,
     std::list<Point_2> &seg_output)
  {
    typename std::set<Hot_Pixel<Rep_> *,
    hot_pixel_dir_cmp<Rep_> >::iterator hot_pixel_iter =
    inp_hot_pixels_intersected_set.begin();
    ++hot_pixel_iter;

    while(hot_pixel_iter != inp_hot_pixels_intersected_set.end()) {
      seg_output.push_back(Point_2((*hot_pixel_iter)->get_x(),
                           (*hot_pixel_iter)->get_y()));
      ++hot_pixel_iter;
    }

  }

template<class Rep_>
void Snap_rounding_2<Rep_>::reroute_isr(std::set<Hot_Pixel<Rep_> *,
   hot_pixel_dir_cmp<Rep_> > &inp_hot_pixels_intersected_set,
   std::list<Point_2> &seg_output,int number_of_intersections,bool first_time)
  {
    typename std::set<Hot_Pixel<Rep_> *,hot_pixel_dir_cmp<Rep_> >::
      iterator hot_pixel_iter,next_hot_pixel_iter,before_last_hot_pixel_iter;
    Segment_data<Rep_> seg;
    std::set<Hot_Pixel<Rep_> *,hot_pixel_dir_cmp<Rep_> >
      hot_pixels_intersected_set;

    if(number_of_intersections > 2 || first_time) {
      before_last_hot_pixel_iter = inp_hot_pixels_intersected_set.end();
      --before_last_hot_pixel_iter;

      for(hot_pixel_iter = inp_hot_pixels_intersected_set.begin();
          hot_pixel_iter != before_last_hot_pixel_iter;++hot_pixel_iter) {
        next_hot_pixel_iter = hot_pixel_iter;
        ++next_hot_pixel_iter;
        seg.set_data((*hot_pixel_iter)->get_x(),(*hot_pixel_iter)->get_y(),
            (*next_hot_pixel_iter)->get_x(),(*next_hot_pixel_iter)->get_y());
        seg.determine_direction();
        find_intersected_hot_pixels(seg,hot_pixels_intersected_set,
            number_of_intersections);
        reroute_isr(hot_pixels_intersected_set,seg_output,
            number_of_intersections,false);
      }
    } else {
      // insert second hot pixel
      hot_pixel_iter = inp_hot_pixels_intersected_set.begin();
      ++hot_pixel_iter;
      seg_output.push_back(Point_2((*hot_pixel_iter)->get_x(),
          (*hot_pixel_iter)->get_y()));
    }
  }


template<class Rep_>
void Snap_rounding_2<Rep_>::iterate()
  {
    std::list<Point_2> seg_output;
    std::set<Hot_Pixel<Rep_> *,hot_pixel_dir_cmp<Rep_> >
      hot_pixels_intersected_set;
    typename std::set<Hot_Pixel<Rep_> *,hot_pixel_dir_cmp<Rep_> >::
      iterator hot_pixel_iter;
    int number_of_intersections;
    Hot_Pixel<Rep_> *hp;

    segments_output_list.clear();
    for(typename std::list<Segment_data<Rep_> >::iterator iter =
        seg_list.begin();iter != seg_list.end();++iter) {
      seg_output.clear();
      iter->determine_direction();
      find_intersected_hot_pixels(*iter,hot_pixels_intersected_set,
        number_of_intersections);

      // hot_pixels_intersected_set must have at least two hot pixels when the
      // segment is not in entirely inside a hot pixel enter first hot pixel
      hot_pixel_iter = hot_pixels_intersected_set.begin();
      if(hot_pixel_iter == hot_pixels_intersected_set.end()) {
        // segment entirely inside a pixel
        hp = new Hot_Pixel<Rep_>(iter->get_x1(),iter->get_y1(),pixel_size);
        seg_output.push_back(Point_2(hp->get_x(),hp->get_y()));
        erase_hp = true;
        delete(hp);
        erase_hp = false;
      } else {
        seg_output.push_back(Point_2((*hot_pixel_iter)->get_x(),
                             (*hot_pixel_iter)->get_y()));
        if(number_of_intersections > 1) {
          // segments that have at most one intersecting hot pixel are
          // done(it was inserted)
          if(wheteher_to_do_isr)
            reroute_isr(hot_pixels_intersected_set,seg_output,
                        number_of_intersections,true);
          else
            reroute_sr(hot_pixels_intersected_set,seg_output);
	}
      }

      segments_output_list.push_back(seg_output);
    }
  }

template<class Rep_>
Snap_rounding_2<Rep_>::Snap_rounding_2(Segment_const_iterator
  begin,Segment_const_iterator end,
  NT inp_pixel_size,bool inp_do_isr,int inp_number_of_kd_trees)
  {
    // initialize approximation angles map    
    erase_hp = false;
    wheteher_to_do_isr = inp_do_isr;
    pixel_size = inp_pixel_size;
    number_of_segments = 0;
    number_of_kd_trees = inp_number_of_kd_trees;
    need_sr = true;
    // copy segments list
    while(begin != end) {
      seg_list.push_back(Segment_data<Rep_>(begin->source().x(),
                         begin->source().y(),begin->target().x(),
                         begin->target().y()));
      seg_2_list.push_back(*begin);
      ++number_of_segments;
      ++begin;
    }
  }

// cctor
template<class Rep_>
Snap_rounding_2<Rep_>::Snap_rounding_2(const Snap_rounding_2 &other)
  {
    erase_hp = false;
    wheteher_to_do_isr = other.wheteher_to_do_isr;
    pixel_size = other.pixel_size;
    number_of_segments = other.number_of_segments;
    number_of_kd_trees = other.number_of_kd_trees;
    need_sr = true;
    seg_list = other.seg_list;
    seg_2_list = other.seg_2_list;
  }

template<class Rep_>
Snap_rounding_2<Rep_>::Snap_rounding_2(
    NT inp_pixel_size,bool inp_do_isr,int inp_number_of_kd_trees)
  {
    // initialize approximation angles map
    need_sr = true;
    erase_hp = false;
    wheteher_to_do_isr = inp_do_isr;
    pixel_size = inp_pixel_size;
    number_of_segments = 0;
    number_of_kd_trees = inp_number_of_kd_trees;
  }

template<class Rep_>
bool Snap_rounding_2<Rep_>::insert(Segment_2 seg)
  {
    need_sr = true;
    seg_list.push_back(Segment_data<Rep_>(
                         seg.source().x(),
                         seg.source().y(),
                         seg.target().x(),
                         seg.target().y()));

    seg_2_list.push_back(seg);
    ++number_of_segments;

    return(true);
  }

template<class Rep_>
bool Snap_rounding_2<Rep_>::push_back(Segment_2 seg)
  {
    return(insert(seg));
  }

template < class Rep_ >
template < class InputIterator >
int
Snap_rounding_2<Rep_>::insert(InputIterator first, InputIterator last)
  {
    need_sr = true;
    int n = 0;
    while(first != last){
      if(insert(*first)){
	n++;
      }
      ++first;
    }
    return n;
  }


template<class Rep_>
const Snap_rounding_2<Rep_>::Polyline_const_iterator
      Snap_rounding_2<Rep_>::polylines_begin()
{
  if(need_sr) {
    need_sr = false;
    find_hot_pixels_and_create_kd_trees();
    iterate();
  }    

  return(segments_output_list.begin());
}

template<class Rep_>
const Snap_rounding_2<Rep_>::Polyline_const_iterator
      Snap_rounding_2<Rep_>::polylines_end()
{
  if(need_sr) {
    need_sr = false;
    find_hot_pixels_and_create_kd_trees();
    iterate();
  }    

  return(segments_output_list.end());
}

template<class Rep_>
bool Snap_rounding_2<Rep_>::remove(Segment_2 seg)
  {
    need_sr = true;
    bool found = false;
    Segment_data<Rep> s;

    for(std::list<Segment_data<Rep> >::iterator i1 = seg_list.begin();
      i1 != seg_list.end();++i1) {
       s = *i1;  
      if(s.equal(seg)) {
        found = true;
        seg_list.erase(i1);
        --number_of_segments;
        break;
      }
    }

    if(found) {
      for(Segment_iterator i2 = seg_2_list.begin();
        i2 != seg_2_list.end();++i2) {
        if(seg == *i2) {
          seg_2_list.erase(i2);
          break;
        }
      }
    }

    return(found);
  }

template<class Rep_>
void Snap_rounding_2<Rep_>::clear()
  { 
    need_sr = true;
    seg_list.clear();
    seg_2_list.clear();
  }

template<class Rep_>
bool Snap_rounding_2<Rep_>::change_number_of_kd_trees(
          int inp_number_of_kd_trees)
  {
    if(inp_number_of_kd_trees > 0) {
      number_of_kd_trees = inp_number_of_kd_trees;
      return(true);
    } else
      return(false);
  }

template<class Rep_>
bool Snap_rounding_2<Rep_>::change_pixel_size(NT inp_pixel_size)
  {
    if(inp_pixel_size > 0) {
      pixel_size = inp_pixel_size;
      return(true);
    } else
      return(false);
  }

template<class Rep_>
void Snap_rounding_2<Rep_>::do_isr(bool inp_do_isr)
  { 
    wheteher_to_do_isr = inp_do_isr;
  }

template<class Rep_>
template<class Out> void Snap_rounding_2<Rep_>::output(Out &o)
  {
    o << number_of_segments << std::endl;
    for(typename std::list<std::list<Point_2> >::iterator iter1 =
        segments_output_list.begin();iter1 != segments_output_list.end();
        ++iter1) {
      for(typename std::list<Point_2>::iterator iter2 = iter1->begin();
          iter2 != iter1->end();++iter2)
        o << iter2->x().to_double() << " " << iter2->y().to_double() << " ";

      o << std::endl;
    }
  }

/*#ifdef ISR_DEBUG
template<class Rep_>
template<class Out> void Snap_rounding_2<Rep_>::output_distances(Out &o)
  {
    double max = 0,max_seg_dis,evarage_dis,cur_dis;
    typename std::list<Segment_data<Rep_> >::iterator
            orig_iter = seg_list.begin();

    for(typename std::list<std::list<std::pair<NT,NT> > >::iterator iter1 =
        segments_output_list.begin();iter1 != segments_output_list.end();
        ++iter1) {
      max_seg_dis = 0;
      for(typename std::list<std::pair<NT,NT> >::iterator
            iter2 = iter1->begin();
          iter2 != iter1->end();++iter2) {
        cur_dis = sqrt(CGAL::squared_distance(Point_2(iter2->first,
          iter2->second),Segment_2(Point_2(orig_iter->get_x1(),
          orig_iter->get_y1()),Point_2(orig_iter->get_x2(),
          orig_iter->get_y2()))).to_double());
        if(cur_dis > max_seg_dis)
          max_seg_dis = cur_dis;
      }

      evarage_dis += max_seg_dis;
      if(max_seg_dis > max)
        max = max_seg_dis;

      ++orig_iter;
    }
    evarage_dis /= seg_list.size();

    std::cerr << "max distance between output and original is " << max << endl;
    std::cerr << "evarage distance between output and original is " <<
                 evarage_dis << endl;
  }


template<class Rep_>
void Snap_rounding_2<Rep_>::window_output(Window_stream &w,bool wait_for_click)
  {
    double x,y;
    bool seg_painted;

    w << CGAL::BLACK;

    for(typename std::set<Hot_Pixel<Rep_> *,hot_pixel_auclidian_cmp<Rep_> >::
        iterator iter = hp_set.begin();
        iter != hp_set.end();++iter)
      (*iter)->draw(w);

    // draw original segments
    w << CGAL::BLACK;
    for(typename std::list<Segment_data<Rep_> >::iterator iter =
        seg_list.begin();iter != seg_list.end();++iter) {
      if(iter->get_x1() == iter->get_x2() && iter->get_y1() == iter->get_y2())
        w << Point_2(iter->get_x1(),iter->get_y1());
      else
        w << Segment_2(Point_2(iter->get_x1(),iter->get_y1()),
                       Point_2(iter->get_x2(),iter->get_y2()));
    }

    // draw isr polylines
    w << CGAL::RED;
    typename std::list<Point_2>::iterator iter2,iter3;

    for(typename std::list<std::list<Point_2> >::iterator iter1 =
        segments_output_list.begin();iter1 != segments_output_list.end();
        ++iter1) {
      if(wait_for_click)
        w.read_mouse(x,y);
      iter2 = iter3 = iter1->begin();
      seg_painted = false;
      for(++iter2;iter2 != iter1->end();++iter2) {
        seg_painted = true;
        w << Segment_2(*iter2,*iter3);
        ++iter3;
      }

      if(!seg_painted) { // segment entirely inside hot pixel
        --iter2;
        w << *iter2;
      }
    }

    int mouse_input;
    while(true) {
      mouse_input = w.read_mouse(x,y);
      if(mouse_input == 1)
        return;
    }
  }
#endif
*/

template<class Rep>
typename Snap_rounding_2<Rep>::Direction Snap_rounding_2<Rep>::seg_dir;

template<class Rep>
bool Snap_rounding_2<Rep>::erase_hp = false;

CGAL_END_NAMESPACE

#endif // CGAL_ISR_2_H
