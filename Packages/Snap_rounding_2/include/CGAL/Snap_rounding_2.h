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

//#ifndef CGAL_ENUM_H
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/enum.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Random.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/intersection_2.h>

//#endif

//#include <CGAL/Sweep_line_tight_2.h>
#include <CGAL/Sweep_line_2.h>

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits.h>
#include <list>
#include <set>
#include <CGAL/leda_real.h>
#include <CGAL/Snap_rounding_kd_2.h>

#include <CGAL/utility.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>

CGAL_BEGIN_NAMESPACE

template<class Rep_>
class Segment_data {

typedef Rep_                                Rep;
typedef typename Rep::FT                    NT;
typedef CGAL::Segment_2<Rep>                Segment_2;
typedef CGAL::Point_2<Rep>                  Point_2;

private:
 Point_2 p;
 Point_2 q;

 Rep_   _gt;

public:
  Segment_data();
  Segment_data(Point_2 p_inp,Point_2 q_inp);
  Point_2 source() const {return(p);}
  Point_2 target() const {return(q);}
  NT get_x1() const;
  NT get_y1() const;
  NT get_x2() const;
  NT get_y2() const;
  inline void set_data(NT inp_x1,NT inp_y1,NT inp_x2,NT inp_y2);
  void determine_direction();
  bool equal(Segment_2 s);
  Segment_data(const Segment_data &other);
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

typedef CGAL::Arr_segment_traits_2<Rep_ >            Traits;// !!!! remove

typedef Rep_                                         Rep;
typedef typename Rep::FT                             NT;
typedef typename Traits::X_curve                     X_curve;
typedef typename Traits::Curve                       Curve;
typedef std::list<X_curve>                           CurveContainer;
typedef typename CurveContainer::iterator            CurveContainerIter;

public:
  friend class Segment_data<Rep>;
  friend class Hot_Pixel<Rep>;
  friend class hot_pixel_dir_cmp<Rep>;

  typedef CGAL::Segment_2<Rep> Segment_2;

  //  typedef typename Rep::Segment_2 Segment_2;!!!! change everywhere

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

  static bool erase_hp;
  static inline bool get_erase_hp() {return(erase_hp);}

  //! A constructor
  Snap_rounding_2(Segment_const_iterator begin,
                  Segment_const_iterator end,
                  NT inp_pixel_size,bool inp_do_isr = true,
                  int inp_number_of_kd_trees = default_number_of_kd_trees);

  //! A constructor
  Snap_rounding_2(NT inp_pixel_size,bool inp_do_isr = true,
                  int inp_number_of_kd_trees = default_number_of_kd_trees);

  //! A copy constructor
  Snap_rounding_2(const Snap_rounding_2& other);

  //! An operator =
  Snap_rounding_2& operator =(const Snap_rounding_2& other);

   //! Returns a constant iterator to the first input segment.
  inline Segment_const_iterator segments_begin() const {
             return(seg_2_list.begin());}

  //! Returns a constant iterator to the after-the-last of the input segments.
  inline Segment_const_iterator segments_end() const {
    return(seg_2_list.end());}

  //! Returns an iterator to the first input segment.
  inline Segment_iterator segments_begin() {return(seg_2_list.begin());}

  //! Returns an iterator to the after-the-last of the input segments.
  inline Segment_iterator segments_end() {return(seg_2_list.end());}

  //! Returns a constant iterator to the output of the first input segment.
  const Polyline_const_iterator polylines_begin();

  /*! Returns a constant iterator to the after-the-last of the output of
   *  the input segments.
   */
  const Polyline_const_iterator polylines_end();

  //! insert a segment
  bool insert(Segment_2 seg);

  //! The STL standard member function for insertion.
  bool push_back(Segment_2 seg);

  //! Insertion of an iterator range.
  template < class InputIterator >
    int insert(InputIterator first, InputIterator last);
  
  //! Remove a segment
  bool remove(Segment_2 seg);

  //! Remove all segments  .           .
  void clear();

  //! Change the number of kd-trees used for segment-hot pixel queries.
  bool change_number_of_kd_trees(int inp_number_of_kd_trees);

  //! Change the pixel size
  bool change_pixel_size(NT inp_pixel_size);

  /*! Determine whether to apply Iterated Snap Rounding (ISR)
      or Snap Rounding (SR).
   */
  void do_isr(bool inp_do_isr);

  template<class Out>
  void output(Out &o);

private:
  enum Direction {UP_RIGHT,UP_LEFT,DOWN_RIGHT,DOWN_LEFT,UP,DOWN,LEFT,
                  RIGHT,POINT_SEG};

  static Direction seg_dir;
  // the next variable is for lazy evaluation:
  // it determines whether an isr/sr work has
  // to be done (at the beginning, after insertion, etc) 
  bool need_sr;

  static const int default_number_of_kd_trees = 1;

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

  void reroute_sr(std::set<Hot_Pixel<Rep> *,hot_pixel_dir_cmp<Rep> >
                  &inp_hot_pixels_intersected_set,std::list<Point_2>
                  &seg_output);
  void reroute_isr(std::set<Hot_Pixel<Rep> *,hot_pixel_dir_cmp<Rep> >
                   &inp_hot_pixels_intersected_set,std::list<Point_2>
                   &seg_output,int number_of_intersections,bool first_time);
  void iterate();
  void copy(const Snap_rounding_2<Rep_>& other);
  static inline Direction get_direction() {return(seg_dir);}
  static inline void set_direction(Direction dir) {seg_dir = dir;}
};

// ctor
template<class Rep_>
Segment_data<Rep_>::Segment_data() {}
template<class Rep_>
Segment_data<Rep_>::Segment_data(Point_2 p_inp,Point_2 q_inp) :
                    p(p_inp), q(q_inp) {}

// cctor
template<class Rep_>
Segment_data<Rep_>::Segment_data(const Segment_data& other)
{
  p = other.p;
  q = other.q;
}

template<class Rep_>
typename Rep_::FT Segment_data<Rep_>::get_x1() const {return(p.x());}

template<class Rep_>
typename Rep_::FT Segment_data<Rep_>::get_y1() const {return(p.y());}

template<class Rep_>
typename Rep_::FT Segment_data<Rep_>::get_x2() const {return(q.x());}

template<class Rep_>
typename Rep_::FT Segment_data<Rep_>::get_y2() const {return(q.y());}

template<class Rep_>
inline void Segment_data<Rep_>::set_data(NT inp_x1,NT inp_y1,NT inp_x2,
            NT inp_y2)
{
  p = Point_2(inp_x1,inp_y1);
  q = Point_2(inp_x2,inp_y2);
}

template<class Rep_>
bool Segment_data<Rep_>::equal(Segment_2 s)
{
  return(s.source() == p && s.target() == q);
}

template<class Rep_>
void Segment_data<Rep_>::determine_direction()
{
  Comparison_result cx = _gt.compare_x_2_object()(p,q);
  Comparison_result cy = _gt.compare_y_2_object()(p,q);

  if(cx == SMALLER) {
   if(cy == SMALLER)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::UP_RIGHT);
    else if(cy == EQUAL)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::RIGHT);
    else
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::DOWN_RIGHT);
  } else if(cx == EQUAL) {
    if(cy == SMALLER)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::UP);
    else if(cy == EQUAL)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::POINT_SEG);
    else
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::DOWN);
  } else {
    if(cy == SMALLER)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::UP_LEFT);
    else if(cy == EQUAL)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::LEFT);
    else
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::DOWN_LEFT);
  }
}

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

    result = CGAL::intersection(seg,*left_seg);//!!!! change: create instance and call its intersection

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

    /*// get subcurves with overlapping ********** 
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
      if((*iter)->intersect(seg))
        hot_pixels_intersected_set.insert(*iter);
    }

    number_of_intersections = hot_pixels_intersected_set.size();
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
      seg_list.push_back(Segment_data<Rep_>(begin->source(),
                                            begin->target()));
      seg_2_list.push_back(*begin);
      ++number_of_segments;
      ++begin;
    }
  }

template<class Rep_>
void Snap_rounding_2<Rep_>::copy(const Snap_rounding_2<Rep_>& other)
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

// cctor
template<class Rep_>
Snap_rounding_2<Rep_>::Snap_rounding_2(const Snap_rounding_2<Rep_>& other)
  {
    copy(other);
  }

// operator =
template<class Rep_>
Snap_rounding_2<Rep_>&
Snap_rounding_2<Rep_>::operator =(const Snap_rounding_2<Rep_>& other)
  {
    copy(other);

    return(*this);
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
					  seg.source(),
                         seg.target()));

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
      need_sr = true;
      return(true);
    } else
      return(false);
  }

template<class Rep_>
bool Snap_rounding_2<Rep_>::change_pixel_size(NT inp_pixel_size)
  {
    if(inp_pixel_size > 0) {
      pixel_size = inp_pixel_size;
      need_sr = true;
      return(true);
    } else
      return(false);
  }

template<class Rep_>
void Snap_rounding_2<Rep_>::do_isr(bool inp_do_isr)
  { 
    wheteher_to_do_isr = inp_do_isr;
    need_sr = true;
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

template<class Rep>
typename Snap_rounding_2<Rep>::Direction Snap_rounding_2<Rep>::seg_dir;

template<class Rep>
bool Snap_rounding_2<Rep>::erase_hp = false;

CGAL_END_NAMESPACE

#endif // CGAL_ISR_2_H
