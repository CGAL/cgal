// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 13
//
// file          : include/CGAL/Arr_segment_exact_cached_traits.h
// package       : arr (1.03)
// author(s)     : Iddo Hanniel
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_ARR_SEGMENT_EXACT_CACHED_TRAITS
#define CGAL_ARR_SEGMENT_EXACT_CACHED_TRAITS

#include <CGAL/Arr_segment_exact_traits.h>

CGAL_BEGIN_NAMESPACE

template <class R>
class Arr_segment_exact_cached_traits : public Arr_segment_exact_traits<R>
{
public:
  typedef Arr_segment_exact_traits<R> Base;
  typedef typename Base::Point Point;
  typedef typename Base::Curve_point_status Curve_point_status;

  typedef typename Base::X_curve Segment;

  struct X_curve : public Segment {
    X_curve() : Segment(), os() {}
    X_curve(const Segment& cv) : Segment(cv), os(cv) {}
    X_curve(const X_curve& cv) : Segment(cv), os(cv.os) {}
    X_curve(const Point& s, const Point& t) : Segment(s,t), os(s,t) {}

    X_curve& operator=(const X_curve& cv) {
      Segment::operator=(cv);
      os=cv.os;
      return *this;
    }

    void set_o_seg(const Segment& s) {os=s;}
    const Segment& o_seg() const {return os;}

  private:
    typename Base::X_curve os; //original segment
  };

  typedef X_curve Curve;

  Arr_segment_exact_cached_traits() : Base() 
  { }
  
  void make_x_monotone(const Curve& cv, std::list<X_curve>& l) {}

  X_curve curve_flip(const X_curve& cv) const {
    X_curve c(cv.target(),cv.source());
    c.set_o_seg(cv.o_seg());
    return c;
  }

  void curve_split(const X_curve& cv, X_curve& c1, X_curve& c2, 
                   const Point& split_pt)
  {
    Base::curve_split(cv,c1,c2,split_pt);
    c1.set_o_seg(cv);
    c2.set_o_seg(cv);
  }

  bool nearest_intersection_to_right(const X_curve& c1,
                                      const X_curve& c2,
                                      const Point& pt,
                                     Point& p1,
                                     Point& p2) const
  {
    if (!do_intersect_to_right(c1,c2,pt)) return false;

    Object res;
    Segment seg;
    res=intersection(c1.o_seg(),c2.o_seg());

    if (assign(seg,res)) {
      //p1, p2 will always be ordered left,right (make seg left to right)
      if (compare_lexicographically_xy(curve_source(seg), 
				       curve_target(seg))==LARGER)
          seg=Base::curve_flip(seg);

      if (compare_lexicographically_xy(curve_target(seg),pt)==LARGER) {
        p2=curve_target(seg);
        if (compare_lexicographically_xy(curve_source(seg),pt)==LARGER)
          p1=curve_source(seg);
        else
          p1=pt;

        return true;
      }
      else {
        return false;
      }
    }
    
    if (assign(p1,res)) {
      if (compare_lexicographically_xy(p1,pt)==LARGER) {
        p2=p1;

        return true;
      }

      return false;
    }

    return false;


    //normalization for speed up - when using leda_rational
    //#ifdef CGAL_LEDA_RATIONAL_H  //normalize if we are with rational numbers
    //return Point(ip.x().normalize(),ip.y().normalize());
    //#else
    //return ip;
    //#endif

  }





};


CGAL_END_NAMESPACE

#endif //CGAL_ARR_SEGMENT_EXACT_CACHED_TRAITS
