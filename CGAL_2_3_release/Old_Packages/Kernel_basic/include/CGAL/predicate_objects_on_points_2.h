// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 
//
// file          : include/CGAL/predicate_objects_on_points_2.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$ 
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_PREDICATE_OBJECTS_ON_POINTS_2_H
#define CGAL_PREDICATE_OBJECTS_ON_POINTS_2_H

#include <CGAL/user_classes.h>

CGAL_BEGIN_NAMESPACE

template <class Point>
class p_Left_of_line_2p
{
public:
   typedef bool    result_type;
   typedef  Arity_tag< 1 >   Arity;

        p_Left_of_line_2p(const Point& a, const Point& b)
         : p_a(a), p_b(b)
        {}

  bool  operator()(const Point& c) const
        { return left_turn( p_a, p_b, c ); }

private:
    Point p_a;
    Point p_b;
};

template <class Point>
class p_Right_of_line_2p
{
public:
   typedef bool    result_type;
   typedef  Arity_tag< 1 >   Arity;

        p_Right_of_line_2p(const Point& a, const Point& b)
         : p_a(a), p_b(b)
        {}

  bool  operator()(const Point& c) const
        { return right_turn( p_a, p_b, c ); }

private:
    Point p_a;
    Point p_b;
};

template <class Point>
class p_Left_of_line_2p_safer
{
public:
  typedef bool    result_type;
   typedef  Arity_tag< 1 >   Arity;

        p_Left_of_line_2p_safer(const Point& a, const Point& b)
         : p_a(a), p_b(b)
        {}

  bool  operator()(const Point& c) const
        {
          if ( (c == p_a) || ( c == p_b ) ) return false;
          return left_turn( p_a, p_b, c );
        }

private:
    Point p_a;
    Point p_b;
};


template <class Point>
struct p_Less_xy
{
  typedef bool    result_type;
  typedef  Arity_tag< 2 >   Arity;

  bool operator()( const Point& p1, const Point& p2) const
       { return lexicographically_xy_smaller( p1, p2); }
};

template <class Point>
struct p_Greater_xy
{
  typedef bool    result_type;
  typedef  Arity_tag< 2 >   Arity;

  bool operator()( const Point& p1, const Point& p2) const
       { return lexicographically_xy_larger( p1, p2); }
};

template <class Point>
struct p_Less_yx
{
  typedef bool    result_type;
  typedef  Arity_tag< 2 >   Arity;

  bool operator()( const Point& p1, const Point& p2) const
       { return lexicographically_yx_smaller( p1, p2); }
};

template <class Point>
struct p_Greater_yx
{
  typedef bool    result_type;
  typedef  Arity_tag< 2 >   Arity;

  bool operator()( const Point& p1, const Point& p2) const
       { return lexicographically_yx_larger( p1, p2); }
};

template <class Point>
class p_Less_dist_to_line_2p
{
public:
  typedef bool    result_type;
  typedef  Arity_tag< 2 >   Arity;

        p_Less_dist_to_line_2p(const Point& a, const Point& b)
         : p_a(a), p_b(b)
        {}

  bool  operator()(const Point& c, const Point& d) const
        {
          Comparison_result 
            res = compare_signed_distance_to_line( p_a, p_b, c, d);
          if ( res == LARGER )
          {
              return false;
          }
          else if ( res == SMALLER )
          {
              return true;
          }
          else
          {
              return lexicographically_xy_smaller( c, d );
          }
        }

private:
  Point           p_a;
  Point           p_b;
};

template <class Point>
class p_Less_dist_to_line_2
{
public:
  typedef bool    result_type;
  typedef  Arity_tag< 4 >   Arity;

  bool  operator()(const Point&a, const Point& b,
                   const Point& c, const Point& d) const
        {
          Comparison_result 
            res = compare_signed_distance_to_line( a, b, c, d);
          if ( res == LARGER )
          {
              return false;
          }
          else if ( res == SMALLER )
          {
              return true;
          }
          else
          {
              return lexicographically_xy_smaller( c, d );
          }
        }

};

template <class Point>
class p_Less_negative_dist_to_line_2p
{
public:
  typedef bool    result_type;
  typedef  Arity_tag< 2 >   Arity;

        p_Less_negative_dist_to_line_2p(const Point& a, const Point& b)
         : p_a(a), p_b(b)
        {}

  bool  operator()(const Point& c, const Point& d) const
        {
          Comparison_result 
            res = compare_signed_distance_to_line( p_a, p_b, c, d);
          if ( res == LARGER )
          {
              return true;
          }
          else if ( res == SMALLER )
          {
              return false;
          }
          else
          {
              return lexicographically_xy_smaller( c, d );
          }
        }

private:
  Point           p_a;
  Point           p_b;
};

template <class Point>
class p_Less_rotate_ccw
{
public:
  typedef bool    result_type;
  typedef  Arity_tag< 3 >   Arity;

  bool  operator()(const Point& r, const Point& p, const Point& q) const
        {
          Orientation ori = orientation(r, p, q);
          if ( ori == LEFTTURN )
          {
              return true;
          }
          else if ( ori == RIGHTTURN )
          {
              return false;
          }
          else
          {
              if (p == r) return false;
              if (q == r) return true;
              if (p == q)         return false;
              return  collinear_are_ordered_along_line( r, q, p);
          }
        }
};

template <class Point>
class p_Less_rotate_ccw_safer
{
public:
  typedef bool    result_type;
  typedef  Arity_tag< 2 >   Arity;

        p_Less_rotate_ccw_safer(const Point& p) : rot_point(p)
        {}

  bool  operator()(const Point& p, const Point& q) const
        {
          if (p == rot_point) return false;
          if (q == rot_point) return true;
          if (p == q)         return false;
          Orientation ori = orientation(rot_point, p, q);
          if ( ori == LEFTTURN )
          {
              return true;
          }
          else if ( ori == RIGHTTURN )
          {
              return false;
          }
          else
          {
              return  collinear_are_ordered_along_line( rot_point, q, p);
          }
        }
  
  void  set_rotation_center( const Point& p)
        { rot_point = p; }


private:
  Point  rot_point;
};

template <class Point>
class p_Less_rotate_ccw_E
{
public:
  typedef bool    result_type;
  typedef  Arity_tag< 2 >   Arity;

        p_Less_rotate_ccw_E(const Point& p) : rot_point(p)
        {}

  bool  operator()(const Point& p, const Point& q) const
        {
          Orientation ori = orientation(rot_point, p, q);
          if ( ori == LEFTTURN )
          {
              return true;
          }
          else if ( ori == RIGHTTURN )
          {
              return false;
          }
          else
          {
              return  has_larger_dist_to_point( rot_point, p, q) ;
          }
        }

  void  set_rotation_center( const Point& p)
        { rot_point = p; }


private:
  Point  rot_point;
};

template <class R>
class r_Right_of_line
{
public:
  typedef bool    result_type;
  typedef  Arity_tag< 1 >   Arity;

  typedef typename R::Point_2  Point;
  typedef typename R::Line_2   Line;

        r_Right_of_line(const Point& a, const Point& b) : l_ab( a, b )
        {}

  bool  operator()(const Point& c) const
        {
          if ( l_ab.is_degenerate() ) return false;
          return (l_ab.oriented_side(c) == ON_NEGATIVE_SIDE);
        }

private:
  Line    l_ab;
};

template <class R>
class r_Left_of_line
{
public:
  typedef bool                 result_type;
  typedef  Arity_tag< 3 >   Arity;

  typedef typename R::Point_2  Point;
  typedef typename R::Line_2   Line;

        r_Left_of_line(): line_constructed(false)
        {}

  bool  operator()(const Point& a, const Point& b, const Point& c) const
        { 
           if (!line_constructed)
           {
              line_constructed = true;
              l_ab = Line(a,b);
              is_deg = a == b;
           }
           return ( !is_deg && (l_ab.oriented_side(c) == ON_POSITIVE_SIDE)); 
        }


private:
  mutable bool line_constructed;
  mutable Line l_ab;
  mutable bool is_deg;
};

template <class Point>
class p_Less_dist_to_point
{
 public:
  typedef bool    result_type;
  typedef  Arity_tag< 3 >   Arity;

       p_Less_dist_to_point( ) 
       {}

  bool operator()( const Point& p0, const Point& p1, const Point& p2) const
       { return has_smaller_dist_to_point(p0, p1, p2); }
};

template <class R>
class r_Less_negative_dist_to_line
{
public:
  typedef bool                 result_type;
  typedef  Arity_tag< 2 >   Arity;

  typedef typename R::Point_2  Point;
  typedef typename R::Line_2   Line;

        r_Less_negative_dist_to_line(const Point& a, const Point& b) : 
           l_ab( a, b )
        {}

  bool  operator()(const Point& c, const Point& d) const
        {
          Comparison_result res = compare_signed_distance_to_line(l_ab, c, d);
          if ( res == SMALLER )
          {
              return false;
          }
          else if ( res == EQUAL )
          {
              return lexicographically_xy_smaller( c, d );
          }
          else
          {
              return true;
          }
        }

private:
  Line    l_ab;
};

template <class R>
class r_Less_dist_to_line
{
public:
  typedef bool    result_type;
  typedef  Arity_tag< 4 >   Arity;

  typedef typename R::Point_2  Point;
  typedef typename R::Line_2   Line;

        r_Less_dist_to_line() : line_constructed( false )
        { }
       
  bool  operator()(const Point& a, const Point& b,
                   const Point& c, const Point& d) const
        {
          if (!line_constructed)
          {
             line_constructed = true;
             l_ab = Line(a,b);
          }
          Comparison_result res = compare_signed_distance_to_line(l_ab, c, d);
          if ( res == SMALLER )
          {
              return true;
          }
          else if ( res == EQUAL )
          {
              return lexicographically_xy_smaller( c, d );
          }
          else
          {
              return false;
          }
        }

private:
  mutable bool line_constructed;
  mutable Line l_ab;
};

template <class R>
class r_Less_in_direction
{
public:
  typedef bool    result_type;
  typedef  Arity_tag< 2 >   Arity;

  typedef typename  R::RT           RT;
  typedef typename  R::Direction_2  Direction;
  typedef typename  R::Point_2      Point;
  typedef typename  R::Line_2       Line;

        r_Less_in_direction( const Direction& dir ) : l(Point(RT(0) , RT(0)),
             Direction(-(dir.dy()), dir.dx() ))
        {}

  bool  operator()(const Point& c, const Point& d) const
        {
          Comparison_result res = compare_signed_distance_to_line(l, c, d);
          if ( res == LARGER )
          {
              return true;
          }
          else if ( res == EQUAL )
          {
              return  lexicographically_xy_smaller( c, d) ;
          }
          else
          {
              return false;
          }
        }

private:
  Line  l;
};

template <class Point>
struct p_Left_turn
{
  typedef bool    result_type;
  typedef  Arity_tag< 3 >   Arity;

  bool  operator()(const Point& p, const Point& q, const Point& r) const
        { return left_turn(p,q,r); }
};

template <class Point>
struct p_Right_turn
{
  typedef bool    result_type;
  typedef  Arity_tag< 3 >   Arity;

  bool  operator()(const Point& p, const Point& q, const Point& r) const
        { return right_turn(p,q,r); }
};

#ifndef CGAL_NO_DEPRECATED_CODE
template <class Point>
struct p_Leftturn
{
  typedef bool    result_type;
  typedef  Arity_tag< 3 >   Arity;

  bool  operator()(const Point& p, const Point& q, const Point& r) const
        { return left_turn(p,q,r); }
};
template <class Point>
struct p_Rightturn
{
  typedef bool    result_type;
  typedef  Arity_tag< 3 >   Arity;

  bool  operator()(const Point& p, const Point& q, const Point& r) const
        { return right_turn(p,q,r); }
};
#endif // CGAL_NO_DEPRECATED_CODE

template <class Point>
struct p_Orientation
{
  typedef Orientation    result_type;

  Orientation  
        operator()(const Point& p, const Point& q, const Point& r) const
        { return orientation(p,q,r); }
  Orientation  
        operator()(const Point& p, const Point& q, const Point& r, 
                   const Point& s) const
        { return orientation(p,q,r,s); }
};

CGAL_END_NAMESPACE

#endif // CGAL_PREDICATE_OBJECTS_ON_POINTS_2_H
