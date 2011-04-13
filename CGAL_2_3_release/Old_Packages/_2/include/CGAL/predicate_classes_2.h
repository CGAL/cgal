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
// release_date  : 
// 
// file          : predicate_classes_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 
#ifndef CGAL_PREDICATES_CLASSES_2_H
#define CGAL_PREDICATES_CLASSES_2_H

CGAL_BEGIN_NAMESPACE

template <class Point>      struct X_equal_2;
template <class Point>      struct Y_equal_2;
template <class Point>      struct Compare_x_2;
template <class Point>      struct Compare_y_2;
template <class Point>      struct Compare_deltax_deltay_2;
template <class Point>      struct Compare_xy_2;
template <class Point>      struct Less_xy_2;
template <class Point>      struct Compare_yx_2;
template <class Point>      struct Less_yx_2;
template <class Point>      struct Are_ordered_along_line_2 ;
template <class Point>      struct Collinear_are_ordered_along_line_2;
template <class Point>      struct Are_strictly_ordered_along_line_2;
template <class Point>      struct Collinear_are_strictly_ordered_along_line_2;
template <class Point>      struct Collinear_2;
template <class Point>      struct Left_turn_2;
template <class Point>      struct Right_turn_2;
#ifndef CGAL_NO_DEPRECATED_CODE
template <class Point>      struct Leftturn_2;
template <class Point>      struct Rightturn_2;
#endif
template <class Point>      struct Orientation_2;
template <class Point>      struct Side_of_oriented_circle_2;
template <class Point>      struct Side_of_bounded_circle_2;
template <class Point>      class  Compare_distance_to_point_2;
template <class Point>      class  Less_distance_to_point_2;
template <class Point, class Line> class  Compare_signed_distance_to_line_2;
template <class Point, class Line> class  Less_signed_distance_to_line_2;
template <class Point>      class Compare_signed_distance_to_implicit_line_2;
template <class Point>      class Less_signed_distance_to_implicit_line_2;
template <class Point, class Line> struct Compare_y_at_point_2;
template <class Line>       struct Compare_y_at_implicit_point_2;
template <class Point, class Line> struct Compare_y_of_lines_at_point_2;
template <class Line>       struct Compare_y_of_lines_at_implicit_point_2;
template <class Point, class Line> struct Compare_x_implicit_point_2;
template <class Line>       struct Compare_x_implicit_points_2;
template <class Line>       struct Compare_x_implicit_points_same_line_2;
template <class Point, class Line> struct Compare_y_implicit_point_2;
template <class Line>       struct Compare_y_implicit_points_2;
template <class Line>       struct Compare_y_implicit_points_same_line_2;

template <class Point>
struct X_equal_2
 : public CGAL_STD::binary_function< Point, Point, bool>
{
  bool operator()( const Point& p, const Point& q) const
       { return x_equal( p, q); }
};

template <class Point>
struct Y_equal_2
 : public CGAL_STD::binary_function< Point, Point, bool>
{
  bool operator()( const Point& p, const Point& q) const
       { return y_equal( p, q); }
};

template <class Point>
struct Compare_x_2
 : public CGAL_STD::binary_function< Point, Point, Comparison_result>
{
  Comparison_result operator()( const Point& p, const Point& q) const
       { return compare_x( p, q); }
};

template <class Point>
struct Compare_y_2
 : public CGAL_STD::binary_function< Point, Point, Comparison_result>
{
  Comparison_result operator()( const Point& p, const Point& q) const
       { return compare_y( p, q); }
};

template <class Point>
struct Compare_deltax_deltay_2
{
  Comparison_result operator()( const Point& p, const Point& q,
                                     const Point& r, const Point& s )
       { return compare_deltax_deltay( p,q,r,s); }
};

template <class Point>
struct Compare_xy_2
 : public CGAL_STD::binary_function< Point, Point, Comparison_result>
{
  Comparison_result operator()( const Point& p, const Point& q) const
      { return compare_lexicographically_xy( p,q); }
};

template <class Point>
struct p_Equal_xy
 : public CGAL_STD::binary_function< Point, Point, bool>
{
  bool operator()( const Point& p, const Point& q) const
       { return equal_xy( p,q); }
};

template <class Point>
struct Less_xy_2
 : public CGAL_STD::binary_function< Point, Point, bool>
{
  bool operator()( const Point& p, const Point& q) const
       { return lexicographically_xy_smaller( p,q); }
};

template <class Point>
struct Compare_yx_2
 : public CGAL_STD::binary_function< Point, Point, Comparison_result>
{
  Comparison_result operator()( const Point& p, const Point& q) const
      { return compare_lexicographically_yx( p,q); }
};

template <class Point>
struct Less_yx_2
 : public CGAL_STD::binary_function< Point, Point, bool>
{
  bool operator()( const Point& p, const Point& q) const
       { return lexicographically_yx_smaller( p,q); }
};


template <class Point>
struct Are_ordered_along_line_2
{
  bool operator()( const Point& p, const Point& q, const Point& r) const
       { return are_ordered_along_line( p,q,r); }
};

template <class Point>
struct Collinear_are_ordered_along_line_2
{
  bool operator()( const Point& p, const Point& q, const Point& r) const
       { return collinear_are_ordered_along_line( p,q,r); }
};

template <class Point>
struct Are_strictly_ordered_along_line_2
{
  bool operator()( const Point& p, const Point& q, const Point& r) const
       { return are_strictly_ordered_along_line( p,q,r); }
};

template <class Point>
struct Collinear_are_strictly_ordered_along_line_2
{
  bool operator()( const Point& p, const Point& q, const Point& r) const
       { return collinear_are_strictly_ordered_along_line( p,q,r); }
};

template <class Point>
struct Collinear_2
{
  bool operator()( const Point& p, const Point& q, const Point& r) const
       { return collinear( p,q,r); }
};


template <class Point>
struct Left_turn_2
{
  bool operator()( const Point& p, const Point& q, const Point& r) const
       { return left_turn( p,q,r); }
};

template <class Point>
struct Right_turn_2
{
  bool operator()( const Point& p, const Point& q, const Point& r) const
       { return right_turn( p,q,r); }
};

#ifndef CGAL_NO_DEPRECATED_CODE
template <class Point>
struct Leftturn_2
{
  bool operator()( const Point& p, const Point& q, const Point& r) const
       { return left_turn( p,q,r); }
};

template <class Point>
struct Rightturn_2
{
  bool operator()( const Point& p, const Point& q, const Point& r) const
       { return right_turn( p,q,r); }
};
#endif

template <class Point>
struct Orientation_2
{
  Orientation operator()( const Point& p, const Point& q, const Point& r) const
       { return orientation( p,q,r); }
};

template <class Point>
struct Side_of_oriented_circle_2
{
  Oriented_side operator()( const Point& p, const Point& q,
                                 const Point& r, const Point& s)
       { return side_of_oriented_circle( p,q,r,s); }
};

template <class Point>
struct Side_of_bounded_circle_2
{
  Oriented_side operator()( const Point& p, const Point& q,
                                 const Point& r, const Point& s)
       { return side_of_bounded_circle( p,q,r,s); }
};



template <class Point>
class Compare_distance_to_point_2
 : public CGAL_STD::binary_function< Point, Point, Comparison_result>
{
 public:
  Compare_distance_to_point_2( const Point& p) : _p(p) {}
  Comparison_result
       operator()( const Point& q, const Point& r)
       { return compare_distance_to_point( _p,q,r); }
 private:
  Point _p;
};

template <class Point>
class Less_distance_to_point_2
{
 public:
  Less_distance_to_point_2( ) {}
  bool operator()( const Point& p, const Point& q, const Point& r)
       { return has_smaller_distance_to_point(p,q,r); }
};

template <class Point, class Line>
class Compare_signed_distance_to_line_2
 : public CGAL_STD::binary_function< Point, Point, Comparison_result>
{
 public:
  Compare_signed_distance_to_line_2( const Line& l) : _l(l) {}
  Comparison_result
       operator()( const Point& p, const Point& q)
       { return compare_signed_distance_to_line( _l, p, q); }
 private:
  Line _l;
};

template <class Point, class Line>
class Less_signed_distance_to_line_2
 : public CGAL_STD::binary_function< Point, Point, bool>
{
 public:
  Less_signed_distance_to_line_2( const Line& l) : _l(l) {}
  bool operator()( const Line& l, const Point& p, const Point& q)
       { return has_smaller_signed_distance_to_line( l, p, q); }
 private:
  Line _l;
};

template <class Point>
class Compare_signed_distance_to_implicit_line_2
 : public CGAL_STD::binary_function< Point, Point, Comparison_result>
{
 public:
  Compare_signed_distance_to_implicit_line_2( const Point& p, const Point& q)
   : _p(p), _q(q)
  {}
  Comparison_result
       operator()( const Point& r, const Point& s)
       { return compare_signed_distance_to_line( _p, _q, r, s); }
 private:
  Point _p;
  Point _q;
};

template <class Point>
class Less_signed_distance_to_implicit_line_2
 : public CGAL_STD::binary_function< Point, Point, bool>
{
 public:
  Less_signed_distance_to_implicit_line_2( const Point& p, const Point& q)
   : _p(p), _q(q)
  {}
  bool operator()( const Point& r, const Point& s)
       { return has_smaller_signed_distance_to_line( _p, _q, r, s); }
 private:
  Point _p;
  Point _q;
};



template <class Point, class Line>
struct Compare_y_at_point_2
 : public CGAL_STD::binary_function< Point, Line, Comparison_result>
{
  Comparison_result operator()( const Point& p, const Line& l)
       { return compare_y_at_x( p, l); }
};

template <class Line>
struct Compare_y_at_implicit_point_2
{
  Comparison_result operator()( const Line& l1, const Line& l2, const Line& l)
       { return compare_y_at_x( l1, l2, l); }
};

template <class Point, class Line>
struct Compare_y_of_lines_at_point_2
{
  Comparison_result operator()( const Point& p, const Line& l1, const Line& l2)
       { return compare_y_at_x( p, l1, l2); }
};

template <class Line>
struct Compare_y_of_lines_at_implicit_point_2
{
  Comparison_result operator()( const Line& l1, const Line& l2,
                                     const Line& h1, const Line& h2)
       { return compare_y_at_x( l1, l2, h1, h2); }
};

template <class Point, class Line>
struct Compare_x_implicit_point_2
{
  Comparison_result operator()( const Point& p,
                                     const Line& l1, const Line& l2)
       { return compare_x( p, l1, l2); }
};

template <class Line>
struct Compare_x_implicit_points_2
{
  Comparison_result operator()( const Line& l1, const Line& l2,
                                     const Line& h1, const Line& h2)
       { return compare_x( l1, l2, h1, h2); }
};

template <class Line>
struct Compare_x_implicit_points_same_line_2
{
  Comparison_result
       operator()( const Line& l, const Line& h1, const Line& h2)
       { return compare_x( l, h1, h2); }
};

template <class Point, class Line>
struct Compare_y_implicit_point_2
{
  Comparison_result operator()( const Point& p,
                                     const Line& l1, const Line& l2)
       { return compare_y( p, l1, l2); }
};

template <class Line>
struct Compare_y_implicit_points_2
{
  Comparison_result operator()( const Line& l1, const Line& l2,
                                     const Line& h1, const Line& h2)
       { return compare_y( l1, l2, h1, h2); }
};

template <class Line>
struct Compare_y_implicit_points_same_line_2
{
  Comparison_result
       operator()( const Line& l, const Line& h1, const Line& h2)
       { return compare_y( l, h1, h2); }
};

CGAL_END_NAMESPACE

#endif // CGAL_PREDICATES_CLASSES_2_H
