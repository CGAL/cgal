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
// file          : include/CGAL/Pm_traits_checker.h
// package       : pm (4.08)
// author(s)     : Oren Nechushtan	<theoren@math.tau.ac.il>
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
/*

  Edited by 

  Iddo Hanniel
  Oren Nechushtan

*/
#ifndef CGAL_PM_TRAITS_CHECKER_H
#define CGAL_PM_TRAITS_CHECKER_H

/* Pm_traits_checker is a traits template
that recieves two template parameters, which are traits themselves.
The predicates will work well iff all of the parameter traits predicates did, 
giving the same answer.
It is assumed through out this traits that the parameter traits do not require any storage.
*/
CGAL_BEGIN_NAMESPACE

// Tr1 and Tr2 are traits with the same predicate names
// C is a convertion between the geometric objects of Tr1 to those of Tr2
template <class Tr1,class Tr2>
class Pm_traits_checker_default_adaptor;

template <class Tr1,class Tr2,class C=Pm_traits_checker_default_adaptor<Tr1,Tr2> >
class Pm_traits_checker : public Tr1
/* The inheritance here is to ensure that constants, enumerations, public classes, etc are inherited to the checker whenever this is required */
{
  Tr1 t;
  Tr2 b;
  C P; //P converts Tr1 objects to Tr2 objects (to enable comparison)
public:

  typedef typename Tr1::X_curve X_curve;
  typedef typename Tr1::Point_2 Point;
  typedef typename Tr1::Point_2 Point_2;
  typedef typename Tr1::Curve_point_status Curve_point_status;

  Pm_traits_checker() : Tr1(),t(),b(),P() {};
  ~Pm_traits_checker() {}

  Point curve_source(const X_curve & cv) const 
  { 
    CGAL_assertion(P(t.curve_source(cv))==b.curve_source(P(cv)));
    return t.curve_source(cv);
  }
  Point curve_target(const X_curve & cv) const 
  {
    CGAL_assertion(P(t.curve_target(cv))==b.curve_target(P(cv)));
    return t.curve_target(cv);
  }
  bool curve_is_vertical(const X_curve & cv) const 
  {
    CGAL_assertion(t.curve_is_vertical(cv)==b.curve_is_vertical(P(cv)));
    return t.curve_is_vertical(cv);
  }	
  bool curve_is_in_x_range(const X_curve & cv, const Point & q) const
    {
    CGAL_assertion(t.curve_is_in_x_range(cv,q)==b.curve_is_in_x_range(P(cv),P(q)));
    return t.curve_is_in_x_range(cv,q);
    }


  Curve_point_status curve_get_point_status(const X_curve &cv, const Point & p) const
  {
    CGAL_assertion((int)(t.curve_get_point_status(cv,p))==(int)(b.curve_get_point_status(P(cv),P(p))));
    return t.curve_get_point_status(cv,p);
  }
  
  Comparison_result 
  curve_compare_at_x(const X_curve &cv1, const X_curve &cv2, const Point &q) 
    const 
  {
    CGAL_assertion(t.curve_compare_at_x(cv1,cv2,q)==b.curve_compare_at_x(P(cv1),P(cv2),P(q)));
    return t.curve_compare_at_x(cv1,cv2,q);
  }
  
  Comparison_result 
  curve_compare_at_x_left(const X_curve &cv1, const X_curve &cv2, 
                          const Point &q) const 
  {
    CGAL_assertion(t.curve_compare_at_x_left(cv1,cv2,q)==b.curve_compare_at_x_left(P(cv1),P(cv2),P(q)));
    return t.curve_compare_at_x_left(cv1,cv2,q);
  }
  
  Comparison_result 
  curve_compare_at_x_right(const X_curve &cv1, const X_curve &cv2, const Point & q) const 
  {
    CGAL_assertion(t.curve_compare_at_x_right(cv1,cv2,q)==b.curve_compare_at_x_right(P(cv1),P(cv2),P(q)));
    return t.curve_compare_at_x_right(cv1,cv2,q);
  }


  bool curve_is_between_cw(const X_curve &cv, 
                           const X_curve &first, 
                           const X_curve &second, 
                           const Point &point)	const
    {
    CGAL_assertion(t.curve_is_between_cw(cv,first,second,point)==b.curve_is_between_cw(P(cv),P(first),P(second),P(point)));
    return t.curve_is_between_cw(cv,first,second,point);
    }

  Comparison_result compare_x(const Point &p1, const Point &p2) const
  {
    CGAL_assertion(t.compare_x(p1,p2)==b.compare_x(P(p1),P(p2)));
    return t.compare_x(p1,p2);
  }
  Comparison_result compare_y(const Point &p1, const Point &p2) const
  {
    CGAL_assertion(t.compare_y(p1,p2)==b.compare_y(P(p1),P(p2)));
    return t.compare_y(p1,p2);
  }
  Point point_to_left(const Point& p) const  
    {
    CGAL_assertion(P(t.point_to_left(p))==b.point_to_left(P(p)));
    return t.point_to_left(p);
    }
  Point point_to_right(const Point& p) const 
    {
    CGAL_assertion(P(t.point_to_right(p))==b.point_to_right(P(p)));
    return t.point_to_right(p);
    }
  bool curve_is_same(const X_curve &cv1, const X_curve &cv2) const
    {
    CGAL_assertion(t.curve_is_same(cv1,cv2)==b.curve_is_same(P(cv1),P(cv2)));
    return t.curve_is_same(cv1,cv2);
    }

};

template <class Tr1,class Tr2>
class Pm_traits_checker_default_adaptor
{
public:
	const typename Tr2::Point operator()(const typename Tr1::Point& p) const {return typename Tr2::Point(p);}
	const typename Tr2::X_curve operator()(const typename Tr1::X_curve& c) const {return typename Tr2::X_curve(c);}
};

CGAL_END_NAMESPACE

#endif // CGAL_PM_SEGMENT_EXACT_TRAITS_CHECKER_H





