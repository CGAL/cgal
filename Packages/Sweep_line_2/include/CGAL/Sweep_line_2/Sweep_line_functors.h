// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-44 $
// release_date  : $CGAL_Date: 2001/03/09 $
//
// file          : include/CGAL/Sweep_line_functors.h
// package       : arr (1.87)
// maintainer    : Tali Zvi <talizvi@post.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_SWEEP_LINE_FUNCTORS_H
#define CGAL_SWEEP_LINE_FUNCTORS_H

CGAL_BEGIN_NAMESPACE

template <class Point, class SweepLineTraits_2>
class Point_less_functor 
{
public:
  typedef SweepLineTraits_2 Traits;
  
  Point_less_functor(Traits *traits) : m_traits(traits) {}
  
  bool operator()(const Point& p1, const Point& p2) const 
  { 
    return (m_traits->compare_xy(p1,p2) == SMALLER);
  }

private:

  /*! a pointer to a trits class */
  Traits *m_traits;
};


template <class SweepLineTraits_2, class Subcurve>
class Status_line_curve_less_functor 
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::X_curve_2 X_curve_2;
  typedef bool (Status_line_curve_less_functor::*func)
                      (const Subcurve*, const Subcurve*) const;

  struct Compare_param {
    Compare_param(Traits *t) : m_compare_func(1), m_traits(t)  {}
    int m_compare_func;
    Traits *m_traits;
  };

  Status_line_curve_less_functor(Compare_param *p) : m_compare_param(p) {
    m_compare[0] = &Status_line_curve_less_functor::compare_at;
    m_compare[1] = &Status_line_curve_less_functor::compare_right;
  }

  bool operator()(const Subcurve* c1, const Subcurve* c2) const {
    return (this->*m_compare[m_compare_param->m_compare_func])(c1, c2);
  }


  bool compare_at(const Subcurve* c1, const Subcurve* c2)  const 
  {
#if 1  // this may ot work with conics. Need to verify
    const Point_2 *p = c1->getReferencePoint();
    Comparison_result r = 
      m_compare_param->m_traits->curve_compare_at_x(c1->getCurve(), 
    				   c2->getCurve(), 
    				   *p);
    if ( r == SMALLER) {
      return true;
    } 
    return false;
#else
    if ( m_compare_param->m_traits->curve_get_point_status(c2->getCurve(),
                                                          c1->getLeftEnd())
	 == LARGER)
      return true;
    return false;
#endif
  }

  bool compare_right(const Subcurve* c1, const Subcurve* c2)  const 
  {

    const Point_2 *p = c1->getReferencePoint();
#if 0
    std::cout << "\t\tComparing between:" << *p << "\n"
	      << "\t\t  " << c1->getCurve() << "\n"
	      << "\t\t  " << c2->getCurve() << "\n";
#endif
    
    const X_curve_2 &cv1 = c1->getCurve();
    const X_curve_2 &cv2 = c2->getCurve();
    Traits *t = m_compare_param->m_traits;
    if ( t->curve_is_vertical(cv1) )
    {
      if (t->curve_is_in_x_range(cv2, c1->getSource()) &&
	  t->curve_get_point_status(cv2, c1->getTopEnd()) == LARGER )
      {
	return true;
      }
      return false;
    }
    if ( t->curve_is_vertical(cv2))
    {
      if (t->curve_is_in_x_range(cv1, c2->getSource()) &&
	  t->curve_get_point_status(cv1, c2->getBottomEnd()) == SMALLER)
      {
	return true;
      }
      return false;
    }

    // non of the curves is vertical... 
    Comparison_result r =  t->curve_compare_at_x (c1->getCurve(), 
						  c2->getCurve(), 
						  *p);

    if (r == EQUAL)
      r = t->curve_compare_at_x_right(c1->getCurve(), 
					     c2->getCurve(), 
					     *p);
    if ( r == SMALLER) {
      return true;
    } 
    if ( r == LARGER ) {
      return false;
    }

    // r = EQUAL
    return ( c1->getId() < c2->getId() );
  }

  void setReference(Point_2 point) {
    m_point = point;
  }

  Compare_param *m_compare_param;

private:
  func m_compare[2];


};

CGAL_END_NAMESPACE

#endif // CGAL_SWEEP_LINE_FUNCTORS_H
