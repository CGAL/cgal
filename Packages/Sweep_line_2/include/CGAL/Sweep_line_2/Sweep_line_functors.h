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

int g_compare_func = 1;

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
#if 0
template <class SweepLineTraits_2, class Subcurve>
class Curve_less_functor 
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::X_curve_2 X_curve_2;
  //typedef Sweep_line_subcurve<Traits> Subcurve;

  Curve_less_functor(Traits *traits) : m_traits(traits) {}
  
  bool operator()(const Subcurve* c1, const Subcurve* c2) const { 
    const Point_2 *p = c1->getReferencePoint();

    if ( !m_traits->curve_is_in_x_range(c1->getCurve(), *p))
      p = &(c1->getLastPoint());
    if ( !m_traits->curve_is_in_x_range(c2->getCurve(), *p))
      p = &(c2->getLastPoint());

   Comparison_result r = 
          m_traits->curve_compare_at_x_right(c1->getCurve(), 
				       c2->getCurve(), 
				       *p);
    if ( r == SMALLER) {
      return true;
    }
    return false;

  }

  void setReference(Point_2 point) {
    m_point = point;
  }

private:

  /*! a pointer to a traits class */
  Traits *m_traits;
};
#endif


template <class SweepLineTraits_2, class Subcurve>
class Status_line_curve_less_functor 
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::X_curve_2 X_curve_2;
  typedef bool (Status_line_curve_less_functor::*func)
                      (const Subcurve*, const Subcurve*) const;

  Status_line_curve_less_functor(Traits *traits) : m_traits(traits) {
    m_compare[0] = &Status_line_curve_less_functor::compare_at;
    m_compare[1] = &Status_line_curve_less_functor::compare_right;
  }
  
  bool operator()(const Subcurve* c1, const Subcurve* c2) const {
    return (this->*m_compare[g_compare_func])(c1, c2);
  }

  func m_compare[2];

  bool compare_at(const Subcurve* c1, const Subcurve* c2)  const 
  {
    const Point_2 *p = c1->getReferencePoint();
    Comparison_result r = 
      m_traits->curve_compare_at_x(c1->getCurve(), 
				   c2->getCurve(), 
				   *p);
    if ( r == SMALLER) {
      return true;
    } 
    return false;
  }

  bool compare_right(const Subcurve* c1, const Subcurve* c2)  const {
    const Point_2 *p = c1->getReferencePoint();
#if 0
    std::cout << "\t\tComparing between:" << *p << "\n"
	      << "\t\t  " << c1->getCurve() << "\n"
	      << "\t\t  " << c2->getCurve() << "\n";
#endif
    
    const X_curve_2 &cv1 = c1->getCurve();
    const X_curve_2 &cv2 = c2->getCurve();
    if ( m_traits->curve_is_vertical(cv1) )
    {
      if (m_traits->curve_is_in_x_range(cv2, c1->getSource()) &&
	  m_traits->curve_get_point_status(cv2, c1->getSource()) == LARGER &&
	  m_traits->curve_is_in_x_range(cv2, c1->getTarget()) &&
	  m_traits->curve_get_point_status(cv2, c1->getTarget()) == LARGER)
      {
	return true;
      }
      return false;
    }
    if ( m_traits->curve_is_vertical(cv2))
    {
      if (m_traits->curve_is_in_x_range(cv1, c2->getSource()) &&
	  m_traits->curve_get_point_status(cv1, c2->getSource()) == SMALLER &&
	  m_traits->curve_is_in_x_range(cv1, c2->getTarget()) &&
	  m_traits->curve_get_point_status(cv1, c2->getTarget()) == SMALLER)
      {
	return true;
      }
      return false;
    }

    // non of the curves is vertical... 
    Comparison_result r = 
      m_traits->curve_compare_at_x_right(c1->getCurve(), 
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

private:

  /*! a pointer to a traits class */
  Traits *m_traits;
};

CGAL_END_NAMESPACE

#endif // CGAL_SWEEP_LINE_FUNCTORS_H
