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
  
  bool operator()(const Point& p1, const Point& p2) const { 

    Comparison_result rx = m_traits->compare_x(p1,p2);
    if (rx == SMALLER)
      return true;
    if (rx == LARGER)
      return false;

    Comparison_result ry = m_traits->compare_y(p1,p2);
    if (ry == SMALLER)
      return true;
    else if (ry == LARGER)
      return false;
  
    return false;  // get here only if p1 == p2.
  }
private:

  /*! a pointer to a trits class */
  Traits *m_traits;
};

template <class SweepLineTraits_2>
class Curve_less_functor 
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::X_curve_2 X_curve_2;
  typedef Sweep_line_subcurve<Traits> Subcurve;

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

template <class SweepLineTraits_2>
class Status_line_curve_less_functor 
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::X_curve_2 X_curve_2;
  typedef Sweep_line_subcurve<Traits> Subcurve;

  Status_line_curve_less_functor(Traits *traits) : m_traits(traits) {}
  
  bool operator()(const Subcurve* c1, const Subcurve* c2) const { 
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
      if ( m_traits->curve_get_point_status(cv2, c1->getSource()) == 
	   Traits::UNDER_CURVE &&
	   m_traits->curve_get_point_status(cv2, c1->getTarget()) == 
	   Traits::UNDER_CURVE) {
	return true;
      }
      return false;
    }
    if ( m_traits->curve_is_vertical(cv2))
    {
      if ( m_traits->curve_get_point_status(cv1, c2->getSource()) == 
	   Traits::ABOVE_CURVE &&
	   m_traits->curve_get_point_status(cv1, c2->getTarget()) == 
	   Traits::ABOVE_CURVE) {
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
    return false;
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
