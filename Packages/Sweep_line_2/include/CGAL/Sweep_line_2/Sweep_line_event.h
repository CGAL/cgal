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
// file          : include/CGAL/Sweep_line_event.h
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
#ifndef CGAL_SWEEP_LINE_EVENT_H
#define CGAL_SWEEP_LINE_EVENT_H

#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <vector>
#include <set>

CGAL_BEGIN_NAMESPACE

template<class SweepLineTraits_2>
class Sweep_line_event
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::X_curve_2 X_curve_2;
  typedef typename Traits::Point_2 Point_2;

#ifdef OLD_IMPL
  typedef Sweep_line_subcurve<Traits> SubCurve;
  typedef Curve_less_functor<Traits> CurveLess;
  typedef typename std::set<SubCurve *, CurveLess> SubcurveContainer;
  typedef typename SubcurveContainer::iterator SubCurveIter;

  Sweep_line_event(const Point_2 &point, Traits *traits) {
    CurveLess lessFunc(traits);
    m_leftCurves = new SubcurveContainer(lessFunc);
    m_rightCurves = new SubcurveContainer(lessFunc);
    m_point = point;
  }
#else
  typedef Sweep_line_subcurve<Traits> SubCurve;
  typedef typename std::list<SubCurve *> SubcurveContainer;
  typedef typename SubcurveContainer::iterator SubCurveIter;

  Sweep_line_event(const Point_2 &point, Traits *traits) {
    m_leftCurves = new SubcurveContainer();
    m_rightCurves = new SubcurveContainer();
    m_point = point;
    m_traits = traits;
  }
#endif

  ~Sweep_line_event() {
    delete m_leftCurves;
    delete m_rightCurves;
  }

  void addCurveToLeft(SubCurve *curve) 
  {
    if (m_leftCurves->empty())
      m_leftCurves->push_back(curve);
    else 
    {
      SubCurveIter iter = m_leftCurves->begin();
      while ( iter != m_leftCurves->end() &&
	      m_traits->curve_compare_at_x_right(curve->getCurve(),
					  (*iter)->getCurve(), 
					   *(curve->getReferencePoint())) 
	      == LARGER)
      {
	++iter;
      }
      if ( iter == m_leftCurves->end() ||
	   !m_traits->curve_is_same((*iter)->getCurve(), curve->getCurve()))
	m_leftCurves->insert(iter, curve);
    }
  }

  void addCurveToRight(SubCurve *curve) 
  {
    if (m_rightCurves->empty())
      m_rightCurves->push_back(curve);
    else 
    {
      SubCurveIter iter = m_rightCurves->begin();
      while ( iter != m_rightCurves->end() &&
	      m_traits->curve_compare_at_x_right(curve->getCurve(),
						 (*iter)->getCurve(), 
						 *(curve->getReferencePoint())) 
	      == LARGER)
      {
	++iter;
      }
      if ( iter == m_rightCurves->end() ||
	   !m_traits->curve_is_same((*iter)->getCurve(), curve->getCurve()))
      {
	m_rightCurves->insert(iter, curve);
      }
    }
  }

  SubCurveIter leftCurvesBegin() {
    return m_leftCurves->begin();
  }

  SubCurveIter leftCurvesEnd() {
    return m_leftCurves->end();
  }

  SubCurveIter rightCurvesBegin() {
    return m_rightCurves->begin();
  }

  SubCurveIter rightCurvesEnd() {
    return m_rightCurves->end();
  }

  int getNumRightCurves() {
    return m_rightCurves->size();
  }

  bool hasLeftCurves() {
    return !m_leftCurves->empty();
  }

  const Point_2 &getPoint() {
    return m_point;
  }

  void Print();
 
private:

  Traits *m_traits;

  SubcurveContainer *m_leftCurves;
  SubcurveContainer *m_rightCurves;

  /*! this is true only if this is a right end of one of the original curves */
  bool m_isEndPoint;

  Point_2 m_point;
  
};


template<class SweepLineTraits_2>
void 
Sweep_line_event<SweepLineTraits_2>::
Print() 
{
  std::cout << "\tLeft curves: \n" ;
  for ( SubCurveIter iter = m_leftCurves->begin() ;
	iter != m_leftCurves->end() ; ++iter )
  {
    const X_curve_2 &c = (*iter)->getCurve();
    std::cout << "\t(" << c << ") \n";
  }
  std::cout << std::endl;
  std::cout << "\tRight curves: \n" ;
  for ( SubCurveIter iter = m_rightCurves->begin() ;
	iter != m_rightCurves->end() ; ++iter )
  {
    const X_curve_2 &c = (*iter)->getCurve();
    std::cout << "\t(" << c << ") \n";
  }
  std::cout << std::endl;
}





/*
  void addCurveToLeft(SubCurve *curve) {
#ifndef VERBOSE
    m_leftCurves->insert(curve); 
#else
    std::pair<SubCurveIter, bool> res = m_leftCurves->insert(curve);
    if (res.second!=true)
      std::cout << "Warning: addCurveToLeft\n";
#endif

  }
  void addCurveToRight(SubCurve *curve) {
#ifndef VERBOSE
    m_rightCurves->insert(curve);
#else
    std::pair<SubCurveIter, bool> res = m_rightCurves->insert(curve);
    if (res.second!=true)
      std::cout << "Warning: addCurveToRight\n";
#endif
  }
*/
CGAL_END_NAMESPACE

#endif // CGAL_SWEEP_LINE_EVENT_H
