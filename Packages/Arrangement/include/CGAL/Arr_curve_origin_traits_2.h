// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>

#ifndef CGAL_CURVE_ORIGIN_TRAITS_2_H
#define CGAL_CURVE_ORIGIN_TRAITS_2_H

#include <list>

CGAL_BEGIN_NAMESPACE

/*! A generic traits class for maintaining an arrangement of X-monotone curves,
 * such that each point to their respective origin curve. This traits class
 * is templated with an base ordinary traits class and derived from it.
 * The X_monotone_curve_2 type is extracted from the base traits class, and
 * redefines it to have a pointer to a Curve_2 type as an extra field.
 *
 * The pointer is updated when the curves are converted from Curve_2 to
 * X_monotone_curve_2, and when the X_monotone_curve_2 curves are split. All
 * other methods are inherited from the base ordinary traits class.
 */
template <class Traits_>
class Arr_curve_origin_traits_2 : public Traits_ {
public:
  typedef Traits_                                       Traits;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           Org_x_monotone_curve_2;

  class X_monotone_curve_2 : public Org_x_monotone_curve_2 {
  private:
    enum {CURVE_2, X_MONOTONE_CURVE_2} m_type;
    union {
      Curve_2 * m_curve;
      X_monotone_curve_2 * m_x_monotone_curve;
    } m_origin;
    
  public:
    /*! Parameterless constructor */
    X_monotone_curve_2() {}
    
    /*! Constructor */
    X_monotone_curve_2(const Org_x_monotone_curve_2 & cv,
                       Curve_2 * origin) :
      Org_x_monotone_curve_2(cv)
    {
      m_type = CURVE_2;
      m_origin.m_curve = origin;
    }

    /*! Constructor */
    X_monotone_curve_2(const Org_x_monotone_curve_2 & cv,
                       X_monotone_curve_2 * origin) :
      Org_x_monotone_curve_2(cv)
    {
      m_type = X_MONOTONE_CURVE_2;
      m_origin.m_x_monotone_curve = origin;
    }

    /*! Constructor */
    X_monotone_curve_2(const X_monotone_curve_2 * cv)
    {
      m_type = cv->m_type;
      if (m_type == CURVE_2)
        m_origin.m_curve = cv->m_origin.m_curve;
      else
        m_origin.m_x_monotone_curve = cv->m_origin.m_x_monotone_curve;
    }

    /*! obtains the curve type */
    bool is_origin_x_monotone() const
    {
      return (m_type == X_MONOTONE_CURVE_2);
    }

    /*! obtains the curve origin */
    Curve_2 * get_origin_curve() const
    {
      return m_origin.m_curve;
    }

    /*! obtains the curve x-monotone origin */
    X_monotone_curve_2 * get_origin_x_monotone_curve() const
    {
      return m_origin.m_x_monotone_curve;
    }
  };

  typedef typename std::list<Curve_2*>                  List_curve_2;
  typedef typename std::list<X_monotone_curve_2*>       List_x_monotone_curve_2;

  // For backward compatibility:
  typedef Point_2                               Point;
  typedef X_monotone_curve_2                    X_curve;
  typedef Curve_2                               Curve;
  
private:
  Traits m_traits;

  mutable List_curve_2 m_curves;

  mutable List_x_monotone_curve_2 m_x_monotone_curves;

public:
  Arr_curve_origin_traits_2() : m_traits() {}

  ~Arr_curve_origin_traits_2()
  {
    typename List_curve_2::const_iterator cit;
    for (cit = m_curves.begin(); cit != m_curves.end(); ++cit)
      delete *cit;
    m_curves.clear();
    typename List_x_monotone_curve_2::const_iterator xcit;
    for (xcit = m_x_monotone_curves.begin(); xcit != m_x_monotone_curves.end();
         ++xcit)
      delete *xcit;
    m_x_monotone_curves.clear();
  }
  
  /*! Cut the given curve into x-monotone subcurves and insert them to the
   * given output iterator. While segments are x_monotone, still need to pass
   * them out.
   * \param cv the curve.
   * \param o the output iterator
   * \return the past-the-end iterator
   */
  template<class OutputIterator>
  OutputIterator curve_make_x_monotone(const Curve_2 & cv,
                                       OutputIterator o) const
  {
    std::list<Org_x_monotone_curve_2>  org_x_curves;
    m_traits.curve_make_x_monotone(cv, std::back_inserter(org_x_curves));
    typename std::list<Org_x_monotone_curve_2>::const_iterator it;
    Curve_2 * stored_cv = new Curve_2(cv);
    m_curves.push_back(stored_cv);
    for (it = org_x_curves.begin(); it != org_x_curves.end(); it++)
      *o++ = X_monotone_curve_2(*it, stored_cv);
    return o;
  } 

  /*! Split a given curve at a given split point into two sub-curves.
   * \param cv the curve to split
   * \param c1 the output first part of the split curve. Its source is the
   * source of the original curve.
   * \param c2 the output second part of the split curve. Its target is the
   * target of the original curve.
   * \param p the split point.
   * \pre p lies on cv but is not one of its end-points.
   */
  void curve_split(const X_monotone_curve_2 & cv, 
                   X_monotone_curve_2 & c1, X_monotone_curve_2 & c2, 
                   const Point_2 & p) const
  {
    Org_x_monotone_curve_2 org_c1, org_c2;
    m_traits.curve_split(cv, org_c1, org_c2, p);
    X_monotone_curve_2 * stored_cv;
    if (cv.is_origin_x_monotone())
      stored_cv = cv.get_origin_x_monotone_curve();
    else {
      stored_cv = new X_monotone_curve_2(cv);
      m_x_monotone_curves.push_back(stored_cv);
    }
    c1 = X_monotone_curve_2(org_c1, stored_cv);
    c2 = X_monotone_curve_2(org_c2, stored_cv);
  }
};

CGAL_END_NAMESPACE

#endif
