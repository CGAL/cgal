// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef GENERAL_POLYGON_2_H
#define GENERAL_POLYGON_2_H

#include <list>

CGAL_BEGIN_NAMESPACE

template <class Arr_traits>
class General_polygon_2
{
public:

  typedef typename Arr_traits::X_monotone_curve_2    X_monotone_curve_2;
  typedef std::list<X_monotone_curve_2>              Containter;
  typedef typename Containter::iterator              Curve_iterator;
  typedef typename Containter::const_iterator        Curve_const_iterator;

protected:
  std::list<X_monotone_curve_2>    m_xcurves;

public:

  General_polygon_2(){}

  template <class CurveIterator>
  General_polygon_2(CurveIterator begin,
                    CurveIterator end) : m_xcurves(begin, end)
  {}

  template <class CurveIterator>
  void insert(CurveIterator begin, CurveIterator end)
  {
    m_xcurves.insert(begin, end, m_xcurves.end());
  }

  bool is_empty() const
  {
    return m_xcurves.empty();
  }

  unsigned int size() const
  {
    return m_xcurves.size();
  }

  Curve_iterator curves_begin()
  {
    return m_xcurves.begin();
  }

  Curve_iterator curves_end()
  {
    return m_xcurves.end();
  }

  Curve_const_iterator curves_begin() const
  {
    return m_xcurves.begin();
  }

  Curve_const_iterator curves_end() const
  {
    return m_xcurves.end();
  }
 
};

CGAL_END_NAMESPACE

#endif
