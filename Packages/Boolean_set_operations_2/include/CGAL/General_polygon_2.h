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

  typedef Arr_traits                                 Traits_2;
  typedef typename Traits_2::X_monotone_curve_2      X_monotone_curve_2;
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
    m_xcurves.insert(m_xcurves.end(), begin, end);
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

  void push_back(const X_monotone_curve_2& cv)
  {
    m_xcurves.push_back(cv);
  }

  void clear()
  {
    m_xcurves.clear();
  }

  Curve_iterator erase(Curve_iterator it)
  {
    return (m_xcurves.erase(it));
  }


  template <class OutputIterator>
  void approximate(OutputIterator oi, unsigned int n) const
  {
    for(Curve_const_iterator itr = m_xcurves.begin();
        itr != m_xcurves.end();
        ++itr)
    {
      itr->approximate(oi, n);
    }
  }
    
 
};


//-----------------------------------------------------------------------//
//                          operator>>
//-----------------------------------------------------------------------//

template <class Traits>
std::istream &operator>>(std::istream &is, General_polygon_2<Traits>& p)
{
  int n; // number of edges
  is >> n;
  typename Traits::X_monotone_curve_2 cv;
 
  if (is) {
      p.clear();
      for (int i=0; i<n; i++) {
        is >> cv;
        p.push_back(cv);
      }
  }
 
  return is;
}

//-----------------------------------------------------------------------//
//                          operator<<
//-----------------------------------------------------------------------//

template <class Traits>
std::ostream
&operator<<(std::ostream &os, const General_polygon_2<Traits>& p)
{
  typename General_polygon_2<Traits>::Curve_const_iterator i;

  switch(os.iword(IO::mode)) {
    case IO::ASCII :
      os << p.size() << ' ';
      for (i = p.curves_begin(); i != p.curves_end(); ++i) {
        os << *i << ' ';
      }
      return os;

    case IO::BINARY :
      os << p.size();
      for (i = p.curves_begin(); i != p.curves_end(); ++i) {
        os << *i;
      }
      return os;

    default:
      os << "Polygon_2(" << std::endl;
      for (i = p.curves_begin(); i != p.curves_end(); ++i) {
        os << "  " << *i << std::endl;
      }
      os << ")" << std::endl;
      return os;
  }
}


CGAL_END_NAMESPACE

#endif
