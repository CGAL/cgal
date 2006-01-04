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
  typedef typename Traits_2::Point_2                 Point_2;
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

  Orientation orientation() const
  {
    Traits_2 tr;
    typename Traits_2::Compare_xy_2 cmp_xy = tr.compare_xy_2_object();
    typename Traits_2::Compare_y_at_x_right_2 cmp_y_at_x_right = 
      tr.compare_y_at_x_right_2_object();

    Point_2 left_most_v = _source(m_xcurves.front());
    Curve_const_iterator from_left_most = m_xcurves.begin();
    Curve_const_iterator into_left_most = m_xcurves.end();
    --into_left_most;

    Curve_const_iterator ci = m_xcurves.begin();

    for(++ci; ci!=m_xcurves.end(); ++ci)
    {
      Comparison_result res_xy = cmp_xy( _source(*ci), left_most_v);
      if(res_xy == LARGER)
        continue;
      if(res_xy == SMALLER)
      {
        left_most_v =  _source(*ci);
        from_left_most = into_left_most = ci;
        --into_left_most;
      }
      else
      {
        // res_xy == EQUAL
        Curve_const_iterator tmp_from_left_most = ci;
        Curve_const_iterator tmp_into_left_most = ci;
        --tmp_into_left_most;

        Comparison_result res_from = cmp_y_at_x_right(*from_left_most,
                                                      *tmp_from_left_most,
                                                      left_most_v);

        Comparison_result res_to = cmp_y_at_x_right(*into_left_most, 
                                                    *tmp_into_left_most,
                                                    left_most_v);
        
        CGAL_assertion(res_from != EQUAL && res_to != EQUAL);
        if(res_from == LARGER && res_to == SMALLER)
        {
          if(cmp_y_at_x_right(*tmp_from_left_most,
                              *into_left_most,
                              left_most_v) == LARGER)
          {
            from_left_most = tmp_from_left_most;
            into_left_most = tmp_into_left_most;
          }
        }
        else
          if(res_from == SMALLER && res_to == LARGER)
          {
            if(cmp_y_at_x_right(*tmp_into_left_most,
                                *from_left_most,
                                left_most_v) == LARGER)
            {
              from_left_most = tmp_from_left_most;
              into_left_most = tmp_into_left_most;
            }
          }
      }
    }// end for
    Comparison_result res = cmp_y_at_x_right(*into_left_most, 
                                             *from_left_most,
                                             left_most_v);
    CGAL_assertion(res != EQUAL);
    if(res == SMALLER)
      return (CLOCKWISE);
    return (COUNTERCLOCKWISE);
  }


  void reverse_orientation()
  {
    m_xcurves.reverse();
    Traits_2 tr;
    typename Traits_2::Construct_opposite_2 ctr_opp =
      tr.construct_opposite_2_object();
    for(Curve_iterator ci = m_xcurves.begin(); ci != m_xcurves.end(); ++ci)
    {
      const X_monotone_curve_2& opp_cv = ctr_opp(*ci);
      *ci = opp_cv;
    }
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

  private:

  Point_2 _source(const X_monotone_curve_2& cv) const
  {
    Traits_2 tr;
    if(tr.compare_endpoints_xy_2_object()(cv) == SMALLER)
      return (tr.construct_min_vertex_2_object()(cv));
    return (tr.construct_max_vertex_2_object()(cv));
  }

  Point_2 _target(const X_monotone_curve_2& cv) const
  {
    Traits_2 tr;
    if(tr.compare_endpoints_xy_2_object()(cv) == LARGER)
      return (tr.construct_min_vertex_2_object()(cv));
    return (tr.construct_max_vertex_2_object()(cv));
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
