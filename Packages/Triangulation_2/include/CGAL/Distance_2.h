// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Distance_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_DISTANCE_2_H
#define CGAL_DISTANCE_2_H

#include<CGAL/number_utils.h>

CGAL_BEGIN_NAMESPACE


template <class I>
class Distance_2{
public:
    typedef typename I::Point Point;
    Distance_2(const I* traits = NULL)
    {}

    Distance_2(const Point& p0,
                    const I* traits = NULL)
        : _p0(p0)
    {}


    Distance_2(const Point& p0,
                    const Point& p1,
                    const I* traits = NULL)
        : _p0(p0), _p1(p1)
    {}


    Distance_2(const Point& p0,
                    const Point& p1,
                    const Point& p2,
                    const I* traits = NULL)
        : _p0(p0), _p1(p1), _p2(p2)
    {}

    void
    set_point(int i, const Point& p)
    {
        CGAL_triangulation_precondition(i == 0 || i == 1 || i == 2);
        switch(i){
        case 0:
            _p0 = p;
            break;
        case 1:
            _p1 = p;
            break;
        default:
            _p2 = p;
        }
    }

    Point
    get_point(int i) const
    {
      CGAL_triangulation_precondition(i == 0 || i == 1 || i == 2);
      switch(i){
      case 0:
        return _p0;
      case 1:
        return _p1;
      }
      return _p2;
    }

    Comparison_result
    compare() const
    {
        return CGAL::compare(
          (typename Point::R::FT)squared_distance(_p0, _p1),
          (typename Point::R::FT)squared_distance(_p0, _p2));

    }

private:
    Point _p0, _p1, _p2;
};

CGAL_END_NAMESPACE

#endif // CGAL_DISTANCE_2_H
