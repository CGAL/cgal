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
// source        : web/Distance_2.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// ============================================================================


#ifndef CGAL_DISTANCE_2_H
#define CGAL_DISTANCE_2_H

template <class I>
class CGAL_Distance_2{
public:
    typedef typename I::Point Point;
    CGAL_Distance_2(const I* traits = NULL)
    {}

    CGAL_Distance_2(const Point& p0,
                    const I* traits = NULL)
        : _p0(p0)
    {}


    CGAL_Distance_2(const Point& p0,
                    const Point& p1,
                    const I* traits = NULL)
        : _p0(p0), _p1(p1)
    {}


    CGAL_Distance_2(const Point& p0,
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

    CGAL_Comparison_result
    compare() const
    {
        return CGAL_compare(
          (typename Point::R::FT)CGAL_squared_distance(_p0, _p1),
          (typename Point::R::FT)CGAL_squared_distance(_p0, _p2));

    }

private:
    Point _p0, _p1, _p2;
};


#endif // CGAL_DISTANCE_2_H
