// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : demo/Robustness/include/CGAL/further_point_generators_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_FURTHER_POINT_GENERATORS_2_H
#define CGAL_FURTHER_POINT_GENERATORS_2_H

#include <CGAL/Cartesian.h>
#include <iterator>
#include <CGAL/Random.h>

namespace CGAL {

class Point_on_horizontal_bar
{
  public:
    typedef Cartesian<double>::Point_2  Point_2;
    typedef std::output_iterator_tag    iterator_category;
    typedef Point_2                     value_type;
    typedef void                        difference_type;
    typedef void                        pointer;
    typedef void                        reference;

    Point_on_horizontal_bar( int lo_x, int hi_x, int Y)
     : a(lo_x), b(hi_x), y(Y)
    { generate_point(); }

    void
    generate_point()
    {
      int x = RS.get_int(a,b);
      p = Point_2( x, y);
    }

    Point_2
    operator*()
    { return p; }

    Point_on_horizontal_bar&
    operator++()
    {
        generate_point();
        return *this;
    }

    Point_on_horizontal_bar
    operator++(int)
    {
        Point_on_horizontal_bar tmp = *this;
        ++(*this);
        return tmp;
    }


  private:
    int a, b;
    Random  RS;
    int  y;
    Point_2 p;

};

class Point_on_vertical_bar
{
  public:
    typedef Cartesian<double>::Point_2  Point_2;
    typedef std::output_iterator_tag    iterator_category;
    typedef Point_2                     value_type;
    typedef void                        difference_type;
    typedef void                        pointer;
    typedef void                        reference;

    Point_on_vertical_bar( int lo_y, int hi_y, int X)
     : a(lo_y), b(hi_y), x(X)
    { generate_point(); }

    void
    generate_point()
    {
      int y = RS.get_int(a,b);
      p = Point_2( x, y);
    }

    Point_2
    operator*()
    { return p; }

    Point_on_vertical_bar&
    operator++()
    {
        generate_point();
        return *this;
    }

    Point_on_vertical_bar
    operator++(int)
    {
        Point_on_vertical_bar tmp = *this;
        ++(*this);
        return tmp;
    }


  private:
    int a, b;
    Random RS;
    int  x;
    Point_2 p;

};

} // namespace CGAL
#endif //  CGAL_FURTHER_POINT_GENERATORS_2_H
