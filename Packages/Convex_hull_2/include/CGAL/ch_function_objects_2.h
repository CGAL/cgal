// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       :
// release_date  :
//
// file          : include/CGAL/ch_function_objects_2.h
// package       : Convex_hull_2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Susan Hert
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CH_FUNCTION_OBJECTS_2_H
#define CH_FUNCTION_OBJECTS_2_H

#include <CGAL/functional_base.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class R>
class r_Less_dist_to_line
{
public:
  typedef bool    result_type;
  typedef  Arity_tag< 4 >   Arity;

  typedef typename R::Point_2  Point;
  typedef typename R::Line_2   Line;

        r_Less_dist_to_line() : line_constructed( false )
        { }

  bool  operator()(const Point& a, const Point& b,
                   const Point& c, const Point& d) const
        {
          if (!line_constructed)
          {
             line_constructed = true;
             l_ab = Line(a,b);
          }
          Comparison_result res = compare_signed_distance_to_line(l_ab, c, d);
          if ( res == SMALLER )
          {
              return true;
          }
          else if ( res == EQUAL )
          {
              return lexicographically_xy_smaller( c, d );
          }
          else
          {
              return false;
          }
        }

private:
  mutable bool line_constructed;
  mutable Line l_ab;
};

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CH_FUNCTION_OBJECTS_2_H	

