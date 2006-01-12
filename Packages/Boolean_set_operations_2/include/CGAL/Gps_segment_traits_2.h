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

#ifndef GPS_SEGMENT_TRAITS_2_H
#define GPS_SEGMENT_TRAITS_2_H


#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2/Polygon_2_curve_iterator.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h>

CGAL_BEGIN_NAMESPACE

template < class Kernel_, 
           class Container_ = std::vector<typename Kernel_::Point_2>,
           class Arr_seg_traits_ = Arr_segment_traits_2<Kernel_> >
class Gps_segment_traits_2 : public Arr_seg_traits_
{
  typedef Arr_seg_traits_                               Base;
  typedef Gps_segment_traits_2<Kernel_>                 Self;

public:

  typedef CGAL::Polygon_2<Kernel_, Container_>          Polygon_2;

  //Polygon_with_holes_2 can be non-simple polygon, with holes that are 
  //entirely inside home , or some vertices of the polygon and its holes
  // may overlap.
  typedef CGAL::Polygon_with_holes_2<Kernel_, Container_>    
                                                        Polygon_with_holes_2;
  
  typedef typename Base::X_monotone_curve_2             X_monotone_curve_2;

  typedef Polygon_2_curve_iterator<X_monotone_curve_2, 
                                   Polygon_2> Curve_const_iterator;

 
  typedef typename Base::Point_2                        Point_2;

  typedef typename Base::Equal_2                        Equal_2;
  typedef typename Base::Construct_min_vertex_2         Construct_min_vertex_2;
  typedef typename Base::Construct_max_vertex_2         Construct_max_vertex_2;
  typedef typename Base::Compare_endpoints_xy_2         Compare_endpoints_xy_2;


  class Construct_polygon_2
  {
  public:

    template<class XCurveIterator>
    void operator()(XCurveIterator begin,
                    XCurveIterator end,
                    Polygon_2& pgn)
    {
      Self self;
      Compare_endpoints_xy_2 cmp_ends = self.compare_endpoints_xy_2_object();
      Construct_min_vertex_2 min_v    = self.construct_min_vertex_2_object();
      Construct_max_vertex_2 max_v    = self.construct_max_vertex_2_object();
      for(XCurveIterator itr = begin; itr != end; ++itr)
      {
        if(cmp_ends(*itr) == SMALLER)
          pgn.push_back(max_v(*itr));
        else
          pgn.push_back(min_v(*itr));
      }
    }
  };

  Construct_polygon_2 construct_polygon_2_object() const
  {
    return Construct_polygon_2();
  }


  class Construct_curves_2
  {
  public:

    std::pair<Curve_const_iterator,
              Curve_const_iterator> operator()(const Polygon_2& pgn)
    {
      Curve_const_iterator c_begin(&pgn, pgn.edges_begin());
      Curve_const_iterator c_end(&pgn, pgn.edges_end());

      return (std::make_pair(c_begin, c_end));
    }
  };

  Construct_curves_2 construct_curves_2_object()
  {
    return Construct_curves_2();
  }

   typedef Is_valid_2<Self>    Is_valid_2;
 
  Is_valid_2 is_valid_2_object()
  {
    return Is_valid_2();
  }
};


CGAL_END_NAMESPACE

#endif
