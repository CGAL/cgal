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
// $URL$
// $Id$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_GPS_SEGMENT_TRAITS_2_H
#define CGAL_GPS_SEGMENT_TRAITS_2_H


#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2/Polygon_2_curve_iterator.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h>
#include <CGAL/Boolean_set_operations_2/Gps_traits_adaptor.h>

CGAL_BEGIN_NAMESPACE

template < class Kernel_, 
           class Container_ = std::vector<typename Kernel_::Point_2>,
           class Arr_seg_traits_ = Arr_segment_traits_2<Kernel_> >
class Gps_segment_traits_2 : public Arr_seg_traits_
{
  typedef Arr_seg_traits_                               Base;
  typedef Gps_segment_traits_2<Kernel_,
                               Container_,
                               Arr_seg_traits_>         Self;

public:

  typedef CGAL::Polygon_2<Kernel_, Container_>          Polygon_2;

  //Polygon_with_holes_2 can be a simple polygon , with holes that are 
  //entirely inside him , or some vertices of the polygon and its holes
  // may overlap.
  typedef CGAL::Polygon_with_holes_2<Kernel_, Container_>    
                                                        Polygon_with_holes_2;
  
  typedef typename Base::X_monotone_curve_2             X_monotone_curve_2;

  typedef Polygon_2_curve_iterator<X_monotone_curve_2, 
                                   Polygon_2> Curve_const_iterator;

 
  typedef typename Base::Point_2                        Point_2;

 

  class Construct_polygon_2
  {
    typedef Gps_segment_traits_2<Kernel_,
                               Container_,
                               Arr_seg_traits_>         Self;
    typedef Gps_traits_adaptor<Base>                    Traits_adaptor;

  public:

    template<class XCurveIterator>
    void operator()(XCurveIterator begin,
                    XCurveIterator end,
                    Polygon_2& pgn)
    {
      Traits_adaptor tr;
      typename Traits_adaptor::Construct_vertex_2 ctr_v =
        tr.construct_vertex_2_object();

      for(XCurveIterator itr = begin; itr != end; ++itr)
      {
        pgn.push_back(ctr_v(*itr, 1));
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


  typedef Gps_traits_adaptor<Base>                             Traits_adaptor;
  typedef Is_valid_2<Base, Polygon_2, Polygon_with_holes_2,
                     Curve_const_iterator, Construct_curves_2,
                     typename Traits_adaptor::Construct_vertex_2,
                     typename Traits_adaptor::Orientation_2>   Is_valid_2;

  Is_valid_2 is_valid_2_object()
  {
    Traits_adaptor tr;

    return (Is_valid_2 (this->construct_curves_2_object(),
                        tr.construct_vertex_2_object(),
                        tr.orientation_2_object()));
  }
};


CGAL_END_NAMESPACE

#endif
