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


#ifndef ARR_BATCHED_POINT_LOCATION_H
#define ARR_BATCHED_POINT_LOCATION_H

#include <CGAL/Sweep_line_2/Sweep_line_2_impl.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Arr_point_location/Arr_batched_point_location_event.h>
#include <CGAL/Arr_point_location/Arr_batched_point_location_visitor.h>
#include <CGAL/Arr_point_location/Arr_batched_point_location_meta_traits.h>
#include <vector>

CGAL_BEGIN_NAMESPACE


template<class _Arrangement, class PointsIterator, class OutputIterator> 
OutputIterator locate(const _Arrangement& arr,
                      PointsIterator points_begin,
                      PointsIterator points_end,
                      OutputIterator out)
{
  typedef typename _Arrangement::Traits_2              Traits;
  typedef typename Traits::X_monotone_curve_2        Base_X_monotone_curve_2;
  typedef typename _Arrangement::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename _Arrangement::Halfedge_const_iterator   Halfedge_const_iterator;

  typedef  Arr_batched_point_location_meta_traits
    <Traits,Halfedge_const_handle>        Meta_traits;
  typedef  typename Meta_traits::X_monotone_curve_2              
    X_monotone_curve_2;

  typedef Arr_batched_point_location_visitor<Meta_traits,
                                             OutputIterator,
                                             _Arrangement>            Visitor;
  typedef Sweep_line_subcurve<Meta_traits>                            Subcurve; 
  typedef Arr_batched_point_location_event<Meta_traits, Subcurve>     Event;

  typedef Sweep_line_2_impl<Meta_traits,
                            Event,
                            Subcurve,
                            Visitor,
                            CGAL_ALLOCATOR(int)>              Sweep_line;
  std::vector<X_monotone_curve_2>      xcurves_vec;
  Traits tr;
  for (Halfedge_const_iterator eit = arr.halfedges_begin();
       eit != arr.halfedges_end();
       ++eit,++eit) 
  {
    if(tr.compare_xy_2_object()((*eit).source().point(),
                     (*eit).target().point()) == LARGER)
      xcurves_vec.push_back(X_monotone_curve_2((*eit).curve(),*eit));
    else
      xcurves_vec.push_back(X_monotone_curve_2((*eit).curve(),(*eit).twin()));
  }
  Visitor visitor_obj(out, arr);
  Sweep_line sweep_line_obj(&visitor_obj);
  sweep_line_obj.init(xcurves_vec.begin(),
                      xcurves_vec.end(),
                      points_begin, 
                      points_end,false);
  sweep_line_obj.sweep();
  return out;
}


CGAL_END_NAMESPACE

#endif
