// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>

#ifdef CGAL_ALPHA_WINDOW_STREAM

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template < class Dt >
Window_stream&
Alpha_shape_2<Dt>::op_window(Window_stream& W) const
{

  typedef  typename Alpha_shape_2<Dt>::Interval_edge_map 
    Interval_edge_map;
  typename Interval_edge_map::const_iterator edge_alpha_it;

  const typename Alpha_shape_2<Dt>::Interval3* pInterval;

  if (get_mode() == Alpha_shape_2<Dt>::REGULARIZED) 
    {

      // it is much faster looking at the sorted intervals 
      // than looking at all sorted faces
      // alpha must be larger than the mid boundary
      // and alpha is smaller than the upper boundary
      for (edge_alpha_it = _interval_edge_map.begin(); 
	   edge_alpha_it != _interval_edge_map.end() &&
	     (*edge_alpha_it).first.first < get_alpha();
	   ++edge_alpha_it) 
	{

	  pInterval = &(*edge_alpha_it).first;

	  CGAL_triangulation_assertion(pInterval->second != Infinity);
	  // since this happens only for convex hull of dimension 1
	  // thus singular

	  if(pInterval->second < get_alpha() &&
	     (pInterval->third >= get_alpha()
	      || pInterval->third == Infinity)) 
	    {
	      // alpha must be larger than the mid boundary
	      // and alpha is smaller than the upper boundary
	      // which might be infinity 
	      // visualize the boundary
	    
 CGAL_triangulation_assertion((classify((*edge_alpha_it).second.first,
					(*edge_alpha_it).second.second) ==
			       Alpha_shape_2<Dt>::REGULAR));
	      // if we used Edelsbrunner and Muecke's definition
	      // regular means incident to a higher-dimensional face
	      // thus we would write to many vertices
	      W << segment((*edge_alpha_it).second.first,
			   (*edge_alpha_it).second.second);

	      // to debug the edge descrition...
// 	      W << Segment((*edge_alpha_it).second.first->vertex(0)->point(),
// 			   (*edge_alpha_it).second.first->vertex(1)->point());
// 	      W << Segment((*edge_alpha_it).second.first->vertex(1)->point(),
// 			   (*edge_alpha_it).second.first->vertex(2)->point());
// 	      W << Segment((*edge_alpha_it).second.first->vertex(2)->point(),
// 			   (*edge_alpha_it).second.first->vertex(0)->point());

	    }
	}
    }
  else 
    { // get_mode() == GENERAL

      // draw the edges
      for (edge_alpha_it = _interval_edge_map.begin(); 
	   edge_alpha_it != _interval_edge_map.end() &&
	     (*edge_alpha_it).first.first < get_alpha();
	   ++edge_alpha_it) 
	{
	
	  pInterval = &(*edge_alpha_it).first;

	  if (pInterval->first == UNDEFINED) 
	    {
	    
	      CGAL_triangulation_assertion(pInterval->second != Infinity);
	      // since this happens only for convex hull of dimension 1
	      // thus singular

	      if(pInterval->second < get_alpha() &&
		 (pInterval->third >= get_alpha()
		  || pInterval->third == Infinity)) 
		{
		  // alpha must be larger than the mid boundary
		  // and alpha is smaller than the upper boundary
		  // which might be infinity 
		  // visualize the boundary
		
 CGAL_triangulation_assertion((classify((*edge_alpha_it).second.first,
					(*edge_alpha_it).second.second) ==
			       Alpha_shape_2<Dt>::REGULAR));
		  W << segment((*edge_alpha_it).second.first,
			       (*edge_alpha_it).second.second);
		}
	    }
	  else 
	    {
	   

	      if(pInterval->third >= get_alpha()
		 || pInterval->third == Infinity) 
		{
		  // if alpha is smaller than the upper boundary
		  // which might be infinity 
		  // visualize the boundary
		
 CGAL_triangulation_assertion(((classify((*edge_alpha_it).second.first,
					 (*edge_alpha_it).second.second) ==
				Alpha_shape_2<Dt>::REGULAR) || 
			       (classify((*edge_alpha_it).second.first,
					 (*edge_alpha_it).second.second) ==
				Alpha_shape_2<Dt>::SINGULAR)));
		  W << segment((*edge_alpha_it).second.first,
			       (*edge_alpha_it).second.second);
		}
	    }
	}
    }

  return W;
}

//-------------------------------------------------------------------

template < class Dt >
Window_stream& 
operator<<(Window_stream& W, const Alpha_shape_2<Dt>& As) 
{
  return As.op_window(W);
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif//CGAL_ALPHA_WINDOW_STREAM
