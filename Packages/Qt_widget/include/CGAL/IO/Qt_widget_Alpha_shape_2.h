#ifndef CGAL_QT_WIDGET_ALPHA_SHAPE_2_H
#define CGAL_QT_WIDGET_ALPHA_SHAPE_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Alpha_shape_2.h>
namespace CGAL{

template< class Dt >
Qt_widget&
operator << ( Qt_widget& ws, const CGAL::Alpha_shape_2<Dt>& As)
{
  typedef typename CGAL::Alpha_shape_2<Dt>::Interval_edge_map 
	  Interval_edge_map;
  typename Interval_edge_map::const_iterator edge_alpha_it;
  const typename CGAL::Alpha_shape_2<Dt>::Interval3* pInterval;
  if (As.get_mode() == Alpha_shape_2<Dt>::REGULARIZED) 
  {
      // it is much faster looking at the sorted intervals 
      // than looking at all sorted faces
      // alpha must be larger than the mid boundary
      // and alpha is smaller than the upper boundary
    /*
      for (edge_alpha_it = As._interval_edge_map.begin(); 
	   edge_alpha_it != As._interval_edge_map.end() &&
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
	    ws << segment((*edge_alpha_it).second.first,
		   (*edge_alpha_it).second.second);

	    // to debug the edge descrition...
	    //ws << Segment((*edge_alpha_it).second.first->vertex(0)->point(),
	    //	    (*edge_alpha_it).second.first->vertex(1)->point());
	    //ws << Segment((*edge_alpha_it).second.first->vertex(1)->point(),
	    //      (*edge_alpha_it).second.first->vertex(2)->point());
	    //ws << Segment((*edge_alpha_it).second.first->vertex(2)->point(),
	    //	    (*edge_alpha_it).second.first->vertex(0)->point());

	    }//endif
	}//endfor

  } else {
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
	    ws << segment((*edge_alpha_it).second.first,
			       (*edge_alpha_it).second.second);
	  }
	} else {
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
	  }//endif
	}//endif
      }//endfor
      */
  }//endif
  return ws;
}

}//end namespace CGAL

#endif
