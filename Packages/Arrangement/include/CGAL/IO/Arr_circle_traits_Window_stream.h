// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 13
//
// file          : include/CGAL/IO/Arr_circle_traits_Window_stream.h
// package       : arr (1.03)
// author(s)     : Iddo Hanniel
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifdef CGAL_ARR_CIRCLES_REAL_TRAITS_H
#ifndef CGAL_ARR_CIRCLES_TRAITS_WINDOW_STREAM_H
#define CGAL_ARR_CIRCLES_TRAITS_WINDOW_STREAM_H

#include <CGAL/IO/Window_stream.h>

CGAL_BEGIN_NAMESPACE

#ifndef ARR_IDDO_DEBUG
//a simple version of the windowstream operator (sufficient for X_curve)
template <class NT>
Window_stream& operator<<(Window_stream& os,
			  //const typename 
			  //Arr_circles_real_traits<NT>::Curve &cv)
                          const Circ_Curve<NT>& cv)
{
  //This is not good enough - it assumes s and t have different x coord, 
  //but for x-monotone arcs it is sufficient (and that's all I need).
  //runs faster than above
  double px= CGAL::to_double((cv.source().x()+cv.target().x())/2);
  double R2= CGAL::to_double(cv.circle().squared_radius());
  double sqr = CGAL::sqrt(R2 - 
			  (CGAL::to_double(px-cv.circle().center().x())*
			   CGAL::to_double(px-cv.circle().center().x())));
  
  double py;
  if ((cv.source().x()-cv.target().x()) * cv.circle().orientation() < 0) 
    //underpart
    py= CGAL::to_double(cv.circle().center().y())-sqr;
  else
    py= CGAL::to_double(cv.circle().center().y())+sqr;
  
  
  os.draw_arc(leda_point(CGAL::to_double(cv.source().x()),
                         CGAL::to_double(cv.source().y())),
              leda_point(px,py),
              leda_point(CGAL::to_double(cv.target().x()),
                         CGAL::to_double(cv.target().y())));
  
  return os;
}

#else //ARR_IDDO_DEBUG defined - use the complicated version for general Curves
template <class NT>
Window_stream& operator<<(Window_stream& os,
                          //const Arr_circles_real_traits<NT>::Curve &cv)
                          const Circ_Curve<NT>& cv)
{
  double px,py; //middle point coordinates
  double R2= CGAL::to_double(cv.circle().squared_radius());

  //checking for X-monotone case
  //the folowing is equivelent to "if (curve is x-monotone)"
  if (cv.is_x_monotone()) {
    px= CGAL::to_double((cv.source().x()+cv.target().x()))/2;
    double sqr = CGAL::sqrt(R2 - 
			    (CGAL::to_double(px-cv.circle().center().x())*
			     CGAL::to_double(px-cv.circle().center().x())));
    if (CGAL::sign(cv.source().x()-cv.target().x()) * 
	cv.circle().orientation() < 0) //under part
      py= CGAL::to_double(cv.circle().center().y())-sqr;
    else
      py= CGAL::to_double(cv.circle().center().y())+sqr;
  }
  else { //if not x-monotone the above is not good enough
    if (cv.source()==cv.target()) { //closed circle
      return os << cv.circle() ;
    }
        
    py=CGAL::to_double(cv.circle().center().y());          
    if (CGAL::compare_y(cv.source(),cv.circle().center()) *
	cv.circle().orientation() > 0) {
      //either s is under center and orient is cw or
      //s is above and orient is ccw
      px=CGAL::to_double(cv.circle().center().x())-CGAL::sqrt(R2);
    }
    else
      if (CGAL::compare_y(cv.source(), cv.circle().center()) * 
	  cv.circle().orientation() < 0) {
	//either s is under center and orient is ccw or
	//s is above and orient is cw
	px=CGAL::to_double(cv.circle().center().x())+CGAL::sqrt(R2);
      }
      else 
	{ //s is one of the endpoints of the circle choos other endpoint
	  if (CGAL::compare_x(cv.source(),cv.circle().center())==SMALLER)
	    px=CGAL::to_double(cv.circle().center().x())+CGAL::sqrt(R2);
	  else
	    px=CGAL::to_double(cv.circle().center().x())-CGAL::sqrt(R2);
	}
  }

  os.draw_arc(leda_point(CGAL::to_double(cv.source().x()),
			 CGAL::to_double(cv.source().y())),
	      leda_point(px,py),
	      leda_point(CGAL::to_double(cv.target().x()),
			 CGAL::to_double(cv.target().y())));


  return os;
}
#endif //ARR_IDDO_DEBUG

CGAL_END_NAMESPACE




#endif //CGAL_ARR_CIRCLES_TRAITS_WINDOW_STREAM_H
#endif //CGAL_ARR_CIRCLES_REAL_TRAITS_H




