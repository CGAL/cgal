// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Radu Ursu


#ifndef CGAL_QT_WIDGET_CONIC_2_H
#define CGAL_QT_WIDGET_CONIC_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Conic_2.h>
#include <CGAL/Simple_cartesian.h>

namespace CGAL{

template< class R >
Qt_widget&
operator << ( Qt_widget& ws, const CGAL::Conic_2<R>& c)
{
  // pixel dimensions of window
  int dim = (std::max)( ws.width(), ws.height());
  // length of a pixel in window-coordinates
  double pixel_x = (ws.x_max() - ws.x_min())/dim;
  double pixel_y = (ws.y_max() - ws.y_min())/dim;
  // pixel coordinates, stored for faster output
  typedef CGAL::Simple_cartesian<double>::Point_2 Point;
  std::vector<Point> vcoordinates;
  // conic coordinates
  double r = CGAL::to_double (c.r()),
         s = CGAL::to_double (c.s()),
         t = CGAL::to_double (c.t()),
         u = CGAL::to_double (c.u()),
         v = CGAL::to_double (c.v()),
         w = CGAL::to_double (c.w());

  // Phase I (drawing in x-direction)
  // solve conic equation for y
  if (s != 0.0)
      for (double x = ws.x_min(); x <= ws.x_max(); x+=pixel_x) {
          double discr = (t*t-4.0*r*s)*(x*x) + (2.0*t*v-4.0*s*u)*x +
                            v*v - 4.0*s*w;
          if (discr >= 0.0) {
              double y1 = (-t*x - v - CGAL::sqrt(discr))/(2.0*s);
              double y2 = (-t*x - v + CGAL::sqrt(discr))/(2.0*s);
              vcoordinates.push_back(Point(x, y1));
              vcoordinates.push_back(Point(x, y2));} }
  else
      for (double x = ws.x_min(); x <= ws.x_max(); x+=pixel_x) {
          double denom = t*x + v;
          if (denom != 0.0) {
              double y = -(r*x*x + u*x + w)/denom;
              vcoordinates.push_back(Point(x, y)); } }

  // Phase II (drawing in y-direction)
  // solve conic equation for x
  if (r != 0.0)
      for (double y = ws.y_min(); y <= ws.y_max(); y+=pixel_y) {
          double discr = (t*t-4.0*r*s)*(y*y) + (2.0*t*u-4.0*r*v)*y +
                            u*u - 4.0*r*w;
          if (discr >= 0.0) {
              double x1 = (-t*y - u - CGAL::sqrt(discr))/(2.0*r);
              double x2 = (-t*y - u + CGAL::sqrt(discr))/(2.0*r);
              vcoordinates.push_back(Point(x1, y));
              vcoordinates.push_back(Point(x2, y));} }
  else
      for (double y = ws.y_min(); y <= ws.y_max(); y+=pixel_y) {
          double denom = t*y + u;
          if (denom != 0.0) {
              double x = -(s*y*y + v*y + w)/denom;
              vcoordinates.push_back(Point(x, y));} }

  typedef typename std::vector<Point>::const_iterator CIT;
  for(CIT it1 = vcoordinates.begin(); it1!= vcoordinates.end(); ++it1)
    ws << *it1;

  return ws;
}

}//end namespace CGAL

#endif
