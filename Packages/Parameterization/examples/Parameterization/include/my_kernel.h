/***************************************************************************

    begin                : jan 02
    copyright            : (C) 2002 by Pierre Alliez
    email                : pierre.alliez@sophia.inria.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef _MY_KERNEL_
#define _MY_KERNEL_

#include <CGAL/Cartesian.h>
#include <algorithm>
#include <vector>
#include <stdio.h>


typedef CGAL::Cartesian<double> Kernel;

class My_kernel : public Kernel
{
  public:
  static double len(const My_kernel::Vector_3 &u) { return std::sqrt(u*u); }
  static double len(const My_kernel::Vector_2 &u) { return std::sqrt(u*u); }

  static double sqdistance(My_kernel::Point_2 &p1,My_kernel::Point_2 &p2)
  {
    double dx = p1.x()-p2.x();
    double dy = p1.y()-p2.y();
    return dx*dx + dy*dy;
  }
  static double distance(My_kernel::Point_2 &p1,My_kernel::Point_2 &p2)
  {
    return std::sqrt(sqdistance(p1,p2));
  }

  // Angle between two vectors (in degrees)
  // we use this formula
  // uv = |u||v| cos(u,v)
  // u  ^ v  = w
  // |w| = |u||v| |sin(u,v)|
  //**************************************************
  static double evaluate_angle(My_kernel::Vector_3 &u,
                               My_kernel::Vector_3 &v)
  {
    static const double PI = 3.14159265359;
    static const double conv = 1.0/PI*180.0;
    return conv*evaluate_angle_rad(u,v);
  }

  // Angle between two vectors (in radians)
  // we use this formula
  // uv = |u||v| cos(u,v)
  // u  ^ v  = w
  // |w| = |u||v| |sin(u,v)|
  //**************************************************
  static double evaluate_angle_rad(My_kernel::Vector_3 &u,
                                   My_kernel::Vector_3 &v)
  {
    static const double PI = 3.14159265359;

	// check
    double product = len(u)*len(v);
    if(product == 0)
      return 0.0;

    // cosine
    double dot = (u*v);
    double cosine = dot / product;

    // sine
    My_kernel::Vector_3 w = CGAL::cross_product(u,v);
    double AbsSine = len(w) / product;

    if(cosine >= 0)
      return std::asin(fix_sine(AbsSine));
    else
      return PI-std::asin(fix_sine(AbsSine));
  }

  // Angle between two vectors (in rad)
  // uv = |u||v| cos(u,v)
  //**************************************************
  static double angle(My_kernel::Vector_2 &u,
                      My_kernel::Vector_2 &v)
  {
    // check
    double product = len(u)*len(v);
    if(product == 0.0)
      return 0.0;
    double cosine = (u*v) / product;
    bool ccw = ((u.x()*v.y()-u.y()*v.x()) > 0.0) ;
    if(ccw)
      return std::acos(fix_sine(cosine));
    else
      return -std::acos(fix_sine(cosine));
  }
  
  //**********************************************
  // fix sine
  //**********************************************
  static double fix_sine(double sine)
  {
    if(sine >= 1)
      return 1;
    else
      if(sine <= -1)
        return -1;
      else
        return sine;
  }
};


#endif // _MY_KERNEL_

