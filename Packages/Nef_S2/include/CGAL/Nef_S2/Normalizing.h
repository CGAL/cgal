// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_NORMALIZING_H
#define CGAL_NORMALIZING_H

#include <CGAL/basic.h>
#include <CGAL/assertions.h>
#include <CGAL/Nef_S2/Sphere_point.h>
#include <CGAL/Nef_S2/Sphere_circle.h>
#include <CGAL/Nef_S2/Sphere_direction.h>

#undef _DEBUG
#define _DEBUG 307
#include <CGAL/Nef_3/debug.h>

template <typename R>
CGAL::Point_3<R> normalized(const CGAL::Point_3<R>& p)
{ 
  typedef typename R::RT     RT;

  RT g = p.hw();
  g = (p.hx() == 0 ? g : gcd(g,p.hx()));
  g = (p.hy() == 0 ? g : gcd(g,p.hy()));
  g = (p.hz() == 0 ? g : gcd(g,p.hz()));
  
  RT x = p.hx()/g;
  RT y = p.hy()/g;
  RT z = p.hz()/g;
  RT w = p.hw()/g;
  
  return CGAL::Point_3<R>(x,y,z,w);
}

template <typename R>
CGAL::Sphere_point<R> normalized(CGAL::Sphere_point<R>& p)
{
  typedef typename R::RT     RT;

  RT g = (p.x()==0) ? ((p.y()==0) ? ((p.z()==0) ? 1: p.z()): p.y()): p.x();
  
  if(p.y() != 0) g = gcd(g,p.y());
  if(p.z() != 0) g = gcd(g,p.z());

  if(g<0) g = -g;

  RT x = p.x()/g;
  RT y = p.y()/g;
  RT z = p.z()/g;

  return CGAL::Sphere_point<R>(x,y,z);
}

template <typename R>
CGAL::Vector_3<R> normalized(const CGAL::Vector_3<R>& p)
{
  typedef typename R::RT     RT;

  RT g = (p.hx()==0) ? ((p.hy()==0) ? ((p.hz()==0) ? 1: p.hz()): p.hy()): p.hx();
  
  if(p.hy() != 0) g = gcd(g,p.hy());
  if(p.hz() != 0) g = gcd(g,p.hz());

  if(g<0) g = -g;

  RT x = p.hx()/g;
  RT y = p.hy()/g;
  RT z = p.hz()/g;

  return CGAL::Vector_3<R>(x,y,z);
}

template <typename R>
CGAL::Sphere_direction<R> normalized(CGAL::Sphere_direction<R>& c)
{
  CGAL::Plane_3<R> h = c.plane();
  CGAL_assertion(!(h.a()==0 && h.b()==0 && h.c()==0 && h.d()==0));
  
  typedef typename R::RT     RT;

  RT x = (h.a()==0) ? ((h.b()==0) ? ((h.c()==0) ? ((h.d()==0) ? 1 
						              : h.d())
                                                : h.c())
                                  : h.b())
                    : h.a();

  
  if(h.b() != 0)
    x = gcd(x,h.b());
  if(h.c() != 0)
    x = gcd(x,h.c());
  if(h.d() !=0)
    x = gcd(x,h.d());
 
  x = CGAL_NTS abs(x);

  RT pa = h.a()/x;
  RT pb = h.b()/x;
  RT pc = h.c()/x;
  RT pd = h.d()/x;
  
  return CGAL::Sphere_direction<R>(CGAL::Plane_3<R>(pa,pb,pc,pd));
}

template <typename R>
CGAL::Plane_3<R> normalized(CGAL::Plane_3<R>& h) { 

  CGAL_assertion(!(h.a()==0 && h.b()==0 && h.c()==0 && h.d()==0));
  
  typedef typename R::RT     RT;

  RT x = (h.a()==0) ? ((h.b()==0) ? ((h.c()==0) ? ((h.d()==0) ? 1 
						              : h.d())
                                                : h.c())
                                  : h.b())
                    : h.a();

  
  TRACE("gcd... i"<<' ');
  
  if(h.b() != 0)
    x = gcd(x,h.b());
  TRACE(x<<' ');
  if(h.c() != 0)
    x = gcd(x,h.c());
  TRACE(x<<' ');
  if(h.d() !=0)
    x = gcd(x,h.d());
  TRACEN(x);
 
  x = CGAL_NTS abs(x);

  RT pa = h.a()/x;
  RT pb = h.b()/x;
  RT pc = h.c()/x;
  RT pd = h.d()/x;

  TRACEN("  after normalizing "  << CGAL::Plane_3<R>(pa,pb,pc,pd));
  return CGAL::Plane_3<R>(pa,pb,pc,pd);
}


template <typename R>
CGAL::Sphere_circle<R> normalized(CGAL::Sphere_circle<R>& c)
{ 
  CGAL::Plane_3<R> h = c.plane();
  CGAL_assertion(!(h.a()==0 && h.b()==0 && h.c()==0 && h.d()==0));
  
  typedef typename R::RT     RT;

  RT x = (h.a()==0) ? ((h.b()==0) ? ((h.c()==0) ? ((h.d()==0) ? 1 
						              : h.d())
                                                : h.c())
                                  : h.b())
                    : h.a();

  
  if(h.b() != 0)
    x = gcd(x,h.b());
  if(h.c() != 0)
    x = gcd(x,h.c());
  if(h.d() !=0)
    x = gcd(x,h.d());
 
  x = CGAL_NTS abs(x);

  RT pa = h.a()/x;
  RT pb = h.b()/x;
  RT pc = h.c()/x;
  RT pd = h.d()/x;
  
  return CGAL::Sphere_circle<R>(CGAL::Plane_3<R>(pa,pb,pc,pd));
}

#endif // CGAL_NORMALIZING_H
