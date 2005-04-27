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

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 307
#include <CGAL/Nef_2/debug.h>

CGAL_BEGIN_NAMESPACE

class Homogeneous_tag;
class Cartesian_tag;
template<typename Tag> class Normalizing;

template<>
class Normalizing<Homogeneous_tag> {
 public:
  template <typename R> static
  CGAL::Point_3<R> normalized(const CGAL::Point_3<R>& p) { 
    
    typedef typename R::RT     RT;
    
    RT g = p.hw();
    g = (p.hx() == 0 ? g : CGAL_NTS gcd(g,p.hx()));
    g = (p.hy() == 0 ? g : CGAL_NTS gcd(g,p.hy()));
    g = (p.hz() == 0 ? g : CGAL_NTS gcd(g,p.hz()));
    
    RT x = p.hx()/g;
    RT y = p.hy()/g;
    RT z = p.hz()/g;
    RT w = p.hw()/g;
    
    return CGAL::Point_3<R>(x,y,z,w);
  }

  template <typename R> static
  CGAL::Sphere_point<R> normalized(CGAL::Sphere_point<R>& p) {

    typedef typename R::RT     RT;
    
    RT g = (p.x()==0) ? ((p.y()==0) ? ((p.z()==0) ? 1: p.z()): p.y()): p.x();
    
    if(p.y() != 0) g = CGAL_NTS gcd(g,p.y());
    if(p.z() != 0) g = CGAL_NTS gcd(g,p.z());
    
    if(g<0) g = -g;

    RT x = p.x()/g;
    RT y = p.y()/g;
    RT z = p.z()/g;
    
    return CGAL::Sphere_point<R>(x,y,z);
  }

  template <typename R> static
  CGAL::Vector_3<R> normalized(const CGAL::Vector_3<R>& p) {

    typedef typename R::RT     RT;
    
    RT g = (p.hx()==0) ? ((p.hy()==0) ? ((p.hz()==0) ? 1: p.hz()): p.hy()): p.hx();
    
    if(p.hy() != 0) g = CGAL_NTS gcd(g,p.hy());
    if(p.hz() != 0) g = CGAL_NTS gcd(g,p.hz());
    
    if(g<0) g = -g;
    
    RT x = p.hx()/g;
    RT y = p.hy()/g;
    RT z = p.hz()/g;
    
    return CGAL::Vector_3<R>(x,y,z);
  }

  template <typename R> static
  CGAL::Sphere_direction<R> normalized(CGAL::Sphere_direction<R>& c) {

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

  template <typename R> static
  CGAL::Plane_3<R> normalized(CGAL::Plane_3<R>& h) { 
  
    CGAL_assertion(!(h.a()==0 && h.b()==0 && h.c()==0 && h.d()==0));
    
    typedef typename R::RT     RT;
    
    RT x = (h.a()==0) ? ((h.b()==0) ? ((h.c()==0) ? ((h.d()==0) ? 1 
						     : h.d())
				       : h.c())
			 : h.b())
      : h.a();
      
    CGAL_NEF_TRACE("gcd... i"<<' ');
    
    if(h.b() != 0)
      x = CGAL_NTS gcd(x,h.b());
    CGAL_NEF_TRACE(x<<' ');
    if(h.c() != 0)
      x = CGAL_NTS gcd(x,h.c());
    CGAL_NEF_TRACE(x<<' ');
    if(h.d() !=0)
      x = CGAL_NTS gcd(x,h.d());
    CGAL_NEF_TRACEN(x);
    
    x = CGAL_NTS abs(x);
    
    RT pa = h.a()/x;
    RT pb = h.b()/x;
    RT pc = h.c()/x;
    RT pd = h.d()/x;
    
    CGAL_NEF_TRACEN("  after normalizing "  << CGAL::Plane_3<R>(pa,pb,pc,pd));
    return CGAL::Plane_3<R>(pa,pb,pc,pd);
  }

  template <typename R> static
  CGAL::Sphere_circle<R> normalized(CGAL::Sphere_circle<R>& c) { 

    CGAL::Plane_3<R> h = c.plane();
    CGAL_assertion(!(h.a()==0 && h.b()==0 && h.c()==0 && h.d()==0));
    
    typedef typename R::RT     RT;
    
    RT x = (h.a()==0) ? ((h.b()==0) ? ((h.c()==0) ? ((h.d()==0) ? 1 
						     : h.d())
				       : h.c())
			 : h.b())
      : h.a();
    
    if(h.b() != 0)
      x = CGAL_NTS gcd(x,h.b());
    if(h.c() != 0)
      x = CGAL_NTS gcd(x,h.c());
    if(h.d() !=0)
      x = CGAL_NTS gcd(x,h.d());
    
    x = CGAL_NTS abs(x);
    
    RT pa = h.a()/x;
    RT pb = h.b()/x;
    RT pc = h.c()/x;
    RT pd = h.d()/x;
    
    return CGAL::Sphere_circle<R>(CGAL::Plane_3<R>(pa,pb,pc,pd));
  }
};

template<>
class Normalizing<Cartesian_tag> {
 public:
  template <typename R> static
  CGAL::Point_3<R> normalized(const CGAL::Point_3<R>& p) { 
    return p;
  }

  template <typename R> static
  CGAL::Vector_3<R> normalized(const CGAL::Vector_3<R>& p) { 
    return p;
  }

  template <typename R> static
  CGAL::Sphere_point<R> normalized(const CGAL::Sphere_point<R>& p) { 

    typedef typename R::RT     RT;
    
    RT g = (p.hx() != 0 ? p.hx() : (p.hy() != 0 ? p.hy() : p.hz()));
    g = CGAL_NTS abs(g);

    RT x = p.hx()/g;
    RT y = p.hy()/g;
    RT z = p.hz()/g;
    
    return CGAL::Sphere_point<R>(CGAL::Point_3<R>(x,y,z,1));
  }

  template <typename R> static
  CGAL::Sphere_direction<R> normalized(CGAL::Sphere_direction<R>& c) {
    return c;
  }

  template <typename R> static
  CGAL::Plane_3<R> normalized(CGAL::Plane_3<R>& h) { 
   CGAL_assertion(!(h.a()==0 && h.b()==0 && h.c()==0 && h.d()==0));
    
    typedef typename R::FT     FT;
    
    FT x = (h.a()==0) ? ((h.b()==0) ? ((h.c()==0) ? ((h.d()==0) ? 1 
						     : h.d())
				       : h.c())
			 : h.b())
      : h.a();
    x = CGAL_NTS abs(x);
    
    FT pa = h.a()/x;
    FT pb = h.b()/x;
    FT pc = h.c()/x;
    FT pd = h.d()/x;
    
    CGAL_NEF_TRACEN("  after normalizing "  << CGAL::Plane_3<R>(pa,pb,pc,pd));
    return CGAL::Plane_3<R>(pa,pb,pc,pd);
  }

  template <typename R> static
  CGAL::Sphere_circle<R> normalized(CGAL::Sphere_circle<R>& c) { 
    return c;
  }
};

template <typename R>
CGAL::Point_3<R> normalized(const CGAL::Point_3<R>& p) { 
  return Normalizing<typename R::Kernel_tag>::normalized(p);
}

template <typename R>
CGAL::Sphere_point<R> normalized(CGAL::Sphere_point<R>& p) {
  return Normalizing<typename R::Kernel_tag>::normalized(p);
}

template <typename R>
CGAL::Vector_3<R> normalized(const CGAL::Vector_3<R>& p) {
  return Normalizing<typename R::Kernel_tag>::normalized(p);
}

template <typename R>
CGAL::Sphere_direction<R> normalized(CGAL::Sphere_direction<R>& c) {
  return Normalizing<typename R::Kernel_tag>::normalized(c);
}

template <typename R>
CGAL::Plane_3<R> normalized(CGAL::Plane_3<R>& h) { 
  return Normalizing<typename R::Kernel_tag>::normalized(h);
}

template <typename R>
CGAL::Sphere_circle<R> normalized(CGAL::Sphere_circle<R>& c) { 
   return Normalizing<typename R::Kernel_tag>::normalized(c);
}

CGAL_END_NAMESPACE
#endif // CGAL_NORMALIZING_H
