// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_NORMALIZING_H
#define CGAL_NORMALIZING_H

#include <CGAL/license/Nef_S2.h>


#include <CGAL/Nef_S2/Sphere_point.h>
#include <CGAL/Nef_S2/Sphere_circle.h>
#include <CGAL/Nef_S2/Sphere_direction.h>
#include <CGAL/Fraction_traits.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 307
#include <CGAL/Nef_2/debug.h>

#ifdef CCGAL_USE_LEDA
#include <CGAL/Cartesian.h>
#include <CGAL/leda_rational.h>
#endif

namespace CGAL {

class Homogeneous_tag;
class Cartesian_tag;
template<typename Tag> class Normalizing;

template<>
class Normalizing<Homogeneous_tag> {
 public:

  template <typename iterator> static
  void normalized(iterator begin, iterator end) { 
    typedef typename std::iterator_traits<iterator>::value_type RT;

    iterator i = begin;
    while(i!=end && *i == RT(0)) ++i;
    if(i==end)
      return;
    
    RT g = *i;
    for(iterator j=i+1; j!=end; ++j)
      g = (*j == 0 ? g : CGAL_NTS gcd(g,*j)); 
    g=CGAL_NTS abs(g);

    for(; i!=end; ++i)
      *i = CGAL::integral_division(*i,g);
  }

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
    
    return typename R::Point_3(x,y,z,w);
  }

  template <typename R> static
  CGAL::Sphere_point<R> normalized(const CGAL::Sphere_point<R>& p) {

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
    
    return typename R::Vector_3(x,y,z);
  }

  template <typename R> static
  CGAL::Sphere_direction<R> normalized(const CGAL::Sphere_direction<R>& c) {

    typename R::Plane_3 h = c.plane();
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
    
    return CGAL::Sphere_direction<R>(typename R::Plane_3(pa,pb,pc,pd));
  }

  template <typename R> static
  CGAL::Plane_3<R> normalized(const CGAL::Plane_3<R>& h) { 
  
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
    
    CGAL_NEF_TRACEN("  after normalizing "  << typename R::Plane_3(pa,pb,pc,pd));
    return typename R::Plane_3(pa,pb,pc,pd);
  }

  template <typename R> static
  CGAL::Sphere_circle<R> normalized(const CGAL::Sphere_circle<R>& c) { 

    typename R::Plane_3 h = c.plane();
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
    
    return CGAL::Sphere_circle<R>(typename R::Plane_3(pa,pb,pc,pd));
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
    
    return CGAL::Sphere_point<R>(x,y,z);
  }

  template <typename R> static
  CGAL::Sphere_direction<R> normalized(const CGAL::Sphere_direction<R>& c) {
    return c;
  }

#ifdef CCGAL_USE_LEDA
// specialization: Plane_3 < Cartesian < leda_rational > >

  
  static Plane_3<CGAL::Cartesian<leda_rational> > 
       normalized(Plane_3<CGAL::Cartesian<leda_rational> >& h) { 

    CGAL_assertion(!(h.a()==0 && h.b()==0 && h.c()==0 && h.d()==0));
    
    typedef leda_rational     FT;
    
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

    pa.normalize();
    pb.normalize();
    pc.normalize();
    pd.normalize();
    
    CGAL_NEF_TRACEN("  after normalizing "  << CGAL::Plane_3<CGAL::Cartesian<leda_rational> >(pa,pb,pc,pd));
    return CGAL::Plane_3<CGAL::Cartesian<leda_rational> >(pa,pb,pc,pd);
  }

#endif
  
  template <typename R> static
  CGAL::Plane_3<R> normalized(const CGAL::Plane_3<R>& h,Tag_true) {
    CGAL_assertion(!(h.a()==0 && h.b()==0 && h.c()==0 && h.d()==0));

    typedef typename R::FT FT;
    typedef Fraction_traits<FT> FracTraits;
    typedef std::vector<typename FracTraits::Numerator_type> NV;

    typename FracTraits::Numerator_type num;
    typename FracTraits::Denominator_type denom;
    typename FracTraits::Decompose decomposer; 
    typename FracTraits::Compose composer; 
    NV vec;
    
    decomposer(h.a(),num,denom);
    vec.push_back(num);
    vec.push_back(denom);
    vec.push_back(denom);
    vec.push_back(denom);
    decomposer(h.b(),num,denom);
    vec[0]*=denom;
    vec[1]*=num;
    vec[2]*=denom;
    vec[3]*=denom;
    decomposer(h.c(),num,denom);
    vec[0]*=denom;
    vec[1]*=denom;
    vec[2]*=num;
    vec[3]*=denom;
    decomposer(h.d(),num,denom);
    vec[0]*=denom;
    vec[1]*=denom;
    vec[2]*=denom;
    vec[3]*=num;
    
    Normalizing<Homogeneous_tag>::
      normalized(vec.begin(),vec.end());
    return typename R::Plane_3(composer(vec[0],1),
			       composer(vec[1],1),
			       composer(vec[2],1),
			       composer(vec[3],1));
  }
  
  
  
  template <typename R> static
  CGAL::Plane_3<R> normalized(const CGAL::Plane_3<R>& h,Tag_false) {   
    CGAL_assertion(!(h.a()==0 && h.b()==0 && h.c()==0 && h.d()==0));
    typedef typename R::FT FT;
    if (h.a()!=0)
      return typename R::Plane_3(FT(1),h.b()/h.a(),h.c()/h.a(),h.d()/h.a());
    if (h.b()!=0)
      return typename R::Plane_3(h.a()/h.b(),FT(1),h.c()/h.b(),h.d()/h.b());
    if (h.c()!=0)
      return typename R::Plane_3(h.a()/h.c(),h.b()/h.c(),FT(1),h.d()/h.c());
    return typename R::Plane_3(h.a()/h.d(),h.b()/h.d(),h.c()/h.d(),FT(1));
  }

  template <typename R> static
  CGAL::Plane_3<R> normalized(const CGAL::Plane_3<R>& h) {   
    return normalized( h,typename Fraction_traits<typename R::FT>::Is_fraction() );
  }

  template <typename R> static
  CGAL::Sphere_circle<R> normalized(CGAL::Sphere_circle<R>& c) { 
    return c;
  }
};

/*
template <typename R, typename iterator>
void normalized(iterator begin, iterator end) {
  Normalizing<typename R::Kernel_tag>::normalized(begin, end);
}
*/

template <typename R>
CGAL::Point_3<R> normalized(const CGAL::Point_3<R>& p) { 
  return Normalizing<typename R::Kernel_tag>::normalized(p);
}

template <typename R>
CGAL::Sphere_point<R> normalized(const CGAL::Sphere_point<R>& p) {
  return Normalizing<typename R::Kernel_tag>::normalized(p);
}

template <typename R>
CGAL::Vector_3<R> normalized(const CGAL::Vector_3<R>& p) {
  return Normalizing<typename R::Kernel_tag>::normalized(p);
}

template <typename R>
CGAL::Sphere_direction<R> normalized(const CGAL::Sphere_direction<R>& c) {
  return Normalizing<typename R::Kernel_tag>::normalized(c);
}

template <typename R>
CGAL::Plane_3<R> normalized(const CGAL::Plane_3<R>& h) { 
  return Normalizing<typename R::Kernel_tag>::normalized(h);
}

template <typename R>
CGAL::Sphere_circle<R> normalized(const CGAL::Sphere_circle<R>& c) { 
   return Normalizing<typename R::Kernel_tag>::normalized(c);
}

} //namespace CGAL
#endif // CGAL_NORMALIZING_H
