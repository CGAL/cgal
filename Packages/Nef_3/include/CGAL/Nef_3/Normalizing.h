// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related docuhmentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : include/CGAL/Nef_3/SNC_constructor.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>
// maintainer    : Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// helper functions for normalizing Plane_3, etc.
// ============================================================================
#ifndef CGAL_NORMALIZING_H
#define CGAL_NORMALIZING_H

#include <CGAL/Nef_3/Infimaximal_box.h>

#include <CGAL/basic.h>
#undef _DEBUG
#define _DEBUG 307
#include <CGAL/Nef_3/debug.h>

template <typename R>
CGAL::Point_3<R> normalized(CGAL::Point_3<R>& p)
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
CGAL::Plane_3<R> normalized(CGAL::Plane_3<R>& h)
{ 

  CGAL_nef3_assertion(!(h.a()==0 && h.b()==0 && h.c()==0 && h.d()==0));
  
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
 
  x = x.abs();

  RT pa = h.a()/x;
  RT pb = h.b()/x;
  RT pc = h.c()/x;
  RT pd = h.d()/x;
  
  TRACEN("  after normalizing "  << CGAL::Plane_3<R>(pa,pb,pc,pd));
  return CGAL::Plane_3<R>(pa,pb,pc,pd);
}


/*
template <typename R>
CGAL::Plane_3<R> normalized_old(CGAL::Plane_3<R>& h)
{ 
  //  TRACEN("  before normalizing "<<h);
  CGAL_nef3_assertion(Infi_box::degree(h.a())==0 && 
		      Infi_box::degree(h.b())==0 && 
		      Infi_box::degree(h.c())==0 && 
		      Infi_box::degree(h.d())<2);
  
  typedef typename R::RT::NT NT;
  typedef typename R::RT     RT;
  NT a(h.a()[0]),b(h.b()[0]),c(h.c()[0]),d(h.d()[0]);
  if(h.d().degree()==1)
    if(d==0)
      d=h.d()[1];
    else
      d=CGAL_NTS gcd(d,h.d()[1]);

  NT x = (a==0) ? ((b==0) ? ((c==0) ? ((d==0) ? 1: d): c): b): a;
  TRACE("gcd... i"<<x<<' ');
  x = ( a != 0 ? a : x);
  TRACE(x<<' ');
  x = ( b != 0 ? CGAL_NTS gcd(x,b) : x );
  TRACE(x<<' ');
  x = ( c != 0 ? CGAL_NTS gcd(x,c) : x );
  TRACE(x<<' ');
  x = ( d != 0 ? CGAL_NTS gcd(x,d) : x );
  TRACEN(x);
  CGAL_nef3_assertion( h == CGAL::Plane_3<R>(a/x,b/x,c/x,d/x));

 
  RT pa = a/x;
  RT pb = b/x;
  RT pc = c/x;
  RT pd = h.d()/x;
  
  return CGAL::Plane_3<R>(pa,pb,pc,pd);
}
*/

#endif // CGAL_NORMALIZING_H
