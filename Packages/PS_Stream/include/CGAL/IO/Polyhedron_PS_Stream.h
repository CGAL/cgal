// Copyright (c) 2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Stephane Postollec

#ifndef CGAL_IO_POLYHEDRON_OS_STREAM_H
#define CGAL_IO_POLYHEDRON_OS_STREAM_H

#ifdef CGAL_POLYHEDRON_3_H

CGAL_BEGIN_NAMESPACE

template <class Traits, class HDS>
PS_Stream_3& operator <<(PS_Stream_3& ps, const Polyhedron_3<Traits,HDS>& P)
{
  typedef Polyhedron_3<Traits,HDS>::Facet_const_iterator FCI;
  typedef
        Polyhedron_3<Traits,HDS>::Halfedge_around_facet_const_circulator HFCC;

  FCI fi = P.facets_begin();
  for(; fi != P.facets_end();++fi)
  {
    Polygon p;
    double scal =ps.compute_plane_equations(*fi).orthogonal_vector() *
                 ps.direction().vector();
   
    std::cout << scal << std::endl;

    if (scal > 0)
      {
	double gris = scal/(ps.norme(
                ps.compute_plane_equations(*fi).orthogonal_vector().x(),
                ps.compute_plane_equations(*fi).orthogonal_vector().y(),
                ps.compute_plane_equations(*fi).orthogonal_vector().z()) *
                            ps.norme(ps.direction().vector().x(),
                                     ps.direction().vector().y(),
                                     ps.direction().vector().z()));

	HFCC hfc = fi->facet_begin();
	HFCC hfc_end = hfc;
	do {
	    p.push_back(ps.transform(ps.trans(),(hfc->vertex())->point()));
	    ++hfc;
	   }while(hfc != hfc_end);
	
        ps.os() << gris <<" " <<"setgray" <<std::endl;
	ps <<p;
	
     }
  }
  return ps;
} 

CGAL_END_NAMESPACE

#endif // CGAL_POLYHEDRON_3_H

#endif // CGAL_IO_POLYHEDRON_OS_STREAM_H
