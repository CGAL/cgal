// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/IO/Polyhedron_PS_Stream.h
// package       : PS_Stream
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stephane Postollec
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

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
