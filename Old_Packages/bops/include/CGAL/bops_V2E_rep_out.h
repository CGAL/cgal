// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 01
//
// file          : include/CGAL/bops_V2E_rep_out.h
// package       : bops (2.2)
// source        : include/CGAL/bops_V2E_rep_out.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL__V2E_REP_OUT_H
#define CGAL__V2E_REP_OUT_H

CGAL_BEGIN_NAMESPACE 

#ifdef CGAL__V2E_REP_H
template<class vertex, class edge>
ostream& operator<<( ostream& o,
                     const _V2E_rep_base_type<vertex, edge>& v2e )
{
  _V2E_rep_base_type<vertex,edge>::header_const_iterator h_it;
  _V2E_rep_base_type<vertex,edge>::vertex_const_iterator v_it;

  for( h_it= v2e.header_begin(); h_it != v2e.header_end(); h_it++) {
    v_it= (*h_it).vertex();
    o << "this-index= " << (*v_it).index()
      << " size= " << (*h_it).size() << endl;

    for( v_it= v2e.vertex_begin(h_it); v_it != v2e.vertex_end(v_it);
         v_it= v2e.vertex_next(v_it) ) {
      o << "       "
        << "index= " << (*v_it).index()
        << ", vertex= " << *(*v_it).vertex()
        << ", next= " << (*v_it).next()
        << ", T_vertex= " << (*v_it).vertex()
        << ", T_edge= " << (*v_it).edge()
        << endl;
    }
  }

  o << endl << " LIST" << endl;
  int i= 0;
  for( v_it= v2e.vertex_begin(); v_it != v2e.vertex_end(); v_it++) {
      o << i++ << " | "
        << "index= " << (*v_it).index()
        << ", vertex= " << *((*v_it).vertex())
        << ", next= " << (*v_it).next() 
        << endl;
  }
  return o;
}
#endif

CGAL_END_NAMESPACE

#endif
