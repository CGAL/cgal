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
// file          : include/CGAL/bops_dcel_face.h
// package       : bops (2.2)
// source        : include/CGAL/bops_dcel_face.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL__DCEL_FACE_H
#define CGAL__DCEL_FACE_H

#include <CGAL/bops_dcel_element.h>

CGAL_BEGIN_NAMESPACE

/*
  FACE in the DCEL:
  -----------------
  face_type:        header (HF), color
  face:             typedef const _Dcel_face_type*   _Dcel_face;
  container:        vector<_Dcel_face_type>
  face_iterator:    typedef vector<_Dcel_face_type>::const_iterator
  conversion:       face and face_iterator are type-identical
*/


template<class I>
class _Dcel_face_type : public _Dcel_element_type {
public:
#ifdef CGAL_CFG_INCOMPLETE_TYPE_BUG_4
  typedef _Dcel_edge_type<I>* edge;
#else
  typedef typename I::const_edges_iterator  edge;
#endif

  _Dcel_face_type(int ind, _Dcel_Color col = _NO_COLOR)
        : _Dcel_element_type(ind, col), _header(NULL) {}
  _Dcel_face_type( _Dcel_Color col = _NO_COLOR )
        : _Dcel_element_type(col), _header(NULL) {}

  edge  header() const { return _header; }

protected:
  edge&  header()  { return _header; }

private:
  edge _header;
  friend class _Dcel_base<I>;
};

CGAL_END_NAMESPACE

#endif /* CGAL__DCEL_FACE_H */
