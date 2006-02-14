// ======================================================================
//
// Copyright (c) 2002 SurfLab of CISE of University of Florida
//
// File          : libs/src/cgalExt/Polyhedron_extension.h
// Description   : 
// Creation_date : 21 Feb 2002
// Author(s)     : Le-Jeng Shiue <sle-jeng@cise.ufl.edu>
//
// ======================================================================

// $Id$

/** @file Polyhedron_extension.h
*/

#ifndef _POLYHEDRON_EXTENSION_H_02212002
#define _POLYHEDRON_EXTENSION_H_02212002

#include <SurfLab/config.h>
#include <CGAL/Polyhedron_3.h>

SURFLAB_BEGIN_NAMESPACE

// ======================================================================
///
template <class _P, class _EF, class _EHE, class _EV>
class Polyhedron_extension {
  typedef typename _P::Traits                 Traits;
  typedef typename _P::Items                  Items;
  typedef typename _P::HalfedgeDS             HalfedgeDS;

  typedef typename Items::Face_wrapper<HalfedgeDS, Traits> Face_wrapper;
  typedef typename Items::Halfedge_wrapper<HalfedgeDS,Traits> Halfedge_wrapper;
  typedef typename Items::Vertex_wrapper<HalfedgeDS, Traits> Vertex_wrapper;
  
  typedef typename Face_wrapper::Face         _Face;
  typedef typename Halfedge_wrapper::Halfedge _Halfedge;
  typedef typename Vertex_wrapper::Vertex     _Vertex;

  typedef _EF                                 ExtFType;
  typedef _EHE                                ExtHEType;
  typedef _EV                                 ExtVType;

  template <class _R>
  struct ExtFace : public CGAL::HalfedgeDS_face_base<_R> { ExtFType  ext; };
  template <class _R>
  struct ExtHalfedge:public CGAL::HalfedgeDS_halfedge_base<_R>{ExtHEType ext;};
  //template <class _R>
  //struct ExtVertex : public CGAL::HalfedgeDS_vertex_base<_R>{ExtVType ext;};

  struct ExtItems : public CGAL::Polyhedron_items_3 {
    template <class _R, class _T>
    struct Face_wrapper { 
      typedef ExtFace<_R>                      Face; 
    };
    template <class _R, class _T>
    struct Halfedge_wrapper { 
      typedef ExtHalfedge<_R>                  Halfedge; 
    };
    //template <class _R, class _T>
    //struct Vertex_wrapper { 
    //typedef ExtVertex                        Vertex;
    //};
  };

 public:
  typedef CGAL::Polyhedron_3<Traits, ExtItems> Polyhedron;
};

SURFLAB_END_NAMESPACE

#endif //_POLYHEDRON_EXTENSION_H_02212002
