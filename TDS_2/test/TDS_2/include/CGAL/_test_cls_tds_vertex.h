// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// source        :
// file          : include/CGAL/_test_cls_tds_vertex.h
// revision      :
// revision_date :
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <cassert>
namespace CGAL {

template <class Tds>
void
_test_cls_tds_vertex( const Tds&)
{
  std::cout << "    vertex" << std::endl;

  typedef typename Tds::Face_handle          Face_handle;
  typedef typename Tds::Vertex_handle        Vertex_handle;

  Tds tds;
  Face_handle fh1 = tds.create_face();

  // Test constructors
  Vertex_handle vh1=tds.create_vertex();

  // Test  face() and set_face()
  vh1->set_face(fh1);
  assert( vh1->face() == fh1 );

  // The functions degree(), is_valid()
  // and deprecated incident_faces(), incident_vertices(),
  // incident_edges() need a vertex in some
  // triangulation,
  // They are tested in _test_tds__2.h
}

} //namespace CGAL
