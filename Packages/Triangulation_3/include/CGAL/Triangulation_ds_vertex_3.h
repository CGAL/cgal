// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/Triangulation_ds_vertex_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

// vertex of a combinatorial triangulation of any dimension <=3

#ifndef CGAL_TRIANGULATION_DS_VERTEX_3_H
#define CGAL_TRIANGULATION_DS_VERTEX_3_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_3.h>

CGAL_BEGIN_NAMESPACE

template < class Vb >
struct Triangulation_ds_vertex_3 
  : public Vb
{
  bool is_valid(bool verbose = false, int level = 0) const;
};


// Should this be moved to TDS::is_valid(Vertex_handle, verbose, level) ?
template < class Vb >
bool
Triangulation_ds_vertex_3<Vb>::is_valid(bool verbose, int level) const
{
  bool result = Vb::is_valid(verbose,level);
  // result = result && cell()->has_vertex(handle());
  result = result && ( &*cell()->vertex(0) == this ||
                       &*cell()->vertex(1) == this ||
                       &*cell()->vertex(2) == this ||
                       &*cell()->vertex(3) == this );
  if ( ! result ) {
    if ( verbose )
      std::cerr << "invalid vertex" << std::endl;
    CGAL_triangulation_assertion(false);
  }
  return result;
}

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_DS_VERTEX_3_H
