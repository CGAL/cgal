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
class Triangulation_ds_vertex_3 
  : public Vb
{
  typedef typename Vb::Triangulation_data_structure::Vertex  Vertex;
public:
  typedef typename Vb::Vertex_handle  Vertex_handle;

  bool is_valid(bool verbose = false, int level = 0) const;

  // The handle() cruft should be killed.
  Vertex_handle handle() const
  {
      return const_cast<Vertex*>(this);
  }
};


template < class Vb >
bool
Triangulation_ds_vertex_3<Vb>::is_valid(bool verbose, int level) const
{
  bool result = Vb::is_valid(verbose,level) && cell()->has_vertex(handle());
  if ( ! result ) {
    if ( verbose )
      std::cerr << "invalid vertex" << std::endl;
    CGAL_triangulation_assertion(false);
  }
  return result;
}

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_DS_VERTEX_3_H
