// Copyright (c) 1999  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>

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
