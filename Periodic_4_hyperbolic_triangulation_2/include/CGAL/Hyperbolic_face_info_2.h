// Copyright (c) 2010  INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/candidate-packages/Triangulation_2/include/CGAL/Delaunay_triangulation_2.h $
// $Id: Delaunay_triangulation_2.h 57509 2010-07-15 09:14:09Z sloriot $
// 
//
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_DELAUNAY_HYPERBOLIC_FACE_INFO_2_H
#define CGAL_DELAUNAY_HYPERBOLIC_FACE_INFO_2_H

#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <stack>
#include <set>

namespace CGAL {
    
class Hyperbolic_face_info_2
{
public:
  Hyperbolic_face_info_2() : _is_finite_invisible(false), _invisible_edge(UCHAR_MAX)
  {
  }
  
  bool is_finite_invisible() const
  {
    return _is_finite_invisible;
  }
  
  void set_finite_invisible(bool is_finite_invisible)
  {
    _is_finite_invisible = is_finite_invisible;
  }
  
  // Supposed to be called before "get_invisible_edge"
  bool has_invisible_edge() const
  {
    return _invisible_edge <= 2;
  }
  
  // Higly recommended to call "has_invisible_edge" before 
  unsigned char get_invisible_edge() const
  {
    assert(_is_finite_invisible);
    assert(_invisible_edge <= 2);
    
    return _invisible_edge;
  }
  
  void set_invisible_edge(unsigned char invisible_edge)
  {
    assert(_is_finite_invisible);
    assert(invisible_edge <= 2); 
    
    _invisible_edge = invisible_edge;
  }
  
private:
  // a face is invisible if its circumscribing circle intersects the circle at infinity
  bool _is_finite_invisible;
  
  // defined only if the face is finite and invisible
  unsigned char _invisible_edge;
};
  
} //namespace CGAL

#endif // CGAL_DELAUNAY_HYPERBOLIC_FACE_INFO_2_H
