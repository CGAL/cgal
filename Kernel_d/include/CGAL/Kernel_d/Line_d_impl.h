// Copyright (c) 2000,2001  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Seel

#ifndef CGAL_LINE_D_C
#define CGAL_LINE_D_C
namespace CGAL {

template <class R> 
Line_d<R> Segment_d<R>::supporting_line() const
{ CGAL_assertion_msg((!is_degenerate()), 
  "Segment_d::supporting_line(): degenerate segment cannot be converted.");
  return Line_d<R>(Base(*this)); 
} 

template <class R>
Line_d<R> Ray_d<R>::supporting_line() const
{ return Line_d<R>(Base(*this)); } 

} //namespace CGAL
#endif //CGAL_LINE_D_C
