// Copyright (c) 2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Stephane Postollec

#ifndef CGAL_IO_POLYGON_PS_STREAM_H
#define CGAL_IO_POLYGON_PS_STREAM_H

#ifdef CGAL_POLYGON_2_H

CGAL_BEGIN_NAMESPACE

template <class _Traits, class _Container>
PS_Stream &operator<<(PS_Stream &ps, const Polygon_2<_Traits,_Container>& p)
{
  typedef Polygon_2<_Traits,_Container>::Vertex_const_iterator VI;
  
  VI i = p.vertices_begin();
  VI end = p.vertices_end();
  ps << point_style(PS_Stream::NONE);
  ps.os() << "/poly {newpath " << std::endl;
  ps << *i;
  ps.os() << "mt" << std::endl;

  do { 
     ps << *i;
     ps.os() << "lt" << std::endl;
     i++;
  } while ( i != end);
 
  ps.os() << "closepath" << std::endl;
  ps.os() << "} def" << std::endl;
  ps.os() << "gsave" << std::endl;
  ps.os() << " poly fill" << std::endl;
  ps.os() << "0 setgray " << std::endl;
  ps.os() << " poly st" << std::endl;
  ps.os() << "grestore" << std::endl;
  return ps;
}

CGAL_END_NAMESPACE

#endif // CGAL_POLYGON_2_H
#endif // CGAL_IO_POLYGON_PS_STREAM_H
