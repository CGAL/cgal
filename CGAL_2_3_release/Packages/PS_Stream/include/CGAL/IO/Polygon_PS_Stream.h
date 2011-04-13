// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/IO/Polygon_PS_Stream.h
// package       : PS_Stream
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stephane Postollec
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

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
