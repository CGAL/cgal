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
// release       : $CGAL_Revision: CGAL-2.3-I-62 $
// release_date  : $CGAL_Date: 2001/05/11 $
//
// file          : include/CGAL/IO/Arr_iostream.h
// package       : Arrangement (1.82)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// author(s)     : Eti Ezra <estere@post.tau.ac.il>
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================

#ifndef CGAL_ARR_IOSTREAM_H
#define CGAL_ARR_IOSTREAM_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#ifndef CGAL_ARRANGEMENT_2_H
#include <CGAL/Arrangement_2.h>
#endif

#ifndef CGAL_IO_ARR_FILE_WRITER_H
#include <CGAL/IO/Arr_file_writer.h>
#endif // CGAL_IO_ARR_FILE_WRITER_H

#ifndef CGAL_IO_WRITE_ARR_H
#include <CGAL/IO/write_arr.h>
#endif // CGAL_IO_WRITE_ARR_H


#include <iostream>

CGAL_BEGIN_NAMESPACE

template <class Dcel, class Traits, class Base_node> inline
::std::ostream& operator << (::std::ostream& o, 
                             const Arrangement_2<Dcel,Traits, 
                             Base_node>& arr) 
{
  //print_OFF(o, arr);
  
  Arr_file_writer< Arrangement_2<Dcel,Traits, Base_node> >  writer(o, arr);
  
  write_arr(arr, writer, o);
  
  return o;
}
 
template <class Dcel, class Traits, class Base_node> inline
::std::istream& operator>>( std::istream& in, 
                            Arrangement_2<Dcel,Traits, Base_node>& arr) {
  // reads a polyhedron from `in' and appends it to P.
  
  arr.read(in);
  
  return in;
}


CGAL_END_NAMESPACE

#endif








