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

#include <CGAL/basic.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/IO/Arr_file_writer.h>
#include <CGAL/IO/write_arr.h>
#include <iostream>

CGAL_BEGIN_NAMESPACE

template <class Dcel, class Traits, class Base_node> inline
::std::ostream & operator << (::std::ostream & o, 
                              const Arrangement_2<Dcel,Traits,Base_node> & arr)
{
  typedef Arrangement_2<Dcel,Traits,Base_node>        Arr_2;
  typedef Arr_file_writer<Arr_2>                      Writer;

  //print_OFF(o, arr);
  
  Writer writer(o, arr);
  write_arr<Arr_2,Writer>(arr, writer, o);
  return o;
}

template <class Dcel, class Traits, class Base_node> inline
::std::istream & operator >> (std::istream & in, 
                              Arrangement_2<Dcel,Traits, Base_node> & arr)
{
  // reads a polyhedron from `in' and appends it to P.
  arr.read(in);
  return in;
}

CGAL_END_NAMESPACE

#endif
