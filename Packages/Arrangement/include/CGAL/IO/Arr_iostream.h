// Copyright (c) 2001  Tel-Aviv University (Israel).
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
// Author(s)     : Eti Ezra <estere@post.tau.ac.il>

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
