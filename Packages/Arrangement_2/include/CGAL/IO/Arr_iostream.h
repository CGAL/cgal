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
// $Revision$
// $Name$
//
// Author(s)     : Eti Ezra <estere@post.tau.ac.il>

#ifndef CGAL_ARR_IOSTREAM_H
#define CGAL_ARR_IOSTREAM_H

#include <CGAL/basic.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/IO/Arrangement_2_writer.h>
#include <CGAL/IO/Arrangement_2_reader.h>
#include <iostream>

CGAL_BEGIN_NAMESPACE

/*! Exporter operator */
template <class Traits, class Dcel> inline
std::ostream & operator << (std::ostream & os, 
                            const Arrangement_2<Traits,Dcel> & arr)
{
  typedef Arrangement_2<Traits,Dcel>                   Arrangement_2;
  typedef Arrangement_2_writer<Arrangement_2>          Arr_writer;
  typedef Arrangement_2_ascii_formatter<Arrangement_2> Ascii_formatter;

  Ascii_formatter ascii(os);
  Arr_writer writer(arr);
  writer(ascii);
    
  return os;
}

/*! Improter operator */
template <class Traits, class Dcel> inline
std::istream & operator >> (std::istream & is, 
                            Arrangement_2<Traits,Dcel> & arr)
{
  typedef Arrangement_2<Traits,Dcel>                   Arrangement_2;
  typedef Arrangement_2_reader<Arrangement_2>          Arr_reader;
  typedef Arrangement_2_ascii_formatter<Arrangement_2> Ascii_formatter;

  Ascii_formatter ascii(is);
  Arr_reader reader(arr);
  reader(ascii);
  
  return is;
}

CGAL_END_NAMESPACE

#endif
