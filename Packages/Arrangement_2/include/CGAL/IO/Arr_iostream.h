// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Michal Meyerovitch <gorgymic@post.tau.ac.il>
//                 Ron Wein           <wein@post.tau.ac.il>
//                 (based on old version by Ester Ezra)

#ifndef CGAL_ARR_IOSTREAM_H
#define CGAL_ARR_IOSTREAM_H

/*! \file
 * Definition of the I/O operators for the Arrangement_2<Traits,Dcel> class.
 */

#include <CGAL/basic.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/IO/Arrangement_2_writer.h>
#include <CGAL/IO/Arrangement_2_reader.h>
#include <iostream>

CGAL_BEGIN_NAMESPACE

/*!
 * Output operator (importer).
 * \param os The output stream.
 * \param arr The arrangement.
 */
template <class Traits, class Dcel>
std::ostream& operator<< (std::ostream & os, 
                          const Arrangement_2<Traits,Dcel>& arr)
{
  typedef Arrangement_2<Traits,Dcel>                    Arrangement_2;
  typedef Arrangement_2_writer<Arrangement_2>           Arr_writer;
  typedef Arrangement_2_ascii_formatter<Arrangement_2>  Ascii_formatter;

  Ascii_formatter ascii_format (os);
  Arr_writer      writer (arr);

  writer (ascii_format);
  return (os);
}

/*!
 * Output operator (exporter).
 * \param is The input stream.
 * \param arr The arrangement.
 */
template <class Traits, class Dcel>
std::istream& operator>> (std::istream& is, 
                          Arrangement_2<Traits,Dcel>& arr)
{
  typedef Arrangement_2<Traits,Dcel>                    Arrangement_2;
  typedef Arrangement_2_reader<Arrangement_2>           Arr_reader;
  typedef Arrangement_2_ascii_formatter<Arrangement_2>  Ascii_formatter;

  Ascii_formatter ascii_format (is);
  Arr_reader      reader(arr);
  
  reader (ascii_format);
  return (is);
}

CGAL_END_NAMESPACE

#endif
