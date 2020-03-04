// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ron Wein           <wein@post.tau.ac.il>
//                 (based on old version by Michal Meyerovitch and Ester Ezra)

#ifndef CGAL_ARR_IOSTREAM_H
#define CGAL_ARR_IOSTREAM_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Definition of the I/O operators for the class-template
 * Arrangement_on_surface_2<GeomTraits,TopTraits>.
 */

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/IO/Arr_text_formatter.h>
#include <CGAL/IO/Arrangement_2_writer.h>
#include <CGAL/IO/Arrangement_2_reader.h>
#include <iostream>

namespace CGAL {

/*!
 * Write an arrangement to an output stream using a given formatter.
 * \param arr The arrangement.
 * \param os The output stream.
 * \param format The formatter.
 */
template <class GeomTraits, class TopTraits, class Formatter>
std::ostream&
    write (const Arrangement_on_surface_2<GeomTraits,TopTraits>& arr,
           std::ostream& os,
           Formatter& format)
{
  typedef Arrangement_on_surface_2<GeomTraits,TopTraits>  Arrangement_2;
  typedef Arrangement_2_writer<Arrangement_2>             Arr_writer;

  Arr_writer      writer (arr);

  format.set_out (os);
  writer (format);
  return (os);
}

/*!
 * Output operator (importer).
 * \param os The output stream.
 * \param arr The arrangement.
 */
template <class GeomTraits, class TopTraits>
std::ostream&
    operator<< (std::ostream& os,
                const Arrangement_on_surface_2<GeomTraits,TopTraits>& arr)
{
  typedef Arrangement_on_surface_2<GeomTraits,TopTraits>  Arrangement_2;
  typedef Arrangement_2_writer<Arrangement_2>             Arr_writer;
  typedef Arr_text_formatter<Arrangement_2>               Text_formatter;

  Text_formatter text_format (os);
  Arr_writer      writer (arr);

  writer (text_format);
  return (os);
}

/*!
 * Read an arrangement from an input stream using a given formatter.
 * \param arr The arrangement.
 * \param os The output stream.
 * \param format The formatter.
 */
template <class GeomTraits, class TopTraits, class Formatter>
std::istream&
    read (Arrangement_on_surface_2<GeomTraits,TopTraits>& arr,
          std::istream& is,
          Formatter& format)
{
  typedef Arrangement_on_surface_2<GeomTraits,TopTraits>  Arrangement_2;
  typedef Arrangement_2_reader<Arrangement_2>             Arr_reader;

  Arr_reader      reader(arr);

  format.set_in (is);
  reader (format);
  return (is);
}

/*!
 * Output operator (exporter).
 * \param is The input stream.
 * \param arr The arrangement.
 */
template <class GeomTraits, class TopTraits>
std::istream&
    operator>> (std::istream& is,
                Arrangement_on_surface_2<GeomTraits,TopTraits>& arr)
{
  typedef Arrangement_on_surface_2<GeomTraits,TopTraits>  Arrangement_2;
  typedef Arrangement_2_reader<Arrangement_2>             Arr_reader;
  typedef Arr_text_formatter<Arrangement_2>               Text_formatter;

  Text_formatter text_format (is);
  Arr_reader      reader(arr);

  reader (text_format);
  return (is);
}

} //namespace CGAL

#endif
