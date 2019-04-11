// Copyright (c) 2019 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_VERBOSITY_H
#define CGAL_KSR_VERBOSITY_H

#include <ostream>

// General verbosity

#ifndef CGAL_KSR_VERBOSE_LEVEL
#define CGAL_KSR_VERBOSE_LEVEL 0
#endif

namespace CGAL
{
namespace KSR
{

struct Null_stream : public std::basic_ostream<char> { };

template <typename Type>
Null_stream& operator<< (Null_stream& null_stream, const Type&)
{
  return null_stream;
}

Null_stream& null_stream()
{
  static Null_stream null_stream;
  return null_stream;
}

} // namespace KSR

} // namespace CGAL

#define CGAL_KSR_NULL_STREAM CGAL::KSR::null_stream()
#define CGAL_KSR_CERR_STREAM std::cerr

#if CGAL_KSR_VERBOSE_LEVEL > 0
#  define CGAL_KSR_CERR_1 CGAL_KSR_CERR_STREAM
#else
#  define CGAL_KSR_CERR_1 CGAL_KSR_NULL_STREAM
#endif

#if CGAL_KSR_VERBOSE_LEVEL > 1
#  define CGAL_KSR_CERR_2 CGAL_KSR_CERR_STREAM
#else
#  define CGAL_KSR_CERR_2 CGAL_KSR_NULL_STREAM
#endif

#if CGAL_KSR_VERBOSE_LEVEL > 2
#  define CGAL_KSR_CERR_3 CGAL_KSR_CERR_STREAM
#else
#  define CGAL_KSR_CERR_3 CGAL_KSR_NULL_STREAM
#endif

#if CGAL_KSR_VERBOSE_LEVEL > 3
#  define CGAL_KSR_CERR_4 CGAL_KSR_CERR_STREAM
#else
#  define CGAL_KSR_CERR_4 CGAL_KSR_NULL_STREAM
#endif


#endif // CGAL_KSR_SILENT
