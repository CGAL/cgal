// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri


#ifndef CGAL_IO_TAGS_H
#define CGAL_IO_TAGS_H

#include <cstddef>
#include <CGAL/config.h>

namespace CGAL {

struct io_Read_write{};
struct io_Extract_insert{};
struct io_Operator{};

template<class T>
struct Io_traits{
    typedef io_Operator Io_tag;
};

template<> struct Io_traits<char>{ typedef io_Read_write Io_tag; };

template<> struct Io_traits<short> { typedef io_Read_write Io_tag; };
template<> struct Io_traits<int> { typedef io_Read_write Io_tag; };
template<> struct Io_traits<long> { typedef io_Read_write Io_tag; };

template<> struct Io_traits<unsigned char>{ typedef io_Read_write Io_tag; };

template<> struct Io_traits<unsigned short> { typedef io_Read_write Io_tag; };
template<> struct Io_traits<unsigned int> { typedef io_Read_write Io_tag; };
template<> struct Io_traits<unsigned long> { typedef io_Read_write Io_tag; };

template<> struct Io_traits<long long> { typedef io_Read_write Io_tag; };
template<> struct Io_traits<unsigned long long> { typedef io_Read_write Io_tag; };

template<> struct Io_traits<float> { typedef io_Read_write Io_tag; };
template<> struct Io_traits<double> { typedef io_Read_write Io_tag; };
template<> struct Io_traits<long double> { typedef io_Read_write Io_tag; };



} //namespace CGAL

#endif // CGAL_IO_TAGS_H
