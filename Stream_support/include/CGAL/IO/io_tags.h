// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
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

#ifndef CGAL_CFG_NO_CPP0X_LONG_LONG
template<> struct Io_traits<long long> { typedef io_Read_write Io_tag; };
template<> struct Io_traits<unsigned long long> { typedef io_Read_write Io_tag; };
#endif

template<> struct Io_traits<float> { typedef io_Read_write Io_tag; };
template<> struct Io_traits<double> { typedef io_Read_write Io_tag; };
template<> struct Io_traits<long double> { typedef io_Read_write Io_tag; };



} //namespace CGAL

#endif // CGAL_IO_TAGS_H
