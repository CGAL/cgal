// Copyright (c) 2019
// GeometryFactory.  All rights reserved. 
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
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_HAS_DATA_H
#define CGAL_HAS_DATA_H
#include <CGAL/Has_member.h>
#include <boost/mpl/if.hpp>
namespace CGAL{
CGAL_GENERATE_MEMBER_DETECTOR(read_data);
CGAL_GENERATE_MEMBER_DETECTOR(write_data);
template <class T>
bool has_extra_data(const T& , typename boost::enable_if<has_read_data<T> >::type* = NULL)
{
  return true;
}

template <class T>
bool has_extra_data(const T& , typename boost::enable_if<boost::mpl::not_<has_read_data<T> > >::type* = NULL)
{
  return false;
}
}
#endif // CGAL_HAS_DATA_H
