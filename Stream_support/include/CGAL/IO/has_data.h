// Copyright (c) 2019
// GeometryFactory.  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
