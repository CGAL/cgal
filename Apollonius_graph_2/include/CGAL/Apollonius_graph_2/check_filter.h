// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_CHECK_FILTER_H
#define CGAL_CHECK_FILTER_H

#include <CGAL/license/Apollonius_graph_2.h>


#undef CGAL_IA_NEW_FILTERS

namespace CGAL {

template < class T>
void must_be_filtered(const T&)
{}

#if defined CGAL_ARITHMETIC_FILTER_H
template < class CT, class ET, class Type, bool Protection, class Cache>
void must_be_filtered(const Filtered_exact<CT, ET, Type, Protection,
		      Cache> &)
{ dont_compile(CT(), ET()); }
#endif

}

#endif
