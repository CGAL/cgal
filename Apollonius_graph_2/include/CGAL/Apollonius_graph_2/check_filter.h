// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_CHECK_FILTER_H
#define CGAL_CHECK_FILTER_H

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
