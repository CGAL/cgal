// Copyright (c) 2010   GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_INTERNAL_INFO_CHECK_H
#define CGAL_INTERNAL_INFO_CHECK_H


#include <boost/mpl/has_xxx.hpp>

namespace CGAL {

namespace internal{

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_typedef_Info,Info,false)

template <class T,bool has_info=Has_typedef_Info<T>::value>
 struct Info_check{
 struct type{};
};

template <class T>
struct Info_check<T,true>{
 typedef typename T::Info type;
};

} } //namespace CGAL::internal

#endif //CGAL_INTERNAL_INFO_CHECK_H
