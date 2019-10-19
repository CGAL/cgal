// Copyright (C) 2018  GeometryFactory Sarl
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//

#ifndef CGAL_CARTESIAN_CONVERTER_FWD_H
#define CGAL_CARTESIAN_CONVERTER_FWD_H

/// \file Cartesian_converter_fwd.h
/// Forward declarations of the `Cartesian_converter` class.

#ifndef DOXYGEN_RUNNING
namespace CGAL {

namespace internal {
template < typename K1, typename K2 >
struct Default_converter;
}//internal

template < class K1, class K2,
           class Converter = typename internal::Default_converter<K1, K2>::Type >
class Cartesian_converter;

} // CGAL
#endif

#endif /* CGAL_CARTESIAN_CONVERTER_FWD_H */


