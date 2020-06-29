// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany),
// INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>



#ifndef CGAL_CONVERT_TO_BFI_H
#define CGAL_CONVERT_TO_BFI_H

#include <CGAL/basic.h>
#include <CGAL/Get_arithmetic_kernel.h>
#include <CGAL/Cache.h>

namespace CGAL {

template <class NTX>
typename Get_arithmetic_kernel<NTX>::Arithmetic_kernel::Bigfloat_interval
convert_to_bfi(const NTX& x) {
    typedef typename Get_arithmetic_kernel<NTX>::Arithmetic_kernel AK;
    typedef typename AK::Bigfloat_interval BFI;
    typedef CGAL::Coercion_traits<NTX,BFI> CT;
    return typename CT::Cast()(x);
}

} //namespace CGAL

#endif // CGAL_CONVERT_TO_BFI_H
