// Copyright (c) 2017
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
// Author(s)     : Sylvain Pion,
//                 Mael Rouxel-Labbé

#ifndef CGAL_EXACT_KERNEL_SELECTOR_H
#define CGAL_EXACT_KERNEL_SELECTOR_H

// This class uses the Kernel tag to automatically choose
// whether Cartesian_converter or Homogeneous_converter should be used

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Simple_homogeneous.h>

#include <CGAL/internal/Exact_type_selector.h>

#include <CGAL/representation_tags.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Homogeneous_converter.h>

namespace CGAL {

template <class CK, class Rep = typename CK::Rep_tag /* Cartesian_tag */>
struct Exact_kernel_selector
{
  typedef typename internal::Exact_field_selector<typename CK::RT>::Type  Exact_nt;
  typedef typename internal::Exact_ring_selector <typename CK::RT>::Type  Exact_rt;
  typedef Simple_cartesian<Exact_nt>                                      Exact_kernel;
  typedef Simple_cartesian<Exact_rt>                                      Exact_kernel_rt;

  typedef Cartesian_converter<CK, Exact_kernel>                           C2E;
  typedef Cartesian_converter<CK, Exact_kernel_rt>                        C2E_rt;
  typedef Cartesian_converter<Exact_kernel, CK>                           E2C;
  typedef Cartesian_converter<Exact_kernel_rt, CK>                        E2C_rt;
};

template <class CK>
struct Exact_kernel_selector<CK, Homogeneous_tag>
{
  typedef typename internal::Exact_ring_selector<typename CK::RT>::Type  Exact_nt;
  typedef Simple_homogeneous<Exact_nt>                                   Exact_kernel;
  typedef Exact_kernel                                                   Exact_kernel_rt;

  typedef Homogeneous_converter<CK, Exact_kernel>                        C2E;
  typedef C2E                                                            C2E_rt;
  typedef Homogeneous_converter<Exact_kernel, CK>                        E2C;
  typedef E2C                                                            E2C_rt;
};

} // namespace CGAL

#endif // CGAL_EXACT_KERNEL_SELECTOR_H
