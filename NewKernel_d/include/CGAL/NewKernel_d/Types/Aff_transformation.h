// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_KD_TYPE_AFF_TRANSFORMATION_H
#define CGAL_KD_TYPE_AFF_TRANSFORMATION_H
#include <CGAL/config.h>
#include <CGAL/NewKernel_d/store_kernel.h>
#include <boost/preprocessor/repetition.hpp>

// Dummy, that's all the Kernel_d concept requires, so a useful class will wait.

namespace CGAL {
template<class R_>
struct Aff_transformation {
  typedef R_ R;
};
namespace CartesianDKernelFunctors {
template<class R_> struct Construct_aff_transformation {
  CGAL_FUNCTOR_INIT_IGNORE(Construct_aff_transformation)
  typedef R_ R;
  typedef typename Get_type<R, Aff_transformation_tag>::type result_type;
  template<class...T>
  result_type operator()(T&&...)const{return result_type();}
};
}
CGAL_KD_DEFAULT_TYPE(Aff_transformation_tag,(CGAL::Aff_transformation<K>),(),());
CGAL_KD_DEFAULT_FUNCTOR(Construct_ttag<Aff_transformation_tag>,(CartesianDKernelFunctors::Construct_aff_transformation<K>),(Aff_transformation_tag),());

}
#endif
