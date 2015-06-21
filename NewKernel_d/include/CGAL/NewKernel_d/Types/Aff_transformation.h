// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
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
#ifdef CGAL_CXX11
  template<class...T>
  result_type operator()(T&&...)const{return result_type();}
#else
  result_type operator()()const{
    return result_type();
  }
#define CGAL_CODE(Z,N,_) template<BOOST_PP_ENUM_PARAMS(N,class U)> \
  result_type operator()(BOOST_PP_ENUM_BINARY_PARAMS(N,U,const& BOOST_PP_INTERCEPT))const{ \
    return result_type(); \
  }
  BOOST_PP_REPEAT_FROM_TO(1, 9, CGAL_CODE, _ )
#undef CGAL_CODE

#endif
};
}
CGAL_KD_DEFAULT_TYPE(Aff_transformation_tag,(CGAL::Aff_transformation<K>),(),());
CGAL_KD_DEFAULT_FUNCTOR(Construct_ttag<Aff_transformation_tag>,(CartesianDKernelFunctors::Construct_aff_transformation<K>),(Aff_transformation_tag),());

}
#endif
