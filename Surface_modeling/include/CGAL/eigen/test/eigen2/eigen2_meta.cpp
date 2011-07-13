// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
// Copyright (C) 2008 Gael Guennebaud <g.gael@free.fr>
//
// Eigen is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// Eigen is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// Eigen. If not, see <http://www.gnu.org/licenses/>.

#include "main.h"

void test_eigen2_meta()
{
  typedef float & FloatRef;
  typedef const float & ConstFloatRef;
  
  VERIFY((ei_meta_if<(3<4),ei_meta_true, ei_meta_false>::ret::ret));
  VERIFY(( ei_is_same_type<float,float>::ret));
  VERIFY((!ei_is_same_type<float,double>::ret));
  VERIFY((!ei_is_same_type<float,float&>::ret));
  VERIFY((!ei_is_same_type<float,const float&>::ret));
  
  VERIFY(( ei_is_same_type<float,ei_cleantype<const float&>::type >::ret));
  VERIFY(( ei_is_same_type<float,ei_cleantype<const float*>::type >::ret));
  VERIFY(( ei_is_same_type<float,ei_cleantype<const float*&>::type >::ret));
  VERIFY(( ei_is_same_type<float,ei_cleantype<float**>::type >::ret));
  VERIFY(( ei_is_same_type<float,ei_cleantype<float**&>::type >::ret));
  VERIFY(( ei_is_same_type<float,ei_cleantype<float* const *&>::type >::ret));
  VERIFY(( ei_is_same_type<float,ei_cleantype<float* const>::type >::ret));

  VERIFY(( ei_is_same_type<float*,ei_unconst<const float*>::type >::ret));
  VERIFY(( ei_is_same_type<float&,ei_unconst<const float&>::type >::ret));
  VERIFY(( ei_is_same_type<float&,ei_unconst<ConstFloatRef>::type >::ret));
  
  VERIFY(( ei_is_same_type<float&,ei_unconst<float&>::type >::ret));
  VERIFY(( ei_is_same_type<float,ei_unref<float&>::type >::ret));
  VERIFY(( ei_is_same_type<const float,ei_unref<const float&>::type >::ret));
  VERIFY(( ei_is_same_type<float,ei_unpointer<float*>::type >::ret));
  VERIFY(( ei_is_same_type<const float,ei_unpointer<const float*>::type >::ret));
  VERIFY(( ei_is_same_type<float,ei_unpointer<float* const >::type >::ret));
  
  VERIFY(ei_meta_sqrt<1>::ret == 1);
  #define VERIFY_META_SQRT(X) VERIFY(ei_meta_sqrt<X>::ret == int(ei_sqrt(double(X))))
  VERIFY_META_SQRT(2);
  VERIFY_META_SQRT(3);
  VERIFY_META_SQRT(4);
  VERIFY_META_SQRT(5);
  VERIFY_META_SQRT(6);
  VERIFY_META_SQRT(8);
  VERIFY_META_SQRT(9);
  VERIFY_META_SQRT(15);
  VERIFY_META_SQRT(16);
  VERIFY_META_SQRT(17);
  VERIFY_META_SQRT(255);
  VERIFY_META_SQRT(256);
  VERIFY_META_SQRT(257);
  VERIFY_META_SQRT(1023);
  VERIFY_META_SQRT(1024);
  VERIFY_META_SQRT(1025);
}
