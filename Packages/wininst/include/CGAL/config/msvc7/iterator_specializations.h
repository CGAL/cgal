// Copyright (c) 1997-2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : 

#ifndef CGAL_ITER_VC7
#define CGAL_ITER_VC7
#include <iterator>
#include <CGAL/config/msvc7/stl_iterator_base.h>
namespace std {
  
  // af: added the following:
  struct iterator_traits<unsigned int> {
    typedef _Int_iterator_tag iterator_category;
  };
  

template<class C__> inline
typename iterator_traits<C__>::iterator_category
_Iter_cat(const C__&)
  {
    typedef typename iterator_traits<C__>::iterator_category c;
    return c();
  }

template <class _Iter> inline 
typename iterator_traits<_Iter>::difference_type*
  _Dist_type(const _Iter&)
  {
    typedef typename iterator_traits<_Iter>::difference_type _diff_type;
    return static_cast<_diff_type*>(0);
  }


template <class _Iter> inline 
typename iterator_traits<_Iter>::value_type*
  _Val_type(const _Iter&)
{
  typedef typename iterator_traits<_Iter>::value_type _value_type;
  return static_cast<_value_type*>(0);
}


#if defined (__INTEL_COMPILER)
template <class _Tp>
struct iterator_traits<const _Tp*> {
  typedef random_access_iterator_tag iterator_category;
  typedef _Tp                         value_type;
  typedef ptrdiff_t                   difference_type;
  typedef const _Tp*                  pointer;
  typedef const _Tp&                  reference;
};

template <class _Tp>
struct iterator_traits<_Tp*> {
  typedef random_access_iterator_tag iterator_category;
  typedef _Tp                         value_type;
  typedef ptrdiff_t                   difference_type;
  typedef _Tp*                        pointer;
  typedef _Tp&                        reference;
};
#endif // defined (__INTEL_COMPILER)

}
#endif // CGAL_ITER_VC7
