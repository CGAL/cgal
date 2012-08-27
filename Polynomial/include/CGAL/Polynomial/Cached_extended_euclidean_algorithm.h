// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
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
//
// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>
//                 
// ============================================================================

#ifndef CGAL_POLYNOMIAL_CACHED_EXTENDED_EUCLIDEAN_ALGORITHM_H
#define CGAL_POLYNOMIAL_CACHED_EXTENDED_EUCLIDEAN_ALGORITHM_H

#include <CGAL/basic.h>
#include <CGAL/Cache.h>
#include <CGAL/extended_euclidean_algorithm.h>

namespace CGAL {
namespace internal{

template <class UFD, int i = 0 > 
struct Cached_extended_euclidean_algorithm{

  struct Extended_euclidean_algorithm{
    typedef std::pair<UFD,UFD> result_type;
    typedef std::pair<UFD,UFD> first_argument_type; 
    result_type operator()(const first_argument_type& pq){
      result_type result; 
      CGAL::extended_euclidean_algorithm(
          pq.first, pq.second, result.first, result.second);
      return result;
    }
  };
  
  typedef std::pair<UFD,UFD> PAIR; 
  typedef Extended_euclidean_algorithm FUNC;
  typedef CGAL::Cache<PAIR,PAIR,FUNC> CACHE;
  
  static CACHE cache;
  
  void operator()(const UFD& p, const UFD& q, UFD& s, UFD& t){
    PAIR pq(p,q);
    PAIR result = cache(pq);
    s = result.first;
    t = result.second;
  }    
};

template <class UFD, int i> 
typename Cached_extended_euclidean_algorithm<UFD,i>::CACHE 
Cached_extended_euclidean_algorithm<UFD,i>::cache;


} // namespace internal
} //namespace CGAL

#endif//CGAL_POLYNOMIAL_CACHED_EXTENDED_EUCLIDEAN_ALGORITHM_H
