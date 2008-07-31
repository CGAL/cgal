// ----------------------------------------------------------------------------
//
// Library       : CGAL
// File          : include/CGAL/Arithmetic_kernel.h
// CGAL_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>
//                 
// ============================================================================

#ifndef CGAL_CACHED_EXTENDED_EUCLIDEAN_ALGORITHM_H
#define CGAL_CACHED_EXTENDED_EUCLIDEAN_ALGORITHM_H

#include <CGAL/basic.h>
#include <CGAL/Cache.h>
#include <CGAL/extended_euclidean_algorithm.h>

CGAL_BEGIN_NAMESPACE
namespace CGALi{

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


} // namespace CGALi
CGAL_END_NAMESPACE

#endif//CGAL_CACHED_EXTENDED_EUCLIDEAN_ALGORITHM_H
