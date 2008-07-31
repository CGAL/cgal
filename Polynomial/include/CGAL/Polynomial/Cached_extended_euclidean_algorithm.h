
#ifndef CGAL_POLYNOMIAL_CACHED_EEA_H
#define CGAL_POLYNOMIAL_CACHED_EEA_H


#include <CGAL/basic.h>
#include <CGAL/extended_euclidean_algorithm.h>
#include <CGAL/Cache.h>

namespace CGAL{ 
namespace CGALi{

template <class UFD> 
struct Extended_euclidean_algorithm{
    typedef std::pair<UFD,UFD> result_type;
    typedef std::pair<UFD,UFD> first_argument_type; 
    result_type operator()(const first_argument_type& pq){
        result_type result; 
        CGAL::extended_euclidean_algorithm(
                pq.first, pq.second,result.first, result.second);
        return result;
    }
};

template <class UFD> 
struct Cached_extended_euclidean_algorithm{
    typedef std::pair<UFD,UFD> PAIR; 
    typedef Extended_euclidean_algorithm<UFD> FUNC;
    typedef CGAL::Cache<PAIR,PAIR,FUNC,CGAL::Identity<PAIR>, CGAL::Less<PAIR> > CACHE;
    
    static CACHE cache;
    
    void operator()(const UFD& p, const UFD& q, UFD& s, UFD& t){
        PAIR pq(p,q);
        PAIR result = cache(pq);
        s = result.first;
        t = result.second;
    }    
};

template <class UFD> 
typename Cached_extended_euclidean_algorithm<UFD>::CACHE 
Cached_extended_euclidean_algorithm<UFD>::cache;

template <class UFD>
void
cached_extended_euclidean_algorithm(
        const UFD& p, 
        const UFD& q, 
        UFD& s, 
        UFD& t){
    // CGAL::extended_euclidean_algorithm(p,q,s,t); return ;
    Cached_extended_euclidean_algorithm<UFD> eea_cache;
    eea_cache(p,q,s,t);
    
}

} // namespace CGALi
} // namespace CGAL

#endif // CGAL_POLYNOMIAL_CACHED_EEA_H
