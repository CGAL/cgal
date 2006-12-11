//Author(s) : Michael Hemmer <mhemmer@uni-mainz.de>

#ifndef CGAL_CHINESE_REMAINDER_H
#define CGAL_CHINESE_REMAINDER_H 1

#include <CGAL/basic.h>
#include <CGAL/extended_euclidean_algorithm.h>

namespace CGAL {

// this is just the version for 'integers'
// NT must be model of RealEmbeddable
// NT must be model of EuclideanRing
template <class NT>
void chinese_remainder(
        NT  m1, NT  u1,
        NT  m2, NT  u2,
        NT& m , NT& u  ){
    typedef Algebraic_structure_traits<NT> AST;
    typename AST::Mod mod;
    //typename AST::Unit_part unit_part;
    typename AST::Integral_division idiv; 
    
    if(u1 < NT(0) ) u1 += m1;
    if(u2 < NT(0) ) u2 += m2;
    
    CGAL_precondition(0  < m1);
    CGAL_precondition(u1 < m1);
    CGAL_precondition(u1 >= NT(0));
        
    CGAL_precondition(0  < m2);
    CGAL_precondition(u2 < m2);
    CGAL_precondition(u2 >= NT(0));
    
    
    NT tmp,c,dummy;
    tmp = CGAL::extended_euclidean_algorithm(m1,m2,c,dummy);
    CGAL_postcondition(tmp == NT(1));
    CGAL_postcondition(m1*c + m2*dummy == NT(1));
    
    
    m = m1*m2;
    NT v = mod(c*(u2-u1),m2);
    u = m1*v + u1;
    
    // u is not unique yet!
    NT m_half = idiv(m-mod(m,NT(2)),NT(2));
    if (u  >  m_half) u -= m ;
    if (u <= -m_half) u += m ; 

}

}///namespace CGAL

#endif //#ifnedef CGAL_CHINESE_REMAINDER_H 1
 
