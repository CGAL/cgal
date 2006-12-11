//Author(s) : Michael Hemmer <mhemmer@uni-mainz.de>

#ifndef CGAL_EXTENDED_EUCLIDEAN_ALGORITHM_H
#define CGAL_EXTENDED_EUCLIDEAN_ALGORITHM_H 1

#include <CGAL/basic.h>

namespace CGAL {

template< class NT > 
NT extended_euclidean_algorithm(const NT& a_, const NT& b_, NT& u, NT& v){

    typedef Algebraic_structure_traits<NT> AST;
    typename AST::Div_mod div_mod;
    typename AST::Unit_part unit_part;
    typename AST::Integral_division idiv;
    
    NT unit_part_a(unit_part(a_));
    NT unit_part_b(unit_part(b_));
    
    NT a(idiv(a_,unit_part_a));
    NT b(idiv(b_,unit_part_b));
  
    
    NT x(0),y(1),last_x(1),last_y(0);
    NT temp, quotient; 
//    typename AST::Div div;
//    typename AST::Mod mod;
  
    //TODO: unroll to avoid swapping 
    while (b != 0){
        temp = b;
        div_mod(a,b,quotient,b);
        a = temp;
        temp = x;
        x = last_x-quotient*x;
        last_x = temp;
        
        temp = y;
        y = last_y-quotient*y;
        last_y = temp;
    }
    u = last_x * unit_part_a;
    v = last_y * unit_part_b;
    
//     std::cout <<"a_: "<<a_ <<std::endl; 
//     std::cout <<"b_: "<<b_ <<std::endl; 
//     std::cout <<"gcd: "<<a <<std::endl; 
//     std::cout <<"u: "<<u <<std::endl; 
//     std::cout <<"v: "<<v <<std::endl; 
//     std::cout <<std::endl; 
    CGAL_precondition(unit_part(a) == NT(1));
    CGAL_precondition(a == a_*u + b_*v);
    return a; 
}

} // namespace CGAL

#endif // CGAL_EXTENDED_EUCLIDEAN_ALGORITHM_H //
