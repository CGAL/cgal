#ifndef CGAL_POLYNOMIAL_MISC_H
#define CGAL_POLYNOMIAL_MISC_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial/fwd.h>

namespace CGAL{
namespace CGALi{

// template meta function Innermost_coefficient
// returns the tpye of the innermost coefficient 
template <class T> struct Innermost_coefficient{ typedef T Type; };
template <class Coefficient> 
struct Innermost_coefficient<Polynomial<Coefficient> >{
    typedef typename Innermost_coefficient<Coefficient>::Type Type; 
};

// template meta function Dimension
// returns the number of variables 
template <class T> struct Dimension{ static const int value = 0;};
template <class Coefficient> 
struct Dimension<Polynomial<Coefficient> > {
    static const int value = Dimension<Coefficient>::value + 1 ; 
};

} // namespace CGALi
} // namespace CGAL

#endif // CGAL_POLYNOMIAL_MISC_H
