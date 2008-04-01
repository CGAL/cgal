#include <CGAL/basic.h>

#ifdef CGAL_USE_GMPXX

#include <iostream>
#include <CGAL/mpz_class.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>


template <class A, class B>
inline
typename CGAL::Algebraic_structure_traits<
typename CGAL::Coercion_traits<A,B>::Type
>::Integral_division::result_type
integral_division(const A& a, const B& b){
    typedef CGAL::Coercion_traits<A,B> CT;
    typedef typename CT::Type Type;
    typename CGAL::Algebraic_structure_traits<
    typename CGAL::Coercion_traits<A,B>::Type
        >::Integral_division integral_division;
    return integral_division(a,b);
}

template <class A, class B>
inline
int
test_coercion(const A& a, const B& b){
    std::cout << "START TEST" << std::endl;
    typedef CGAL::Coercion_traits<A,B> CT;
    typedef typename CT::Type Type;
    typename CT::Cast cast;
    Type x = cast(a);

    typename CGAL::Algebraic_structure_traits<
    typename CGAL::Coercion_traits<A,B>::Type
        >::Integral_division integral_division;
    
    return 1;
    //return integral_division(a,b);
}




int main() {
    
//  mpz_class x(1);
//  test_coercion(x,x);
//  std::cout << test_coercion(x,x) << std::endl;
//  std::cout << test_coercion(x+x,x*x) << std::endl;
//  std::cout << integral_division(x+x,x*x) << std::endl;
//  std::cout << CGAL::integral_division(x+x,x*x) << std::endl;
    
    {
        typedef mpz_class NT;
        typedef CGAL::Euclidean_ring_tag Tag;
        typedef CGAL::Tag_true Is_exact;
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>();
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(-15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6), NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(-15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(-15));
  
        CGAL::test_real_embeddable<NT>();
        
    }
    return 0;
}
#else 
int main()
{
  return 0;
}
#endif //CGAL_USE_GMPXX
