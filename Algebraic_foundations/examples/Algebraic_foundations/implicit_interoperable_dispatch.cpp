#include <CGAL/Coercion_traits.h>
#include <CGAL/Quotient.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/IO/io.h>

// this is the implementation for ExplicitInteroperable types
template <typename A, typename B>
typename CGAL::Coercion_traits<A,B>::Type
binary_function_(const A& a , const B& b, CGAL::Tag_false){
    std::cout << "Call for ExplicitInteroperable types: " << std::endl;
    typedef CGAL::Coercion_traits<A,B> CT;
    typename CT::Cast cast;
    return cast(a)*cast(b);
}

// this is the implementation for ImplicitInteroperable types
template <typename A, typename B>
typename CGAL::Coercion_traits<A,B>::Type
binary_function_(const A& a , const B& b, CGAL::Tag_true){
    std::cout << "Call for ImpicitInteroperable types: " << std::endl;
    return a*b;
}

// this function selects the correct implementation
template <typename A, typename B>
typename CGAL::Coercion_traits<A,B>::Type
binary_func(const A& a , const B& b){
    typedef CGAL::Coercion_traits<A,B> CT;
    typedef typename CT::Are_implicit_interoperable Are_implicit_interoperable;
    return binary_function_(a,b,Are_implicit_interoperable());
}

int main(){
    CGAL::set_pretty_mode(std::cout);

    // Function call for ImplicitInteroperable types
    std::cout<< binary_func(double(3), int(5)) << std::endl;

    // Function call for ExplicitInteroperable types
    CGAL::Quotient<int>           rational(1,3);    // == 1/3
    CGAL::Sqrt_extension<int,int> extension(1,2,3); // == 1+2*sqrt(3)
    CGAL::Sqrt_extension<CGAL::Quotient<int>,int> result = binary_func(rational, extension);
    std::cout<< result << std::endl;

    return 0;
}
