#include <CGAL/Coercion_traits.h>
#include <CGAL/IO/io.h>

// this is an implementation for ExplicitInteroperable types
// the result type is determined via Coercion_traits<A,B>
template <typename A, typename B>
typename CGAL::Coercion_traits<A,B>::Type
binary_func(const A& a , const B& b){
    typedef CGAL::Coercion_traits<A,B> CT;

    // check for explicit interoperability
    CGAL_static_assertion((CT::Are_explicit_interoperable::value));

    // CT::Cast is used to to convert both types into the coercion type
    typename CT::Cast cast;
    // all operations are performed in the coercion type
    return cast(a)*cast(b);
}

int main(){
    // Function call for the interoperable types
    std::cout<< binary_func(double(3), int(5)) << std::endl;
    // Note that Coercion_traits is symmetric
    std::cout<< binary_func(int(3), double(5)) << std::endl;
    return 0;
}
