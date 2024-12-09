#include <CGAL/use.h>
#include <CGAL/int.h>
#include <CGAL/Scalar_factor_traits.h>
#include <cassert>


int main(){
    typedef CGAL::Scalar_factor_traits<int> SFT;
    CGAL_USE_TYPE(SFT);
    static_assert(::std::is_same<int, SFT::Type>::value);
    static_assert(::std::is_same<int, SFT::Scalar>::value);

    typedef SFT::Scalar_factor Scalar_factor;
    {
        typedef Scalar_factor::result_type result_type;
        CGAL_USE_TYPE(result_type);
        static_assert(::std::is_same<int, result_type>::value);

        typedef Scalar_factor::argument_type argument_type;
        CGAL_USE_TYPE(argument_type);
        static_assert(::std::is_same<int, argument_type>::value);
    }
    typedef SFT::Scalar_div Scalar_div;
    {
        typedef Scalar_div::result_type result_type;
        CGAL_USE_TYPE(result_type);
        static_assert(::std::is_same<void, result_type>::value);

        typedef Scalar_div::first_argument_type first_argument_type;
        CGAL_USE_TYPE(first_argument_type);
        static_assert(::std::is_same<int&, first_argument_type>::value);
        typedef Scalar_div::second_argument_type second_argument_type;
        CGAL_USE_TYPE(second_argument_type);
        static_assert(::std::is_same<int, second_argument_type>::value);
    }

    int i;
    i = 0 ;  CGAL::remove_scalar_factor(i); assert( 0 == i);
    i = 1 ;  CGAL::remove_scalar_factor(i); assert( 1 == i);
    i = 2 ;  CGAL::remove_scalar_factor(i); assert( 1 == i);
    i =-2 ;  CGAL::remove_scalar_factor(i); assert(-1 == i);
}
