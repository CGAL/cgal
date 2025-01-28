#include <CGAL/use.h>
#include <CGAL/assertions.h>
#include <CGAL/Coercion_traits.h>
#include <cassert>
int main(){
    {
    typedef CGAL::Coercion_traits<int,int> CT;
    CGAL_USE_TYPE(CT);
    static_assert( std::is_same<CT::Type,int>::value);
    static_assert( std::is_same<CT::Are_implicit_interoperable,CGAL::Tag_true>::value);
    static_assert( std::is_same<CT::Are_explicit_interoperable,CGAL::Tag_true>::value);
    assert( 5 == CT::Cast()(5));
    }
    {
    typedef CGAL::Coercion_traits<CGAL::Tag_true,CGAL::Tag_false> CT;
    CGAL_USE_TYPE(CT);
//    static_assert( std::is_same<CT::Type,CGAL::Null_type>::value);
    static_assert(std::is_same<CT::Are_implicit_interoperable,CGAL::Tag_false>::value);
    static_assert(std::is_same<CT::Are_explicit_interoperable,CGAL::Tag_false>::value);
    static_assert(std::is_same<CT::Cast,CGAL::Null_functor>::value);
    }
}
