#include <CGAL/use.h>
#include <CGAL/assertions.h>
#include <CGAL/Coercion_traits.h>
#include <cassert>
int main(){
    {
    typedef CGAL::Coercion_traits<int,int> CT;
    CGAL_USE_TYPE(CT);
    CGAL_static_assertion(( boost::is_same<CT::Type,int>::value));
    CGAL_static_assertion(
            ( boost::is_same<CT::Are_implicit_interoperable,CGAL::Tag_true>::value));
    CGAL_static_assertion(
            ( boost::is_same<CT::Are_explicit_interoperable,CGAL::Tag_true>::value));
    assert( 5 == CT::Cast()(5));
    }
    {
    typedef CGAL::Coercion_traits<CGAL::Tag_true,CGAL::Tag_false> CT;
    CGAL_USE_TYPE(CT);
//    CGAL_static_assertion(( boost::is_same<CT::Type,CGAL::Null_type>::value));
    CGAL_static_assertion(
            ( boost::is_same<CT::Are_implicit_interoperable,CGAL::Tag_false>::value));
    CGAL_static_assertion(
            ( boost::is_same<CT::Are_explicit_interoperable,CGAL::Tag_false>::value));
    CGAL_static_assertion(
            ( boost::is_same<CT::Cast,CGAL::Null_functor>::value));
    }
}
