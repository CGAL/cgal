#include <CGAL/basic.h>
#include <CGAL/Coercion_traits.h>
#include <CGAL/Testsuite/assert.h>
int main(){
    {
    typedef CGAL::Coercion_traits<int,int> CT;
    BOOST_STATIC_ASSERT(( boost::is_same<CT::Type,int>::value));
    BOOST_STATIC_ASSERT(
            ( boost::is_same<CT::Are_implicit_interoperable,CGAL::Tag_true>::value));
    BOOST_STATIC_ASSERT(
            ( boost::is_same<CT::Are_explicit_interoperable,CGAL::Tag_true>::value));
    CGAL_test_assert( 5 == CT::Cast()(5));
    }   
    {
    typedef CGAL::Coercion_traits<CGAL::Tag_true,CGAL::Tag_false> CT;
//    BOOST_STATIC_ASSERT(( boost::is_same<CT::Type,CGAL::Null_type>::value));
    BOOST_STATIC_ASSERT(
            ( boost::is_same<CT::Are_implicit_interoperable,CGAL::Tag_false>::value));
    BOOST_STATIC_ASSERT(
            ( boost::is_same<CT::Are_explicit_interoperable,CGAL::Tag_false>::value));
    BOOST_STATIC_ASSERT(
            ( boost::is_same<CT::Cast,CGAL::Null_functor>::value));
    }      
}
