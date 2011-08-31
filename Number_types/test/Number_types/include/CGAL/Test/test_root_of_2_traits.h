#include <CGAL/basic.h>

template <class T, class RootOf1, class RootOf2>
void test_root_of_traits(){
    // pure type checking  
    typedef CGAL::Root_of_traits<T> RoT;
    typedef typename RoT::Root_of_1 Root_of_1;
    typedef typename RoT::Root_of_2 Root_of_2;
    
    CGAL_static_assertion((::boost::is_same<RootOf1,Root_of_1>::value));
    CGAL_static_assertion((::boost::is_same<RootOf2,Root_of_2>::value));
    
    typedef typename RoT::Make_root_of_2 Make_root_of_2;
    typedef typename Make_root_of_2::result_type result_type;
    CGAL_static_assertion((::boost::is_same<Root_of_2,result_type>::value));

    const Make_root_of_2& make_root_of_2 = Make_root_of_2();
    Root_of_2 r  = make_root_of_2(T(0),T(-1),T(2));        //-sqrt(2)
    Root_of_2 rl = make_root_of_2(T(1),T(0),T(-2),true);  //-sqrt(2);
    Root_of_2 rr = make_root_of_2(T(1),T(0),T(-2),false); //+sqrt(2)
    assert(r == rl);
    assert(rl != rr);
    
    assert( r * Root_of_1(2) == CGAL::make_root_of_2(T(0),T(-2),T(2)));
    assert( r * T(2) == CGAL::make_root_of_2(T(0),T(-2),T(2)));
}
