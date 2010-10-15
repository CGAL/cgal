// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Michael Hemmer <mhemmer@uni-mainz.de>
// 

// Test program for the CGAL::Root_of_traits 



#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Test/_test_real_embeddable.h>
#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Root_of_2.h>
#include <CGAL/IO/io_tags.h>
#include <CGAL/Lazy_exact_nt.h>  



template <class T, class RootOf1, class RootOf2>
void test_root_of_traits(){
    // pure type checking  
    typedef CGAL::Root_of_traits<T> RoT;
    typedef typename RoT::Root_of_1 Root_of_1;
    typedef typename RoT::Root_of_2 Root_of_2;
    
    BOOST_STATIC_ASSERT((::boost::is_same<RootOf1,Root_of_1>::value));
    BOOST_STATIC_ASSERT((::boost::is_same<RootOf2,Root_of_2>::value));
    
    typedef typename RoT::Make_root_of_2 Make_root_of_2;
    typedef typename Make_root_of_2::result_type result_type;
    BOOST_STATIC_ASSERT((::boost::is_same<Root_of_2,result_type>::value));

    Root_of_2 r = CGAL::make_root_of_2(T(0),T(-1),T(2)); //-sqrt(2)
    Root_of_2 rl = CGAL::make_root_of_2(T(1),T(0),T(-2),true); //-sqrt(2);
    Root_of_2 rr = CGAL::make_root_of_2(T(1),T(0),T(-2),false); //+sqrt(2)
    assert(r == rl);
    assert(rl != rr);
    
    assert( r * Root_of_1(2) == CGAL::make_root_of_2(T(0),T(-2),T(2)));
    assert( r * T(2) == CGAL::make_root_of_2(T(0),T(-2),T(2)));
}

int main(){
    test_root_of_traits< double , double , double >();
    
#ifdef CGAL_USE_GMP
    //TODO: switch to Gmpq
    {
        typedef CGAL::Gmpz RT;
        typedef CGAL::Gmpq Root_of_1;
        typedef CGAL::Root_of_2<CGAL::Gmpz> Root_of_2;
        
        test_root_of_traits<RT,Root_of_1,Root_of_2>();
        test_root_of_traits<Root_of_1,Root_of_1,Root_of_2>();
    }{
        typedef CGAL::Lazy_exact_nt<CGAL::Gmpq> RT;
        typedef CGAL::Lazy_exact_nt<CGAL::Gmpq> Root_of_1;
        typedef CGAL::Lazy_exact_nt<CGAL::Root_of_2<CGAL::Gmpz> > Root_of_2;
        test_root_of_traits<RT,Root_of_1,Root_of_2>(); 
        test_root_of_traits<Root_of_1,Root_of_1,Root_of_2>();
    }
#endif
    return 0;
}
