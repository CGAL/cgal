#include <CGAL/basic.h>

namespace CGAL{
namespace Test{

template <class T, class RootOf1, class RootOf2>
void test_root_of_traits(){
    // pure type checking  
    typedef CGAL::Root_of_traits<T> RoT;
    typedef typename RoT::Root_of_1 Root_of_1;
    typedef typename RoT::Root_of_2 Root_of_2;
    
    CGAL_static_assertion((::boost::is_same<RootOf1,Root_of_1>::value));
    CGAL_static_assertion((::boost::is_same<RootOf2,Root_of_2>::value));
    
    typedef typename RoT::Make_root_of_2 Make_root_of_2;
    typedef typename RoT::Make_sqrt      Make_sqrt;
    typedef typename RoT::Inverse        Inverse;  
    typedef typename RoT::Square         Square; 

    const Make_root_of_2& make_root_of_2 = Make_root_of_2();
    const Make_sqrt&      make_sqrt      = Make_sqrt();
    const Inverse&        inverse        = Inverse();  
    const Square&         square         = Square(); 
    
    CGAL_static_assertion((::boost::is_same<Root_of_2,typename Make_root_of_2::result_type>::value));
    CGAL_static_assertion((::boost::is_same<Root_of_2,typename Make_sqrt::result_type>::value));
    CGAL_static_assertion((::boost::is_same<Root_of_2,typename Inverse::result_type>::value));
    CGAL_static_assertion((::boost::is_same<Root_of_2,typename Square::result_type>::value));


    {
      Root_of_2 r  = make_root_of_2(T(0),T(-1),T(2)); //-sqrt(2)
      Root_of_2 rl = make_root_of_2(T(1),T(0),T(-2),true); //-sqrt(2);
      Root_of_2 rr = make_root_of_2(T(1),T(0),T(-2),false); //+sqrt(2)
      assert(r == rl);
      assert(rl != rr);
      
      assert( r * Root_of_1(2) == make_root_of_2(T(0),T(-2),T(2)));
      assert( r * T(2) == make_root_of_2(T(0),T(-2),T(2)));
    }{
      Root_of_2 r  = CGAL::make_root_of_2(T(0),T(-1),T(2)); //-sqrt(2)
      Root_of_2 rl = CGAL::make_root_of_2(T(1),T(0),T(-2),true); //-sqrt(2);
      Root_of_2 rr = CGAL::make_root_of_2(T(1),T(0),T(-2),false); //+sqrt(2)
      assert(r == rl);
      assert(rl != rr);
      
      assert( r * Root_of_1(2) == CGAL::make_root_of_2(T(0),T(-2),T(2)));
      assert( r * T(2) == CGAL::make_root_of_2(T(0),T(-2),T(2)));
    }

   
    {
      Root_of_2 r  = make_sqrt(T(2)); //sqrt(2)
      Root_of_2 rr = make_root_of_2(T(1),T(0),T(-2),false); //+sqrt(2)
      assert(r == rr);
    }{
      Root_of_2 r  = CGAL::make_sqrt(T(2)); //sqrt(2)
      Root_of_2 rr = CGAL::make_root_of_2(T(1),T(0),T(-2),false); //+sqrt(2)
      assert(r == rr);
    }

    {
      Root_of_2 r  = inverse(CGAL::make_sqrt(T(2))); 
      Root_of_2 rr = 1/CGAL::make_sqrt(T(2));
      assert(r == rr);
    }{
        Root_of_2 r  = CGAL::inverse(CGAL::make_sqrt(T(2)));
        Root_of_2 rr = 1/CGAL::make_sqrt(T(2)); 
        assert(r == rr);
    }

    {
      Root_of_2 r  = square(CGAL::make_sqrt(T(2)));
      Root_of_2 rr = CGAL::make_sqrt(T(2))*CGAL::make_sqrt(T(2));
      assert(r == rr);
    }{
      Root_of_2 r  = CGAL::square(CGAL::make_sqrt(T(2)));
      Root_of_2 rr = CGAL::make_sqrt(T(2))*CGAL::make_sqrt(T(2));
      assert(r == rr);
    }

    bool is_not_exact = !CGAL::Algebraic_structure_traits<T>::Is_exact::value; 
    {
      std::vector<Root_of_2> roots;
      CGAL::compute_roots_of_2(T(1),T(0),T(-2),std::back_inserter(roots));
      assert(roots.size()==2);  
      assert(roots[0]==-CGAL::make_sqrt(T(2)) || is_not_exact );
      assert(roots[1]== CGAL::make_sqrt(T(2)) || is_not_exact );
    }  
    {
      Root_of_2 roots[2]= {Root_of_2(1),Root_of_2(1)};
      CGAL::compute_roots_of_2(T(13),T(4),T(-23),roots); 
      assert(roots[0]==CGAL::make_root_of_2(T(13),T(4),T(-23),true)  || is_not_exact );
      assert(roots[1]==CGAL::make_root_of_2(T(13),T(4),T(-23),false) || is_not_exact );
    }   
    {
      std::vector<Root_of_2> roots;
      CGAL::compute_roots_of_2(T(1),T(-6),T(9),std::back_inserter(roots));
      assert(roots.size()==1);
      assert(roots[0]==Root_of_2(3) || is_not_exact );
    }  
    {
      std::vector<Root_of_2> roots;
      CGAL::compute_roots_of_2(T(1),T(0),T(2),std::back_inserter(roots));
      assert(roots.size()==0);
    } 
    {
      std::vector<Root_of_2> roots;
      CGAL::compute_roots_of_2(T(0),T(2),T(3),std::back_inserter(roots));
      assert(roots.size()==1);
      assert(roots[0]==-Root_of_2(3)/Root_of_2(2) || is_not_exact );
    } 
    
}

} // namespace Test 
} // namespace CGAL 
