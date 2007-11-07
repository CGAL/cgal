// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>

/*! \file CGAL/Modular.C
  test for number type modul 
*/

#include <CGAL/basic.h>
#include <CGAL/Testsuite/assert.h>
#include <CGAL/Modular.h>

#include <CGAL/_test_algebraic_structure.h>

int main()
{   
    typedef CGAL::Modular NT;
    typedef CGAL::Field_tag Tag;
    typedef CGAL::Tag_true Is_exact;
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>();
    
    int old_prime = NT::get_current_prime(); 
    CGAL_test_assert(old_prime == NT::set_current_prime(7));
    CGAL_test_assert(7 == NT::get_current_prime()); 
    
    NT x(4),y(5),z(12),t;
    
    // operator ==
    CGAL_test_assert(!(x==y));
    CGAL_test_assert(y==z);
    // operator !=
    CGAL_test_assert(x!=y);
    CGAL_test_assert(!(z!=y));

    // constructor
    CGAL_test_assert(NT(2)==NT(2-5*NT::get_current_prime()));
    CGAL_test_assert(NT(2)==NT(2+5*NT::get_current_prime()));

    // operator unary +
    CGAL_test_assert((+x)==x);
    // operator unary -
    CGAL_test_assert(-x==x*NT(-1));
    
    // operator binary +
    CGAL_test_assert((x+y)==NT(2));
    // operator binary -
    CGAL_test_assert((x-y)==NT(6));
    // operator *
    CGAL_test_assert((x*y)==NT(6));
    // operator /
    CGAL_test_assert((x/y)==NT(5));

    
    // operator +=
    t=x; CGAL_test_assert((x+y)==(t+=y));
    // operator -=
    t=x; CGAL_test_assert((x-y)==(t-=y));
    // operator *=
    t=x; CGAL_test_assert((x*y)==(t*=y));
    // operator /=
    t=x; CGAL_test_assert((x/y)==(t/=y));

    // left/right Hand
    // operator ==
    CGAL_test_assert(x==4);
    CGAL_test_assert(5==y);
    // operator !=
    CGAL_test_assert(x!=5);
    CGAL_test_assert(4!=y);
    // operator +
    t=x; CGAL_test_assert((x+5)==(5+x));
    // operator -
    t=x; CGAL_test_assert((x-5)==(4-y));
    // operator *
    t=x; CGAL_test_assert((x*5)==(5*x));
    // operator =
    t=x; CGAL_test_assert((x/5)==(4/y));

    //cout << x << endl;
        
}


