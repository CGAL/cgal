// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>

/*! \file CGAL/Modular.C
  test for number type modul 
*/

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Modular.h>

#include <CGAL/Test/_test_algebraic_structure.h>

int main()
{   
    typedef CGAL::Modular NT;
    typedef CGAL::Field_tag Tag;
    typedef CGAL::Tag_true Is_exact;
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>();
    
    int old_prime = NT::get_current_prime(); 
    assert(old_prime == NT::set_current_prime(7));
    assert(7 == NT::get_current_prime()); 
    
    NT x(4),y(5),z(12),t;
    
    // operator ==
    assert(!(x==y));
    assert(y==z);
    // operator !=
    assert(x!=y);
    assert(!(z!=y));

    // constructor
    assert(NT(2)==NT(2-5*NT::get_current_prime()));
    assert(NT(2)==NT(2+5*NT::get_current_prime()));

    // operator unary +
    assert((+x)==x);
    // operator unary -
    assert(-x==x*NT(-1));
    
    // operator binary +
    assert((x+y)==NT(2));
    // operator binary -
    assert((x-y)==NT(6));
    // operator *
    assert((x*y)==NT(6));
    // operator /
    assert((x/y)==NT(5));

    
    // operator +=
    t=x; assert((x+y)==(t+=y));
    // operator -=
    t=x; assert((x-y)==(t-=y));
    // operator *=
    t=x; assert((x*y)==(t*=y));
    // operator /=
    t=x; assert((x/y)==(t/=y));

    // left/right Hand
    // operator ==
    assert(x==4);
    assert(5==y);
    // operator !=
    assert(x!=5);
    assert(4!=y);
    // operator +
    t=x; assert((x+5)==(5+x));
    // operator -
    t=x; assert((x-5)==(4-y));
    // operator *
    t=x; assert((x*5)==(5*x));
    // operator =
    t=x; assert((x/5)==(4/y));

    //cout << x << endl;
        
}


