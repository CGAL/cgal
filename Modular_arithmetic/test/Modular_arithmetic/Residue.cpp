// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>

/*! \file CGAL/Residue.C
  test for number type modul
*/

#include <CGAL/config.h>
#include <cassert>
#include <CGAL/Residue.h>
#include <CGAL/FPU.h>
#include <CGAL/Modular_traits.h>

#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/number_utils.h>

int main()
{
    // Enforce IEEE double precision and rounding mode to nearest
    CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST);

    typedef CGAL::Residue NT;
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

    assert(7 == NT::set_current_prime(old_prime));
    typedef long Integer;

    Integer       int_x(7);
    Integer       prime(NT::get_current_prime());
    CGAL::Residue mod_x(7);
    for(int i = 0; i < 10000; i++){
        assert(mod_x == CGAL::modular_image(int_x));
        int_x *= int_x; int_x = CGAL::mod(int_x, prime);
        mod_x *= mod_x;
        assert(mod_x == CGAL::modular_image(int_x));
        int_x += int_x; int_x = CGAL::mod(int_x, prime);
        mod_x += mod_x;
        assert(mod_x == CGAL::modular_image(int_x));

        assert(mod_x == CGAL::modular_image(int_x));
        int_x *= int_x; int_x = CGAL::mod(int_x, prime);
        mod_x *= mod_x;

        assert(mod_x == CGAL::modular_image(int_x));
        int_x -= int_x; int_x = CGAL::mod(int_x, prime);
        mod_x -= (CGAL::Residue&)mod_x;
    }
    {
        CGAL::Residue::set_current_prime(67111043);
        CGAL::Residue x(-33546401);
        CGAL::Residue y(23950928);
        CGAL::Residue q = CGAL::integral_division(x,y);
        assert(x == q*y);
    }
}


