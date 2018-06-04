#include <iostream>
#include <CGAL/CORE_Expr.h>
#include <CGAL/internal/Exact_complex.h>
#include <CGAL/Point_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Cartesian.h>

using namespace CGAL;
using namespace std;

int main(void) {

    typedef CORE::Expr                  NT;
    typedef Exact_complex<NT>           cplx;
    typedef Point_2< Cartesian<NT> >    Point;
    typedef Circle_2< Cartesian<NT> >   Circle;

    cplx n;
    n.set_imag(-1);
    n.set_real(3);
    cout << "n = " << n << endl;
    cout << "sq_mod(n) = " << n.square_modulus() << endl;
    cout << "mod(n)    = " << n.modulus() << endl;
    cout << endl;
    cplx m(2, 3);
    cout << "m = " << m << endl;
    cout << "sq_mod(m) = " << m.square_modulus() << endl;
    cout << "mod(m)    = " << m.modulus() << endl;
    cout << endl;

    cout << "m + n = " << m + n << endl;
    cout << "m - n = " << m - n << endl;
    cout << "m * n = " << m * n << endl;
    cout << "m / n = " << m / n << endl;
    cout << "n / m = " << n / m << endl;

    cplx mon = m/n;
    cplx nom = n/m;

    cplx res1(NT(3)/NT(10), NT(11)/NT(10));
    cplx res2(NT(3)/NT(13), NT(-11)/NT(13));

    CGAL_assertion(mon == res1);
    CGAL_assertion(nom == res2);

    cout << "n < m: " << (n < m ? "true" : "false") << endl;

    Point pt(NT(2)/NT(10), NT(6)/NT(10));
    cplx p(pt.x(), pt.y());
    cout << "p = " << p << endl;
    cout << "recip = " << p.reciprocal() << endl;
    cout << "inver = " << p.invert_in_unit_circle() << endl;

    Circle c(Point(NT(2)/NT(10), NT(3)/NT(4)), NT(1)/NT(2));
    cplx ipc = p.invert_in_circle(c);
    cout << "inverted in circle: " << ipc << endl;

    return 0;
}