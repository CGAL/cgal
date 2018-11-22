#include <iostream>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Point_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Cartesian.h>
// #include <CGAL/Hyperbolic_octagon_translation_matrix.h>
#include <CGAL/internal/Exact_complex.h>
#include <CGAL/Hyperbolic_octagon_translation.h>
#include <vector>

using namespace CGAL;
using namespace std;

int main(void) {

    typedef CORE::Expr                                          NT;
    typedef Exact_complex<NT>                                   ECplx;
    typedef Hyperbolic_octagon_translation<NT>                  Word;

    Word w;
    // cout << "empty word: " << w << ", matrix: " << w.matrix() << endl;

    Word a(0), ab(0, 5), abc(0, 5, 2), abcd(0, 5, 2, 7), dcb(7, 2, 5), dc(7, 2), d(7);
    // cout << "a    = " << a << ",    matrix: " << a.matrix()     << endl;
    // cout << "ab   = " << ab << ",   matrix: " << ab.matrix()    << endl;
    // cout << "abc  = " << abc << ",  matrix: " << abc.matrix()   << endl;
    // cout << "abcd = " << abcd << ", matrix: " << abcd.matrix()  << endl;
    // cout << "dcb  = " << dcb << ",  matrix: " << dcb.matrix()   << endl;
    // cout << "dc   = " << dc << ",   matrix: " << dc.matrix()    << endl;
    // cout << "d    = " << d << ",    matrix: " << d.matrix()     << endl;    

    return 0;
}