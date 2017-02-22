#include <iostream>
#include <CGAL/CORE_Expr.h>
#include <CGAL/exact_complex.h>
#include <CGAL/Point_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Hyperbolic_octagon_translation_matrix.h>
#include <vector>

using namespace CGAL;
using namespace std;

int main(void) {

    typedef CORE::Expr                                  NT;
    typedef Hyperbolic_octagon_translation_matrix<NT>   Matrix;
    typedef Point_2< Cartesian<NT> >                    Point;
    // typedef Circle_2< Cartesian<NT> >   Circle;

    Matrix m;
    cout << "Identity matrix: " << m << endl;

    vector<Matrix> gens;
    get_generators(gens);
    for (int i = 0; i < gens.size(); i++) {
        cout << "g[" << i << "] = " << gens[i] << endl;
    }

    CGAL_assertion(gens[0]*gens[4] == m);
    CGAL_assertion(gens[1]*gens[5] == m);
    CGAL_assertion(gens[2]*gens[6] == m);
    CGAL_assertion(gens[3]*gens[7] == m);

    Point o(NT(0), NT(0));
    vector<Point> imp;
    for (int i = 0; i < gens.size(); i++) {
        imp.push_back(gens[i].apply(o));
        cout << "imp[" << i << "] = " << imp[i] << endl;
    }

    CGAL_assertion(imp[0] == Point(CGAL::sqrt(NT(2))/(CGAL::sqrt(NT(1)+CGAL::sqrt(NT(2)))), NT(0)));
    CGAL_assertion(imp[1] == Point(NT(1)/(CGAL::sqrt(NT(1)+CGAL::sqrt(NT(2)))), NT(1)/(CGAL::sqrt(NT(1)+CGAL::sqrt(NT(2))))));

    cout << "test concluded successfully!" << endl;
    return 0;
}