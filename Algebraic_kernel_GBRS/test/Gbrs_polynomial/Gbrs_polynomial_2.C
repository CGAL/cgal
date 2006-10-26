#include <CGAL/Gbrs_polynomial_2.h>

typedef CGAL::Rational_polynomial_2	Polynomial;

int main () {
	Polynomial p (2, 1);	// d_x=2, d_y=1
	Polynomial q (7, 8);
	Polynomial zero (1, 4);	// zero(x) = 0

	p.set_coef (2, 1, 2);	// 2*x^2*y
	p.set_coef (1, 0, 3);	// 3*x
	p.set_coef (0, 0, 7);	// 7
	// p = 2*x^2*y + 3*x + 7
	q = -p;
	p.set_coef (0, 0, 5);
	std::cout << "p = " << p << "\nq = " << q << "\n0 = " << zero << "\n";
	std::cout << "p==q : " << (p==q) << "\np!=q : " << (p!=q) << "\n";

	Polynomial r;
	Polynomial s (1, 0);
	s.set_coef (1, 0, 1);	// s = 1 * x^1
	r = p*s;
	std::cout << "s = " << s << "\nr = p*s = " << r <<
		"\np*=(s*2) = " << (p*=(s*2)) << std::endl;

	return 0;
}

