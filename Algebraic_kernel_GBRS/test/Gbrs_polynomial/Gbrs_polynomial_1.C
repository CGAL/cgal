#include <CGAL/Gbrs_algebraic_kernel.h>
#include <gmp.h>
#include <vector>

typedef CGAL::GBRS_algebraic_kernel<CGAL::Gmpz>	AlgKernel;
typedef AlgKernel::Coefficient			Coefficient;
typedef AlgKernel::Algebraic_real_1		Algebraic;
typedef AlgKernel::Polynomial_1			Polynomial;

int main () {
	AlgKernel ker;
	// construct the polynomial x^3-2x
	std::vector<Coefficient> coefs;
	coefs.push_back (Coefficient (0));	// x^0
	coefs.push_back (Coefficient (-2));	// x^1
	coefs.push_back (Coefficient (0));	// x^2
	coefs.push_back (Coefficient (1));	// x^3
	coefs.push_back (Coefficient (0));	// x^4
	coefs.push_back (Coefficient (0));	// x^5
	// the container has now all the coefficients in increasing monomial
	// order: <0, -2, 0, 1, 0, 0>
	Polynomial p = ker.construct_polynomial_1_object()
		(coefs.begin (), coefs.end ());
	std::cout << "p(x) = " << p << std::endl;
	std::cout << "p(x) * (int)3 = " << p*3 << std::endl;
	mpz_t x;
	mpz_init (x);
	mpz_set_ui (x, 3);
	std::cout << "p(x) * (mpz_t)3 = " << p*x << std::endl;
	mpz_clear (x);

	// find the roots of p
	std::vector<Algebraic> roots (3);
	std::vector<Algebraic>::iterator r_end =
		ker.construct_solve_1_object ()
			(p, roots.begin (), false);
	std::vector<Algebraic>::iterator it;
	std::cout << "roots:" << std::endl;
	for (it = roots.begin (); it != r_end; ++it)
		std::cout << *it << std::endl;

	return 0;
}

