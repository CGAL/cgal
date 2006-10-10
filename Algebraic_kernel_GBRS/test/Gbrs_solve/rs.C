#include <CGAL/Gbrs_algebraic_kernel.h>
#include <vector>

typedef CGAL::GBRS_algebraic_kernel<CGAL::Gmpq>	AlgKernel;
typedef AlgKernel::Coefficient			Coefficient;
typedef AlgKernel::Algebraic_real_1		Algebraic;
typedef AlgKernel::Polynomial_1			Polynomial;

int main () {
	AlgKernel ker;
	// construct the polynomial x^3-2x
	std::vector<Coefficient> coefs;
	coefs.push_back (Coefficient (1));
	coefs.push_back (Coefficient (0));
	coefs.push_back (Coefficient (-2));
	coefs.push_back (Coefficient (0));
	coefs.push_back (Coefficient (0));
	coefs.push_back (Coefficient (0));
	Polynomial p = ker.construct_polynomial_1_object()
		(coefs.begin (), coefs.end ());
	std::cout << "p(x) = " << p << std::endl;

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

