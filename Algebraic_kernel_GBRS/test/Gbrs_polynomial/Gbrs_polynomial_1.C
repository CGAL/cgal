#include <CGAL/Gbrs_algebraic_kernel.h>
#include <gmp.h>
#include <vector>

typedef CGAL::GBRS_algebraic_kernel<CGAL::Gmpz>	AlgKernel;
typedef AlgKernel::Coefficient			Coefficient;
typedef AlgKernel::Algebraic_real_1		Algebraic;
typedef AlgKernel::Polynomial_1			Polynomial;

std::ostream &print_sign (std::ostream &o, CGAL::Sign s) {
	switch (s) {
		case CGAL::POSITIVE: o << "positive";
				     return o;
		case CGAL::NEGATIVE: o << "negative";
				     return o;
		case CGAL::ZERO: o << "zero";
				 return o;
		default: CGAL_assertion_msg (false, "UNDECIDED not implemented");
	}
	return o;	// never reached
}

int main () {
	AlgKernel ker;
	// construct the polynomial x^3-2x
	std::vector<Coefficient> coefsp;
	coefsp.push_back (Coefficient (0));	// x^0
	coefsp.push_back (Coefficient (-2));	// x^1
	coefsp.push_back (Coefficient (0));	// x^2
	coefsp.push_back (Coefficient (1));	// x^3
	coefsp.push_back (Coefficient (0));	// x^4
	coefsp.push_back (Coefficient (0));	// x^5
	// the container has now all the coefficients in increasing monomial
	// order: <0, -2, 0, 1, 0, 0>
	Polynomial p = ker.construct_polynomial_1_object()
		(coefsp.begin (), coefsp.end ());
	std::cout << "p(x) = " << p << std::endl;
	/*
	// we test some overloaded operators:
	std::cout << "p(x) * (int)3 = " << p*3 << std::endl;
	mpz_t x;
	mpz_init (x);
	mpz_set_ui (x, 3);
	std::cout << "p(x) * (mpz_t)3 = " << p*x << std::endl;
	mpz_clear (x);
	*/

	// find the roots of p
	std::vector<Algebraic> rootsp (3);
	std::vector<Algebraic>::iterator p_end =
		ker.construct_solve_1_object () (p, rootsp.begin (), false);
	std::vector<Algebraic>::iterator itp;
	std::cout << "roots of p:" << std::endl;
	for (itp = rootsp.begin (); itp != p_end; ++itp) {
		std::cout << *itp << " is root of " << itp->pol();
		std::cout << " (root " << itp->nr() << ", multiplicity " <<
			itp->mult() << ")" << std::endl;
	}

	// now we create the polynomial q = x-1
	std::vector<Coefficient> coefsq;
	coefsq.push_back (Coefficient(-1));
	coefsq.push_back (Coefficient(1));
	Polynomial q = ker.construct_polynomial_1_object()
		(coefsq.begin (), coefsq.end ());
	std::cout << "\nq(x) = " << q << std::endl;
	// we solve it
	std::vector<Algebraic> rootsq (1);
	std::vector<Algebraic>::iterator q_end =
		ker.construct_solve_1_object () (q, rootsq.begin (), false);
	std::vector<Algebraic>::iterator itq;
	std::cout << "roots of q:" << std::endl;
	for (itq = rootsq.begin (); itq != q_end; ++itq) {
		std::cout << *itq << " is root of " << itq->pol();
		std::cout << " (root " << itq->nr() << ", multiplicity " <<
			itq->mult() << ")" << std::endl;
	}

	// now, we evaluate some roots
	std::cout << "\nsign of p evaluated at the root 0 of q: ";
	print_sign (std::cout, ker.construct_signat_1_object()(p, rootsq[0]));
	std::cout << std::endl;

	std::cout << "sign of q evaluated at the root 2 of p: ";
	print_sign (std::cout, ker.construct_signat_1_object()(q, rootsp[2]));
	std::cout << std::endl;

	std::cout << "\np(3) = " << p.eval(3) << std::endl;
	std::cout << "sign of p evaluated at (coefficient)3: ";
	print_sign (std::cout,
			ker.construct_signat_1_object()(p, Coefficient(3)));
	std::cout << std::endl;

	return 0;
}

