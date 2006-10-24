#include <CGAL/Gbrs_algebraic_kernel.h>
#include <gmp.h>
#include <vector>

typedef CGAL::GBRS_algebraic_kernel<CGAL::Gmpz>	AlgKernel;
typedef AlgKernel::Coefficient			Coefficient;
typedef AlgKernel::Algebraic_real_1		Algebraic;
typedef AlgKernel::Polynomial_1			Polynomial;

inline std::ostream &show_alg (std::ostream &o, Algebraic &a) {
	return (o << a << " is root of " << a.pol() << " (root " << a.nr() <<
			", multiplicity " << a.mult() << ", precision " <<
			a.rsprec() << ")");
}

std::ostream &print_sign (std::ostream &o, CGAL::Sign s) {
	switch (s) {
		case CGAL::POSITIVE: return (o << "positive");
		case CGAL::NEGATIVE: return (o << "negative");
		case CGAL::ZERO: return (o << "zero");
		default: CGAL_assertion_msg(false,"UNDECIDED not implemented");
	}
	return o;	// never reached
}

std::ostream &print_comp (std::ostream &o, CGAL::Comparison_result r) {
	switch (r) {
		case CGAL::SMALLER: return (o << "smaller");
		case CGAL::LARGER: return (o << "larger");
		default: return (o << "equal");
	}
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
		ker.construct_solve_1_object() (p, rootsp.begin (), false);
	std::vector<Algebraic>::iterator itp;
	std::cout << "roots of p:" << std::endl;
	for (itp = rootsp.begin (); itp != p_end; ++itp)
		show_alg (std::cout, *itp) << std::endl;
	// now we create the polynomial q = x-1
	std::vector<Coefficient> coefsq;
	coefsq.push_back (Coefficient(-1414));
	coefsq.push_back (Coefficient(1000));
	Polynomial q = ker.construct_polynomial_1_object()
		(coefsq.begin (), coefsq.end ());
	std::cout << "\nq(x) = " << q << std::endl;
	// we solve it
	std::vector<Algebraic> rootsq (1);
	std::vector<Algebraic>::iterator q_end =
		ker.construct_solve_1_object() (q, rootsq.begin (), false);
	std::vector<Algebraic>::iterator itq;
	std::cout << "roots of q:" << std::endl;
	for (itq = rootsq.begin (); itq != q_end; ++itq)
		show_alg (std::cout, *itq) << std::endl;
	// now, we evaluate some roots
	std::cout << "\nsign of p evaluated at the root 0 of q: ";
	print_sign (std::cout, ker.construct_signat_1_object()(p, rootsq[0]));
	std::cout << std::endl;

	std::cout << "sign of q evaluated at the root 2 of p: ";
	print_sign (std::cout, ker.construct_signat_1_object()(q, rootsp[2]));
	std::cout << std::endl;

	std::cout << "sign of p evaluated at the root 0 of p: ";
	print_sign (std::cout, ker.construct_signat_1_object()(p, rootsp[0]));
	std::cout << std::endl;

	std::cout << "\np(3) = " << p.eval(3) <<
		"; sign of p evaluated at (coefficient)3: ";
	print_sign (std::cout,
			ker.construct_signat_1_object()(p, Coefficient(3)));
	std::cout << std::endl;

	std::cout << "\ncomparison between root 0 of p and root 0 of q: ";
	print_comp (std::cout,
			ker.construct_compare_1_object()(rootsp[0], rootsq[0]))
			<< std::endl;
	std::cout << "comparison between root 0 of q and root 0 of p: ";
	print_comp (std::cout,
			ker.construct_compare_1_object()(rootsq[0], rootsp[0]))
		<< std::endl;

	// using a small default precision (<=10), easily we can see how roots
	// are refined when needed in comparisons
	std::cout << "\ncomparison between root 2 of p and root 0 of q: ";
	print_comp (std::cout,
			ker.construct_compare_1_object()(rootsp[2], rootsq[0]))
		<< std::endl;
	show_alg (std::cout, rootsp[2]) << std::endl;
	show_alg (std::cout, rootsq[0]) << std::endl;

	return 0;
}

