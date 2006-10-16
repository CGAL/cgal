#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
#include <mpfr.h>
#include <CGAL/Algebraic_1.h>

bool test_int () {
	int one = 1;
	int two = 2;

	CGAL::Algebraic_1 interval1 (1);
	CGAL::Algebraic_1 interval2 (0.5, 1.5);

	interval1.is_point ();
	interval2.is_point ();

	interval1 == one;

	try { (interval1 <= one);
	} catch (CGAL::comparison_overlap_exn &o) { true;}

	try { (interval1 < one);
	} catch (CGAL::comparison_overlap_exn &o) { true;}

	(interval1 == two);

	try { (interval2 == one);
	} catch (CGAL::comparison_overlap_exn &o) { true;}

	try { (interval2 == two);
	} catch (CGAL::comparison_overlap_exn &o) { true;}

	return true;
}

bool test_intervals () {
	CGAL::Algebraic_1 i1 (7);
	CGAL::Algebraic_1 i2 (9);
	(i1==i2);
	(i1!=i2);
	(i1<i2);
	(i1>i2);
	(i1<=i2);
	(i1>=i2);
	i2 = 1 + i1;
	(i2 == 8);
	(i2 != 8);
	(i2 < 8);
	(i2 > 8);
	(i2 <= 8);
	(i2 >= 8);
	mpfr_t x;
	mpfr_init (x);
	mpfr_set_si (x, 1, GMP_RNDN);
	i2 = x + i1;
	mpz_t y;
	mpz_init_set_si (y, 1);
	i2 = y + i1;
	CGAL::Gmpz z(1);
	i2 = z + i1;
	//--------------------------------------------------
	// std::cout << "i1 = " << i1 << std::endl;
	// std::cout << "i2 = " << i2 << std::endl;
	// std::cout << "i1==i2: " << (i1==i2) << std::endl;
	// std::cout << "i1!=i2: " << (i1!=i2) << std::endl;
	// std::cout << "i1<i2: " << (i1<i2) << std::endl;
	// std::cout << "i1>i2: " << (i1>i2) << std::endl;
	// std::cout << "i1<=i2: " << (i1<=i2) << std::endl;
	// std::cout << "i1>=i2: " << (i1>=i2) << std::endl;
	// i2 = 1 + i1;
	// std::cout << "(int)1 + i1 =\t" << i2 << std::endl;
	// std::cout << "1+i1 == 8? : " << (i2 == 8) << std::endl;
	// std::cout << "1+i1 != 8? : " << (i2 != 8) << std::endl;
	// std::cout << "1+i1 < 8? : " << (i2 < 8) << std::endl;
	// std::cout << "1+i1 > 8? : " << (i2 > 8) << std::endl;
	// std::cout << "1+i1 <= 8? : " << (i2 <= 8) << std::endl;
	// std::cout << "1+i1 >= 8? : " << (i2 >= 8) << std::endl;
	// mpfr_t x;
	// mpfr_init (x);
	// mpfr_set_si (x, 1, GMP_RNDN);
	// i2 = x + i1;
	// std::cout << "(mpfr_t)1 + i1 =\t" << i2 << std::endl;
	// mpz_t y;
	// mpz_init_set_si (y, 1);
	// i2 = y + i1;
	// std::cout << "(mpz_t)1 + i1 =\t" << i2 << std::endl;
	// CGAL::Gmpz z(1);
	// i2 = z + i1;
	// std::cout << "(Gmpz)1 + i1 =\t" << i2 << std::endl;
	//-------------------------------------------------- 
	return true;
}

// you can use CGAL::Gmpz and CGAL::Gmpq as T
template <class T>
bool test_class () {
	T one (1);
	T two (2);

	CGAL::Algebraic_1 interval1 (1);
	CGAL::Algebraic_1 interval2 (0.5, 1.5);

	interval1.is_point ();
	interval2.is_point ();

	(interval1 == one);
	(interval1 != one);

	try { (interval1 <= one);
	} catch (CGAL::comparison_overlap_exn &o) { true;}

	try { (interval1 < one);
	} catch (CGAL::comparison_overlap_exn &o) { true;}

	(interval1 == two);

	try { (interval2 == one);
	} catch (CGAL::comparison_overlap_exn &o) { true;}

	try { (interval2 == two);
	} catch (CGAL::comparison_overlap_exn &o) { true;}

	return true;
}

bool test_mpfr () {
	mpfr_t one, two;
	mpfr_inits (one, two, NULL);
	mpfr_set_si (one, 1, GMP_RNDN);
	mpfr_set_si (two, 2, GMP_RNDN);

	CGAL::Algebraic_1 interval1 (1);
	CGAL::Algebraic_1 interval2 (0.5, 1.5);

	interval1.is_point ();
	interval2.is_point ();

	(interval1 == one);
	(one == interval1);

	try { (interval1 <= one);
	} catch (CGAL::comparison_overlap_exn &o) { true;}

	try { (interval1 < one);
	} catch (CGAL::comparison_overlap_exn &o) { true;}

	(interval1 == two);

	try { (interval2 == one);
	} catch (CGAL::comparison_overlap_exn &o) { true;}

	try { (interval2 == two);
	} catch (CGAL::comparison_overlap_exn &o) { true;}

	return true;
}

int main () {
	if ((test_int ()) &&
			(test_intervals ()) &&
			(test_class<CGAL::Gmpz> ()) &&
			(test_class<CGAL::Gmpq> ()) &&
			(test_mpfr ()))
		return 0;
	return 1;
}
