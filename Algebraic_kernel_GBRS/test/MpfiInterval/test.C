#include <CGAL/Gmpz.h>
#include <mpfr.h>
#include <CGAL/MpfiInterval.h>

void test_intervals () {
	CGAL::MpfiInterval i1 (7);
	CGAL::MpfiInterval i2 (9);
	std::cout << "i1 = " << i1 << std::endl;
	std::cout << "i2 = " << i2 << std::endl;
	std::cout << "i1==i2: " << (i1==i2) << std::endl;
	std::cout << "i1!=i2: " << (i1!=i2) << std::endl;
	std::cout << "i1<i2: " << (i1<i2) << std::endl;
	std::cout << "i1>i2: " << (i1>i2) << std::endl;
	std::cout << "i1<=i2: " << (i1<=i2) << std::endl;
	std::cout << "i1>=i2: " << (i1>=i2) << std::endl;
	i2 = 1 + i1;
	std::cout << "(int)1 + i1 =\t" << i2 << std::endl;
	std::cout << "1+i1 == 8? : " << (i2 == 8) << std::endl;
	std::cout << "1+i1 != 8? : " << (i2 != 8) << std::endl;
	std::cout << "1+i1 < 8? : " << (i2 < 8) << std::endl;
	std::cout << "1+i1 > 8? : " << (i2 > 8) << std::endl;
	std::cout << "1+i1 <= 8? : " << (i2 <= 8) << std::endl;
	std::cout << "1+i1 >= 8? : " << (i2 >= 8) << std::endl;
	mpfr_t x;
	mpfr_init (x);
	mpfr_set_si (x, 1, GMP_RNDN);
	i2 = x + i1;
	std::cout << "(mpfr_t)1 + i1 =\t" << i2 << std::endl;
	mpz_t y;
	mpz_init_set_si (y, 1);
	i2 = y + i1;
	std::cout << "(mpz_t)1 + i1 =\t" << i2 << std::endl;
	CGAL::Gmpz z(1);
	i2 = z + i1;
	std::cout << "(Gmpz)1 + i1 =\t" << i2 << std::endl;
}

int main () {
	test_intervals ();
	return 0;
}
