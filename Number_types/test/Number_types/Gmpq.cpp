#include <CGAL/basic.h>

#ifdef CGAL_USE_GMP

#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>

#ifdef CGAL_USE_GMPXX
#include <CGAL/gmpxx.h>
#endif

#include <CGAL/Interval_nt.h>
#include <cassert>

#include <sstream>
#include <vector>

typedef CGAL::Gmpz Gmpz;
typedef CGAL::Gmpq Gmpq;

template < typename NT >
void test_overflow_to_interval(const NT&)
{
  std::cout << "Tests if to_interval() handles overflow nicely or not."
            << std::endl;

  NT val = 2;
  for (int i=0; i<20; ++i) {
    val = val*val;
    CGAL::Interval_nt<> inter = CGAL_NTS to_interval(val);
    CGAL::Interval_nt<> minter = CGAL_NTS to_interval(-val);
    assert(CGAL::is_valid(inter));
    assert(CGAL_NTS is_finite(inter.inf()));
    assert(CGAL::is_valid(minter));
    assert(CGAL_NTS is_finite(minter.sup()));
  }
}

void test_overflow_to_double()
{
  std::cout << "Tests if to_double(Gmpq) overflows or not." << std::endl;

  Gmpq val = Gmpq(1)/2;
  for (int i=0; i<3000; ++i) {
    // std::cout << CGAL_NTS to_double(val) << std::endl;
    // std::cout << val.numerator() << " , " << val.denominator() << std::endl;
    val = val * (1LL<<16);
    val = val / (1<<16);
  }
  assert(CGAL_NTS to_double(val) == 0.5);
}

typedef std::pair<std::string, Gmpq> Test_pair;
typedef std::vector<Test_pair> Test_set;
Test_set make_derived_tests (const Test_pair& pair) {
  // generate exponents, signs, and appended letter
  Test_set result;
  std::string input = pair.first;
  Gmpq should_be = pair.second;
  std::string minus = "-";
  for (std::string l = ""; l != "ff"; l += "f") {
    // append char 'f' to test whether it is correctly not eaten
    for (std::string sign = ""; sign != "done";) {
      // prepend signs "", "+", "-"	
      Gmpq s = should_be;
      if (sign == "-") s = -s;
      result.push_back(Test_pair(sign+input+l , s));
      result.push_back(Test_pair(sign+input+"e-2"+l , s/100));
      result.push_back(Test_pair(sign+input+"e1"+l, s*10));
      result.push_back(Test_pair(sign+input+"e+3"+l, s*1000));
      result.push_back(Test_pair(sign+input+"e0"+l , s));
      if (sign == "-") sign = "done";
      if (sign == "+") sign = "-";
      if (sign == "") sign = "+";
    }
  }
  return result;
}

void test_input_from_float()
{
  std::cout << "Tests Gmpq input from floats." << std::endl;

  Test_set test_set;
  // nonnegative integers
  test_set.push_back (Test_pair (std::string ("123"), Gmpq(123,1)));
  test_set.push_back (Test_pair (std::string ("0"), Gmpq(0,1)));

  // nonnegative floats, with or without digits before/after .
  test_set.push_back (Test_pair (std::string ("0.0"), Gmpq(0,1)));
  test_set.push_back (Test_pair (std::string ("0.123"), Gmpq(123,1000)));
  test_set.push_back (Test_pair (std::string ("00.1"), Gmpq(1,10)));
  test_set.push_back (Test_pair (std::string ("01.1"), Gmpq(11,10)));
  test_set.push_back (Test_pair (std::string ("0."), Gmpq(0,1)));
  test_set.push_back (Test_pair (std::string (".0"), Gmpq(0,1)));
  test_set.push_back (Test_pair (std::string ("12.34"), Gmpq(1234,100)));
  test_set.push_back (Test_pair (std::string ("0.56"), Gmpq(56,100)));
  test_set.push_back (Test_pair (std::string (".78"), Gmpq(78,100)));
  test_set.push_back (Test_pair (std::string ("90."), Gmpq(90,1)));
  test_set.push_back (Test_pair (std::string ("90.0"), Gmpq(90,1)));
  test_set.push_back (Test_pair (std::string ("7.001"), Gmpq(7001,1000)));

  // exponents and signs are automatically generated in
  // make_derived_tests

  // now the actual test
  std::cout << " Running " << test_set.size()  << " master tests..."
	    << std::endl;
  for (unsigned int i=0; i<test_set.size(); ++i) {
    Test_set derived = make_derived_tests (test_set[i]);
    std::cout << " Running " <<derived.size() << " derived tests..."
	      << std::endl;
    for (unsigned int j=0; j <derived.size(); ++j) {
      std::stringstream is(derived[j].first);
      Gmpq should_be (derived[j].second);
      Gmpq q; is >> q;
      assert(!is.fail());
      if (!is.eof()) assert(is.peek()=='f');
      if (q != should_be)
	std::cout << " Error: " << derived[j].first + "  read as "
		  << q << std::endl;
      assert(q == should_be);
    }
  }
}

template<class NumberType>
int test_operators(){
        typedef CGAL::Gmpq      Gmpq;
        typedef NumberType      NT;
        Gmpq a(5.5);
        NT b(2);
        if((-a)==(-5.5)&&(a+b)==7.5&&(a-b)==3.5&&(a*b)==11&&(a/b)==2.75)
                return 0;
        else{
                std::cerr<<"failed!"<<std::endl;
                exit(-1);
        }
}

#define TEST(_string,_code) \
        std::cerr<<"testing "<<_string<<": "<<std::flush; \
        _code; \
        std::cerr<<"OK"<<std::endl;

int main() {

  // Added by Bernd Gaertner
  test_input_from_float();

  // Added by Daniel Russel
  std::cout << "Testing IO" << std::endl;
  {
    const char *buf="12345678";
    Gmpz z;
    std::istringstream iss(buf);
    iss >> z;
    assert(!iss.fail());
    std::cout << z << "== 12345678" << std::endl;
    assert(z== Gmpz(12345678));
  }
  {
    const char *buf="-65758345";
    Gmpz z;
    std::istringstream iss(buf);
    iss >> z;
    assert(!iss.fail());
    std::cout << z << "== -65758345" << std::endl;
    assert(z== Gmpz(-65758345));
  }
  {
    const char *buf="  -  65758345";
    Gmpz z;
    std::istringstream iss(buf);
    iss >> z;
    assert(!iss.fail());
    std::cout << z << " == -65758345" << std::endl;
    assert(z== Gmpz(-65758345));
  }
  {
    const char *buf="12345678a";
    Gmpz z;
    std::istringstream iss(buf);
    char c;
    iss >> z >> c;
    assert(!iss.fail());
    std::cout << z << " == 12345678" << std::endl;
    assert(z== Gmpz(12345678));
    assert(c=='a');
  }

  {
    const char *buf="asadf";
    Gmpz z(12);
    std::istringstream iss(buf);
    iss >> z;
    std::cout << z << " fails with 12" << std::endl;
    assert(iss.fail());
    assert(z== Gmpz(12));
  }

  {
    const char *buf="-asadf";
    Gmpz z(12);
    std::istringstream iss(buf);
    iss >> z;
    std::cout << z << " fails with 12" << std::endl;
    assert(iss.fail());
    assert(z== Gmpz(12));
  }
  {
    const char *buf="";
    Gmpz z(12);
    std::istringstream iss(buf);
    iss >> z;
    std::cout << z << " fails with 12" << std::endl;
    assert(iss.fail());
    assert(z== Gmpz(12));
  }
  {
    const char *buf="100/1";
    Gmpq iot;
    std::istringstream iss(buf);
    iss >> iot;
    std::cout << iot << "== 100 " << std::endl;
    assert(!iss.fail());
    assert(iot== Gmpq(100.0));
  }

  Gmpq q;
  Gmpq q1(12);
  Gmpq q1bis(12u); // test construction from unsigned int
  // For the moment Gmpq is not interoperable with
  // unsigned int: the following assertions would fail to compile.
  //   assert(q1bis == 12);
  //   assert(q1bis == 12u);
  //   assert(q1bis >= 12u);
  Gmpq q2(3.1415);
  Gmpz z1(1), z2(2);
  Gmpq q3(z1);
  Gmpq q4(z1, z2);
  assert(q4.numerator() == Gmpz(1));
  assert(q4.denominator() == Gmpz(2));
  Gmpq q5(q4);

  Gmpq qi1(0,1); // int int
  Gmpq qi2(3, -3);
  assert(qi2.numerator() == Gmpz(-1)); // because Gmpq is canonicalized
  Gmpq qi3(-3, 3);
  assert(qi2 == qi3);

  Gmpq q6((signed long)1, (unsigned long)2);
  Gmpq q7((unsigned long)1, (unsigned long)2);

  Gmpq qs("-100/-111", 2);
  assert(qs.numerator() == Gmpz(4));
  assert(qs.denominator() == Gmpz(7));

  Gmpq qs2("100/1");
  assert(qs2.numerator() == Gmpz(100));
  assert(qs2.denominator() == Gmpz(1));

  test_overflow_to_double();
  test_overflow_to_interval(Gmpz());
  test_overflow_to_interval(Gmpq());

#ifdef CGAL_USE_GMPXX
  test_overflow_to_interval(mpz_class());
  test_overflow_to_interval(mpq_class());
#endif

  // test operators Gmpq/Gmpfr (added by Luis)
  TEST("operators Gmpq",test_operators<CGAL::Gmpq>();)
  TEST("operators int",test_operators<int>();)
  TEST("operators long",test_operators<long>();)
  TEST("operators double",test_operators<double>();)
  TEST("operators Gmpz",test_operators<CGAL::Gmpz>();)
  TEST("operators Gmpfr",test_operators<CGAL::Gmpfr>();)

  return 0;
}
#undef TEST
#else
int main()
{
  return 0;
}
#endif //CGAL_USE_GMP
