// STL.
#include <array>
#include <limits>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>

// Boost.
#include <boost/type_index.hpp>

// Kernels.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

// CGAL.
#include <CGAL/Random.h>
#include <CGAL/Quotient.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Lazy_exact_nt.h>

#ifdef CGAL_USE_BOOST_MP

#if false // https://github.com/boostorg/multiprecision/issues/370

void test_boost_eval_lehmer() {

  const boost::multiprecision::cpp_int a("500000000000000052504760255204421627079393309355027816932345132815919505535709229444276879024105562954502314530690391078574434507015318513443905076213688875017942541908041275407131568575177172639474548726709751235383681696449966404295647940685784470144122251803020020951078103818191513659921807053133698549053838430992170843235673537548059987539601671975279280846041564435631581262016246808786828637048154067265620710396778995313534536353760281048487250661054626168637371167135426013683337484254647996964562455566714879467026196409913165805735073230830136024016362543811769017875638974011487488573436");
  const boost::multiprecision::cpp_int b("1500000000000000157514280765613264881238179928065083450797035398447758516607127688332830637072316688863506943592071173235723303521045955540331715228641066625053827625724123826221394705725531517918423646180129253706151045089349899212886943822057353410432366755409060062853234311454574540979765421159386595647161515292215193506006556519037965168192736708179557957863203557666055574947146355487693991882510747766220045897624670399027877365714431356466054500731862264092476764347207739651025585146903094168986610767496468412336047796468657032646893153521091155634158263410282629846280069312485301157888001");
  const boost::multiprecision::cpp_int r = boost::multiprecision::gcd(a, b);
}

#endif // issue 370

void test_minimal_boost_gcd() {

  boost::multiprecision::cpp_int u = 1;
  for (unsigned i = 1; i <= 50; ++i) {
    u *= i;
  }
  std::cout << "u: " << u << std::endl;

  boost::multiprecision::cpp_int v = 1;
  for (unsigned i = 1; i <= 100; ++i) {
    v *= i;
  }
  std::cout << "v: " << v << std::endl;

  // const auto r = boost::multiprecision::gcd(u, v); // fail
  const boost::multiprecision::cpp_int r = boost::multiprecision::gcd(u, v); // pass
  std::cout << "r: " << r << std::endl;

  u = u / r;
  v = v / r;

  std::cout << "new u: " << u << std::endl;
  std::cout << "new v: " << v << std::endl;
}

#if false // _MM_ROUND_UP is not available on all platforms

void test_minimal_nextafter() {

  _MM_SET_ROUNDING_MODE(_MM_ROUND_UP); // fail
  // _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST); // pass

  const boost::multiprecision::cpp_int x("1312729512902970206056841780066779136");

  double i = x.template convert_to<double>();
  double s = i;

  const double inf = std::numeric_limits<double>::infinity();
  assert(i != inf && s != inf);
  const int cmp = x.compare(i);
  if (cmp > 0) {
    s = nextafter(s, +inf);
    assert(x.compare(s) < 0);
  } else if (cmp < 0) {
    i = nextafter(i, -inf);
    assert(x.compare(i) > 0);
  }
}

#endif // test_minimal_nextafter

void test_to_interval_boost() {

  using NT = boost::multiprecision::cpp_int;
  using Quotient = CGAL::Quotient<NT>;
  using Traits = CGAL::Real_embeddable_traits<Quotient>;
  using Interval = typename Traits::To_interval;

  NT n, d;
  Quotient x;
  double i, s;

  n = NT("-15284404573383541");
  d = NT("4503599627370496");

  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: -3.3938195750112902793" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: -3.3938195750112898352" << std::endl;
  std::cout << std::endl;

  // Results for current tight using cpp_rational.
  assert(i == -3.3938195750112902793);
  assert(s == -3.3938195750112898352);
}

// In assert, we use values from impl2.
void test_to_interval_tight_rational_1() {

  // In green, we compare to impl2.
  #define TESTCASE_RAT_10 // pass all three
  #define TESTCASE_RAT_11 // impl1: fails tightness (sup is larger, s = 9.3488310472396616291)
  #define TESTCASE_RAT_12 // pass all three
  #define TESTCASE_RAT_13 // pass all three
  #define TESTCASE_RAT_14 // pass all three
  #define TESTCASE_RAT_15 // pass all three
  #define TESTCASE_RAT_16 // pass all three
  #define TESTCASE_RAT_17 // pass all three
  #define TESTCASE_RAT_18 // pass all three
  #define TESTCASE_RAT_19 // pass all three

  using NT = boost::multiprecision::cpp_int;
  using Quotient = CGAL::Quotient<NT>;
  using Traits = CGAL::Real_embeddable_traits<Quotient>;
  using Interval = typename Traits::To_interval;

  // std::cout << std::endl;
  // std::cout << boost::typeindex::type_id<Quotient>() << std::endl;
  // std::cout << std::endl;
  // std::cout << boost::typeindex::type_id<typename Traits::Type>() << std::endl;
  // std::cout << std::endl;

  NT n, d;
  Quotient x;
  double i, s;

  std::cout << std::endl;
  std::cout << "- T1 testing tight interval for rationals ..." << std::endl;
  std::cout << std::endl;

  #ifdef TESTCASE_RAT_10 // small numbers

  std::cout << "=============" << std::endl;
  std::cout << "CASE0 RESULT:" << std::endl;
  std::cout << "=============" << std::endl;

  n = NT("39792587355159975");
  d = NT("140737488355328");
  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 282.7433388230813307" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 282.74333882308138755" << std::endl;
  std::cout << std::endl;

  assert(i == 282.7433388230813307);
  assert(s == 282.74333882308138755);

  #endif

  #ifdef TESTCASE_RAT_11 // large numbers

  std::cout << "=============" << std::endl;
  std::cout << "CASE1 RESULT:" << std::endl;
  std::cout << "=============" << std::endl;

  n = NT("772537196160711547532081795586792063331305895970601529435744397743492241616327030886637827664482971614281724796166908515292029740442872965475211471498392497954317530347232852540146110053764627070672243390766540271554856759037331142360111552286202392826786995364211101723592791550906796165626083442695020580821188398298798456115881346136681033873");
  d = NT("82634630175374856683315372867724319098240552701588533218371381248009342768269285501674184091886435054368116496214846441734481770666205690731018817430937185570378353100803926136323598244976110318516454816403989543192819758059431171537258117598056453283568595627159988837663160716950017789671313834717457946818990093589809113731838629064768225280");
  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 9.3488310472396563" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 9.3488310472396580764" << std::endl;
  std::cout << std::endl;

  // Results for current tight using cpp_rational.
  assert(i == 9.3488310472396563);
  assert(s == 9.3488310472396580764);

  #endif

  #ifdef TESTCASE_RAT_12

  std::cout << "=============" << std::endl;
  std::cout << "CASE2 RESULT:" << std::endl;
  std::cout << "=============" << std::endl;

  n = NT("772537196160711547532081795586792063331305895970601529435744397743492241616327030886637827664482971614281724796166908515292029740442872965475211471498392497954317530347232852540146110053764627070672243390766540271554856759037331142360111552286202392826786995364211101723592791550906796165626083442695020580821188398298798456115881346136681033873");
  d = NT("1");
  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 1.7976931348623157081e+308" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: inf" << std::endl;
  std::cout << std::endl;

  assert(i == (std::numeric_limits<double>::max)());
  assert(s == std::numeric_limits<double>::infinity());

  #endif

  #ifdef TESTCASE_RAT_13

  std::cout << "=============" << std::endl;
  std::cout << "CASE3 RESULT:" << std::endl;
  std::cout << "=============" << std::endl;

  n = NT("1");
  d = NT("772537196160711547532081795586792063331305895970601529435744397743492241616327030886637827664482971614281724796166908515292029740442872965475211471498392497954317530347232852540146110053764627070672243390766540271554856759037331142360111552286202392826786995364211101723592791550906796165626083442695020580821188398298798456115881346136681033873");
  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 0.0 or higher" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 0.0 or higher" << std::endl;
  std::cout << std::endl;

  assert(i >= 0.0 && i <= (std::numeric_limits<double>::min)() * 2.0);
  assert(s >= 0.0 && s <= (std::numeric_limits<double>::min)() * 2.0);
  assert(i <= s);

  #endif

  #ifdef TESTCASE_RAT_14

  std::cout << "=============" << std::endl;
  std::cout << "CASE4 RESULT:" << std::endl;
  std::cout << "=============" << std::endl;

  n = NT("10");
  d = NT("10");
  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 1" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 1" << std::endl;
  std::cout << std::endl;

  assert(i == 1.0);
  assert(s == 1.0);

  #endif

  #ifdef TESTCASE_RAT_15

  std::cout << "=============" << std::endl;
  std::cout << "CASE5 RESULT:" << std::endl;
  std::cout << "=============" << std::endl;

  n = NT("1");
  d = NT("6");
  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 0.16666666666666665741" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 0.16666666666666668517" << std::endl;
  std::cout << std::endl;

  assert(i == 0.16666666666666665741);
  assert(s == 0.16666666666666668517);

  #endif

  #ifdef TESTCASE_RAT_16

  std::cout << "=============" << std::endl;
  std::cout << "CASE6 RESULT:" << std::endl;
  std::cout << "=============" << std::endl;

  n = +NT("6");
  d = +NT("3");
  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: " << 6.0 / 3.0 << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: " << 6.0 / 3.0 << std::endl;
  std::cout << std::endl;

  assert(i == 2.0);
  assert(s == 2.0);

  #endif

  #ifdef TESTCASE_RAT_17

  std::cout << "=============" << std::endl;
  std::cout << "CASE7 RESULT:" << std::endl;
  std::cout << "=============" << std::endl;

  n = -NT("1");
  d = -NT("2");
  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: " << 1.0 / 2.0 << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: " << 1.0 / 2.0 << std::endl;
  std::cout << std::endl;

  assert(i == 1.0 / 2.0);
  assert(s == 1.0 / 2.0);

  #endif

  #ifdef TESTCASE_RAT_18

  std::cout << "=============" << std::endl;
  std::cout << "CASE8 RESULT:" << std::endl;
  std::cout << "=============" << std::endl;

  n = -NT("1");
  d = +NT("3");
  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: -0.33333333333333337034" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: -0.33333333333333331483" << std::endl;
  std::cout << std::endl;

  assert(i == -0.33333333333333337034);
  assert(s == -0.33333333333333331483);

  #endif

  #ifdef TESTCASE_RAT_19 // small numbers, num > 0 and den < 0

  std::cout << "=============" << std::endl;
  std::cout << "CASE9 RESULT:" << std::endl;
  std::cout << "=============" << std::endl;

  n = +NT("39792587355159975");
  d = -NT("140737488355328");
  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: -282.74333882308138755" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: -282.7433388230813307" << std::endl;
  std::cout << std::endl;

  assert(i == -282.74333882308138755);
  assert(s == -282.7433388230813307);

  #endif

  std::cout << "---SUCCESS ALL---" << std::endl;
  std::cout << std::endl;
}

// In assert, we use values from impl2.
void test_to_interval_tight_rational_2() {

  // In green, we compare to impl2.
  #define TESTCASE_RAT_20  // pass all three
  #define TESTCASE_RAT_21  // impl1: fails tightness (sup is larger, s = 0.43464565325999998668)
  #define TESTCASE_RAT_22  // pass all three
  #define TESTCASE_RAT_23  // pass all three
  #define TESTCASE_RAT_24  // pass all three
  #define TESTCASE_RAT_25  // pass all three
  #define TESTCASE_RAT_26  // pass all three
  #define TESTCASE_RAT_27  // pass all three
  #define TESTCASE_RAT_28  // pass all three
  #define TESTCASE_RAT_29  // pass all three
  #define TESTCASE_RAT_210 // impl1, impl3: fails tightness (sup is larger, s = 5.9425938166208590782e+26)
  #define TESTCASE_RAT_211 // impl1, impl3: fails tightness (inf is smaller, i = 3602879701896396.5, sup is larger, s = 3602879701896398)

  using NT = boost::multiprecision::cpp_int;
  using Quotient = CGAL::Quotient<NT>;
  using Traits = CGAL::Real_embeddable_traits<Quotient>;
  using Interval = typename Traits::To_interval;

  NT n, d;
  Quotient x;
  double i, s;

  std::cout << std::endl;
  std::cout << "- T2 testing tight interval for rationals ..." << std::endl;
  std::cout << std::endl;

  #ifdef TESTCASE_RAT_20

  std::cout << "TEST 0" << std::endl; // num, case 2

  n = NT("-15284404573383541");
  d = NT("4503599627370496");

  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: -3.3938195750112902793" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: -3.3938195750112898352" << std::endl;
  std::cout << std::endl;

  assert(i == -3.3938195750112902793);
  assert(s == -3.3938195750112898352);

  #endif

  #ifdef TESTCASE_RAT_21

  std::cout << "TEST 1" << std::endl; // num, case 4

  const double nn = 0.43464565326;
  x = Quotient(nn);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 0.43464565325999998668" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 0.43464565325999998668" << std::endl;
  std::cout << std::endl;

  assert(i == 0.43464565325999998668);
  assert(s == 0.43464565325999998668);
  assert(i == s);

  #endif

  #ifdef TESTCASE_RAT_22

  std::cout << "TEST 2" << std::endl; // num, case 4

  n = NT("1");
  d = NT("2");

  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 0.5" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 0.5" << std::endl;
  std::cout << std::endl;

  assert(i == 0.5);
  assert(s == 0.5);

  #endif

  #ifdef TESTCASE_RAT_23

  std::cout << "TEST 3" << std::endl; // shift = 0

  n = NT("7725371961607115");
  d = NT("1");

  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 7725371961607115" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 7725371961607115" << std::endl;
  std::cout << std::endl;

  assert(i == 7725371961607115);
  assert(s == 7725371961607115);

  #endif

  #ifdef TESTCASE_RAT_24

  std::cout << "TEST 4" << std::endl;

  n = NT("772537196160711547532081795586792063331305895970601529435744397743492241616327030886637827664482971614281724796166908515292029740442872965475211471498392497954317530347232852540146110053764627070672243390766540271554856759037331142360111552286202392826786995364211101723592791550906796165626083442695020580821188398298798456115881346136681033873");
  d = NT("1");

  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 1.7976931348623157081e+308" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: inf" << std::endl;
  std::cout << std::endl;

  assert(i == (std::numeric_limits<double>::max)());
  assert(s == std::numeric_limits<double>::infinity());

  #endif

  #ifdef TESTCASE_RAT_25

  std::cout << "TEST 5" << std::endl;

  n = NT("1");
  d = NT("772537196160711547532081795586792063331305895970601529435744397743492241616327030886637827664482971614281724796166908515292029740442872965475211471498392497954317530347232852540146110053764627070672243390766540271554856759037331142360111552286202392826786995364211101723592791550906796165626083442695020580821188398298798456115881346136681033873");

  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 0.0 or higher" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 0.0 or higher" << std::endl;
  std::cout << std::endl;

  assert(i >= 0.0 && i <= (std::numeric_limits<double>::min)() * 2.0);
  assert(s >= 0.0 && s <= (std::numeric_limits<double>::min)() * 2.0);
  assert(i <= s);

  #endif

  #ifdef TESTCASE_RAT_26

  std::cout << "TEST 6" << std::endl; // case shift > 0 && p_bits = 51, subcase1

  n = NT("1");
  d = NT("10");

  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 0.099999999999999991673" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 0.10000000000000000555" << std::endl;
  std::cout << std::endl;

  assert(i == 0.099999999999999991673);
  assert(s == 0.10000000000000000555);

  #endif

  #ifdef TESTCASE_RAT_27

  std::cout << "TEST 7" << std::endl; // non representable double

  n = NT("1");
  d = NT("3");

  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 0.33333333333333331483" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 0.33333333333333337034" << std::endl;
  std::cout << std::endl;

  assert(i == 0.33333333333333331483);
  assert(s == 0.33333333333333337034);

  #endif

  #ifdef TESTCASE_RAT_28

  std::cout << "TEST 8" << std::endl; // fails assertion (q_bits == num_dbl_digits || r != 0)

  n = NT("21");
  d = NT("3");

  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 7" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 7" << std::endl;
  std::cout << std::endl;

  assert(i == 7);
  assert(s == 7);

  #endif

  #ifdef TESTCASE_RAT_29

  std::cout << "TEST 9" << std::endl; // case shift > 0 && p_bits = 51, subcase3

  n = NT("17");
  d = NT("3");

  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 5.6666666666666660745" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 5.6666666666666669627" << std::endl;
  std::cout << std::endl;

  assert(i == 5.6666666666666660745);
  assert(s == 5.6666666666666669627);

  #endif

  #ifdef TESTCASE_RAT_210

  std::cout << "TEST 10" << std::endl; // case shift < 0 && p_bits = 51

  n = NT("7725371961607115475320817955");
  d = NT("13");

  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 5.9425938166208577038e+26" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 5.942593816620858391e+26" << std::endl;
  std::cout << std::endl;

  assert(i == 5.9425938166208577038e+26);
  assert(s == 5.942593816620858391e+26);

  #endif

  #ifdef TESTCASE_RAT_211

  std::cout << "TEST 11" << std::endl; // case shift = 0 && p_bits = 51 && cmp = 0, subcase3

  n = NT("36028797018963975");
  d = NT("10");

  x = Quotient(n, d);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 3602879701896397.5" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 3602879701896397.5" << std::endl;
  std::cout << std::endl;

  assert(i == 3602879701896397.5);
  assert(s == 3602879701896397.5);

  #endif

  std::cout << "---SUCCESS ALL---" << std::endl;
  std::cout << std::endl;
}

void test_to_interval_tight_integer() {

  #define TESTCASE_INT_0
  #define TESTCASE_INT_1
  #define TESTCASE_INT_2
  #define TESTCASE_INT_3

  using NT = boost::multiprecision::cpp_int;
  using Traits = CGAL::Real_embeddable_traits<NT>;
  using Interval = typename Traits::To_interval;

  // std::cout << std::endl;
  // std::cout << boost::typeindex::type_id<NT>() << std::endl;
  // std::cout << std::endl;
  // std::cout << boost::typeindex::type_id<typename Traits::Type>() << std::endl;
  // std::cout << std::endl;

  NT x;
  double i, s;

  std::cout << std::endl;
  std::cout << "- T testing tight interval for integers ..." << std::endl;
  std::cout << std::endl;

  #ifdef TESTCASE_INT_0

  std::cout << "=============" << std::endl;
  std::cout << "CASE0 RESULT:" << std::endl;
  std::cout << "=============" << std::endl;

  x = NT("0");
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 0" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 0" << std::endl;
  std::cout << std::endl;

  assert(i == 0.0);
  assert(s == 0.0);

  #endif

  #ifdef TESTCASE_INT_1

  std::cout << "=============" << std::endl;
  std::cout << "CASE1 RESULT:" << std::endl;
  std::cout << "=============" << std::endl;

  x = NT("5");
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 5" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: 5" << std::endl;
  std::cout << std::endl;

  assert(i == 5.0);
  assert(s == 5.0);

  #endif

  #ifdef TESTCASE_INT_2

  std::cout << "=============" << std::endl;
  std::cout << "CASE2 RESULT:" << std::endl;
  std::cout << "=============" << std::endl;

  x = NT(-12);
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: -12" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: -12" << std::endl;
  std::cout << std::endl;

  assert(i == -12.0);
  assert(s == -12.0);

  #endif

  #ifdef TESTCASE_INT_3

  std::cout << "=============" << std::endl;
  std::cout << "CASE3 RESULT:" << std::endl;
  std::cout << "=============" << std::endl;

  x = NT("772537196160711547532081795586792063331305895970601529435744397743492241616327030886637827664482971614281724796166908515292029740442872965475211471498392497954317530347232852540146110053764627070672243390766540271554856759037331142360111552286202392826786995364211101723592791550906796165626083442695020580821188398298798456115881346136681033873");
  std::tie(i, s) = Interval()(x);

  std::cout << std::endl;
  std::cout << "inf: " << i << std::endl;
  std::cout << "ref: 1.7976931348623157081e+308" << std::endl;
  std::cout << "sup: " << s << std::endl;
  std::cout << "ref: inf" << std::endl;
  std::cout << std::endl;

  assert(i == (std::numeric_limits<double>::max)());
  assert(s == std::numeric_limits<double>::infinity());

  #endif

  std::cout << "---SUCCESS ALL---" << std::endl;
  std::cout << std::endl;
}

void test_shift_positive() {
  {
    double d = (1LL << 53) - 1;
    auto shift = std::numeric_limits<double>::max_exponent - (std::numeric_limits<double>::digits);
    auto r = CGAL::Boost_MP_internal::shift_positive_interval({d,d},shift);
    d = ldexp(d,shift);
    assert(!isinf(d));
    assert(r.first == d && d == r.second);
  }
  {
    double d = (1LL << 52);
    auto shift = std::numeric_limits<double>::max_exponent - std::numeric_limits<double>::digits + 1;
    auto r = CGAL::Boost_MP_internal::shift_positive_interval({d,d},shift);
    d = ldexp(d,shift);
    assert(d >= (std::numeric_limits<double>::max)());
    assert(r.first <= d && d <= r.second);
  }
  {
    double d = (1LL << 53) - 1;
    auto shift = std::numeric_limits<double>::min_exponent - std::numeric_limits<double>::digits - 1;
    auto r = CGAL::Boost_MP_internal::shift_positive_interval({d,d},shift);
    d = ldexp(d,shift);
    assert(d <= (std::numeric_limits<double>::min)());
    assert(r.first <= d && d <= r.second);
  }
  {
    double d = (1LL << 53) - 2;
    auto shift = std::numeric_limits<double>::min_exponent - std::numeric_limits<double>::digits - 1;
    auto r = CGAL::Boost_MP_internal::shift_positive_interval({d,d},shift);
    d = ldexp(d,shift);
    assert(d < (std::numeric_limits<double>::min)());
    assert(r.first <= d && d <= r.second);
  }
  {
    double d = (1LL << 52);
    auto shift = std::numeric_limits<double>::min_exponent - std::numeric_limits<double>::digits;
    auto r = CGAL::Boost_MP_internal::shift_positive_interval({d,d},shift);
    d = ldexp(d,shift);
    assert(d == (std::numeric_limits<double>::min)());
    assert(r.first == d && d == r.second);
  }
}
#endif // CGAL_USE_BOOST_MP

int main() {
#ifdef CGAL_USE_BOOST_MP
  test_shift_positive();
#endif

  // Make sure we have the same seed.
  CGAL::get_default_random() = CGAL::Random(0);
  std::cout.precision(20);
  #ifdef CGAL_USE_BOOST_MP

  #if false
  test_boost_eval_lehmer();
  test_minimal_nextafter();
  #endif
  test_minimal_boost_gcd();
  test_to_interval_boost();

  // Test rational types.
  test_to_interval_tight_rational_1();
  test_to_interval_tight_rational_2();

  // Test integer types.
  test_to_interval_tight_integer();

  #endif // CGAL_USE_BOOST_MP
  return EXIT_SUCCESS;
}
