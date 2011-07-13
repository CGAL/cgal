// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2006-2008 Benoit Jacob <jacob.benoit.1@gmail.com>
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// Eigen is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// Eigen is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// Eigen. If not, see <http://www.gnu.org/licenses/>.

#include <cstdlib>
#include <cerrno>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>

#ifdef NDEBUG
#undef NDEBUG
#endif

// bounds integer values for AltiVec
#ifdef __ALTIVEC__
#define EIGEN_MAKING_DOCS
#endif

#ifndef EIGEN_TEST_FUNC
#error EIGEN_TEST_FUNC must be defined
#endif

#define DEFAULT_REPEAT 10

#ifdef __ICC
// disable warning #279: controlling expression is constant
#pragma warning disable 279
#endif

namespace Eigen
{
  static std::vector<std::string> g_test_stack;
  static int g_repeat;
  static unsigned int g_seed;
  static bool g_has_set_repeat, g_has_set_seed;
}

#define EI_PP_MAKE_STRING2(S) #S
#define EI_PP_MAKE_STRING(S) EI_PP_MAKE_STRING2(S)

#define EIGEN_DEFAULT_IO_FORMAT IOFormat(4, 0, "  ", "\n", "", "", "", "")

#ifndef EIGEN_NO_ASSERTION_CHECKING

  namespace Eigen
  {
    static const bool should_raise_an_assert = false;

    // Used to avoid to raise two exceptions at a time in which
    // case the exception is not properly caught.
    // This may happen when a second exceptions is triggered in a destructor.
    static bool no_more_assert = false;
    static bool report_on_cerr_on_assert_failure = true;

    struct eigen_assert_exception
    {
      eigen_assert_exception(void) {}
      ~eigen_assert_exception() { Eigen::no_more_assert = false; }
    };
  }
  // If EIGEN_DEBUG_ASSERTS is defined and if no assertion is triggered while
  // one should have been, then the list of excecuted assertions is printed out.
  //
  // EIGEN_DEBUG_ASSERTS is not enabled by default as it
  // significantly increases the compilation time
  // and might even introduce side effects that would hide
  // some memory errors.
  #ifdef EIGEN_DEBUG_ASSERTS

    namespace Eigen
    {
      namespace internal
      {
        static bool push_assert = false;
      }
      static std::vector<std::string> eigen_assert_list;
    }
    #define eigen_assert(a)                       \
      if( (!(a)) && (!no_more_assert) )     \
      { \
        if(report_on_cerr_on_assert_failure) \
          std::cerr <<  #a << " " __FILE__ << "(" << __LINE__ << ")\n"; \
        Eigen::no_more_assert = true;       \
        throw Eigen::eigen_assert_exception(); \
      }                                     \
      else if (Eigen::internal::push_assert)       \
      {                                     \
        eigen_assert_list.push_back(std::string(EI_PP_MAKE_STRING(__FILE__)" ("EI_PP_MAKE_STRING(__LINE__)") : "#a) ); \
      }

    #define VERIFY_RAISES_ASSERT(a)                                                   \
      {                                                                               \
        Eigen::no_more_assert = false;                                                \
        Eigen::eigen_assert_list.clear();                                                \
        Eigen::internal::push_assert = true;                                                 \
        Eigen::report_on_cerr_on_assert_failure = false;                              \
        try {                                                                         \
          a;                                                                          \
          std::cerr << "One of the following asserts should have been triggered:\n";  \
          for (uint ai=0 ; ai<eigen_assert_list.size() ; ++ai)                           \
            std::cerr << "  " << eigen_assert_list[ai] << "\n";                          \
          VERIFY(Eigen::should_raise_an_assert && # a);                               \
        } catch (Eigen::eigen_assert_exception) {                                        \
          Eigen::internal::push_assert = false; VERIFY(true);                                \
        }                                                                             \
        Eigen::report_on_cerr_on_assert_failure = true;                               \
        Eigen::internal::push_assert = false;                                                \
      }

  #else // EIGEN_DEBUG_ASSERTS
    // see bug 89. The copy_bool here is working around a bug in gcc <= 4.3
    #define eigen_assert(a) \
      if( (!Eigen::internal::copy_bool(a)) && (!no_more_assert) )\
      {                                       \
        Eigen::no_more_assert = true;         \
        if(report_on_cerr_on_assert_failure)  \
          eigen_plain_assert(a);              \
        else                                  \
          throw Eigen::eigen_assert_exception(); \
      }
    #define VERIFY_RAISES_ASSERT(a) {                             \
        Eigen::no_more_assert = false;                            \
        Eigen::report_on_cerr_on_assert_failure = false;          \
        try {                                                     \
          a;                                                      \
          VERIFY(Eigen::should_raise_an_assert && # a);           \
        }                                                         \
        catch (Eigen::eigen_assert_exception&) { VERIFY(true); }     \
        Eigen::report_on_cerr_on_assert_failure = true;           \
      }

  #endif // EIGEN_DEBUG_ASSERTS

  #define EIGEN_USE_CUSTOM_ASSERT

#else // EIGEN_NO_ASSERTION_CHECKING

  #define VERIFY_RAISES_ASSERT(a) {}

#endif // EIGEN_NO_ASSERTION_CHECKING


#define EIGEN_INTERNAL_DEBUGGING
#include <Eigen/QR> // required for createRandomPIMatrixOfRank

static void verify_impl(bool condition, const char *testname, const char *file, int line, const char *condition_as_string)
{
  if (!condition)
  {
    std::cerr << "Test " << testname << " failed in " << file << " (" << line << ")" \
      << std::endl << "    " << condition_as_string << std::endl << std::endl; \
    abort();
  }
}

#define VERIFY(a) verify_impl(a, g_test_stack.back().c_str(), __FILE__, __LINE__, EI_PP_MAKE_STRING(a))

#define VERIFY_IS_EQUAL(a, b) VERIFY(test_is_equal(a, b))
#define VERIFY_IS_APPROX(a, b) VERIFY(test_isApprox(a, b))
#define VERIFY_IS_NOT_APPROX(a, b) VERIFY(!test_isApprox(a, b))
#define VERIFY_IS_MUCH_SMALLER_THAN(a, b) VERIFY(test_isMuchSmallerThan(a, b))
#define VERIFY_IS_NOT_MUCH_SMALLER_THAN(a, b) VERIFY(!test_isMuchSmallerThan(a, b))
#define VERIFY_IS_APPROX_OR_LESS_THAN(a, b) VERIFY(test_isApproxOrLessThan(a, b))
#define VERIFY_IS_NOT_APPROX_OR_LESS_THAN(a, b) VERIFY(!test_isApproxOrLessThan(a, b))

#define VERIFY_IS_UNITARY(a) VERIFY(test_isUnitary(a))

#define CALL_SUBTEST(FUNC) do { \
    g_test_stack.push_back(EI_PP_MAKE_STRING(FUNC)); \
    FUNC; \
    g_test_stack.pop_back(); \
  } while (0)

#ifdef EIGEN_TEST_PART_1
#define CALL_SUBTEST_1(FUNC) CALL_SUBTEST(FUNC)
#else
#define CALL_SUBTEST_1(FUNC)
#endif

#ifdef EIGEN_TEST_PART_2
#define CALL_SUBTEST_2(FUNC) CALL_SUBTEST(FUNC)
#else
#define CALL_SUBTEST_2(FUNC)
#endif

#ifdef EIGEN_TEST_PART_3
#define CALL_SUBTEST_3(FUNC) CALL_SUBTEST(FUNC)
#else
#define CALL_SUBTEST_3(FUNC)
#endif

#ifdef EIGEN_TEST_PART_4
#define CALL_SUBTEST_4(FUNC) CALL_SUBTEST(FUNC)
#else
#define CALL_SUBTEST_4(FUNC)
#endif

#ifdef EIGEN_TEST_PART_5
#define CALL_SUBTEST_5(FUNC) CALL_SUBTEST(FUNC)
#else
#define CALL_SUBTEST_5(FUNC)
#endif

#ifdef EIGEN_TEST_PART_6
#define CALL_SUBTEST_6(FUNC) CALL_SUBTEST(FUNC)
#else
#define CALL_SUBTEST_6(FUNC)
#endif

#ifdef EIGEN_TEST_PART_7
#define CALL_SUBTEST_7(FUNC) CALL_SUBTEST(FUNC)
#else
#define CALL_SUBTEST_7(FUNC)
#endif

#ifdef EIGEN_TEST_PART_8
#define CALL_SUBTEST_8(FUNC) CALL_SUBTEST(FUNC)
#else
#define CALL_SUBTEST_8(FUNC)
#endif

#ifdef EIGEN_TEST_PART_9
#define CALL_SUBTEST_9(FUNC) CALL_SUBTEST(FUNC)
#else
#define CALL_SUBTEST_9(FUNC)
#endif

#ifdef EIGEN_TEST_PART_10
#define CALL_SUBTEST_10(FUNC) CALL_SUBTEST(FUNC)
#else
#define CALL_SUBTEST_10(FUNC)
#endif

#ifdef EIGEN_TEST_PART_11
#define CALL_SUBTEST_11(FUNC) CALL_SUBTEST(FUNC)
#else
#define CALL_SUBTEST_11(FUNC)
#endif

#ifdef EIGEN_TEST_PART_12
#define CALL_SUBTEST_12(FUNC) CALL_SUBTEST(FUNC)
#else
#define CALL_SUBTEST_12(FUNC)
#endif

#ifdef EIGEN_TEST_PART_13
#define CALL_SUBTEST_13(FUNC) CALL_SUBTEST(FUNC)
#else
#define CALL_SUBTEST_13(FUNC)
#endif

#ifdef EIGEN_TEST_PART_14
#define CALL_SUBTEST_14(FUNC) CALL_SUBTEST(FUNC)
#else
#define CALL_SUBTEST_14(FUNC)
#endif

#ifdef EIGEN_TEST_PART_15
#define CALL_SUBTEST_15(FUNC) CALL_SUBTEST(FUNC)
#else
#define CALL_SUBTEST_15(FUNC)
#endif

#ifdef EIGEN_TEST_PART_16
#define CALL_SUBTEST_16(FUNC) CALL_SUBTEST(FUNC)
#else
#define CALL_SUBTEST_16(FUNC)
#endif

namespace Eigen {

template<typename T> inline typename NumTraits<T>::Real test_precision() { return NumTraits<T>::dummy_precision(); }
template<> inline float test_precision<float>() { return 1e-3f; }
template<> inline double test_precision<double>() { return 1e-6; }
template<> inline float test_precision<std::complex<float> >() { return test_precision<float>(); }
template<> inline double test_precision<std::complex<double> >() { return test_precision<double>(); }
template<> inline long double test_precision<long double>() { return 1e-6; }

inline bool test_isApprox(const int& a, const int& b)
{ return internal::isApprox(a, b, test_precision<int>()); }
inline bool test_isMuchSmallerThan(const int& a, const int& b)
{ return internal::isMuchSmallerThan(a, b, test_precision<int>()); }
inline bool test_isApproxOrLessThan(const int& a, const int& b)
{ return internal::isApproxOrLessThan(a, b, test_precision<int>()); }

inline bool test_isApprox(const float& a, const float& b)
{ return internal::isApprox(a, b, test_precision<float>()); }
inline bool test_isMuchSmallerThan(const float& a, const float& b)
{ return internal::isMuchSmallerThan(a, b, test_precision<float>()); }
inline bool test_isApproxOrLessThan(const float& a, const float& b)
{ return internal::isApproxOrLessThan(a, b, test_precision<float>()); }
inline bool test_isApprox(const double& a, const double& b)
{ return internal::isApprox(a, b, test_precision<double>()); }

inline bool test_isMuchSmallerThan(const double& a, const double& b)
{ return internal::isMuchSmallerThan(a, b, test_precision<double>()); }
inline bool test_isApproxOrLessThan(const double& a, const double& b)
{ return internal::isApproxOrLessThan(a, b, test_precision<double>()); }

inline bool test_isApprox(const std::complex<float>& a, const std::complex<float>& b)
{ return internal::isApprox(a, b, test_precision<std::complex<float> >()); }
inline bool test_isMuchSmallerThan(const std::complex<float>& a, const std::complex<float>& b)
{ return internal::isMuchSmallerThan(a, b, test_precision<std::complex<float> >()); }

inline bool test_isApprox(const std::complex<double>& a, const std::complex<double>& b)
{ return internal::isApprox(a, b, test_precision<std::complex<double> >()); }
inline bool test_isMuchSmallerThan(const std::complex<double>& a, const std::complex<double>& b)
{ return internal::isMuchSmallerThan(a, b, test_precision<std::complex<double> >()); }

inline bool test_isApprox(const long double& a, const long double& b)
{
    bool ret = internal::isApprox(a, b, test_precision<long double>());
    if (!ret) std::cerr
        << std::endl << "    actual   = " << a
        << std::endl << "    expected = " << b << std::endl << std::endl;
    return ret;
}

inline bool test_isMuchSmallerThan(const long double& a, const long double& b)
{ return internal::isMuchSmallerThan(a, b, test_precision<long double>()); }
inline bool test_isApproxOrLessThan(const long double& a, const long double& b)
{ return internal::isApproxOrLessThan(a, b, test_precision<long double>()); }

template<typename Type1, typename Type2>
inline bool test_isApprox(const Type1& a, const Type2& b)
{
  return a.isApprox(b, test_precision<typename Type1::Scalar>());
}

// The idea behind this function is to compare the two scalars a and b where
// the scalar ref is a hint about the expected order of magnitude of a and b.
// Therefore, if for some reason a and b are very small compared to ref,
// we won't issue a false negative.
// This test could be: abs(a-b) <= eps * ref
// However, it seems that simply comparing a+ref and b+ref is more sensitive to true error.
template<typename Scalar,typename ScalarRef>
inline bool test_isApproxWithRef(const Scalar& a, const Scalar& b, const ScalarRef& ref)
{
  return test_isApprox(a+ref, b+ref);
}

template<typename Derived1, typename Derived2>
inline bool test_isMuchSmallerThan(const MatrixBase<Derived1>& m1,
                                   const MatrixBase<Derived2>& m2)
{
  return m1.isMuchSmallerThan(m2, test_precision<typename internal::traits<Derived1>::Scalar>());
}

template<typename Derived>
inline bool test_isMuchSmallerThan(const MatrixBase<Derived>& m,
                                   const typename NumTraits<typename internal::traits<Derived>::Scalar>::Real& s)
{
  return m.isMuchSmallerThan(s, test_precision<typename internal::traits<Derived>::Scalar>());
}

template<typename Derived>
inline bool test_isUnitary(const MatrixBase<Derived>& m)
{
  return m.isUnitary(test_precision<typename internal::traits<Derived>::Scalar>());
}

template<typename T, typename U>
bool test_is_equal(const T& actual, const U& expected)
{
    if (actual==expected)
        return true;
    // false:
    std::cerr
        << std::endl << "    actual   = " << actual
        << std::endl << "    expected = " << expected << std::endl << std::endl;
    return false;
}

/** Creates a random Partial Isometry matrix of given rank.
  *
  * A partial isometry is a matrix all of whose singular values are either 0 or 1.
  * This is very useful to test rank-revealing algorithms.
  */
template<typename MatrixType>
void createRandomPIMatrixOfRank(typename MatrixType::Index desired_rank, typename MatrixType::Index rows, typename MatrixType::Index cols, MatrixType& m)
{
  typedef typename internal::traits<MatrixType>::Index Index;
  typedef typename internal::traits<MatrixType>::Scalar Scalar;
  enum { Rows = MatrixType::RowsAtCompileTime, Cols = MatrixType::ColsAtCompileTime };

  typedef Matrix<Scalar, Dynamic, 1> VectorType;
  typedef Matrix<Scalar, Rows, Rows> MatrixAType;
  typedef Matrix<Scalar, Cols, Cols> MatrixBType;

  if(desired_rank == 0)
  {
    m.setZero(rows,cols);
    return;
  }

  if(desired_rank == 1)
  {
    // here we normalize the vectors to get a partial isometry
    m = VectorType::Random(rows).normalized() * VectorType::Random(cols).normalized().transpose();
    return;
  }

  MatrixAType a = MatrixAType::Random(rows,rows);
  MatrixType d = MatrixType::Identity(rows,cols);
  MatrixBType  b = MatrixBType::Random(cols,cols);

  // set the diagonal such that only desired_rank non-zero entries reamain
  const Index diag_size = std::min(d.rows(),d.cols());
  if(diag_size != desired_rank)
    d.diagonal().segment(desired_rank, diag_size-desired_rank) = VectorType::Zero(diag_size-desired_rank);

  HouseholderQR<MatrixAType> qra(a);
  HouseholderQR<MatrixBType> qrb(b);
  m = qra.householderQ() * d * qrb.householderQ();
}

} // end namespace Eigen

template<typename T> struct GetDifferentType;

template<> struct GetDifferentType<float> { typedef double type; };
template<> struct GetDifferentType<double> { typedef float type; };
template<typename T> struct GetDifferentType<std::complex<T> >
{ typedef std::complex<typename GetDifferentType<T>::type> type; };

template<typename T> std::string type_name() { return "other"; }
template<> std::string type_name<float>() { return "float"; }
template<> std::string type_name<double>() { return "double"; }
template<> std::string type_name<int>() { return "int"; }
template<> std::string type_name<std::complex<float> >() { return "complex<float>"; }
template<> std::string type_name<std::complex<double> >() { return "complex<double>"; }
template<> std::string type_name<std::complex<int> >() { return "complex<int>"; }

// forward declaration of the main test function
void EIGEN_CAT(test_,EIGEN_TEST_FUNC)();

using namespace Eigen;

void set_repeat_from_string(const char *str)
{
  errno = 0;
  g_repeat = int(strtoul(str, 0, 10));
  if(errno || g_repeat <= 0)
  {
    std::cout << "Invalid repeat value " << str << std::endl;
    exit(EXIT_FAILURE);
  }
  g_has_set_repeat = true;
}

void set_seed_from_string(const char *str)
{
  errno = 0;
  g_seed = strtoul(str, 0, 10);
  if(errno || g_seed == 0)
  {
    std::cout << "Invalid seed value " << str << std::endl;
    exit(EXIT_FAILURE);
  }
  g_has_set_seed = true;
}

int main(int argc, char *argv[])
{
    g_has_set_repeat = false;
    g_has_set_seed = false;
    bool need_help = false;

    for(int i = 1; i < argc; i++)
    {
      if(argv[i][0] == 'r')
      {
        if(g_has_set_repeat)
        {
          std::cout << "Argument " << argv[i] << " conflicting with a former argument" << std::endl;
          return 1;
        }
        set_repeat_from_string(argv[i]+1);
      }
      else if(argv[i][0] == 's')
      {
        if(g_has_set_seed)
        {
          std::cout << "Argument " << argv[i] << " conflicting with a former argument" << std::endl;
          return 1;
        }
         set_seed_from_string(argv[i]+1);
      }
      else
      {
        need_help = true;
      }
    }

    if(need_help)
    {
      std::cout << "This test application takes the following optional arguments:" << std::endl;
      std::cout << "  rN     Repeat each test N times (default: " << DEFAULT_REPEAT << ")" << std::endl;
      std::cout << "  sN     Use N as seed for random numbers (default: based on current time)" << std::endl;
      std::cout << std::endl;
      std::cout << "If defined, the environment variables EIGEN_REPEAT and EIGEN_SEED" << std::endl;
      std::cout << "will be used as default values for these parameters." << std::endl;
      return 1;
    }

    char *env_EIGEN_REPEAT = getenv("EIGEN_REPEAT");
    if(!g_has_set_repeat && env_EIGEN_REPEAT)
      set_repeat_from_string(env_EIGEN_REPEAT);
    char *env_EIGEN_SEED = getenv("EIGEN_SEED");
    if(!g_has_set_seed && env_EIGEN_SEED)
      set_seed_from_string(env_EIGEN_SEED);

    if(!g_has_set_seed) g_seed = (unsigned int) time(NULL);
    if(!g_has_set_repeat) g_repeat = DEFAULT_REPEAT;

    std::cout << "Initializing random number generator with seed " << g_seed << std::endl;
    srand(g_seed);
    std::cout << "Repeating each test " << g_repeat << " times" << std::endl;

    Eigen::g_test_stack.push_back(EI_PP_MAKE_STRING(EIGEN_TEST_FUNC));

    EIGEN_CAT(test_,EIGEN_TEST_FUNC)();
    return 0;
}
