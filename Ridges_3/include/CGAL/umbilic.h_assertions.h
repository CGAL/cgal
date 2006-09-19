// Replace this remark by a proper copyright notice.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Geert-Jan Giezeman, Sven Schönherr

// Generated from script create_assertions.sh

// macro definitions
// =================
// assertions
// ----------


#if defined(CGAL_UMBILIC.H_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || defined(NDEBUG)
#  define CGAL_umbilic.h_assertion(EX) (static_cast<void>(0))
#  define CGAL_umbilic.h_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_umbilic.h_assertion_code(CODE)
#else
#  define CGAL_umbilic.h_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_umbilic.h_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_umbilic.h_assertion_code(CODE) CODE
#  define CGAL_umbilic.h_assertions 1
#endif // CGAL_UMBILIC.H_NO_ASSERTIONS

#if defined(CGAL_UMBILIC.H_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_UMBILIC.H_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_umbilic.h_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_umbilic.h_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_umbilic.h_exactness_assertion_code(CODE)
#else
#  define CGAL_umbilic.h_exactness_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_umbilic.h_exactness_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_umbilic.h_exactness_assertion_code(CODE) CODE
#  define CGAL_umbilic.h_exactness_assertions 1
#endif // CGAL_UMBILIC.H_NO_ASSERTIONS

#if defined(CGAL_UMBILIC.H_NO_ASSERTIONS) \
  || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_UMBILIC.H_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_umbilic.h_expensive_assertion(EX) (static_cast<void>(0))
#  define CGAL_umbilic.h_expensive_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_umbilic.h_expensive_assertion_code(CODE)
#else
#  define CGAL_umbilic.h_expensive_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_umbilic.h_expensive_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_umbilic.h_expensive_assertion_code(CODE) CODE
#  define CGAL_umbilic.h_expensive_assertions 1
#endif // CGAL_UMBILIC.H_NO_ASSERTIONS

#if defined(CGAL_UMBILIC.H_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_UMBILIC.H_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_UMBILIC.H_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_umbilic.h_expensive_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_umbilic.h_expensive_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_umbilic.h_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_umbilic.h_expensive_exactness_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_umbilic.h_expensive_exactness_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_umbilic.h_expensive_exactness_assertion_code(CODE) CODE
#  define CGAL_umbilic.h_expensive_exactness_assertions 1
#endif // CGAL_UMBILIC.H_NO_ASSERTIONS


// preconditions
// -------------

#if defined(CGAL_UMBILIC.H_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_umbilic.h_precondition(EX) (static_cast<void>(0))
#  define CGAL_umbilic.h_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_umbilic.h_precondition_code(CODE)
#else
#  define CGAL_umbilic.h_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_umbilic.h_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_umbilic.h_precondition_code(CODE) CODE
#  define CGAL_umbilic.h_preconditions 1
#endif // CGAL_UMBILIC.H_NO_PRECONDITIONS

#if defined(CGAL_UMBILIC.H_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_UMBILIC.H_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_umbilic.h_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_umbilic.h_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_umbilic.h_exactness_precondition_code(CODE)
#else
#  define CGAL_umbilic.h_exactness_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_umbilic.h_exactness_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_umbilic.h_exactness_precondition_code(CODE) CODE
#  define CGAL_umbilic.h_exactness_preconditions 1
#endif // CGAL_UMBILIC.H_NO_PRECONDITIONS

#if defined(CGAL_UMBILIC.H_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_UMBILIC.H_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_umbilic.h_expensive_precondition(EX) (static_cast<void>(0))
#  define CGAL_umbilic.h_expensive_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_umbilic.h_expensive_precondition_code(CODE)
#else
#  define CGAL_umbilic.h_expensive_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_umbilic.h_expensive_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_umbilic.h_expensive_precondition_code(CODE) CODE
#  define CGAL_umbilic.h_expensive_preconditions 1
#endif // CGAL_UMBILIC.H_NO_PRECONDITIONS

#if defined(CGAL_UMBILIC.H_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_UMBILIC.H_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_UMBILIC.H_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_umbilic.h_expensive_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_umbilic.h_expensive_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_umbilic.h_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_umbilic.h_expensive_exactness_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_umbilic.h_expensive_exactness_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_umbilic.h_expensive_exactness_precondition_code(CODE) CODE
#  define CGAL_umbilic.h_expensive_exactness_preconditions 1
#endif // CGAL_UMBILIC.H_NO_PRECONDITIONS


// postconditions
// --------------

#if defined(CGAL_UMBILIC.H_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_umbilic.h_postcondition(EX) (static_cast<void>(0))
#  define CGAL_umbilic.h_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_umbilic.h_postcondition_code(CODE)
#else
#  define CGAL_umbilic.h_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_umbilic.h_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_umbilic.h_postcondition_code(CODE) CODE
#  define CGAL_umbilic.h_postconditions 1
#endif // CGAL_UMBILIC.H_NO_POSTCONDITIONS

#if defined(CGAL_UMBILIC.H_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_UMBILIC.H_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_umbilic.h_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_umbilic.h_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_umbilic.h_exactness_postcondition_code(CODE)
#else
#  define CGAL_umbilic.h_exactness_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_umbilic.h_exactness_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_umbilic.h_exactness_postcondition_code(CODE) CODE
#  define CGAL_umbilic.h_exactness_postconditions 1
#endif // CGAL_UMBILIC.H_NO_POSTCONDITIONS

#if defined(CGAL_UMBILIC.H_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_UMBILIC.H_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_umbilic.h_expensive_postcondition(EX) (static_cast<void>(0))
#  define CGAL_umbilic.h_expensive_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_umbilic.h_expensive_postcondition_code(CODE)
#else
#  define CGAL_umbilic.h_expensive_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_umbilic.h_expensive_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_umbilic.h_expensive_postcondition_code(CODE) CODE
#  define CGAL_umbilic.h_expensive_postconditions 1
#endif // CGAL_UMBILIC.H_NO_POSTCONDITIONS

#if defined(CGAL_UMBILIC.H_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_UMBILIC.H_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_UMBILIC.H_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_umbilic.h_expensive_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_umbilic.h_expensive_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_umbilic.h_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_umbilic.h_expensive_exactness_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_umbilic.h_expensive_exactness_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_umbilic.h_expensive_exactness_postcondition_code(CODE) CODE
#  define CGAL_umbilic.h_expensive_exactness_postconditions 1
#endif // CGAL_UMBILIC.H_NO_POSTCONDITIONS


// warnings
// --------

#if defined(CGAL_UMBILIC.H_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || defined(NDEBUG)
#  define CGAL_umbilic.h_warning(EX) (static_cast<void>(0))
#  define CGAL_umbilic.h_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_umbilic.h_warning_code(CODE)
#else
#  define CGAL_umbilic.h_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_umbilic.h_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_umbilic.h_warning_code(CODE) CODE
#  define CGAL_umbilic.h_warnings 1
#endif // CGAL_UMBILIC.H_NO_WARNINGS

#if defined(CGAL_UMBILIC.H_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_UMBILIC.H_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_umbilic.h_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_umbilic.h_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_umbilic.h_exactness_warning_code(CODE)
#else
#  define CGAL_umbilic.h_exactness_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_umbilic.h_exactness_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_umbilic.h_exactness_warning_code(CODE) CODE
#  define CGAL_umbilic.h_exactness_warnings 1
#endif // CGAL_UMBILIC.H_NO_WARNINGS

#if defined(CGAL_UMBILIC.H_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_UMBILIC.H_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_umbilic.h_expensive_warning(EX) (static_cast<void>(0))
#  define CGAL_umbilic.h_expensive_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_umbilic.h_expensive_warning_code(CODE)
#else
#  define CGAL_umbilic.h_expensive_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_umbilic.h_expensive_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_umbilic.h_expensive_warning_code(CODE) CODE
#  define CGAL_umbilic.h_expensive_warnings 1
#endif // CGAL_UMBILIC.H_NO_WARNINGS

#if defined(CGAL_UMBILIC.H_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_UMBILIC.H_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_UMBILIC.H_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_umbilic.h_expensive_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_umbilic.h_expensive_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_umbilic.h_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_umbilic.h_expensive_exactness_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_umbilic.h_expensive_exactness_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_umbilic.h_expensive_exactness_warning_code(CODE) CODE
#  define CGAL_umbilic.h_expensive_exactness_warnings 1
#endif // CGAL_UMBILIC.H_NO_WARNINGS


