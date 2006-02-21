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


#if defined(CGAL_JET_FITTING_3_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || defined(NDEBUG)
#  define CGAL_jet_fitting_3_assertion(EX) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_assertion_code(CODE)
#else
#  define CGAL_jet_fitting_3_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_jet_fitting_3_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_jet_fitting_3_assertion_code(CODE) CODE
#  define CGAL_jet_fitting_3_assertions 1
#endif // CGAL_JET_FITTING_3_NO_ASSERTIONS

#if defined(CGAL_JET_FITTING_3_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_JET_FITTING_3_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_jet_fitting_3_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_exactness_assertion_code(CODE)
#else
#  define CGAL_jet_fitting_3_exactness_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_jet_fitting_3_exactness_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_jet_fitting_3_exactness_assertion_code(CODE) CODE
#  define CGAL_jet_fitting_3_exactness_assertions 1
#endif // CGAL_JET_FITTING_3_NO_ASSERTIONS

#if defined(CGAL_JET_FITTING_3_NO_ASSERTIONS) \
  || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_JET_FITTING_3_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_jet_fitting_3_expensive_assertion(EX) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_expensive_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_expensive_assertion_code(CODE)
#else
#  define CGAL_jet_fitting_3_expensive_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_jet_fitting_3_expensive_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_jet_fitting_3_expensive_assertion_code(CODE) CODE
#  define CGAL_jet_fitting_3_expensive_assertions 1
#endif // CGAL_JET_FITTING_3_NO_ASSERTIONS

#if defined(CGAL_JET_FITTING_3_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_JET_FITTING_3_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_JET_FITTING_3_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_jet_fitting_3_expensive_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_expensive_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_jet_fitting_3_expensive_exactness_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_jet_fitting_3_expensive_exactness_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_jet_fitting_3_expensive_exactness_assertion_code(CODE) CODE
#  define CGAL_jet_fitting_3_expensive_exactness_assertions 1
#endif // CGAL_JET_FITTING_3_NO_ASSERTIONS


// preconditions
// -------------

#if defined(CGAL_JET_FITTING_3_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_jet_fitting_3_precondition(EX) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_precondition_code(CODE)
#else
#  define CGAL_jet_fitting_3_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_jet_fitting_3_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_jet_fitting_3_precondition_code(CODE) CODE
#  define CGAL_jet_fitting_3_preconditions 1
#endif // CGAL_JET_FITTING_3_NO_PRECONDITIONS

#if defined(CGAL_JET_FITTING_3_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_JET_FITTING_3_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_jet_fitting_3_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_exactness_precondition_code(CODE)
#else
#  define CGAL_jet_fitting_3_exactness_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_jet_fitting_3_exactness_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_jet_fitting_3_exactness_precondition_code(CODE) CODE
#  define CGAL_jet_fitting_3_exactness_preconditions 1
#endif // CGAL_JET_FITTING_3_NO_PRECONDITIONS

#if defined(CGAL_JET_FITTING_3_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_JET_FITTING_3_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_jet_fitting_3_expensive_precondition(EX) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_expensive_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_expensive_precondition_code(CODE)
#else
#  define CGAL_jet_fitting_3_expensive_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_jet_fitting_3_expensive_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_jet_fitting_3_expensive_precondition_code(CODE) CODE
#  define CGAL_jet_fitting_3_expensive_preconditions 1
#endif // CGAL_JET_FITTING_3_NO_PRECONDITIONS

#if defined(CGAL_JET_FITTING_3_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_JET_FITTING_3_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_JET_FITTING_3_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_jet_fitting_3_expensive_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_expensive_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_jet_fitting_3_expensive_exactness_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_jet_fitting_3_expensive_exactness_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_jet_fitting_3_expensive_exactness_precondition_code(CODE) CODE
#  define CGAL_jet_fitting_3_expensive_exactness_preconditions 1
#endif // CGAL_JET_FITTING_3_NO_PRECONDITIONS


// postconditions
// --------------

#if defined(CGAL_JET_FITTING_3_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_jet_fitting_3_postcondition(EX) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_postcondition_code(CODE)
#else
#  define CGAL_jet_fitting_3_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_jet_fitting_3_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_jet_fitting_3_postcondition_code(CODE) CODE
#  define CGAL_jet_fitting_3_postconditions 1
#endif // CGAL_JET_FITTING_3_NO_POSTCONDITIONS

#if defined(CGAL_JET_FITTING_3_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_JET_FITTING_3_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_jet_fitting_3_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_exactness_postcondition_code(CODE)
#else
#  define CGAL_jet_fitting_3_exactness_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_jet_fitting_3_exactness_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_jet_fitting_3_exactness_postcondition_code(CODE) CODE
#  define CGAL_jet_fitting_3_exactness_postconditions 1
#endif // CGAL_JET_FITTING_3_NO_POSTCONDITIONS

#if defined(CGAL_JET_FITTING_3_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_JET_FITTING_3_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_jet_fitting_3_expensive_postcondition(EX) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_expensive_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_expensive_postcondition_code(CODE)
#else
#  define CGAL_jet_fitting_3_expensive_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_jet_fitting_3_expensive_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_jet_fitting_3_expensive_postcondition_code(CODE) CODE
#  define CGAL_jet_fitting_3_expensive_postconditions 1
#endif // CGAL_JET_FITTING_3_NO_POSTCONDITIONS

#if defined(CGAL_JET_FITTING_3_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_JET_FITTING_3_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_JET_FITTING_3_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_jet_fitting_3_expensive_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_expensive_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_jet_fitting_3_expensive_exactness_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_jet_fitting_3_expensive_exactness_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_jet_fitting_3_expensive_exactness_postcondition_code(CODE) CODE
#  define CGAL_jet_fitting_3_expensive_exactness_postconditions 1
#endif // CGAL_JET_FITTING_3_NO_POSTCONDITIONS


// warnings
// --------

#if defined(CGAL_JET_FITTING_3_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || defined(NDEBUG)
#  define CGAL_jet_fitting_3_warning(EX) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_warning_code(CODE)
#else
#  define CGAL_jet_fitting_3_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_jet_fitting_3_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_jet_fitting_3_warning_code(CODE) CODE
#  define CGAL_jet_fitting_3_warnings 1
#endif // CGAL_JET_FITTING_3_NO_WARNINGS

#if defined(CGAL_JET_FITTING_3_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_JET_FITTING_3_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_jet_fitting_3_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_exactness_warning_code(CODE)
#else
#  define CGAL_jet_fitting_3_exactness_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_jet_fitting_3_exactness_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_jet_fitting_3_exactness_warning_code(CODE) CODE
#  define CGAL_jet_fitting_3_exactness_warnings 1
#endif // CGAL_JET_FITTING_3_NO_WARNINGS

#if defined(CGAL_JET_FITTING_3_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_JET_FITTING_3_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_jet_fitting_3_expensive_warning(EX) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_expensive_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_expensive_warning_code(CODE)
#else
#  define CGAL_jet_fitting_3_expensive_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_jet_fitting_3_expensive_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_jet_fitting_3_expensive_warning_code(CODE) CODE
#  define CGAL_jet_fitting_3_expensive_warnings 1
#endif // CGAL_JET_FITTING_3_NO_WARNINGS

#if defined(CGAL_JET_FITTING_3_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_JET_FITTING_3_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_JET_FITTING_3_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_jet_fitting_3_expensive_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_expensive_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_jet_fitting_3_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_jet_fitting_3_expensive_exactness_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_jet_fitting_3_expensive_exactness_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_jet_fitting_3_expensive_exactness_warning_code(CODE) CODE
#  define CGAL_jet_fitting_3_expensive_exactness_warnings 1
#endif // CGAL_JET_FITTING_3_NO_WARNINGS


