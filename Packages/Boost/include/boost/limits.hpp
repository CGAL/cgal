
//  (C) Copyright Boost.org 1999. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
//
// use this header as a workaround for missing <limits>

//  See http://www.boost.org/libs/utility/limits.html for documentation.

#ifndef BOOST_LIMITS
#define BOOST_LIMITS

#include <boost/config.hpp>

#ifdef BOOST_NO_LIMITS
# include <boost/detail/limits.hpp>
#else
# include <limits>
#endif

#if (defined(BOOST_HAS_LONG_LONG) && defined(BOOST_NO_LONG_LONG_NUMERIC_LIMITS)) \
      || (defined(BOOST_HAS_MS_INT64) && defined(BOOST_NO_MS_INT64_NUMERIC_LIMITS))
// Add missing specializations for numeric_limits:
#ifdef BOOST_HAS_MS_INT64
#  define BOOST_LLT __int64
#else
#  define BOOST_LLT long long
#endif

namespace std
{
  template<>
  class numeric_limits<BOOST_LLT> 
  {
   public:

      BOOST_STATIC_CONSTANT(bool, is_specialized = true);
#ifdef BOOST_HAS_MS_INT64
      static BOOST_LLT min(){ return 0x8000000000000000i64; }
      static BOOST_LLT max(){ return 0x7FFFFFFFFFFFFFFFi64; }
#elif defined(LLONG_MAX)
      static BOOST_LLT min(){ return LLONG_MIN; }
      static BOOST_LLT max(){ return LLONG_MAX; }
#elif defined(LONGLONG_MAX)
      static BOOST_LLT min(){ return LONGLONG_MIN; }
      static BOOST_LLT max(){ return LONGLONG_MAX; }
#else
      static BOOST_LLT min(){ return 1LL << (sizeof(BOOST_LLT) * CHAR_BIT - 1); }
      static BOOST_LLT max(){ return ~min(); }
#endif
      BOOST_STATIC_CONSTANT(int, digits = sizeof(BOOST_LLT) * CHAR_BIT -1);
      BOOST_STATIC_CONSTANT(int, digits10 = (CHAR_BIT * sizeof (BOOST_LLT) - 1) * 301L / 1000);
      BOOST_STATIC_CONSTANT(bool, is_signed = true);
      BOOST_STATIC_CONSTANT(bool, is_integer = true);
      BOOST_STATIC_CONSTANT(bool, is_exact = true);
      BOOST_STATIC_CONSTANT(int, radix = 2);
      static BOOST_LLT epsilon() throw() { return 0; };
      static BOOST_LLT round_error() throw() { return 0; };

      BOOST_STATIC_CONSTANT(int, min_exponent = 0);
      BOOST_STATIC_CONSTANT(int, min_exponent10 = 0);
      BOOST_STATIC_CONSTANT(int, max_exponent = 0);
      BOOST_STATIC_CONSTANT(int, max_exponent10 = 0);

      BOOST_STATIC_CONSTANT(bool, has_infinity = false);
      BOOST_STATIC_CONSTANT(bool, has_quiet_NaN = false);
      BOOST_STATIC_CONSTANT(bool, has_signaling_NaN = false);
      BOOST_STATIC_CONSTANT(bool, has_denorm = false);
      BOOST_STATIC_CONSTANT(bool, has_denorm_loss = false);
      static BOOST_LLT infinity() throw() { return 0; };
      static BOOST_LLT quiet_NaN() throw() { return 0; };
      static BOOST_LLT signaling_NaN() throw() { return 0; };
      static BOOST_LLT denorm_min() throw() { return 0; };

      BOOST_STATIC_CONSTANT(bool, is_iec559 = false);
      BOOST_STATIC_CONSTANT(bool, is_bounded = false);
      BOOST_STATIC_CONSTANT(bool, is_modulo = false);

      BOOST_STATIC_CONSTANT(bool, traps = false);
      BOOST_STATIC_CONSTANT(bool, tinyness_before = false);
      BOOST_STATIC_CONSTANT(float_round_style, round_style = round_toward_zero);
      
  };

  template<>
  class numeric_limits<unsigned BOOST_LLT> 
  {
   public:

      BOOST_STATIC_CONSTANT(bool, is_specialized = true);
#ifdef BOOST_HAS_MS_INT64
      static unsigned BOOST_LLT min(){ return 0ui64; }
      static unsigned BOOST_LLT max(){ return 0xFFFFFFFFFFFFFFFFui64; }
#elif defined(ULLONG_MAX) && defined(ULLONG_MIN)
      static unsigned BOOST_LLT min(){ return ULLONG_MIN; }
      static unsigned BOOST_LLT max(){ return ULLONG_MAX; }
#elif defined(ULONGLONG_MAX) && defined(ULONGLONG_MIN)
      static unsigned BOOST_LLT min(){ return ULONGLONG_MIN; }
      static unsigned BOOST_LLT max(){ return ULONGLONG_MAX; }
#else
      static unsigned BOOST_LLT min(){ return 0uLL; }
      static unsigned BOOST_LLT max(){ return ~0uLL; }
#endif
      BOOST_STATIC_CONSTANT(int, digits = sizeof(BOOST_LLT) * CHAR_BIT);
      BOOST_STATIC_CONSTANT(int, digits10 = (CHAR_BIT * sizeof (BOOST_LLT)) * 301L / 1000);
      BOOST_STATIC_CONSTANT(bool, is_signed = false);
      BOOST_STATIC_CONSTANT(bool, is_integer = true);
      BOOST_STATIC_CONSTANT(bool, is_exact = true);
      BOOST_STATIC_CONSTANT(int, radix = 2);
      static unsigned BOOST_LLT epsilon() throw() { return 0; };
      static unsigned BOOST_LLT round_error() throw() { return 0; };

      BOOST_STATIC_CONSTANT(int, min_exponent = 0);
      BOOST_STATIC_CONSTANT(int, min_exponent10 = 0);
      BOOST_STATIC_CONSTANT(int, max_exponent = 0);
      BOOST_STATIC_CONSTANT(int, max_exponent10 = 0);

      BOOST_STATIC_CONSTANT(bool, has_infinity = false);
      BOOST_STATIC_CONSTANT(bool, has_quiet_NaN = false);
      BOOST_STATIC_CONSTANT(bool, has_signaling_NaN = false);
      BOOST_STATIC_CONSTANT(bool, has_denorm = false);
      BOOST_STATIC_CONSTANT(bool, has_denorm_loss = false);
      static unsigned BOOST_LLT infinity() throw() { return 0; };
      static unsigned BOOST_LLT quiet_NaN() throw() { return 0; };
      static unsigned BOOST_LLT signaling_NaN() throw() { return 0; };
      static unsigned BOOST_LLT denorm_min() throw() { return 0; };

      BOOST_STATIC_CONSTANT(bool, is_iec559 = false);
      BOOST_STATIC_CONSTANT(bool, is_bounded = false);
      BOOST_STATIC_CONSTANT(bool, is_modulo = false);

      BOOST_STATIC_CONSTANT(bool, traps = false);
      BOOST_STATIC_CONSTANT(bool, tinyness_before = false);
      BOOST_STATIC_CONSTANT(float_round_style, round_style = round_toward_zero);
      
  };
}
#endif 

#endif
