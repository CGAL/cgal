//  boost timer.hpp header file  ---------------------------------------------//

//  (C) Copyright Beman Dawes 1994-99. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

//  See http://www.boost.org/libs/timer for documentation.

//  Revision History
//  01 Apr 01  Modified to use new <boost/limits.hpp> header. (JMaddock)
//  12 Jan 01  Change to inline implementation to allow use without library
//             builds. See docs for more rationale. (Beman Dawes) 
//  25 Sep 99  elapsed_max() and elapsed_min() added (John Maddock)
//  16 Jul 99  Second beta
//   6 Jul 99  Initial boost version

#ifndef BOOST_TIMER_HPP
#define BOOST_TIMER_HPP

#include <sys/times.h>
#include <limits.h>

#ifndef CLK_TCK
#  define CLK_TCK 100
#endif


//  timer  -------------------------------------------------------------------//

//  A timer object measures elapsed time.

//  It is recommended that implementations measure wall clock rather than CPU
//  time since the intended use is performance measurement on systems where
//  total elapsed time is more important than just process or CPU time.

//  Warnings: The maximum measurable elapsed time may well be only 596.5+ hours
//  due to implementation limitations.  The accuracy of timings depends on the
//  accuracy of timing information provided by the underlying platform, and
//  this varies a great deal from platform to platform.

class timer
{
 public:
         timer() { _start_time = times (&_tms); } // postcondition: elapsed()==0
//         timer( const timer& src );      // post: elapsed()==src.elapsed()
//        ~timer(){}
//  timer& operator=( const timer& src );  // post: elapsed()==src.elapsed()
  void   restart() { _start_time = times (&_tms); } // post: elapsed()==0
  double elapsed(struct tms *now) const                  // return elapsed time in seconds
    { 
      clock_t now_clock = times (now);
      now->tms_utime -= _tms.tms_utime;
      now->tms_stime -= _tms.tms_stime;
      now->tms_cutime -= _tms.tms_cutime;
      now->tms_cstime -= _tms.tms_cstime;
      return  double (now_clock - _start_time) / double (CLK_TCK);
    }
  static double convert (clock_t clk) {
    return double (clk) / double (CLK_TCK);
  }
 private:
  std::clock_t _start_time;
  struct tms _tms;
}; // timer


#endif  // BOOST_TIMER_HPP
