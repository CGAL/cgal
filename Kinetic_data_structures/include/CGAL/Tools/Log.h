// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_LOG_H_
#define CGAL_LOG_H_
#include <CGAL/basic.h>
#include <CGAL/Tools/utility_macros.h>
#include <iostream>
#include <fstream>
#include <ios>

#if defined(BOOST_MSVC)
// Disable the warning about dll-interface needed for std::ofstream members
// of Log::State.
#  pragma warning(push)
#  pragma warning(disable:4251)
#endif

namespace CGAL {

class Log
{
public: 
  enum Level {NONE=0, SOME=2, LOTS=3};
  
  enum Target {COUT, FILE, DEVNULL};
private:
  struct CGAL_EXPORT State {
    Target target_;
    Level level_;
    std::ofstream fstream_;
    std::ofstream null_;
    std::ofstream maple_;
    bool maple_is_open_;
    bool output_maple_;
    State(){
      level_= NONE;
      target_= COUT;
      null_.open("/dev/null");
      //maple_.open("maple.log");
      maple_is_open_=false;
      output_maple_=true;
    }
  };
  CGAL_EXPORT static State state_;
public:
  // The different types of logs supported
  /*  MAPLE is a log which should be able to be fed directly in to
      maple and preferably will produce obviously good or bad outwhen
      when evaluated.
  */

  static Level level() {return state_.level_;}
  static void set_level(Level l) {state_.level_=l;}
 

  static std::ostream &stream(Level l) {
    if (is_output(l)) {
      if (state_.target_== COUT) {
	return std::cout;
      }
      else {
	return state_.fstream_;
      }
    }
    else {
      return state_.null_;
    }
  }

  static bool is_output(Level l) {
    return l <= level();
  }
  static Target target() {return state_.target_;}
  static CGAL_SET(Target, target, state_.target_=k);
  static CGAL_SET(std::string, filename, state_.fstream_.open(k.c_str()));

  static bool is_output_maple(){return state_.output_maple_;}
  
  static void set_output_maple(bool b) {
    state_.output_maple_=b;
  }
  std::ofstream &maple_stream() {
    if (!state_.maple_is_open_) {
      state_.maple_is_open_=true;
      state_.maple_.open("maple.log");
    }
    return state_.maple_;
  }
private:
  Log() {
   
  }
};



#ifndef CGAL_DISABLE_LOGGING
#define CGAL_LOG(level, expr) if (CGAL::Log::is_output(level))\
    { CGAL::Log::stream(level) << expr << std::flush;};
#define CGAL_LOG_WRITE(level, expr) if (CGAL::Log::is_output(level))\
{std::ostream &LOG_STREAM= CGAL::Log::stream(level); expr;}
#define CGAL_ERROR(expr) std::cerr << expr << std::endl;
#define CGAL_ERROR_WRITE(expr) {std::ostream &LOG_STREAM= std::cerr; expr; std::cerr << std::endl;}
#define CGAL_SET_LOG_LEVEL(level) CGAL::Log::set_level(level);
#define CGAL_GET_LOG_LEVEL CGAL::Log::level();

template <class T>
inline int CGAL_assertion_strip_unsigned(const T&t) {
  return static_cast<int>(t);
}


/*inline int CGAL_assertion_strip_unsigned(size_t t) {
  return static_cast<int>(t);
  }*/

#define CGAL_assert_equal(a,b) do {if (a != b) { CGAL_ERROR("" #a " = " << a); CGAL_ERROR("" #b " = " << b); CGAL_assertion(a ==b);} } while (0)
#define CGAL_check_bounds(a,b,e) do {if (CGAL::CGAL_assertion_strip_unsigned(a) < CGAL::CGAL_assertion_strip_unsigned(b) || CGAL::CGAL_assertion_strip_unsigned(a) >=CGAL::CGAL_assertion_strip_unsigned(e)){ CGAL_ERROR("" #a " = " << a); CGAL_ERROR("[" #b "..." #e ") = [" << b << "..." << e << ")"); CGAL_error();} } while (0)

#else
#define CGAL_LOG(l,e)
#define CGAL_LOG_WRITE(l,e)
#define CGAL_ERROR(e)
#define CGAL_ERROR_WRITE(e)
#define CGAL_SET_LOG_LEVEL(l)
#define CGAL_assert_equal(a,b) 
#define CGAL_check_bounds(a,b,c)
#endif

struct Set_log_state{
  Set_log_state(Log::Level l) {
    old_level_= CGAL_GET_LOG_LEVEL;
    CGAL_SET_LOG_LEVEL(l);
  }
  ~Set_log_state() {
    CGAL_SET_LOG_LEVEL(old_level_);
  }

  Log::Level old_level_;
};

} //namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif


