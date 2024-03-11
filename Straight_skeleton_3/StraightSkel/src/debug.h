/**
 * @file   debug.h
 * @author Gernot Walzl
 * @date   2011-12-22
 *
 * based on http://www.digitalpeer.com/id/debug
 */

#ifndef DEBUG_H
#define DEBUG_H

#include <cstdio>
#include <cassert>
#include <iostream>
#include "config.h"
#include "util/StackTrace.h"


#ifdef DEBUG

#define DBGOUT std::cout

/* just a helper for code location */
#define LOC DBGOUT << "[DEBUG] " << __FILE__ << ":" << __LINE__ << ": ";

/* macro using var args */
#define DEBUG_PRINT(...) LOC fprintf(stdout, __VA_ARGS__); DBGOUT << std::endl;

/* macro for general debug print statements. */
#define DEBUG_VAL(var) LOC DBGOUT << var << std::endl;

/* macro that dumps a variable name and its actual value */
#define DEBUG_VAR(var) LOC DBGOUT << (#var) << "=" << var << std::endl;

/* macro that warns, if a smart pointer is invalid */
#define DEBUG_SPTR(sptr) if (!sptr) { LOC DBGOUT << "shared pointer is invalid: " << (#sptr) << std::endl; util::StackTrace::print(DBGOUT); }
#define DEBUG_WPTR(wptr) if (wptr.expired()) { LOC DBGOUT << "weak pointer is expired: " << (#wptr) << std::endl; }

#else

/* when debug isn't defined all the macro calls do absolutely nothing */
#define DEBUG_PRINT(...)
#define DEBUG_VAL(var)
#define DEBUG_VAR(var)
#define DEBUG_SPTR(sptr)
#define DEBUG_WPTR(wptr)

#endif

#endif /* DEBUG_H */

