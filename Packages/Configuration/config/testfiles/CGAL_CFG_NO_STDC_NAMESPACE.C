// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.0-I-3 $
// release_date  : $CGAL_Date: 1999/03/08 $
//
// file          : config/testfiles/CGAL_CFG_NO_STDC_NAMESPACE.C
// package       : Configuration (1.26)
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_STDC_NAMESPACE.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| The flag CGAL_CFG_NO_STDC_NAMESPACE is set, if a compiler does not
//| put the parts of the standard library inherited from the standard
//| C library in namespace std. (only tests for the symbols used in CGAL)

// (what about this one?) #include <ciso646>

#include <csetjmp>
using std::longjmp;
using std::jmp_buf;

#include <cstdarg>
using std::va_list;

#include <ctime>
using std::time_t;   
using std::clock_t;  
using std::tm;       
using std::asctime;  
using std::ctime;    
using std::clock;    
using std::difftime; 
using std::gmtime;   
using std::localtime;
using std::mktime;   
using std::time;     
using std::strftime; 

#include <csignal>
using std::sig_atomic_t;
using std::raise;
using std::signal;

#include <cstdlib>
using std::abort;
using std::atexit;
using std::exit;
using std::getenv;
using std::system;
using std::calloc;
using std::malloc;
using std::free;
using std::realloc;
using std::atof;
using std::atoi;
using std::atol;
using std::mblen;
using std::mbstowcs;
using std::mbtowc;
using std::strtod;
using std::strtol;
using std::strtoul;
using std::wcstombs;
using std::wctomb;
using std::bsearch;
using std::qsort;
using std::div_t;
using std::ldiv_t;
using std::abs;
using std::div;
using std::labs;
using std::ldiv;
using std::rand;
using std::srand;

#include <cctype>
using std::isalnum;
using std::isalpha;
using std::iscntrl;
using std::isdigit;
using std::isgraph;
using std::islower;
using std::isprint;
using std::ispunct;
using std::isspace;
using std::isupper;
using std::isxdigit;

#if 0
#include <cwctype>
using std::wctrans_t;
using std::wctype_t;
using std::wint_t;
using std::iswalnum;
using std::iswalpha;
using std::iswcntrl;
using std::iswctype;
using std::iswdigit;
using std::iswgraph;
using std::iswlower;
using std::iswprint;
using std::iswpunct;
using std::iswspace;
using std::iswupper;
using std::iswxdigit;
using std::towctrans;
using std::towlower;
using std::towupper;
using std::wctrans;
using std::wctype;
#endif

#include <cstring>
using std::memchr;
using std::memcmp;
using std::memcpy;
using std::memmove;
using std::memset;
using std::strcat;
using std::strchr;
using std::strcmp;
using std::strcoll;
using std::strcpy;
using std::strcspn;
using std::strerror;
// using std::strlen;
using std::strncat;
using std::strncmp;
using std::strncpy;
using std::strpbrk;
using std::strrchr;
using std::strspn;
using std::strstr;
using std::strtok;
using std::strxfrm;
#if 0
#include <cwchar>
using std::wint_t;
/* These do not exist in egcs-1.1.1
using std::mbstate_t;
using std::btowc;
using std::fwide;
using std::fwprintf;
using std::fwscanf;
using std::mbrlen;
using std::mbrtowc;
using std::mbsinit;
using std::mbsrtowcs;
using std::swprintf;
using std::swscanf;
using std::vfwprintf;
using std::vswprintf;
using std::vwprintf;
using std::wcrtomb;
using std::wcsrtombs;
using std::wcsstr;
using std::wctob;
using std::wmemchr;
using std::wmemcmp;
using std::wmemcpy;
using std::wmemmove;
using std::wmemset;
using std::wprintf;
using std::wscanf;
*/
using std::fgetwc;
using std::fgetws;
using std::fputwc;
using std::fputws;
using std::getwc;
using std::getwchar;
using std::putwc;
using std::putwchar;
using std::ungetwc;
using std::wcscat;
using std::wcschr;
using std::wcscmp;
using std::wcscoll;
using std::wcscpy;
using std::wcscspn;
using std::wcsftime;
using std::wcslen;
using std::wcsncat;
using std::wcsncmp;
using std::wcsncpy;
using std::wcspbrk;
using std::wcsrchr;
using std::wcsspn;
using std::wcstod;
using std::wcstok;
using std::wcstol;
using std::wcstoul;
using std::wcsxfrm;
#endif

#include <cstdio>
using std::FILE;
using std::fpos_t;
using std::clearerr;
using std::fclose;
using std::feof;
using std::ferror;
using std::fflush;
using std::fgetc;
using std::fgetpos;
using std::fgets;
using std::fopen;
using std::fprintf;
using std::fputc;
using std::fputs;
using std::fread;
using std::freopen;
using std::fscanf;
using std::fseek;
using std::fsetpos;
using std::ftell;
using std::fwrite;
using std::getc;
using std::getchar;
using std::gets;
using std::perror;
using std::printf;
using std::putc;
using std::putchar;
using std::puts;
using std::remove;
using std::rename;
using std::rewind;
using std::scanf;
using std::setbuf;
using std::setvbuf;
using std::sprintf;
using std::sscanf;
using std::tmpfile;
using std::tmpnam;
using std::ungetc;
using std::vfprintf;
using std::vprintf;
using std::vsprintf;

#include <climits>
/* Values defined as macros
using std::CHAR_BIT;
using std::CHAR_MIN;
using std::CHAR_MAX;
using std::INT_MIN;
using std::INT_MAX;
using std::LONG_MIN;
using std::LONG_MAX;
using std::MB_LEN_MAX;
using std::SCHAR_MIN;
using std::SCHAR_MAX;
using std::SHRT_MIN;
using std::SHRT_MAX;
using std::UCHAR_MAX;
using std::UINT_MAX;
using std::ULONG_MAX;
using std::USHRT_MAX;
*/

// (what about this one?) #include <clocale>
/* These do not exist in egcs-1.1.1
using std::lconv;
using std::localeconv;
using std::setlocale;
*/

#include <cfloat> 
/* Values defined as macros
using std::DBL_DIG;
using std::DBL_EPSILON;
using std::DBL_MANT_DIG;
using std::DBL_MAX;
using std::DBL_MAX_10_EXP;
using std::DBL_MAX_EXP;
using std::DBL_MIN;
using std::DBL_MIN_10_EXP;
using std::DBL_MIN_EXP;
using std::FLT_DIG;
using std::FLT_EPSILON;
using std::FLT_MANT_DIG;
using std::FLT_MAX;
using std::FLT_MAX_10_EXP;
using std::FLT_MAX_EXP;
using std::FLT_MIN;
using std::FLT_MIN_10_EXP;
using std::FLT_MIN_EXP;
using std::FLT_RADIX;
using std::FLT_ROUNDS;
using std::LDBL_DIG;
using std::LDBL_EPSILON;
using std::LDBL_MANT_DIG;
using std::LDBL_MAX;
using std::LDBL_MAX_10_EXP;
using std::LDBL_MAX_EXP;
using std::LDBL_MIN;
using std::LDBL_MIN_10_EXP;
using std::LDBL_MIN_EXP;
*/

#include <cmath>
using std::acos;
using std::asin;
using std::atan2;
using std::atan;
using std::ceil;
using std::cos;
using std::cosh;
using std::exp;
using std::fabs;
using std::floor;
using std::fmod;
using std::frexp;
using std::ldexp;
using std::log10;
using std::log;
using std::modf;
using std::pow;
using std::sin;
using std::sinh;
using std::sqrt;
using std::tan;
using std::tanh;

#include <cstddef>
using std::ptrdiff_t;
using std::size_t;

int main()
{
  std::strlen("a");
  return 0;
}

// EOF //

