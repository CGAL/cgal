#ifndef NUMBER_TYPE_H
#define NUMBER_TYPE_H

#include <CGAL/config.h>

#include "short_names.h"

#include <CGAL/basic.h>

#if defined(USE_CONIC_TRAITS)
#include <CGAL/leda_real.h>
#else

#if defined(USE_LEDA_KERNEL)
#include <CGAL/leda_rational.h>
#else

#if defined(USE_MY_KERNEL)
#include <CGAL/leda_rational.h>
#else

#if defined(USE_LAZY_RAT)
#include <CGAL/leda_rational.h>
#include <CGAL/Lazy_exact_nt.h>
#else

#if defined(USE_LAZY_QUOTIENT)
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Lazy_exact_nt.h>
#else

#include <CGAL/leda_rational.h>
#endif
#endif
#endif
#endif
#endif

#if defined(USE_CONIC_TRAITS)
typedef leda_real                                       NT;
#define NUMBER_TYPE "Leda Real"
#else

#if defined(USE_LEDA_KERNEL)
typedef leda_rational                                   NT;
#define NUMBER_TYPE "Leda Rat"
#else

#if defined(USE_MY_KERNEL)
typedef leda_rational                                   NT;
#define NUMBER_TYPE "Leda Rat"
#else

#if defined(USE_LAZY_RAT)
typedef leda_rational                                   NT;
typedef CGAL::Lazy_exact_nt<NT>                         WNT;
#define NUMBER_TYPE "Lazy Leda Rat"
#else

#if defined(USE_LAZY_QUOTIENT)
typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Lazy_exact_nt<NT>                         WNT;
#define NUMBER_TYPE "Lazy Quotient Float"
#else

typedef leda_rational                                   NT;
typedef NT                                              WNT;
#define NUMBER_TYPE "Leda Rat"
#endif
#endif
#endif
#endif
#endif

#endif
