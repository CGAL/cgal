#ifndef NUMBER_TYPE_H
#define NUMBER_TYPE_H

#include <CGAL/config.h>

#include "short_names.h"

#include <CGAL/basic.h>

#if defined(USE_CONIC_TRAITS)
#include <CGAL/leda_real.h>
#elif defined(USE_LEDA_KERNEL)
#include <CGAL/leda_rational.h>
#elif defined(USE_MY_KERNEL)
#include <CGAL/leda_rational.h>
#elif defined(USE_LAZY_RAT)
#include <CGAL/leda_rational.h>
#include <CGAL/Lazy_exact_nt.h>
#elif defined(USE_GMPQ)
#include <CGAL/Gmpq.h>
#elif defined(USE_LAZY_GMPQ)
#include <CGAL/Gmpq.h>
#include <CGAL/Lazy_exact_nt.h>
#elif defined(USE_MP_FLOAT)
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#elif defined(USE_LAZY_QUOTIENT)
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Lazy_exact_nt.h>
#else
#include <CGAL/leda_rational.h>
#endif

#if defined(USE_CONIC_TRAITS)
typedef leda_real                                       NT;
typedef NT                                              WNT;
#define NUMBER_TYPE "Leda Real"
#elif defined(USE_LEDA_KERNEL)
typedef leda_rational                                   NT;
#define NUMBER_TYPE "Leda Rat"
#elif defined(USE_MY_KERNEL)
typedef leda_rational                                   NT;
#define NUMBER_TYPE "Leda Rat"
#elif defined(USE_LAZY_RAT)
typedef leda_rational                                   NT;
typedef CGAL::Lazy_exact_nt<NT>                         WNT;
#define NUMBER_TYPE "Lazy Leda Rat"
#elif defined(USE_GMPQ)
typedef CGAL::Gmpq                                      NT;
typedef NT                                              WNT;
#define NUMBER_TYPE "Gmpq"
#elif defined(USE_LAZY_GMPQ)
typedef CGAL::Gmpq                                      NT;
typedef CGAL::Lazy_exact_nt<NT>                         WNT;
#define NUMBER_TYPE "Lazy Gmpq"
#elif defined(USE_MP_FLOAT)
typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef NT                                              WNT;
#define NUMBER_TYPE "Quotient MP Float"
#elif defined(USE_LAZY_QUOTIENT)
typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Lazy_exact_nt<NT>                         WNT;
#define NUMBER_TYPE "Lazy Quotient MP Float"
#else
typedef leda_rational                                   NT;
typedef NT                                              WNT;
#define NUMBER_TYPE "Leda Rat"
#endif

#endif
