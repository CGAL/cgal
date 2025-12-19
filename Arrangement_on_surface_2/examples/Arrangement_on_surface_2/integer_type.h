#ifndef INTEGER_TYPE_H
#define INTEGER_TYPE_H

#include <CGAL/basic.h>

#if CGAL_USE_GMP && CGAL_USE_MPFI
#include <CGAL/Gmpz.h>
using Integer = CGAL::Gmpz;
#elif CGAL_USE_CORE
#include <CGAL/CORE_BigInt.h>
using Integer = CORE::BigInt;
#else
#include <CGAL/leda_integer.h>
using Integer = LEDA::integer;
#endif

#endif
