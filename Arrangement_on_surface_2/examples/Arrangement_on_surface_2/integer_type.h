#ifndef INTEGER_TYPE_H
#define INTEGER_TYPE_H

#include <CGAL/basic.h>

#if CGAL_USE_GMP && CGAL_USE_MPFI
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz      Integer;
#elif CGAL_USE_CORE
#include <CGAL/CORE_BigInt.h>
typedef CORE::BigInt    Integer;
#else
#include <CGAL/leda_integer.h>
typedef LEDA::integer   Integer;
#endif

#endif
