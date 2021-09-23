#ifndef ARR_RAT_FUNCTIONS_H
#define ARR_RAT_FUNCTIONS_H

#include <CGAL/basic.h>
#include <CGAL/CORE_BigInt.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Arr_rational_function_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CORE::BigInt                               Number_type;
typedef CGAL::Algebraic_kernel_d_1<Number_type>    AK1;
typedef CGAL::Arr_rational_function_traits_2<AK1>  Traits;

typedef Traits::Polynomial_1                       Polynomial;
typedef Traits::Algebraic_real_1                   Alg_real;
typedef Traits::Bound                              Bound;

typedef CGAL::Arrangement_2<Traits>                Arrangement;

#endif
