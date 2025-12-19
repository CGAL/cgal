#ifndef ARR_RAT_FUNCTIONS_H
#define ARR_RAT_FUNCTIONS_H

#include <CGAL/basic.h>
#include <CGAL/CORE_BigInt.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Arr_rational_function_traits_2.h>
#include <CGAL/Arrangement_2.h>

using Number_type = CORE::BigInt;
using AK1 = CGAL::Algebraic_kernel_d_1<Number_type>;
using Traits = CGAL::Arr_rational_function_traits_2<AK1>;

using Polynomial = Traits::Polynomial_1;
using Alg_real = Traits::Algebraic_real_1;
using Bound = Traits::Bound;

using Arrangement = CGAL::Arrangement_2<Traits>;

#endif
