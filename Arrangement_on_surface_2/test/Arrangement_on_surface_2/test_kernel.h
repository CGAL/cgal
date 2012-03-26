#ifndef CGAL_TEST_KERNEL_H
#define CGAL_TEST_KERNEL_H

#include "test_number_type.h"

// ============================================================================
// Kernel includes:
// ============================================================================
#if TEST_KERNEL == SIMPLE_CARTESIAN_KERNEL
#include <CGAL/Simple_cartesian.h>

#elif TEST_KERNEL == CARTESIAN_KERNEL
#include <CGAL/Cartesian.h>

#elif TEST_KERNEL == UNIVARIATE_ALGEBRAIC_KERNEL
#include <CGAL/Algebraic_kernel_d_1.h>

#else
#error No kernel (KERNEL) specified!
#endif

// ============================================================================
// Kernel typedef:
// ============================================================================
#if TEST_KERNEL == SIMPLE_CARTESIAN_KERNEL
typedef CGAL::Simple_cartesian<Number_type>             Kernel;
#define KERNEL_TYPE "Simple Cartesian"

#elif TEST_KERNEL == CARTESIAN_KERNEL
typedef CGAL::Cartesian<Number_type>                    Kernel;
#define KERNEL_TYPE "Cartesian"

#elif TEST_KERNEL == UNIVARIATE_ALGEBRAIC_KERNEL
typedef CGAL::Algebraic_kernel_d_1<Number_type>         Kernel;
#define KERNEL_TYPE "Univariate Algebraic"

#else
#error No kernel (KERNEL) specified!
#endif

#endif
