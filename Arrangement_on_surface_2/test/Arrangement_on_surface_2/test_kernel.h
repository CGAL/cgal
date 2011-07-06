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

#else
#error No kernel (KERNEL) specified!
#endif

#endif
