#ifndef KERNEL_TYPE_HPP
#define KERNEL_TYPE_HPP

#include <CGAL/basic.h>

#include "bench_config.hpp"

// Kernel:
#if BENCH_KERNEL == EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#elif BENCH_KERNEL == LEDA_KERNEL
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>

#elif BENCH_KERNEL == MY_KERNEL
#include <CGAL/Arr_segment_traits_leda_kernel_2.h>

#elif BENCH_KERNEL == CARTESIAN_KERNEL
#include <CGAL/Cartesian.h>

#elif BENCH_KERNEL == SIMPLE_CARTESIAN_KERNEL
#include <CGAL/Simple_cartesian.h>

#elif BENCH_KERNEL == LAZY_CARTESIAN_KERNEL
#include <CGAL/Cartesian.h>
#include <CGAL/Lazy_kernel.h>

#elif BENCH_KERNEL == LAZY_SIMPLE_CARTESIAN_KERNEL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Lazy_kernel.h>

#else
#error No kernel (KERNEL) specified!
#endif

// Typedefs:
#if BENCH_KERNEL == EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL
typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
#define KERNEL_TYPE "Exact"

#elif BENCH_KERNEL == LEDA_KERNEL
typedef CGAL::leda_rat_kernel_traits                            Kernel;
#define KERNEL_TYPE "Leda"

#elif BENCH_KERNEL == MY_KERNEL
typedef CGAL::Arr_segment_traits_leda_kernel_2                  Kernel;
#define KERNEL_TYPE "Partial Leda"

#elif BENCH_KERNEL == CARTESIAN_KERNEL
typedef CGAL::Cartesian<Number_type>                            Kernel;
#define KERNEL_TYPE "Cartesian"

#elif BENCH_KERNEL == SIMPLE_CARTESIAN_KERNEL
typedef CGAL::Simple_cartesian<Number_type>                     Kernel;
#define KERNEL_TYPE "Simple Cartesian"

#elif BENCH_KERNEL == LAZY_CARTESIAN_KERNEL
typedef CGAL::Lazy_kernel<CGAL::Cartesian<Number_type> >        Kernel;
#define KERNEL_TYPE "Lazy Cartesian"

#elif BENCH_KERNEL == LAZY_SIMPLE_CARTESIAN_KERNEL
typedef CGAL::Lazy_kernel<CGAL::Simple_cartesian<Number_type> > Kernel;
#define KERNEL_TYPE "Lazy Simple Cartesian"

#else
#error No kernel (KERNEL) specified!
#endif

#endif
