#ifndef BENCH_CONFIG_H
#define BENCH_CONFIG_H

#define CARTESIAN_KERNEL                0
#define SIMPLE_CARTESIAN_KERNEL         1
#define LEDA_KERNEL                     2
#define MY_KERNEL                       3

#define SEGMENT_TRAITS                  0
#define SEGMENT_CACHED_TRAITS           1
#define LEDA_SEGMENT_TRAITS             2
#define POLYLINE_TRAITS                 3
#define CACHED_POLYLINE_TRAITS          4
#define CONIC_TRAITS                    5
#define EXACUS_CONIC_TRAITS             6
#define CK_CONIC_TRAITS                 7
#define CORE_CONIC_TRAITS               8
#define CK_CIRCLE_TRAITS                9

#define DOUBLE_NT                       0
#define MP_FLOAT_NT                     1
#define GMPZ_NT                         2
#define LEDA_RAT_NT                     3
#define QUOTIENT_MP_FLOAT_NT            4
#define QUOTIENT_CGAL_GMPZ_NT           5
#define GMPQ_NT                         6
#define CGAL_GMPQ_NT                    7
#define LAZY_LEDA_RAT_NT                8
#define LAZY_CGAL_GMPQ_NT               9
#define LAZY_QUOTIENT_MP_FLOAT_NT       10
#define LEDA_REAL_NT                    11
#define CORE_EXPR_NT                    12
#define NIX_LEDA_FIELD_WITH_SQRT_NT     13
#define NIX_CORE_FIELD_WITH_SQRT_NT     14
#define LAZY_GMPZ_NT                    15

// Default value based on dependencies:
#if BENCH_KERNEL == LEDA_KERNEL
#if !defined(BENCH_NT)
#define BENCH_NT LEDA_RAT_NT
#endif
#if !defined(BENCH_TRAITS)
#define BENCH_TRAITS SEGMENT_TRAITS
#endif
#endif

#if BENCH_KERNEL == MY_KERNEL
#if !defined(BENCH_NT)
#define BENCH_NT LEDA_RAT_NT
#endif
#if !defined(BENCH_TRAITS)
#define BENCH_TRAITS LEDA_SEGMENT_TRAITS
#endif
#endif

#if BENCH_TRAITS == LEDA_SEGMENT_TRAITS
#if !defined(BENCH_NT)
#define BENCH_NT LEDA_RAT_NT
#endif
#if !defined(BENCH_KERNEL)
#define BENCH_KERNEL LEDA_KERNEL
#endif
#endif

#if BENCH_TRAITS == CONIC_TRAITS
#if !defined(BENCH_KERNEL)
#define BENCH_KERNEL = CARTSIAN_KERNEL
#endif
#if !defined(BENCH_NT)
#define BENCH_NT LEDA_REAL_NT
#endif
#endif

#if BENCH_TRAITS == EXACUS_CONIC_TRAITS
#if !defined(BENCH_NT)
#define BENCH_NT NIX_LEDA_FIELD_WITH_SQRT_NT
#endif
#endif

#if BENCH_TRAITS == CORE_CONIC_TRAITS
#if !defined(BENCH_NT)
#define BENCH_NT CORE_EXPR_NT
#endif
#endif

// Default values:
#if !defined(BENCH_KERNEL)
#define BENCH_KERNEL CARTESIAN_KERNEL
#endif

#if !defined(BENCH_NT)
#define BENCH_NT LEDA_RAT_NT
#endif

#if !defined(BENCH_TRAITS)
#define BENCH_TRAITS SEGMENT_TRAITS
#endif

// Illegal combinations:
#if BENCH_KERNEL == LEDA_KERNEL && BENCH_NT != LEDA_RAT_NT
#error Leda kernel implies rational number type!
#endif

#if BENCH_KERNEL == MY_KERNEL && BENCH_NT != LEDA_RAT_NT
#error My kernel implies rational number type!
#endif

#if BENCH_KERNEL == MY_KERNEL && BENCH_TRAITS != LEDA_SEGMENT_TRAITS
#error My kernel implies leda segment traits!
#endif

#if BENCH_TRAITS == CONIC_TRAITS && BENCH_NT != LEDA_REAL_NT
#error "Conic traits implies real number type!"
#endif

#if BENCH_TRAITS == CONIC_TRAITS && BENCH_KERNEL == LEDA_KERNEL
#error "Conic traits implies non leda kernel!"
#endif

#if BENCH_TRAITS == CONIC_TRAITS && BENCH_KERNEL == MY_KERNEL
#error "Conic traits implies non my kernel!"
#endif

#if BENCH_TRAITS == LEDA_SEGMENT_TRAITS && BENCH_NT != LEDA_RAT_NT
#error "Leda segment traits implies rational number type!"
#endif

#if BENCH_TRAITS == LEDA_SEGMENT_TRAITS && BENCH_KERNEL != LEDA_KERNEL && \
    BENCH_KERNEL != MY_KERNEL
#error "Leda segment traits implies leda kernel or my kernel!"
#endif

#if (BENCH_NT == NIX_LEDA_FIELD_WITH_SQRT_NT || \
     BENCH_NT == NIX_CORE_FIELD_WITH_SQRT_NT) && \
    BENCH_TRAITS != EXACUS_CONIC_TRAITS
#error "NIX number type implies EXACUS conic traits!"
#endif

#if BENCH_TRAITS == CORE_CONIC_TRAITS && BENCH_NT != CORE_EXPR_NT
#error "Core conic traits implies Core Expr number type!"
#endif

#endif
