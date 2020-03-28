
#ifndef CGAL_SDG_VERBOSE
#define CGAL_SDG_DEBUG(a)
#else
#define CGAL_SDG_DEBUG(a) { a }
#endif

#include <iostream>
#include <fstream>
#include <cassert>

#include <CGAL/Exact_rational.h>

// choose number type
typedef CGAL::Exact_rational exact_ring_t;
typedef CGAL::Exact_rational exact_field_t;

#include <CGAL/Simple_cartesian.h>

#include <CGAL/Segment_Delaunay_graph_Linf_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h>

typedef CGAL::Simple_cartesian<double>          CK;
typedef CGAL::Simple_cartesian<exact_ring_t>    EK_ring;
typedef CGAL::Simple_cartesian<exact_field_t>   EK_field;

typedef CGAL::Field_with_sqrt_tag  CK_MTag;
typedef CGAL::Field_tag       EK_MTag;
typedef CGAL::Integral_domain_without_division_tag        EK_MTag_wi;

typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<CK,CK_MTag,
                                                       EK_field,EK_MTag>
Gt;

typedef
CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2<CK,
                                                                     CK_MTag,
                                                                     EK_ring,
                                                                     EK_MTag_wi>
Gt_wi;

typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt>      SDG2;
typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt_wi>   SDG2_wi;

#include "test_types.h"
#include "sdg_main.h"

// Comment for cgal_create_makefile : main()
