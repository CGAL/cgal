
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
typedef CGAL::Exact_rational ring_number_t;
typedef CGAL::Exact_rational field_number_t;

#include <CGAL/Simple_cartesian.h>

#include <CGAL/Segment_Delaunay_graph_Linf_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_traits_2.h>

typedef CGAL::Simple_cartesian<ring_number_t> K_ring;
typedef CGAL::Simple_cartesian<field_number_t> K_field;

typedef CGAL::Field_tag  MTag;
typedef CGAL::Integral_domain_without_division_tag   MTag_wi;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<K_field,MTag> Gt;

typedef
CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2<K_ring,MTag_wi>
Gt_wi;

typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt>      SDG2;
typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt_wi>   SDG2_wi;

#include "test_types.h"
#include "sdg_main.h"

// Comment for cgal_create_makefile : main()
