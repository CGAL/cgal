#include <CGAL/basic.h>

#ifndef CGAL_SDG_VERBOSE
#define CGAL_SDG_DEBUG(a)
#else
#define CGAL_SDG_DEBUG(a) { a }
#endif

#include <iostream>
#include <fstream>
#include <cassert>

// choose number type
#ifdef CGAL_USE_GMP

#  include <CGAL/Gmpq.h>
typedef CGAL::Gmpq                     exact_ring_t;
typedef CGAL::Gmpq                     exact_field_t;

namespace CGAL {
// needed for the drawing methods
Gmpq sqrt(const Gmpq& x) {
  return Gmpq(  sqrt( to_double(x) )  );
}

} //namespace CGAL
#else

#  include <CGAL/MP_Float.h>
#  include <CGAL/Quotient.h>
typedef CGAL::MP_Float                 exact_ring_t;
typedef CGAL::Quotient<exact_ring_t>   exact_field_t;

#endif

typedef exact_ring_t   ring_number_t;
typedef exact_field_t  field_number_t;

#include <CGAL/Simple_cartesian.h>

#include <CGAL/Segment_Delaunay_graph_Linf_hierarchy_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_traits_2.h>

struct K_ring  : public CGAL::Simple_cartesian<ring_number_t> {};
struct K_field : public CGAL::Simple_cartesian<field_number_t> {};

typedef CGAL::Field_tag  MTag;
typedef CGAL::Integral_domain_without_division_tag   MTag_wi;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<K_field,MTag> Gt;

typedef
CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2<K_ring,MTag_wi>
Gt_wi;

typedef CGAL::Segment_Delaunay_graph_Linf_hierarchy_2<Gt>      SDG2;
typedef CGAL::Segment_Delaunay_graph_Linf_hierarchy_2<Gt_wi>   SDG2_wi;

#include "test_types_hierarchy.h"
#include "sdg_hierarchy_main.h"

// Comment for cgal_create_makefile : main()
