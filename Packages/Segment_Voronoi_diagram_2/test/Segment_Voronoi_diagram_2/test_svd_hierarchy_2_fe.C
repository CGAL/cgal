#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <cassert>

// choose number type
#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpq.h>
typedef CGAL::Gmpq                     exact_ring_t;
typedef CGAL::Gmpq                     exact_field_t;
#else
#  include <CGAL/MP_Float.h>
#  include <CGAL/Quotient.h>
typedef CGAL::MP_Float                 exact_ring_t;
typedef CGAL::Quotient<exact_ring_t>   exact_field_t;
#endif

#include <CGAL/Filtered_exact.h>

typedef CGAL::Filtered_exact<double,exact_ring_t>   ring_number_t;
typedef CGAL::Filtered_exact<double,exact_field_t>  field_number_t;

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Number_type_traits.h>

#include <CGAL/Segment_Voronoi_diagram_hierarchy_2.h>
#include <CGAL/Segment_Voronoi_diagram_traits_2.h>

struct K_ring  : public CGAL::Simple_cartesian<ring_number_t> {};
struct K_field : public CGAL::Simple_cartesian<field_number_t> {};

typedef CGAL::Field_tag  MTag;
typedef CGAL::Ring_tag   MTag_wi;

typedef CGAL::Segment_Voronoi_diagram_traits_2<K_field,MTag> Gt;

typedef
CGAL::Segment_Voronoi_diagram_traits_without_intersections_2<K_ring,
							     MTag_wi>
Gt_wi;

typedef CGAL::Segment_Voronoi_diagram_hierarchy_2<Gt>      SVD2;
typedef CGAL::Segment_Voronoi_diagram_hierarchy_2<Gt_wi>   SVD2_wi;

#include "test_types_hierarchy.h"
#include "svd_hierarchy_main.h"


