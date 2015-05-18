#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <cassert>

// choose number type

// choose number type
#ifdef CGAL_USE_LEDA

#  include <CGAL/leda_rational.h>
typedef leda_rational                     exact_ring_t;
typedef leda_rational                     exact_field_t;

namespace CGAL {
// needed for the drawing methods
leda_rational sqrt(const leda_rational& x) {
  return leda_rational(  CGAL::sqrt( CGAL::to_double(x) )  );
}
}
#elif defined( CGAL_USE_GMP )

#  include <CGAL/Gmpq.h>
typedef CGAL::Gmpq                     exact_ring_t;
typedef CGAL::Gmpq                     exact_field_t;

namespace CGAL {
// needed for the drawing methods
Gmpq sqrt(const Gmpq& x) {
  return Gmpq(  CGAL::sqrt( CGAL::to_double(x) )  );
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

#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>

typedef CGAL::Simple_cartesian<ring_number_t> K_ring;
typedef CGAL::Simple_cartesian<field_number_t> K_field;

typedef CGAL::Field_tag  MTag;
typedef CGAL::Integral_domain_without_division_tag   MTag_wi;

typedef CGAL::Segment_Delaunay_graph_traits_2<K_field,MTag> Gt;

typedef
CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K_ring,MTag_wi>
Gt_wi;

typedef CGAL::Segment_Delaunay_graph_2<Gt>      SDG2;
typedef CGAL::Segment_Delaunay_graph_2<Gt_wi>   SDG2_wi;

#include "test_types.h"
#include "sdg_main.h"

// Comment for cgal_create_makefile : main()

