#ifndef CGAL_AOS3_RATIONAL_CROSS_SECTION_H
#define CGAL_AOS3_RATIONAL_CROSS_SECTION_H
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_cross_section.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
class Rational_cross_section{
  typedef Rational_cross_section CGAL_AOS3_TARG This;
  CGAL_AOS3_TRAITS;
  typedef Combinatorial_cross_section CGAL_AOS3_TARG CCS;
  const CCS &ccs_;
  Traits tr_;
  CGAL_AOS3_TYPENAME Traits::FT z_;
public:
   Rational_cross_section(const CCS &ccs, const Traits &tr): ccs_(ccs), tr_(tr){}

  void set_z(CGAL_AOS3_TYPENAME Traits::FT z) {
    z_=z;
  }
  const CGAL_AOS3_TYPENAME Traits::FT & z() const {
    return z_;
  }
  

  CGAL_AOS3_TYPENAME Traits::Point_2 
  compute_rule_rule_intersection(CGAL_AOS3_TYPENAME Traits::Sphere_3_key ra,
				 CGAL_AOS3_TYPENAME Traits::Sphere_3_key rb) const {
 
  return CGAL_AOS3_TYPENAME Traits::Point_2(tr_.sphere_3(ra).center()[plane_coordinate(0).index()], 
		    tr_.sphere_3(rb).center()[plane_coordinate(1).index()]);
}


  CGAL_AOS3_TYPENAME Traits::Sphere_point_3 
  sphere_point(CGAL_AOS3_TYPENAME CCS::Point pt) const;


  CGAL_AOS3_TYPENAME Traits::Point_2 
  rule_rule_intersection(CGAL_AOS3_TYPENAME Traits::Sphere_3_key ra,
			 CGAL_AOS3_TYPENAME Traits::Sphere_3_key rb) const;


  CGAL_AOS3_TYPENAME Traits::Line_3
  positive_line(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
		Coordinate_index i) const ;


 CGAL_AOS3_TYPENAME Traits::Line_3
  negative_line(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
		Coordinate_index i) const ;


  CGAL_AOS3_TYPENAME Traits::Circle_2
  circle(CGAL_AOS3_TYPENAME Traits::Sphere_3_key a) const;
    
  bool intersects(CGAL_AOS3_TYPENAME Traits::Sphere_3_key a) const ;


  CGAL_AOS3_TYPENAME Traits::Point_2
  center_point(CGAL_AOS3_TYPENAME Traits::Sphere_3_key a, 
	       CGAL_AOS3_TYPENAME Traits::Sphere_3_key b) const ;


  void audit() const;

};

CGAL_AOS3_END_INTERNAL_NAMESPACE


#ifdef CGAL_AOS3_USE_TEMPLATES
#include <CGAL/Arrangement_of_spheres_3/Rational_cross_section_impl.h>
#endif
#endif
