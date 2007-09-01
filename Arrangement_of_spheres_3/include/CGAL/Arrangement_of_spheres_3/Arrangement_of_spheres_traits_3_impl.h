#include <CGAL/Arrangement_of_spheres_traits_3.h>
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Coordinate_index.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_key.h>
CGAL_BEGIN_NAMESPACE

							

// predicates -------------------------------------------------------------------------------

CGAL_AOS3_TEMPLATE
bool Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::intersects(Sphere_3_key a,
								Sphere_3_key b) const {
  typedef CGAL_AOS3_TYPENAME Table::Equal_centers_exception SE;
  try {
    Point_3 eqpt= table_->equipower_point(a,b);
    std::cout << eqpt  << std::endl;
    CGAL::Bounded_side bs=table_->sphere(a).bounded_side(eqpt);
    //CGAL::Bounded_side obs=table_->sphere(b).bounded_side(eqpt);
    return bs != CGAL::ON_UNBOUNDED_SIDE;
  } catch (SE q) {
    std::cout << "Equal centers for " << a << " and " << b << std::endl;
    return false;
  }
  /*}
    return ret;*/
}



CGAL_AOS3_TEMPLATE
bool Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::intersects(Sphere_3_key a,
								Sphere_3_key b,
								Sphere_3_key c) const {
  //CGAL_precondition(intersects(a,b));
  //CGAL_precondition(intersects(b,c));
  //CGAL_precondition(intersects(c,a));
  // NOTE could be better, I think
  Plane_3 pab = table_->equipower_plane(a,b);
  Plane_3 pac = table_->equipower_plane(a,c);
  CGAL::Object o= geom_traits_object().intersect_3_object()(pab, pac);
  Line_3 l;
  if (CGAL::assign(l,o)){
    Sphere_point_3 sp(table_->sphere(a), l);
    return sp.is_valid();
  } else {
    std::cout << "Degeneracy for " << a << " " << b << " " << c
	      << " with planes " << pac << " and " << pab << std::endl;
    return pab== pac || pab.opposite() == pac;
  }
}





CGAL_AOS3_TEMPLATE
CGAL::Bounded_side 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::rules_bounded_side_of_sphere(const Sphere_point_3 &t,
                                                                             Sphere_3_key x,
                                                                             Sphere_3_key y,
                                                                             Sphere_3_key sphere) const {
  CGAL_AOS3_CANONICAL_PREDICATE(rule_plane(x, CGAL_AOS3_INTERNAL_NS::plane_coordinate(0)),
                                rule_plane(y, CGAL_AOS3_INTERNAL_NS::plane_coordinate(1)),
                                table_->sphere(sphere),
                                t,
                                return ON_BOUNDARY,
                                return ON_BOUNDED_SIDE,
                                return ON_UNBOUNDED_SIDE,
                                return ON_UNBOUNDED_SIDE,
                                return ON_UNBOUNDED_SIDE);
}



/*CGAL_AOS3_TEMPLATE
CGAL::Bounded_side 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::extremum_bounded_side_of_sphere(const Event_point_3 &t,
										Sphere_3_key pt,
										Rule_direction rd,
										Sphere_3_key sphere) const{
  Bounded_side ibs=point_bounded_side_of_sphere(sphere_events(pt).first, sphere);
  std::cout << "ibs is " << ibs << std::endl;
  if (!intersects(pt, sphere) ) {
    return ibs;
  } else {
    Event_pair ep = sphere_intersect_extremum_events(pt, rd.constant_coordinate(), sphere);
    std::cout << "EPs are " << ep.first << " and " << ep.second << std::endl;
    int cross=0;
    if (ep.first.is_valid()) {
      Comparison_result cr= compare_point_to_rule_c(ep.first, pt, 
						    CGAL_AOS3_INTERNAL_NS::other_plane_coordinate(rd.constant_coordinate()));
      if ((rd.is_positive() && cr == LARGER || rd.is_negative() && cr==SMALLER)
	  && compare_points_c(ep.first, t, CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) == CGAL::SMALLER) {
	++cross;
      }
    }
    if (ep.second.is_valid()) {
      Comparison_result cr= compare_point_to_rule_c(ep.second, pt, other_plane_coordinate(rd.constant_coordinate()));
      if ((rd.is_positive() && cr == LARGER || rd.is_negative() && cr==SMALLER)
	  && compare_points_c(ep.second, t, CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) == CGAL::SMALLER) {
	++cross;
      }
    }
    std::cout << "cross is " << cross << std::endl;
    if (cross%2==1) {
      return Bounded_side(-ibs);
    } else {
      return ibs;
    }
  }
  }*/
 





  
  









CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_to_rule_c(const Sphere_point_3 &pt,
									Sphere_3_key rule,
									Coordinate_index C) const{
  /*std::cout << "Comparing point " << pt << " to rule " << table_->center(rule) 
    << " on coordinate " << C << std::endl;*/
  return pt.compare(table_->center(rule)[C.index()], C);
}




CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_to_rule_c(const Center_point_3 &pt,
								  Sphere_3_key rule,
								  Coordinate_index C) const{
  return compare_sphere_centers_c(pt.key(), rule, C);
}








/*CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_point_to_equipower_point_c(const Sphere_point_3 &sp,
										   Sphere_3_key a, 
										   Sphere_3_key b,
										   Coordinate_index C) const{
  // NOTE redo this with computing it directly
  Point_3 pt= table_->equipower_point(a,b);
  return sp.compare(pt, C);
}

CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_center_to_equipower_point_c(Sphere_3_key k,
										   Sphere_3_key a, 
										   Sphere_3_key b,
										   Coordinate_index C) const{
  // NOTE redo this with computing it directly
  Point_3 pt= table_->equipower_point(a,b);
  return compare(table_->center(k)[C.index()], pt[C.index()]);
}
*/



										   





CGAL_AOS3_TEMPLATE
CGAL::Comparison_result
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_sphere_centers_c(Sphere_3_key a, Sphere_3_key b,
									 Coordinate_index C) const {
  return CGAL::compare(table_->center(a)[C.index()], table_->center(b)[C.index()]);
}






// Internal predicates



























// Constructions----------------------------------------------------------------------------------------






CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Event_pair 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::sphere_events(Sphere_3_key ind) const {
  FT v[3];
  v[CGAL_AOS3_INTERNAL_NS::plane_coordinate(0).index()]=0;
  v[CGAL_AOS3_INTERNAL_NS::plane_coordinate(1).index()]=0;
  v[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()]=1;
  return Event_pair(Sphere_point_3(table_->sphere(ind), 
				   Line_3(table_->center(ind), Vector_3(v[0], v[1], v[2]))),
		    Sphere_point_3(table_->sphere(ind), 
				   Line_3(table_->center(ind), Vector_3(-v[0], -v[1], -v[2]))));
}





CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Event_pair 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::intersection_2_events(Sphere_3_key a, Sphere_3_key b) const {
  std::cout << "Intersection events for " << a << " and " << b << std::endl;
  if (!intersects(a,b)) {
    return Event_pair();
  }
  Plane_3 eqp= table_->equipower_plane(a,b);
  Plane_3 sp= table_->separating_plane(a,b);
  //CGAL_assertion(eqp.has_on(eqpoint));
  CGAL::Object o= geom_traits_object().intersect_3_object()(eqp, sp);
  std::cout << "eqp= " << eqp << std::endl << "sp= " << sp << std::endl;
  //<< "eqpoint= " << eqpoint << std::endl;
  Line_3 l;
  if (!CGAL::assign(l,o)){
    std::cout << "Two spheres with same center" << std::endl;
    Event_point_3 ep(table_->equipower_point(a,b));
    return Event_pair(ep, Event_point_3());
  }
  Line_3 lf, lb;
  
  if (l.to_vector()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()] > 0){
    lf=l;
    lb=l.opposite();
  } else {
    lf=l.opposite();
    lb=l;
  }
  if (l.to_vector()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()]==0) {
    std::cerr << "Degeneracy, picking one for start of circle " 
	      << a << " " << b << std::endl;
  }
  
  std::cout << table_->sphere(a) << std::endl 
	    << "lf= " << lf << std::endl 
	    << "lb= " << lb << std::endl;
  Sphere_point_3 spa(table_->sphere(a), lf);
  Sphere_point_3 spb(table_->sphere(a), lb);
  std::cout << spa << std::endl << spb << std::endl;
  Event_pair ret(spa, spb);
  //CGAL_assertion(compare_depths(ret.first, ret.second) != CGAL::LARGER);
  return ret;
}
 












CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Sphere_3_key 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::new_sphere_3(const Sphere_3 &s) {
  return table_->new_sphere(s);
}







CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Event_pair 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::intersection_3_events(Sphere_3_key a,
								      Sphere_3_key b, 
								      Sphere_3_key c) const {
  CGAL_precondition(a != b && a != c && c!= a);
  if (!intersects(a,b, c)) {
    return Event_pair();
  }
  Plane_3 eqpab= table_->equipower_plane(a,b);
  Plane_3 eqpac= table_->equipower_plane(a,c);

  CGAL::Object o= geom_traits_object().intersect_3_object()(eqpab, eqpac);
  Line_3 l;
  if (!CGAL::assign(l,o)){
    // they intersect on the same circle. Ick.
    CGAL_assertion(0);
  }
  Line_3 lf, lb;
  CGAL_assertion(l.to_vector()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()] != 0);
  if (l.to_vector()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()] > 0){
    lf=l;
    lb=l.opposite();
  } else {
    lf=l.opposite();
    lb=l;
  }
  if (l.to_vector()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()] ==0) {
    std::cerr << "Degeneracy, picking one for start of threesome " 
	      << a << " " << b << " " << c << std::endl;
  }
  
  Event_pair ret(Sphere_point_3(table_->sphere(a), lf),
		 Sphere_point_3(table_->sphere(a), lb));
  //CGAL_assertion(compare_depths(ret.first, ret.second) != CGAL::LARGER);
  return ret;
}






CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Plane_3 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::rule_plane(Sphere_3_key r, Coordinate_index C) const {
  FT n[3]={0,0,0};
  n[C.index()]= 1;
  return Plane_3(table_->center(r), Vector_3(n[0], n[1], n[2]));
}




CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Plane_3 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::const_c_plane(const Sphere_point_3 &pt, Coordinate_index C) const {
  FT v[3]={0,0,0};
  v[other_plane_coordinate(C).index()]=1;

 return Plane_3(pt.line().point(), pt.line().point() + pt.line().to_vector(),
		pt.line().point()+ Vector_3(v[0], v[1], v[2]));
}



CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Event_pair 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::sphere_intersect_extremum_events(Sphere_3_key a, 
										 Coordinate_index C,
										 Sphere_3_key b) const {
  CGAL_precondition(intersects(a,b));
  Plane_3 eqpab= table_->equipower_plane(a,b);
  Plane_3 rp= rule_plane(a,C);
  CGAL::Object o= geom_traits_object().intersect_3_object()(eqpab, rp);
  Line_3 l;
  if (!CGAL::assign(l,o)){
    return Event_pair();
  } else {
    Sphere_point_3 sa(table_->sphere(a), l);
    if (!sa.is_valid()) return Event_pair();
    Event_point_3 pa(sa);
    Event_point_3 pb(Sphere_point_3(table_->sphere(a), l.opposite()));
    if (pb < pa) std::swap(pa, pb);
    return Event_pair(pa, pb);
  }
}


CGAL_AOS3_TEMPLATE
void 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::advance_circle_cross_rule_event(Sphere_3_key a, 
										Sphere_3_key b,
										Sphere_3_key rs,
										Coordinate_index C) {
  CGAL_AOS3_INTERNAL_NS::Sphere_key_pair k(a,b);
  CGAL_AOS3_TYPENAME Table::Pair_data::KC_pair kc(rs, C);
  CGAL_assertion(table_->pair_data(k).cxr_events_.find(kc) != table_->pair_data(k).cxr_events_.end());
  CGAL_AOS3_TYPENAME Table::Event_pair_data &ed = table_->pair_data(k).cxr_events_[kc];
  CGAL_assertion(ed.index_ != -1);
  ++ed.index_;
}






CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Event_point_3 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::circle_cross_rule_event(Sphere_3_key a, 
									Sphere_3_key b,
									Sphere_3_key rs,
									Coordinate_index C) const {
  CGAL_AOS3_INTERNAL_NS::Sphere_key_pair k(a,b);
  CGAL_AOS3_TYPENAME Table::Event_pair_data &ed
    = table_->pair_data(k).cxr_events_[CGAL_AOS3_TYPENAME Table::Pair_data::KC_pair(rs, C)];

  if (ed.index_==-1) {
    std::cout << "Initializing AAR for " << a << " " << b << " " << rs << " " << C << std::endl;
    Event_pair ep= circle_cross_rule_event_internal(a,b,rs,C);
    ed.events_[0]=ep.first;
    ed.events_[1]=ep.second;
    CGAL_assertion(!ed.events_[0].is_valid()
		   || oriented_side_of_separating_plane(ed.events_[0], a,b) == CGAL::ON_POSITIVE_SIDE);
    CGAL_assertion(!ed.events_[1].is_valid()
		   || oriented_side_of_separating_plane(ed.events_[1], a,b) == CGAL::ON_POSITIVE_SIDE);
    ed.index_=0;
  }
  if (ed.index_ < 2) {
    std::cout << "Returning " << ed.index_ << " of AAR" << std::endl;
    CGAL_assertion(!ed.events_[ed.index_].is_valid()
		   || oriented_side_of_separating_plane(ed.events_[ed.index_], a,b) == CGAL::ON_POSITIVE_SIDE);
    return ed.events_[ed.index_];
  }else return Event_point_3();
}







CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Event_pair 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::circle_cross_rule_event_internal(Sphere_3_key a, 
										 Sphere_3_key b,
										 Sphere_3_key rs,
										 Coordinate_index C) const {
  Plane_3 eqp= table_->equipower_plane(a,b);
  /*FT na[3]={0,0,0};
    na[C.index()]=1;*/
  Plane_3 rp = rule_plane(rs, C); //(table_->center(rs), Vector_3(na[0], na[1], na[2]));
  CGAL::Object o=geom_traits_object().intersect_3_object()(eqp, rp);
  Line_3 l;
  if (CGAL::assign(l,o)) {
    Sphere_point_3 sp(table_->sphere(a), l);
    //Event_point_3 sp();
    if (!sp.is_valid()) {
      return Event_pair();
    } else {
      Event_point_3 nep(sp);
      Event_point_3 spp(Sphere_point_3(table_->sphere(a), l.opposite()));
      std::cout << "AAR events are " << nep << " and " << spp << std::endl;
      std::cout << "l is " << l << std::endl;
      Event_pair ep= (nep <spp)? std::make_pair(nep, spp) : std::make_pair(spp, nep);
      Oriented_side os0=oriented_side_of_separating_plane(ep.first, a,b);
      Oriented_side os1=oriented_side_of_separating_plane(ep.second, a,b);
      if (os0==ON_NEGATIVE_SIDE || a>b && os0==ON_ORIENTED_BOUNDARY){
	ep.first=ep.second;
	ep.second=Event_point_3();
	if (os1== CGAL::ON_NEGATIVE_SIDE || a>b && os1 == CGAL::ON_ORIENTED_BOUNDARY){
	  ep.first= Event_point_3();
	}
      } else if (os1== CGAL::ON_NEGATIVE_SIDE || a>b && os1 == CGAL::ON_ORIENTED_BOUNDARY){
	  ep.second= Event_point_3();
      }
      return ep;
    }
  } else {
    return Event_pair();
  }
}

CGAL_AOS3_TEMPLATE
void 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::advance_sphere_intersect_rule_rule_event(Sphere_3_key s, 
											 Sphere_3_key x,
											 Sphere_3_key y) {
  CGAL_AOS3_INTERNAL_NS::Sphere_key_triple k(s,x,y);
  CGAL_AOS3_TYPENAME Table::Event_pair_data &ed= table_->triple_data(k).srr_events_;
  CGAL_assertion(ed.index_ != -1);
  ++ed.index_;
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Event_point_3 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::sphere_intersect_rule_rule_event(Sphere_3_key s, 
										 Sphere_3_key x,
										 Sphere_3_key y) const {
  CGAL_AOS3_INTERNAL_NS::Sphere_key_triple k(s,x,y);
  CGAL_AOS3_TYPENAME Table::Event_pair_data &ed= table_->triple_data(k).srr_events_;
  if (ed.index_==-1) {
    Event_pair ep= sphere_intersect_rule_rule_events(s,x,y);
    ed.events_[0]=ep.first;
    ed.events_[1]=ep.second;
    ed.index_=0;
  }
  if (ed.index_ < 2) return ed.events_[ed.index_];
  else return Event_point_3();
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Event_pair 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::sphere_intersect_rule_rule_events(Sphere_3_key s, 
										  Sphere_3_key x,
										  Sphere_3_key y) const {
  CGAL_precondition(s != Sphere_3_key());
  CGAL_precondition(x != Sphere_3_key());
  CGAL_precondition(y != Sphere_3_key());

  std::cout << "SRR event for " << s << " " << x << " " << y << std::endl;

  Plane_3 x_plane= rule_plane(x, CGAL_AOS3_INTERNAL_NS::plane_coordinate(0));
  Plane_3 y_plane= rule_plane(y, CGAL_AOS3_INTERNAL_NS::plane_coordinate(1));
  CGAL::Object o= geom_traits_object().intersect_3_object()(x_plane, y_plane);
  Line_3 l;
  if (!CGAL::assign(l,o)){
    return Event_pair();
  } else {
    std::cout << "Line is " << l << std::endl;
    std::cout << "Sphere is " << table_->sphere(s) << std::endl;
    Sphere_point_3 sa(table_->sphere(s), l);
    if (!sa.is_valid()) return Event_pair();
    Event_point_3 pa(sa);
    Event_point_3 pb(Sphere_point_3(table_->sphere(s), l.opposite()));
    if (pb < pa) std::swap(pa, pb);
    std::cout << "Ret is " << pa << " " << pb << std::endl;
    return Event_pair(pa, pb);
  }
}


CGAL_END_NAMESPACE
