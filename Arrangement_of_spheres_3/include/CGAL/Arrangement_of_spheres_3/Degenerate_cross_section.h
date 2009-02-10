#ifndef CGAL_AOS3_DEGENERATE_CROSS_SECTION_RULES_H
#define CGAL_AOS3_DEGENERATE_CROSS_SECTION_RULES_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_cross_section.h>
#include <CGAL/Arrangement_of_spheres_3/Rational_cross_section.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
class Degenerate_cross_section {
  CGAL_AOS3_TRAITS;
  typedef Degenerate_cross_section CGAL_AOS3_TARG This;
  typedef Combinatorial_cross_section CGAL_AOS3_TARG CCS;
  typedef Rational_cross_section CGAL_AOS3_TARG RCS;
  typedef Irrational_cross_section_insertion CGAL_AOS3_TARG ICSI;
  typedef Irrational_cross_section_removal CGAL_AOS3_TARG ICSD;
  typedef Irrational_cross_section_location CGAL_AOS3_TARG ICSL;
  typedef Irrational_cross_section_rules CGAL_AOS3_TARG ICSR;
  typedef CGAL_AOS3_TYPENAME Traits::Sphere_3_key Sphere_3_key;
  typedef CGAL_AOS3_TYPENAME Traits::Event_point_3 Event_point_3;
  typedef CGAL_AOS3_TYPENAME Traits::Sphere_point_3 Sphere_point_3;
  typedef CGAL_AOS3_TYPENAME Traits::FT NT;
  typedef CGAL_AOS3_TYPENAME std::vector<Sphere_3_key> Rule_list;
 
  //typedef CGAL_AOS3_TYPENAME Traits::Degeneracy_exception Degeneracy;

  
  class Less_center {
    Traits tr_;
    Coordinate_index i_;
  public:
    Less_center(Traits tr, Coordinate_index i): tr_(tr), i_(i){}
    bool operator()(Sphere_3_key a, Sphere_3_key b) const {
      return tr_.compare_sphere_centers_c(a,b,i_) == SMALLER;
    }
  };

  class Equal_center {
    Traits tr_;
    Coordinate_index i_;
  public:
    Equal_center(Traits tr, Coordinate_index i): tr_(tr), i_(i){}
    bool operator()(Sphere_3_key a, Sphere_3_key b) const {
      return tr_.compare_sphere_centers_c(a,b,i_) == EQUAL;
    }
  };

  class Less_center {
    Traits tr_;
    Coordinate_index ci_;
    
  public:
    Compare_to_center(Traits tr, Coordinate_index i): tr_(tr), ci_(i){}
    bool operator()(Sphere_3_key a, const Sphere_points_3 &sp) const {
      return tr_.compare_to_rule_c(sp, a, ci_) == CGAL::LARGER;
    }
    bool operator()( const Sphere_points_3 &sp, Sphere_3_key a) const {
      return tr_.compare_to_rule_c(sp, a, ci_) == CGAL::SMALLER;
    }
  };
  
public:
 
  typedef CGAL_AOS3_TYPENAME CCS::Halfedge_handle Halfedge_handle;
  typedef CGAL_AOS3_TYPENAME CCS::Vertex_handle Vertex_handle;
  typedef CGAL_AOS3_TYPENAME CCS::Face_handle Face_handle;
  typedef CGAL_AOS3_TYPENAME CCS::Point Point;
  typedef CGAL_AOS3_TYPENAME CCS::Curve Curve;

  
  Degenerate_cross_section(CCS &cs, const Traits &tr): tr_(tr),
						      cs_(cs),
						       rcs_(cs, tr){
    for (int i=0; i< 2; ++i) {
      sorted_lists_[i].insert(sorted_lists_[i].end(),
			      tr.sphere_3_keys_begin(),
			      tr.sphere_3_keys_end());
      std::sort(sorted_lists_[i].begin(),
		sorted_lists_[i].end(), Less_center(tr_,i));
      sorted_lists_[i].erase(std::unique(sorted_lists_[i].begin(),
					 sorted_lists_[i].end(),
					 Equal_center(tr_,i)));
    }
  }

  struct Box {
    int &operator[](int i){
      return d_[i];
    }
    int d_[2];
  };
  //typedef std::pair<int,int> Box;

  struct Degeneracy {
    Event_point_3 location_;
    typedef std::map<Vertex_handle, Box,
		     CGAL_AOS3_TYPENAME CCS::Handle_compare> Vertices;
    //Vertices vertices_;
    Vertex_handle a_vertex_;
    //std::vector<Box> boxes_;
    std::set<Curve> curves_;
    typedef std::map<Point, Box> Points;
    Points points_;
    std::vector<Halfedge_handle> outer_boundary_;
    std::vector<Halfedge_handle> inner_boundary_;
    Halfedge_handle outer_rule_points_[4];
    Halfedge_handle inner_rule_points_[4];
  };


  bool box_ok(Box box, Point pt) {
    for (int i=0; i< 2; ++i) {
      if (box[i]==-1) continue;
      Sphere_point_3 sp= rcs_.sphere_point(pt);
      if (tr_.compare_to_rule_c(sp, sorted_lists_[i][box[i]],
				  plane_coordinate(i)) != LARGER 
	  || tr_.compare_to_rule_c(sp, sorted_lists_[i][box[i]+1],
				     plane_coordinate(i)) != SMALLER)
	return false;
    }
    return true;
  }

  int find_box(const Event_point_3 &ep, Point pt, int i) {
    CGAL_AOS3_TYPENAME Rule_list::iterator it0
      = std::lower_bound(sorted_lists_[i].begin(),
			 sorted_lists_[i].end(),
			 pt,
			 Less_center(tr_, cs_, 
				     plane_coordinate(i), ep));
    int offset= it0- sorted_lists_[i].begin();
    // where would equal be?
    if (tr_.compare_to_rule_c(ep, *it0, ci) == CGAL::EQUAL) {
      if (pt.is_rule_rule()) return -1;
      else if (pt.is_sphere_rule() && pt.constant_coordinate() == plane_coordinate(i)) {
	return -1;
      } else if (pt.is_sphere_sphere()) {
	Comparison_result cr= tr_.compare_to_circle_circle_after_c(Center_point_3(ep, sorted_lits_[i][offset], plane_coordinate(i)));
	if (cr== SMALLER) return i;
	else return i-1;
      } else {
	CGAL_assertion(pt.is_sphere_rule());
	Comparison_result cr= tr_.compare_to_circle_circle_after_c(Center_point_3(ep, sorted_lits_[i][offset], plane_coordinate(i)));
	if (cr== SMALLER) return i;
	else return i-1;
      }
    }

    return offset;
  }

  Box find_box(const Event_point_3 &ep, Point pt) {
    //Sphere_point_3 sp= rcs_.sphere_point(pt);
    int b0= find_box(ep, pt, 0);
    int b1= find_box(ep, pt, 1);
    return Box(b0,b1);
  }
   
    
  void operator()(std::vector<Event_point_3> &event_locations) {


    //! \todo change arrangement initializer to have a layer where you
    //! just pass in a set of arcs
    //! \todo write an expand vertex function
    //! \todo add a star hole function (more or less stitch in)


    std::vector<Degeneracy> degeneracies(event_locations.size());
    CGAL_LOG(Log::SOME, "Degeneracy found at \n");
    for (unsigned int i=0; i< event_locations.size(); ++i) {
      CGAL_LOG(Log::SOME, event_locations[i] << std::endl);
      degeneracies[i].location_=event_locations[i];
    }
    
    
    
    //std::vector<std::vector<Vertex_handle> > vertices(event_locations.size());
    for (unsigned int i=0; i< event_locations.size(); ++i) {
      // find all faces which are compatible
      // for each face, search for compatible edges
      // if there are two, insert vertex not connected by an edge,
      // connect and throw new faces on back
   
      CGAL_assertion(0);
     
      // this needs to include intersection pairs
      //degeneracies[i].a_vertex_=find_points(degeneracies[i].location_,
      // degeneracies[i].points_);
    }
    
 
    for (unsigned int i=0; i< degeneracies.size(); ++i) {
      // * figure out outside edges--and what rules come in
      // * find all arcs involved
      // * collapse to a point
      // build list of points

      sc_.collapse_to_vertex(*degeneracies[i].a_vertex_,
			     Is_degenerate(degeneracies[i].location_),
			     Accum_curves_and_points(degeneracies[i].curves_));
    } 

    for (unsigned int i=0; i< degeneracies.size(); ++i) {
      //degeneracies[i].boxes_.resize(degeneracies[i].vertices_.size());
      for (CGAL_AOS3_TYPENAME Degeneracy::Points::iterator
	     it = degeneracies[i].points_.begin();
	   it != degeneracies[i].points_.end(); ++it) {
	it->second=find_box(degeneracies[i].location_,
	  it->first);
      }
    }


    NT lb=to_interval(event_locations.front()).first;
    NT ub=to_interval(event_locations.front()).second+0.0009765625;
    NT mp;
    while (true) {
      bool ok=true;
      do {
	mp=(ub+lb)/2;
	if (event_locations.front() >= mp) {
	  lb=mp;
	} else break;
      } while (true);
      rcs_.set_z(mp);
      for (unsigned int i=0; i< degeneracies.size(); ++i) {
	for (CGAL_AOS3_TYPENAME Degeneracy::Vertices::iterator
	       it = degeneracies[i].vertices_.begin();
	     it != degeneracies[i].vertices_.end(); ++it) {
	  CGAL_assertion(0);
	  ok= box_ok(it->second, 
		     it->first->point(), 
		     mp);
	  if (!ok) break;
	}
	if (!ok) break;
      }
      if (ok && cs_.visitor().simulator()->next_event_time() > mp) break;
      else {
	ub=mp;
      }
    }
    
    CGAL_LOG(Log::SOME, "Time is " << mp << std::endl);
   


    for (unsigned int i=0; i< degeneracies.size(); ++i) {
      // build arrangement with rules
      //  -- remove any edge with two non-degen vertices
      //  -- split any vertex which is non-degen and has more than one edge
      //  -- create anchors
      // build list of outgoing rules for each vertex
      // add any needed outgoing rules to vertex
      //  -- no problem with degeneracies :-)
    }
    
    for (unsigned int i=0; i< degeneracies.size(); ++i) {
      // delete extra rules
      // figure out which rules are there, delete extra vertices
    }
    for (unsigned int i=0; i< degeneracies.size(); ++i) {
      // stitch in graphs 
      //    just match an outward directed curve
    }
  }
private:
  Traits tr_;
  CCS &cs_;
  RCS rcs_;
  Rule_list sorted_lists_[2];
};

CGAL_AOS3_END_INTERNAL_NAMESPACE

#ifdef CGAL_AOS3_USE_TEMPLATES
#include "Degenerate_cross_section_impl.h"
#endif
#endif
