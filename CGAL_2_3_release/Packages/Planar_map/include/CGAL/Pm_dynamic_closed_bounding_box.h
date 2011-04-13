// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 13
//
// file          : include/CGAL/Pm_dynamic_closed_bounding_box.h
// package       : pm (4.20)
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//                 
//
// coordinator   : Tel-Aviv University (Dan Halperin)
//
// Chapter       : 
// email         : cgal@cs.uu.nl
//
// ======================================================================

#ifndef CGAL_PM_DYNAMIC_CLOSED_BOUNDING_BOX_H
#define CGAL_PM_DYNAMIC_CLOSED_BOUNDING_BOX_H

#ifndef CGAL_PM_BOUNDING_BOX_BASE_H
#include <CGAL/Pm_bounding_box_base.h>
#endif

#ifndef CGAL_CIRCULATOR_H
#include <CGAL/circulator.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Planar_map_>
class Pm_dynamic_closed_bounding_box : 
  public Pm_bounding_box_base<Planar_map_> {
public:
  typedef Planar_map_ Planar_map;
  typedef Pm_dynamic_closed_bounding_box<Planar_map> Self;
  typedef typename Planar_map::Traits Traits;
  typedef typename Planar_map::Traits_wrap Traits_wrap;
  typedef typename Planar_map::Halfedge_handle Halfedge_handle;
  typedef typename Planar_map::Face_handle Face_handle;
  typedef typename Planar_map::Vertex_handle Vertex_handle;
  typedef typename Planar_map::Vertex_const_handle Vertex_const_handle;
  typedef typename Planar_map::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Planar_map::Face_const_handle Face_const_handle;
  typedef typename Planar_map::Vertex_iterator Vertex_iterator;
  typedef typename Planar_map::Halfedge_iterator Halfedge_iterator;
  typedef typename Planar_map::Face_iterator Face_iterator;
  typedef typename Planar_map::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Planar_map::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Planar_map::Face_const_iterator Face_const_iterator;
  typedef std::vector<Halfedge_handle> Halfedge_container;

  typedef typename Planar_map::Locate_type Locate_type;
  typedef typename Traits::X_curve X_curve;
  typedef typename Traits::Point Point;
  typedef typename Traits::Ray Ray;
  typedef typename Traits::Direction Direction;
  typedef typename Traits::Bounding_box Bounding_box;
  typedef typename Traits::Boundary_type Boundary_type;
  typedef typename Traits::Point_container Point_container;
  typedef typename Traits::X_curve_container X_curve_container;
  typedef typename Point_container::iterator Point_iterator;
  typedef typename X_curve_container::iterator X_curve_iterator;
  typedef typename Planar_map::Point_location_base Point_location_base;
  /*
    //  typedef std::list<X_curve> X_curve_container;
    //  typedef Topological_map<_Dcel> TPM;
  */
  typedef typename Planar_map::Halfedge_around_vertex_circulator 
                                             Halfedge_around_vertex_circulator;
  typedef typename Planar_map::Holes_iterator                   Holes_iterator;
  typedef typename Planar_map::Holes_const_iterator       Holes_const_iterator;
  typedef typename Planar_map::Ccb_halfedge_const_circulator 
                                                 Ccb_halfedge_const_circulator;
  typedef typename Planar_map::Ccb_halfedge_circulator Ccb_halfedge_circulator;
  /*
    typedef typename Base::Size Size;
  */
  
  Pm_dynamic_closed_bounding_box() : intersection_on_boundary(false){}
  ~Pm_dynamic_closed_bounding_box(){}
  
  
  void init(Planar_map& pmp, Traits& tr) {
    pm = &pmp;
    traits = (Traits_wrap*)(&tr);
  }
  bool insert(const Point& p) {
    // Returns true if bounding box remained unchanged.
    Bounding_box bold = traits->get_bounding_box(), 
                 b    = traits->increase_bounding_box(p);
    if (b!=bold)
      {
        Traits* PL_traits=(Traits*)get_point_location()->get_traits();
        if (traits!=PL_traits) PL_traits->set_bounding_box(b);
        rebuild_bounding_box();
        return false;
      }
    return true;
  }
  bool insert(const Point_iterator& begin,const Point_iterator& end
#ifdef _MSC_VER
              ,Point* dummy=0
#endif
              )
  { // workaround for MSVC6.0
    // Returns true if bounding box remained unchanged.
    
    Point_iterator it=begin;
    Bounding_box bold = traits->get_bounding_box();
    while(it!=end)
      {
#ifdef CGAL_PM_DEBUG
        CGAL_assertion(debug_invariant());
#endif
        traits->increase_bounding_box(*it);
        ++it;
      }
    Bounding_box b=traits->get_bounding_box();
    if (b!=bold)
      {
        Traits* PL_traits=(Traits*)get_point_location()->get_traits();
        if (traits!=PL_traits) PL_traits->set_bounding_box(b);
#ifdef CGAL_PM_DEBUG
        CGAL_assertion(debug_invariant());
#endif
        rebuild_bounding_box();
        return false;
      }
    return true;
  }
  bool insert(const X_curve& cv) {
    // we check that the finite endpoints of the curve is in the BBox
    // and if not enlarge the BBox to contain them.
    // if the curve is infinite (in both sides) we have to check that
    // it passes through the BBox.
    // Last, we change the "infinite" endpoints to lie on the BBox.
    
    Point p1, p2;
    X_curve c;
    
#ifdef CGAL_PM_DEBUG
    CGAL_assertion(debug_invariant());
#endif
    
    Bounding_box bold = traits->get_bounding_box(), 
                 b    = traits->increase_bounding_box(cv);
    
    if (b!=bold)
      {
        Traits* PL_traits=(Traits*)get_point_location()->get_traits();
        if (traits!=PL_traits) PL_traits->set_bounding_box(b);
        //optimize here ....
        rebuild_bounding_box();
        return false;
      }
    return true;
    
    /*
      while (!traits.is_curve_inside_bounding_box(cv,bt))
      {
      switch (bt)
      {
      case Traits::TOP:    
      find_adjacent_points(top_left, p1, p2);
      c = traits.create_segment(p1,p2);
      delete_vertex(top_left, c);
      break;
      case Traits::RIGHT: 
      find_adjacent_points(top_right, p1, p2);
      c = traits.create_segment(p1,p2);
      // debug
      ppp = top_right->point();
      //end debug
      delete_vertex(top_right, c);
      
      break;
      case Traits::BOTTOM: 
      find_adjacent_points(bottom_right, p1, p2);
      c = traits.create_segment(p1,p2);
      delete_vertex(bottom_right, c);
      break;
      
      case Traits::LEFT:
      find_adjacent_points(bottom_left, p1, p2);
      c = traits.create_segment(p1,p2);
      delete_vertex(bottom_left, c);
      break;
      }
      
      reinsert_vertex_of_BBox(bt);
      
      switch (bt)
      {
      case Traits::TOP:    
      find_adjacent_points(top_right, p1, p2);
      c = traits.create_segment(p1,p2);
      // debug
      ppp = top_right->point();
      //end debug
      delete_vertex(top_right, c);
      
      reinsert_vertex_of_BBox(Traits::RIGHT);
      break;
      case Traits::RIGHT:
      find_adjacent_points(bottom_right, p1, p2);
      c = traits.create_segment(p1,p2);
      // debug
      ppp = bottom_right->point();
      //end debug
      delete_vertex(bottom_right, c);
      
      reinsert_vertex_of_BBox(Traits::BOTTOM);
      break;
      case Traits::BOTTOM: 
      find_adjacent_points(bottom_left, p1, p2);
      c = traits.create_segment(p1,p2);
      // debug
      ppp = bottom_left->point();
      //end debug
      delete_vertex(bottom_left, c);
      
      reinsert_vertex_of_BBox(Traits::LEFT);
      break;
      case Traits::LEFT:   
      find_adjacent_points(top_left, p1, p2);
      c = traits.create_segment(p1,p2);
      // debug
      ppp = top_left->point();
      //end debug
      delete_vertex(top_left, c);
      
      reinsert_vertex_of_BBox(Traits::TOP);
      break;
      }
      }
    */
  }
  
  bool insert(const X_curve_iterator& begin,const X_curve_iterator& end
#ifdef _MSC_VER
              ,X_curve* dummy=0
#endif
              )
  { // workaround for MSVC6.0
    
#ifdef CGAL_PM_DEBUG
    CGAL_assertion(debug_invariant());
#endif
    
    X_curve_iterator it=begin;
    Bounding_box bold = traits->get_bounding_box();
    while(it!=end)
      {
        traits->increase_bounding_box(*it);
        ++it;
      }
    Bounding_box b = traits->get_bounding_box();
    if (b!=bold)
      {
        Traits* PL_traits=(Traits*)get_point_location()->get_traits();
        if (traits!=PL_traits) PL_traits->set_bounding_box(b);
        rebuild_bounding_box();
        return false;
      }
    return true;
  }	
  bool insert(const X_curve& cv,const Ray& r) { 
    /* 
       returns true if the bounding box remains unchanged, false otherwise. 
    */
    const Bounding_box before = traits->get_bounding_box();
    const Bounding_box& after=traits->increase_bounding_box(cv,r,before);
    CGAL_postcondition(after==traits->increase_bounding_box(cv,r,after));
    // makes sure that the intersection of the ray and the curve 
    // (if such intersection exists) is inside the bounding box.
    return before==after;
  }
  
  /* The point location query function may updates the resulting 
     halfedge handle and locate type as expected from the bounding box */
  /* Returns true if bounding box remained unchanged */
  bool locate(const Point& p, Locate_type& lt,Halfedge_handle& h){
    //		apply_bounding_box(pm->unbounded_face(),h,lt);
    return true;
  };
  
  /* Returns true if bounding box remained unchanged */
  bool vertical_ray_shoot(const Point& p, Locate_type& lt, bool up,
                          Halfedge_handle& h){
    CGAL_precondition(lt==Planar_map::VERTEX||lt==Planar_map::EDGE||
                      lt==Planar_map::UNBOUNDED_FACE);
    if (lt==Planar_map::UNBOUNDED_FACE)
      {
        Halfedge_iterator it=pm->halfedges_begin();
        Ray r(p,Direction(0,up ? 1 : -1));
        while(it!=pm->halfedges_end()) 
          {
            const X_curve& cv=it->curve();
            if (!traits->curve_is_in_x_range(cv,p) && !insert(cv,r))
              {
                rebuild_bounding_box();
                return false;
              }
            else {++it;++it;}
          }
        // make sure that the intersection of the curve with the vertical ray 
	// is inside the bounding box
        // (if such intersection exists).
      }
    else if (lt==Planar_map::VERTEX)
      {
        if (!insert(h->target()->point()))
          {
            rebuild_bounding_box();
            return false;
          }
      }
    return true;
    //		apply_bounding_box(pm->unbounded_face(),h,lt);
  };
  
  void split_edge(const X_curve &cv,
                  Halfedge_handle e1,
                  Halfedge_handle e2
                  //additions by iddo for arrangement
                  ,const X_curve& cv1, const X_curve& cv2
                  //end additions
                  ) {}
  void split_boundary_edge(const Halfedge_handle &h,
                           Halfedge_handle h1,
                           Halfedge_handle h2,
                           const Point& p) 
  {
    X_curve cv1,cv2;
    traits->split_curve(h->curve(),cv1,cv2,p);
    h1 = pm->split_edge(h,cv1,cv2);
    h2 = h1->next_halfedge();
#ifdef CGAL_PM_DEBUG
    CGAL_postcondition(h1->target()==h2->source());
    CGAL_postcondition(h1->source()==h->source());
    CGAL_postcondition(h2->target()==h->target());
#endif
  }
  
  void merge_edge(const X_curve &cv1,
                  const X_curve &cv2,
                  Halfedge_handle e
                  //additions by iddo for arrangement
                  ,const X_curve& cv
                  //end additions
                  ) {}
  
  void remove_edge(Halfedge_handle e) {}
  void clear(){traits->set_bounding_box();}
  inline bool is_empty() const 
  {   
    return is_empty(traits->get_bounding_box());
  }
protected:
  inline bool is_empty(const Bounding_box& b) const {
    return b.identical(Traits::unbounded_box());
  }
  
protected:	
  void apply_bounding_box(
                          const Face_handle& unbounded, 
                          const Halfedge_handle& he,Locate_type& lt) const
    // Reevaluates the locate type according to the bounding box.
  {
    switch(lt)
      {
      case Planar_map::UNBOUNDED_FACE:
        break;
      case Planar_map::FACE:
        {
          Ccb_halfedge_const_circulator begin(he->twin());
          Ccb_halfedge_const_circulator circ=begin;
          while (++circ!=begin) 
            { 
	      if (circ->face()==unbounded) 
		lt=Planar_map::UNBOUNDED_FACE; 
	      break;
	    }
        }
        break;
      case Planar_map::EDGE:
        if (he->twin()->face()==unbounded ||
            he->twin()->face()==unbounded) lt=Planar_map::UNBOUNDED_EDGE;
        break;
      case Planar_map::VERTEX:
        {
          Halfedge_around_vertex_circulator 
            begin=he->target()->incident_halfedges(),
            circ=begin;
          while (++circ!=begin) 
            { 
	      if (circ->face()==unbounded) 
		lt=Planar_map::UNBOUNDED_VERTEX; 
	      break;
	    }
        }
        break;
      }
  }
  
  //bool cooriented(const Ccb_halfedge_circulator& h,const X_curve& cv) const
  bool cooriented(const Halfedge_const_handle& h,const X_curve& cv) const
    // Returns weather the halfedge and the curve are cooriented or not.
  {
    CGAL_precondition(traits->curve_is_same(h->curve(),cv));
    const Point& s=h->source()->point(),&t=h->target()->point();
    if (traits->point_is_same_x(s,t))
      {
        CGAL_precondition( 
	  traits->curve_get_status(cv)==Traits::CURVE_VERTICAL_UP || 
	  traits->curve_get_status(cv)==Traits::CURVE_VERTICAL_DOWN);
        return (traits->point_is_higher(t,s) ==
		(traits->curve_get_status(cv)==Traits::CURVE_VERTICAL_UP));
      }
    else
      {
        CGAL_precondition(
          traits->curve_get_status(cv)==Traits::CURVE_RIGHT ||
	  traits->curve_get_status(cv)==Traits::CURVE_LEFT);
        return (traits->point_is_right(t,s) ==
		(traits->curve_get_status(cv)==Traits::CURVE_RIGHT));
      }
  }
  
  Point_location_base* get_point_location() const { 
    // Returns the point location while discarding the const. 
    return (Point_location_base*)pm->get_point_location();	
  }
  
  bool target_is_intersection_point(const Halfedge_handle& h,Point& p) const
    // Returns weather the halfedge's target is a real intersection point 
    // on the boundary.
    // In that case, p is set to that point.
  {
    Halfedge_around_vertex_circulator avc_b=h;
    Halfedge_around_vertex_circulator avc_ntnt=avc_b;
    ++avc_ntnt;
    p = traits->curve_target(avc_ntnt->curve());
#ifdef CGAL_PM_DEBUG
    Halfedge_around_vertex_circulator avc_nt=avc_ntnt;
    CGAL_assertion(avc_nt!=avc_b); 
    // Isolated corners are not allowed.
#endif
    ++avc_ntnt;
    return avc_ntnt!=avc_b; 
    //	Found an intersection point on the boundary, i.e. a point 
  }
  
  void find_reference_point()
    /* Finds an intersection point of the boundary with the 
       Planar map curves, if such exists. */
  {
    CGAL_precondition(
                      pm->unbounded_face()->holes_begin()!=
                      pm->unbounded_face()->holes_end());
    
#ifdef CGAL_PM_DEBUG
    Holes_iterator 
      ho_nit=++pm->unbounded_face()->holes_begin(),
      ho_e=pm->unbounded_face()->holes_end();
    CGAL_precondition(ho_nit==ho_e); 
    // The unbounded face of the Planar map should have at most one hole.
#endif 
    Halfedge_container c;
    Holes_iterator it=pm->unbounded_face()->holes_begin();
    Ccb_halfedge_circulator
      ccb_b = *ccb_hi;
    halfedge_on_boundary=ccb_b;
    Point p;
    do
      {
        if (target_is_intersection_point(halfedge_on_boundary,p))
          {
            intersection_on_boundary=true;
            return; 
          }
        //	Found an intersection point on the boundary, i.e. a point 
        //	that belongs to one of the infinite curves in the planar map.
      }
    while(halfedge_on_boundary++!=ccb_b);
    intersection_on_boundary=false;
  }
  
  
  void rebuild_bounding_box(){
    if (initialized)
      // bounding box is already initialized.
      {
        // Reevaluates the planar map, by recalculating the bounding box.
        find_reference_point();
        // move this to modifying functions.
#ifdef CGAL_PM_DEBUG
        if (intersection_on_boundary)
          {
            Halfedge_around_vertex_circulator avc_b=halfedge_on_boundary;
            Halfedge_around_vertex_circulator avc_ntnt=avc_b;
            ++avc_ntnt;
            ++avc_ntnt;
            CGAL_assertion(avc_ntnt!=avc_b); 
				// Invariant: intersection_on_boundary -> 
				// halfedge_on_boundary->target() is one.
          }
#endif
        
        X_curve_container new_boundary;
        traits->get_x_curve_boundary(new_boundary);
        // Get new boundary from traits.
        
        typedef typename X_curve_container::reverse_iterator 
          reverse_iterator;
        // The original orientation of the boundary is ccw, while we expect cw.
        reverse_iterator new_boundary_rbegin=new_boundary.rbegin();
        /*
          reverse_iterator new_boundary_real_rend=--(new_boundary.rend());
        */
        
        CGAL_postcondition(new_boundary.rbegin()!=new_boundary.rend());
        
        typedef Forward_circulator_from_iterator<reverse_iterator>
          reverse_circulator;
        // Making a circulator from the container. (simpler to work with).
        /*
          reverse_circulator 
          boundary_circulator(new_boundary_rbegin,new_boundary_real_rend);
        */
        // note that the second iterator is not new_boundary.rend();
        
        reverse_circulator::iterator 
          boundary_begin=new_boundary_rbegin,
          boundary_circ=boundary_begin;
        
        //		X_curve curr_cv;
        //		Point curr_endp;
        // Notify the point location about the scheduled update.
        Ccb_halfedge_circulator begin=halfedge_on_boundary,circ=begin;
        Point last_end,curr_endp;
        // last and current points along new boundary.
        
        get_point_location()->update(halfedge_on_boundary);
        
        if (intersection_on_boundary)
	// Synchronize between the new boundary and the old boundary.
          {
            while (
                   !traits->is_point_on_curve(
                      *boundary_circ,
		      halfedge_on_boundary->target()->point())
                   ) boundary_circ++;
	    // new_boundary_it contains an intersection point.
            while(
                  target_is_intersection_point(circ,curr_endp) &&
                  traits->is_point_on_curve(*boundary_circ,curr_endp)
                  )
              {
                halfedge_on_boundary=circ;
                ++circ;
              }
	    // Synchronize boundary curve iterator and halfedge circulator.
	    // Curve turns counter clockwise, while halfedge turns clockwise.
            
            X_curve curr_cv=*boundary_circ;
	    // last curves along old and new boundaries.
            bool last_was_intersecting=false;
	    // last edge's target had an intersection point. 
            Halfedge_around_vertex_circulator inter_edge=circ;
	    // auxiliary halfedge for the intersection point calculation.
            ++inter_edge;
            curr_endp=traits->curve_target(inter_edge->curve());
            ++circ; // circ is oriented clockwise.
            do{
              if (target_is_intersection_point(circ,curr_endp)) 
                {
                  inter_edge=circ;
                  ++inter_edge;
                  curr_endp=traits->curve_target(inter_edge->curve());
                  while (!traits->is_point_on_curve(*boundary_circ,curr_endp))
                    {
                      circ->set_curve(curr_cv);
                      circ->twin()->set_curve(curr_cv);
                      last_was_intersecting=false;
                      ++boundary_circ;
                    }
                  // update the curve according to the newly split curve.
                  {
                    X_curve slower,faster;
                    traits->split_curve(curr_cv,slower,faster,curr_endp);
                    circ->set_curve(slower);
                    circ->twin()->set_curve(slower);
                    curr_cv=faster;
                    last_was_intersecting=true;
                  }
                }
              ++circ; // circ is oriented clockwise.
              
	      // don't forget to remove the old halfedges!!!
              // optimization: need not remove all corver halfedges, 
	      // only in case of a lack, or an over use.
              
            } while (circ!=begin);
            
/*
  const X_curve& cv=vnext->curve();
  Halfedge_handle h=vnext;
  Vertex_handle 
  s=h->source(),
  t=h->target(); 
  // optimization: never handle both h and h->twin().
  if (h->face()==h->twin()->face() && 
  traits->point_is_left_low(s->point(),t->point())) 
  {
  #ifdef CGAL_PM_DEBUG
  CGAL_precondition( 
  traits->point_is_same(
  traits->curve_source(h->curve()),h->source()->point()) &&
  traits->point_is_same(
  traits->curve_target(h->curve()),h->target()->point()) ||
  traits->point_is_same(
  traits->curve_source(h->curve()),h->target()->point()) &&
  traits->point_is_same(
  traits->curve_target(h->curve()),h->source()->point()) );
  #endif
  --circ;
  // turn counter clockwise around face.
  continue;
  }
  if (cooriented(h,cv))
  {
  t->set_point(traits->curve_target(cv));
  s->set_point(traits->curve_source(cv));
  }
  else
  {
  t->set_point(traits->curve_source(cv));
  s->set_point(traits->curve_target(cv));
  }
  CGAL_postcondition(!traits->point_is_same(
  h->source()->point(),h->target()->point()));
  #ifdef CGAL_PM_DEBUG
  CGAL_precondition( 
  traits->point_is_same(
  traits->curve_source(h->curve()),h->source()->point()) &&
  traits->point_is_same(
  traits->curve_target(h->curve()),h->target()->point()) ||
  traits->point_is_same(
  traits->curve_source(h->curve()),h->target()->point()) &&
  traits->point_is_same(
  traits->curve_target(h->curve()),h->source()->point()) );
  #endif
  }
  else
  {
  
  }
  ++circ;
  } while (circ!=begin);
  //		++it;
*/
	    
          }
        else
				// Recreate the whole boundary.
          {
            do
              {
                halfedge_on_boundary->set_curve(*boundary_circ);
                halfedge_on_boundary->twin()->set_curve(*boundary_circ);
                ++boundary_circ;
                ++halfedge_on_boundary;
                // optimization. In general the number of curves in the 
                // boundary should not be a constant.
              }
            while(boundary_circ!=boundary_begin);
          }
      }
    else
      // bounding box is empty.
      {	
        X_curve_container cv_aux; 
        traits->get_x_curve_boundary(cv_aux);
        // returns the curve boundary of the bounding box.
        // e.g. the 4 (infinite) segments of the bounding box
        //		in this order: bottom, right, top, left.
        
        const typename X_curve_container::iterator 
          &bcv = cv_aux.begin(), &ecv = cv_aux.end();
        typename X_curve_container::iterator 
          currcv=bcv, nextcv=bcv;
        ++nextcv;
        const Halfedge_handle& first = 
          pm->insert_in_face_interior(*bcv, pm->unbounded_face());
#ifdef CGAL_PM_DEBUG
        CGAL_postcondition(cooriented(first,*bcv));
#endif
        Halfedge_handle curre=first; 
        currcv=nextcv;
        ++nextcv;
        while(nextcv!=ecv)
          {
            curre = pm->insert_from_vertex(*currcv, curre->target(), true);
            currcv=nextcv;
            ++nextcv;
          }
        curre = pm->insert_at_vertices(*currcv, curre->target(), 
				       first->source());
        
        CGAL_assertion(curre->next_halfedge()==first);
      }
  }
  
#ifdef CGAL_PM_DEBUG
  void debug() const {}
  bool debug_invariant() const{
    CGAL_assertion(traits);
    CGAL_assertion(&traits->get_bounding_box());
    return true;
  }
#endif
	
protected:
	Planar_map* pm;
	Traits_wrap* traits;
	bool initialized;
	Ccb_halfedge_circulator halfedge_on_boundary;
	// Halfedge on the outer boundary that points at an 
	// intersection point on the boundary, if such exists.
	bool intersection_on_boundary;
};

CGAL_END_NAMESPACE

#endif //CGAL_PM_DYNAMIC_CLOSED_BOUNDING_BOX_H
/*
Special situations:
|--------
|  x	|
|  |	|
|  |	|
|-\ /---|

 vrs
 */
 
 /*  Souldn't 		
	Planar_map::Halfedge_handle insert_from_vertex(const X_curve& cv, 
	Vertex_handle v1, bool source) 
	be
	Planar_map::Halfedge_handle insert_from_vertex(const X_curve& cv, 
	Vertex_const_handle v1, bool source) 
	
	 same for ...
	 Halfedge_handle insert_at_vertices(const X_curve& cv, 
	 Vertex_handle v1, Vertex_handle v2) 
	 
	  */
	  
	  /* optimization: It is possible to keep a pointer to a bounding box,
             instead of a virtual base class. */
