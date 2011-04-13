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
// file          : include/CGAL/Pm_dynamic_open_bounding_box.h
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

#ifndef CGAL_PM_DYNAMIC_OPEN_BOUNDING_BOX_H
#define CGAL_PM_DYNAMIC_OPEN_BOUNDING_BOX_H

#ifndef CGAL_PLANAR_MAP_2_H
#include <CGAL/Planar_map_2.h>
#endif

#ifndef CGAL_PM_BOUNDING_BOX_BASE_H
#include <CGAL/Pm_bounding_box_base.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Planar_map_>
class Pm_dynamic_open_bounding_box : public Pm_bounding_box_base<Planar_map_> {
public:
  typedef Planar_map_ Planar_map;
  typedef Pm_bounding_box_base<Planar_map_> Base;
  typedef Pm_dynamic_open_bounding_box<Planar_map> Self;
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
  typedef typename Planar_map::Locate_type Locate_type;
  typedef std::vector<Halfedge_handle> Halfedge_container;

  typedef typename Traits::X_curve X_curve;
  typedef typename Traits::Point Point;
  typedef typename Traits::Ray Ray;
  typedef typename Traits::Direction Direction;
  typedef typename Traits::Bounding_box Bounding_box;
  //  typedef typename Traits::Boundary_type Boundary_type;
  typedef typename Traits::Point_container Point_container;
  typedef typename Point_container::iterator Point_iterator;
  typedef typename Traits::X_curve_container X_curve_container;
  typedef typename X_curve_container::iterator X_curve_iterator;
  typedef typename Planar_map::Point_location_base Point_location_base;
  /*
    typedef std::list<X_curve> X_curve_container;
    typedef Topological_map<_Dcel> TPM;
  */
  typedef typename Planar_map::Halfedge_around_vertex_circulator 
     Halfedge_around_vertex_circulator;
  typedef typename Planar_map::Holes_iterator Holes_iterator;
  typedef typename Planar_map::Holes_const_iterator Holes_const_iterator;
  typedef typename Planar_map::Ccb_halfedge_const_circulator 
    Ccb_halfedge_const_circulator;
  typedef typename Planar_map::Ccb_halfedge_circulator Ccb_halfedge_circulator;
	
  struct Token : public Base::Token {
  // The following class is used to synchronize between two instances of 
  // Bounding_box and Point_location strategies, during 
  // rebuild of the bounding box in the first and an update call in the second.

    typedef typename Traits::Bounding_box Bounding_box;

    Token(const Bounding_box& t) : b(t) {};

    virtual void rebuild_bounding_box(const Point_location_base* p) const 
    {((Traits*)p->get_traits())->set_bounding_box(b);}
  private:
    const Bounding_box& b;
  };

  Pm_dynamic_open_bounding_box(){}
  ~Pm_dynamic_open_bounding_box(){}
 
  void init(Planar_map& pmp, Traits& tr) {
    pm = &pmp;
    traits = (Traits_wrap*)(&tr);
  }
  
  bool insert(const Point& p) {
#ifdef CGAL_PM_DEBUG
    CGAL_assertion(debug_invariant());
#endif

    // Returns true if bounding box remained unchanged.
    Bounding_box bold = traits->get_bounding_box(), 
    b = is_empty(bold) ? traits->get_bounding_box(p) : 
                         traits->get_bounding_box(p,bold);
    if (!traits->bounding_box_is_same(b))
      {
        rebuild_bounding_box(b);
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
    
#ifdef CGAL_PM_DEBUG
    CGAL_assertion(debug_invariant());
#endif

    Point_iterator it=begin;
    Bounding_box bold = traits->get_bounding_box(),b=bold;
    if (it!=end)
      {
        b = is_empty(b) ? traits->get_bounding_box(*it) : 
                          traits->get_bounding_box(*it,b);
        ++it;
      }
    while(it!=end)
      {
#ifdef CGAL_PMBB_DEBUG
        std::cout << "\nPm_dynamic_open_bounding_box<Planar_map>::insert(" 
                  << "begin,end)" << std::endl;
#endif
        
#ifdef CGAL_PM_DEBUG
        CGAL_assertion(debug_invariant());
#endif
        b=traits->get_bounding_box(*it,b);
        ++it;
      }
    if (!traits->bounding_box_is_same(b))
      {
#ifdef CGAL_PM_DEBUG
        CGAL_assertion(debug_invariant());
#endif
        rebuild_bounding_box(b);
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

  Bounding_box 
    bold = traits->get_bounding_box(), 
    b = is_empty(bold) ? traits->get_bounding_box(cv) : 
                         traits->get_bounding_box(cv,bold);

    if (!traits->bounding_box_is_same(b))
    {
      rebuild_bounding_box(b);
      return false;
    }
  return true;
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

    //    Boundary_type bt;
    //		Point p1, p2;
    //		X_curve c;
    X_curve_iterator it=begin;
    Bounding_box bold = traits->get_bounding_box(),b=bold;
    if (it!=end)
      {
        b = is_empty(b) ? traits->get_bounding_box(*it) : 
                             traits->get_bounding_box(*it,b);
        ++it;
      }
    while(it!=end)
      {
        b = traits->get_bounding_box(*it,b);
        ++it;
      }
    if (!traits->bounding_box_is_same(b))
      {
        rebuild_bounding_box(b);
        return false;
      }
    return true;
  }
  
  bool insert(const X_curve& cv,const Ray& r) { 
/* returns true if the bounding box remains unchanged, false otherwise. */

#ifdef CGAL_PM_DEBUG
    CGAL_assertion(debug_invariant());
#endif

    const Bounding_box& bold = traits->get_bounding_box(), 
      b = is_empty(bold) ? traits->get_bounding_box(cv,r) : 
                           traits->get_bounding_box(cv,r,bold);
#ifdef CGAL_PMBB_DEBUG
    CGAL_postcondition(b==traits->get_bounding_box(cv,r,b));
#endif
    // makes sure that the intersection of the ray and the curve 
    // (if such intersection exists) is inside the bounding box.
    if (!traits->bounding_box_is_same(b))
      {
        rebuild_bounding_box(b);
        return false;
      }
    return true;
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
              return false;
            else {++it;++it;}
          }
        // make sure that the intersection of the curve with the vertical ray 
	// is inside the bounding box
        // (if such intersection exists).
      }
    else if (lt==Planar_map::VERTEX)
      {
        if (!insert(h->target()->point()))
          return false;
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
  virtual void split_boundary_edge(const Halfedge_handle &h,
                                   Halfedge_handle h1,
                                   Halfedge_handle h2,
                                   const Point& p) {}
  
  void merge_edge(const X_curve &cv1,
                  const X_curve &cv2,
                  Halfedge_handle e
                  //additions by iddo for arrangement
                  ,const X_curve& cv
                  //end additions
                  ) {}
  
  void remove_edge(Halfedge_handle e) {}
  void clear()
  {
    traits->set_bounding_box();
  }
  inline bool is_empty() const {
    return is_empty(traits->get_bounding_box());
  }
 protected:
  inline bool is_empty(const Bounding_box& b) const {
    return b.identical(Traits::unbounded_box());
  }
 public:
  
#ifdef CGAL_PM_DEBUG
  void debug() const {}
  bool debug_invariant() const{
    CGAL_assertion(traits);
    CGAL_assertion(&traits->get_bounding_box());
    return true;
  }
#endif
  
protected:	
  void apply_bounding_box(const Face_handle& unbounded, 
                          const Halfedge_handle& he,Locate_type& lt) const
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
              if (circ->face()==unbounded) lt=Planar_map::UNBOUNDED_FACE; 
              break;
            }
        }
        break;
      case Planar_map::EDGE:
        if (he->face()==unbounded ||
            he->twin()->face()==unbounded) lt=Planar_map::UNBOUNDED_EDGE;
        break;
      case Planar_map::VERTEX:
        {
          Halfedge_around_vertex_circulator 
            begin=he->target()->incident_halfedges(),
            circ=begin;
          while (++circ!=begin) 
            { 
              if (circ->face()==unbounded) lt=Planar_map::UNBOUNDED_VERTEX; 
              break;}
        }
        break;
      }
  }
  bool cooriented(const Ccb_halfedge_circulator& h,const X_curve& cv) const
  {
#ifdef CGAL_PMBB_DEBUG
    if (!traits->curve_is_same(h->curve(),cv))
      std::cerr << "\nbool cooriented(h," << cv << "); h->curve()=" 
                << h->curve() << std::endl;
#endif
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
#ifdef CGAL_PMBB_DEBUG
        if (
            traits->curve_get_status(cv)!=Traits::CURVE_RIGHT &&
            traits->curve_get_status(cv)!=Traits::CURVE_LEFT)
        {
          std::cout << "\nCGAL_precondition(cv status is RIGHT||LEFT)" 
                    << "cv=" << cv << ",traits->curve_get_status(cv)="
                    << traits->curve_get_status(cv);
        }
#endif
        CGAL_precondition(
          traits->curve_get_status(cv)==Traits::CURVE_RIGHT ||
          traits->curve_get_status(cv)==Traits::CURVE_LEFT);
        return (traits->point_is_right(t,s)==
                (traits->curve_get_status(cv)==Traits::CURVE_RIGHT));
      }
  }
  // discards the const 
  Point_location_base* get_point_location() { 
    return (Point_location_base*)pm->get_point_location();	
  }

  void rebuild_bounding_box(const Bounding_box& b){
    /* The main function in the bounding box. When called, the planar map 
       representation is rebuild.
    Precondition:
       The bounding box contains in its interior all the interesting points.
    */

#ifdef CGAL_PMBB_DEBUG
    std::cout << "\nPm_dynamic_open_bounding_box::rebuild_bounding_box()";
#endif    

    Halfedge_container c;
    Holes_iterator hit=pm->unbounded_face()->holes_begin(),
      hoe=pm->unbounded_face()->holes_end();
    while (hit!=hoe) {
      Ccb_halfedge_circulator begin = *hit,
        circ=begin;
      do {
        if (traits->curve_is_unbounded(circ->curve()))
          {
            // Each halfedge contributes an unbounded curve, due to the 
            // dynamic nature of the bounding box and the invariant that
            // the curves are all bounded inside.
            // We optimize here taking only one representative for each edge.
            
            Halfedge_handle second=circ;
            ++second;
            if (second==circ->twin()) c.push_back(circ);
          }
        ++circ;
      }
      while(circ!=begin);
      ++hit;
    }
    typename Halfedge_container::iterator it=c.begin(),cend=c.end();
		get_point_location()->update(it,cend,Token(b));

#ifdef CGAL_PMBB_DEBUG
    CGAL_assertion(traits->get_bounding_box().identical(b));
#endif

    while (it!=cend)
      {
        Halfedge_handle h=*it;
        const X_curve& cv=h->curve();
        Vertex_handle s=h->source(),t=h->target(); 
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
        if (!(traits->point_is_same(
            traits->curve_source(h->curve()),
            h->source()->point()) &&
          traits->point_is_same(
            traits->curve_target(h->curve()),
            h->target()->point()) ||
          traits->point_is_same(
            traits->curve_source(h->curve()),
            h->target()->point()) &&
          traits->point_is_same(
                                traits->curve_target(h->curve()),
				h->source()->point())))
          {
            std::cout << "\ncurve_source(h->curve())=" 
		      << traits->curve_source(h->curve());
            std::cout << " curve_target(h->curve())=" 
		      << traits->curve_target(h->curve());
            std::cout << " h->source()->point()=" << h->source()->point();
            std::cout << " h->target()->point()=" << h->target()->point();
          }
        CGAL_postcondition(
          traits->point_is_same(
            traits->curve_source(h->curve()),
            h->source()->point()) &&
          traits->point_is_same(
            traits->curve_target(h->curve()),
            h->target()->point()) ||
          traits->point_is_same(
            traits->curve_source(h->curve()),
            h->target()->point()) &&
          traits->point_is_same(
            traits->curve_target(h->curve()),h->source()->point()) );
#endif
        ++it;
      }
    c.clear();
  }
  
protected:
  Planar_map* pm;
  Traits_wrap* traits;
  
};

CGAL_END_NAMESPACE

#endif //CGAL_PM_DYNAMIC_OPEN_BOUNDING_BOX_H
/*
  Special situations:
  |--------
  |  x	|
  |  |	|
  |  |	|
  |-\ /---|
    vrs
*/
