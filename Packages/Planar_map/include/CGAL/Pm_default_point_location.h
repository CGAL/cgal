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
// file          : include/CGAL/Pm_default_point_location.h
// package       : pm (4.08)
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//                 Iddo Hanniel <hanniel@math.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_PM_DEFAULT_POINT_LOCATION_H
#define CGAL_PM_DEFAULT_POINT_LOCATION_H

#ifndef CGAL_PM_POINT_LOCATION_BASE_H
#include <CGAL/Pm_point_location_base.h>
#endif

#ifndef CGAL_TRAPEZOIDAL_DECOMPOSITION_2_H
#include <CGAL/Trapezoidal_decomposition_2.h>
#endif

#ifndef CGAL_TD_TRAITS_H
#include <CGAL/Td_traits.h>
#endif

CGAL_BEGIN_NAMESPACE

////////////////////////////////////////////////////////////////////
//  DEFAULT PLANAR MAP STRATEGY
///////////////////////////////////////////////////////////////////

template <class Planar_map_>
class PL_X_curve_plus: public Planar_map_::Traits::X_curve
{
public:
  typedef Planar_map_ Planar_map;
  typedef typename Planar_map::Traits Traits;
  typedef typename Traits::X_curve curve; 
  typedef typename Traits::Point Point;

  typedef typename Planar_map::Locate_type Locate_type;
  typedef typename Planar_map::Halfedge_handle Halfedge_handle;
  typedef typename Planar_map::Halfedge_iterator Halfedge_iterator;
  
  PL_X_curve_plus() : curve(),parent() {};
  PL_X_curve_plus(const curve &cv,const Halfedge_handle& p) : 
    curve(cv), parent(p) {}
  PL_X_curve_plus(const Halfedge_handle& p) : curve(p->curve()),parent(p){}
  // used when no Halfedge_handle is supplied.
  PL_X_curve_plus(const curve &cv) : curve(cv),parent() {};
  PL_X_curve_plus(const PL_X_curve_plus &cv) : curve(cv),parent(cv.parent){}
  ~PL_X_curve_plus(){}
  PL_X_curve_plus& operator=(const PL_X_curve_plus &cv)
  {
    curve::operator=(cv);
    parent=cv.get_parent();
    return *this;
  }
  bool operator==(const PL_X_curve_plus &cv) const
  {
    return curve::operator==(cv);
  }
  Halfedge_handle get_parent() const
  {
    return parent;
  }
protected:
	Halfedge_handle parent;
};



//------------------------------------------------------------------


template <class Planar_map_>
class Pm_default_point_location : public Pm_point_location_base<Planar_map_> {
public:
  typedef Planar_map_ Planar_map;
  typedef Pm_point_location_base<Planar_map_> Base;
  typedef Pm_default_point_location<Planar_map> Self;
  typedef typename Planar_map::Traits Pm_traits;
  typedef typename Planar_map::Traits_wrap Pm_traits_wrap;
  typedef typename Planar_map::Locate_type Locate_type;
  typedef typename Planar_map::Halfedge_handle Halfedge_handle;
  typedef typename Planar_map::Halfedge_iterator Halfedge_iterator;
  typedef typename Planar_map::Ccb_halfedge_circulator Ccb_halfedge_circulator;
  typedef PL_X_curve_plus<Planar_map> X_curve_plus;
  typedef Td_traits<Pm_traits_wrap,X_curve_plus> Td_traits;
  typedef Trapezoidal_decomposition_2<Td_traits> Trapezoidal_decomposition;
  typedef Pm_bounding_box_base<Planar_map> Bounding_box;
  typedef typename Pm_traits::Point Point;
  typedef typename Pm_traits::X_curve X_curve;
  typedef typename Pm_point_location_base<Planar_map>::
    Halfedge_handle_container Halfedge_handle_container;
  typedef typename Pm_point_location_base<Planar_map>::Halfedge_handle_iterator
    Halfedge_handle_iterator;
  typedef typename Base::Token Token;
  
protected:
  typedef Trapezoidal_decomposition TD;
  typedef Planar_map PM;
  typedef const Self* cPLp;
  
	
public:
  Pm_default_point_location(bool rebuild=true) : traits(NULL){
    td.set_needs_update(rebuild);
  }
  ~Pm_default_point_location() 
  {
    if (traits) delete traits;
  }
  /* 
     Postcondition:
     h->curve() with a reference back to h
     is inserted into TD.
  */
  void insert(Halfedge_handle h,
              const X_curve& cv) //additions by iddo for arrangement
  {
    td.insert(X_curve_plus(cv,h));
  }
  
  Halfedge_handle locate(const Point& p, Locate_type& lt) const;
  Halfedge_handle locate(const Point& p, Locate_type& lt);
  
  Halfedge_handle vertical_ray_shoot(const Point& p, Locate_type& lt, bool up)
    const;
  Halfedge_handle vertical_ray_shoot(const Point& p, Locate_type& lt, bool up);
  
  //the function is called after the combinatoric split
  //we want O(1) is it possible?? no!
  void split_edge(const X_curve &cv,
                  Halfedge_handle e1,
                  Halfedge_handle e2
                  //additions by iddo for arrangement
                  ,const X_curve& cv1, const X_curve& cv2
                  //end additions
                  )
  {
    //	  td.split_edge(X_curve_plus(cv),X_curve_plus(e1),X_curve_plus(e2));
    td.split_edge(X_curve_plus(cv),X_curve_plus(cv1,e1),X_curve_plus(cv2,e2));
  }
  
  /*
    called after combinatoric merge
    e is the new edge cv1,cv2 are the original curves
  */
  void merge_edge(const X_curve &cv1,
                  const X_curve &cv2,
                  Halfedge_handle e 
                  //additions by iddo for arrangement
                  ,const X_curve& cv
                  //end additions
                  )
  {
    td.merge_edge(
                  X_curve_plus(cv1),
                  X_curve_plus(cv2),
                  //		X_curve_plus(e)
                  X_curve_plus(cv,e)
                  );
  }
  
  //called before combinatoric deletion
  void remove_edge(Halfedge_handle e)
  {
    td.remove(X_curve_plus(e));
  }
  void remove_edge(const Halfedge_handle_iterator& begin,
		const Halfedge_handle_iterator& end)
  {
    std::vector<X_curve_plus> c;
    Halfedge_handle_iterator it=begin;
    while (it!=end) { c.push_back((*it)->curve());++it;}
    td.remove(c.begin(),c.end());
  }
  
  inline void clear() { td.clear(); } 

  inline void update(const Halfedge_handle_iterator& begin,
                     const Halfedge_handle_iterator& end,
                     const Token& token)
    // Shuffle curves, remove them from the map and reinsert them afterwards.
  {

#ifdef CGAL_PMBB_DEBUG
		std::cout << "\nPL::update called with traits = "; 
		traits->debug();
#endif

		if (begin!=end)
		{

			Halfedge_handle_container c;
			Halfedge_handle_iterator it=begin;
			while (it!=end)
	      {
					c.push_back(Halfedge_handle(*it));
					++it;
				}
			remove_edge(begin,end);
			Halfedge_handle_iterator cend=c.end();
			it=c.begin();
			token.rebuild_bounding_box(this);
			// 
			while(it!=cend)
	      {
					insert(*it,(*it)->curve());
					++it;
				}
		}
		else
		{
			token.rebuild_bounding_box(this);
		}
#ifdef CGAL_PMBB_DEBUG
		std::cout << "\nPL::update exited with traits = "; 
		traits->debug();
#endif
	}  
  
  const TD* get_trapezoidal_decomposition() const {return &td;}
  
  inline const Pm_traits* get_traits() const {return traits;}

  void init(Planar_map& pmp, Pm_traits& tr) {
    pm = &pmp;
    if (traits) delete traits;
    traits = new Td_traits(tr);
    td.init_traits(traits);
  }

#ifdef CGAL_PM_DEBUG
public:
  void debug()
  {

#ifdef CGAL_TD_DEBUG
    td.debug();
#endif

  }
#endif
  
protected:
  TD td;

private:
  Planar_map* pm;
  Td_traits* traits;

  bool halfedge_represents_point(const Halfedge_handle& h,const Point& p) const
  {
    const Point 
      &source=h->source()->point(),
      &target=h->target()->point();
    return !(!traits->point_is_same_x(target,p)||
             traits->point_is_same_x(source,p)&&
             (traits->point_is_higher(p,target)&&
              traits->point_is_lower(target,source)||
              traits->point_is_lower(p,target)&&
              traits->point_is_higher(target,source)
              ));
  }
  /* 
     description:
     returns if the point is located in the 
     open halfplane on the right side of the
     input curve
     (considering the halfedge orientation)
     preconditions:
     the input curve is not vertical,
     p is in its x range but not on its closure
  */
  bool halfedge_represents_point_inside_face(const Halfedge_handle& h,
                                             const Point& p) const 
  {
    return (traits->curve_get_point_status(h->curve(),p)
	    ==Pm_traits::ABOVE_CURVE)
      ==traits->point_is_left(h->source()->point(),h->target()->point());
  }
  Halfedge_handle halfedge_representing_unbounded_face() const
  {
    CGAL_assertion(pm);
    if (pm->unbounded_face()->holes_begin()!=pm->unbounded_face()->holes_end())
 {
   // case PM is not empty
   //return *(pm->unbounded_face()->holes_begin());
   typename Planar_map::Holes_iterator hot=pm->unbounded_face()->holes_begin();
   return (*hot);
 }
    else
      // case PM is empty
      return pm->halfedges_end();
  }
  
  //use the latter
  //to workaround internal compiler error in egcs1.1
  //	Locate_type convert(const Point& p,const TD::Locate_type& lt,
  // Halfedge_handle& h,bool up=true) const	
  Locate_type convert(const Point& p,const int& lt,Halfedge_handle& h,
		      bool up=true) const
  {
    switch(lt)
      {
        // h->target() should represent p
      case TD::POINT:
        if (!halfedge_represents_point(h,p)) h=h->twin();
        return PM::VERTEX;
      case TD::CURVE:
        /* special case:
           h->source()->point(),p,h->target()->point() have same x
           coardinates.
           return value should be h(no h->twin()).
        */
        // orientation of h
        if (up==traits->point_is_left(h->source()->point(),
                                      h->target()->point())) 
          h=h->twin();
        return PM::EDGE;
      case TD::TRAPEZOID:
        if (!halfedge_represents_point_inside_face(h,p)) h=h->twin();
        CGAL_postcondition(halfedge_represents_point_inside_face(h,p));
        if (!h->face()->is_unbounded())
          return PM::FACE;
        // otherwise pass to UNBOUNDED_TRAPEZOID case
      case TD::UNBOUNDED_TRAPEZOID:
        h=halfedge_representing_unbounded_face();
        CGAL_postcondition(h->face()->is_unbounded());
        return PM::UNBOUNDED_FACE;
      default:
        CGAL_assertion(lt==TD::POINT||lt==TD::CURVE||lt==TD::TRAPEZOID||
                       lt==TD::UNBOUNDED_TRAPEZOID);
      }
    return Locate_type();
  }
  const Bounding_box* get_bounding_box() const {return pm->get_bounding_box();}
};


CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Pm_default_point_location.C>
#endif

#endif


/*
X_curve msvc6 workaround:

typedef typename Traits::X_curve X_curve;
 causes 
error C2086: '<Unknown>' : redefinition
*/
	 
/*
egcs workaround:

To solve the internal compiler errors egcs had when dealing
with nested templated classes we moved them outside:

template <class Planar_map_>
class PL_X_curve_plus: public Planar_map_::Traits::X_curve	 
*/
