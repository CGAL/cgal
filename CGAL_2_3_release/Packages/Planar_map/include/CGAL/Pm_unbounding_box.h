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
// file          : include/CGAL/Pm_unbounding_box.h
// package       : pm (4.20)
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//                 
//
// maintainer(s) : Oren Nechushtan <theoren@math.tau.ac.il>
//                 
//
// coordinator   : Tel-Aviv University (Dan Halperin)
//
// Chapter       : 
// email         : cgal@cs.uu.nl
//
// ======================================================================

#ifndef CGAL_PM_UNBOUNDING_BOX_H
#define CGAL_PM_UNBOUNDING_BOX_H

#ifndef CGAL_PLANAR_MAP_2_H
#include <CGAL/Planar_map_2.h>
#endif

#ifndef CGAL_PM_BOUNDING_BOX_BASE_H
#include <CGAL/Pm_bounding_box_base.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Planar_map_>
class Pm_unbounding_box : public Pm_bounding_box_base<Planar_map_> {
public:
	typedef Planar_map_ Planar_map;
	typedef Pm_unbounding_box<Planar_map> Self;
	//	typedef Planar_map_2<Dcel,Traits,Self> Planar_map;
	/*
	typedef Planar_map_2<Dcel,Traits> Base;
	typedef Planar_map_Bbox_2<Dcel,Traits> Self;
	typedef Planar_map_traits_wrap<Traits> Traits_wrap;
	*/
	typedef typename Planar_map::Traits Traits;
	typedef typename Traits::X_curve X_curve;
	typedef typename Traits::Point Point;	

        typedef typename std::vector<Point>::iterator Point_iterator;
	typedef typename std::vector<X_curve>::iterator X_curve_iterator;

	/*	typedef typename Traits::Bounding_box Bounding_box;
	typedef typename Traits::Boundary_type Boundary_type;
	typedef typename Traits::Point_boundary_container Point_boundary_container;
	typedef typename Traits::X_curve_boundary_container X_curve_boundary_container;
	*/
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
	
	/*
	typedef typename Base::Halfedge_handle Halfedge_handle;
	typedef typename Base::Face_handle Face_handle;
	typedef typename Base::Vertex_handle Vertex_handle;
	typedef typename Base::Vertex_const_handle Vertex_const_handle;
	typedef typename Base::Halfedge_const_handle Halfedge_const_handle;
	typedef typename Base::Face_const_handle Face_const_handle;
	typedef typename Base::Vertex_iterator Vertex_iterator;
	typedef typename Base::Halfedge_iterator Halfedge_iterator;
	typedef typename Base::Face_iterator Face_iterator;
	typedef typename Base::Vertex_const_iterator Vertex_const_iterator;
	typedef typename Base::Halfedge_const_iterator Halfedge_const_iterator;
	typedef typename Base::Face_const_iterator Face_const_iterator;
	typedef typename Base::Locate_type Locate_type;
	typedef Pm_point_location_base<Base> Point_location_base;
	//  typedef std::list<X_curve> X_curve_container;
	//  typedef Topological_map<_Dcel> TPM;
	typedef typename Base::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
	typedef typename Base::Holes_iterator Holes_iterator;
	typedef typename Base::Holes_const_iterator Holes_const_iterator;
	typedef typename Base::Ccb_halfedge_const_circulator Ccb_halfedge_const_circulator;
	typedef typename Base::Ccb_halfedge_circulator Ccb_halfedge_circulator;
	typedef typename Base::Size Size;
	*/	
	
	Pm_unbounding_box(){}
	~Pm_unbounding_box(){}
	
	
	void init(Planar_map& pmp, Traits& tr) {
	/*
    pm = &pmp;
    traits = (Traits_wrap*)(&tr);
		*/
	}
	
	bool insert(const Point& p) {return true;}
	bool insert(const Point_iterator& begin,const Point_iterator& end
#ifndef _MSC_VER
  )
#else
		,Point* dummy=0) 	// workaround for MSVC6.0
#endif
          {return true;}
	bool insert(const X_curve& cv) {return true;}
	bool insert(const X_curve_iterator& begin,const X_curve_iterator& end
#ifndef _MSC_VER
  )
#else
		,X_curve* dummy=0) 	// workaround for MSVC6.0
#endif          

          {return true;}
	/* The point location query function may updates the resulting 
	halfedge handle and locate type as expected from the bounding box */
	bool locate(const Point& p, Locate_type& lt,Halfedge_handle& h){return true;}	
	bool vertical_ray_shoot(const Point& p, Locate_type& lt, bool up,
		Halfedge_handle& h){return true;}
	
	void split_edge(const X_curve &cv,
		Halfedge_handle e1,
		Halfedge_handle e2,
		const X_curve& cv1, 
		const X_curve& cv2
		) {}

	void split_boundary_edge(const Halfedge_handle &h,
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
	inline bool is_empty() const {return false;}

#ifdef CGAL_PM_DEBUG
	void debug() const {}
#endif
	
};

CGAL_END_NAMESPACE

#endif //CGAL_PM_UNBOUNDING_BOX_H
