// Copyright (c) 1999-2016   INRIA Nancy - Grand Est (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Iordan Iordanov  <Iordan.Iordanov@loria.fr>
//                 


#ifndef CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_VERTEX_BASE_2_H
#define CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_VERTEX_BASE_2_H

#include <CGAL/basic.h>
#include <CGAL/Dummy_tds_2.h>

namespace CGAL {

template< class GT, class Vb = CGAL::Triangulation_vertex_base_2<GT> >
class Periodic_4_hyperbolic_triangulation_vertex_base_2 : public Vb {
public:
  typedef Vb                                         Base;
	typedef typename Vb::Triangulation_data_structure  Triangulation_data_structure;
  typedef Triangulation_data_structure               TDS;
	typedef typename TDS::Vertex_handle                Vertex_handle;
	typedef typename TDS::Face_handle 		             Face_handle;
	typedef typename GT::Point_2			                 Point;
	typedef typename GT::Hyperbolic_translation 			                 Hyperbolic_translation;


	template <typename TDS2>
  	struct Rebind_TDS {
      typedef typename Vb::template Rebind_TDS<TDS2>::Other                  Vb2;
    	typedef Periodic_4_hyperbolic_triangulation_vertex_base_2<GT, Vb2>     Other;
  	};

private:
	Face_handle _face;
	Point 		_p;
	Hyperbolic_translation 		_tr; 	// Used only during insert();
	bool 		_stored_translation;

public:
	Periodic_4_hyperbolic_triangulation_vertex_base_2() : 
  Base(), _face(), _stored_translation(false)
	{}

	Periodic_4_hyperbolic_triangulation_vertex_base_2(const Point & p) : 
	Base(p), _face(), _p(p), _stored_translation(false) 
	{}

  Periodic_4_hyperbolic_triangulation_vertex_base_2(const Point & p, Face_handle fh) : 
  Base(p, fh), _face(fh), _p(p), _stored_translation(false) 
  {}

	Periodic_4_hyperbolic_triangulation_vertex_base_2(const Face_handle& fh) :
	Base(), _face(fh), _stored_translation(false) 
	{}

	const Face_handle& face() {
		return _face;
	} 

	void set_face(const Face_handle& fh) {
		_face = fh;
	}

	void set_point(const Point & p) { _p = p; }
  	
  	const Point&  point() const { return _p; }

  	// the non const version of point() is undocumented
  	// but needed to make the point iterator work
  	// using Lutz projection scheme
  	Point&        point() { return _p; }


  	// For use by the Compact_container.
  	void * for_compact_container() const { 
  		return _face.for_compact_container(); 
  	}
  	
  	void * & for_compact_container() { 
  		return _face.for_compact_container(); 
  	}

  	void store_translation(Hyperbolic_translation tr) {
  		if (!_stored_translation) {
  			_tr = tr;
  			_stored_translation = true;
  		}
  	}

  	Hyperbolic_translation get_translation() {
  		return _tr;
  	}

  	void remove_translation() {
  		_tr = Hyperbolic_translation();
  		_stored_translation = false;
  	}

  	bool stored_translation() {
  		return _stored_translation;
  	}

};

template < class TDS >
inline 
std::istream&
operator>>(std::istream &is, Periodic_4_hyperbolic_triangulation_vertex_base_2<TDS> &) {
	return is;
}

template < class TDS >
inline
std::ostream&
operator<<(std::ostream &os, Periodic_4_hyperbolic_triangulation_vertex_base_2<TDS> &) {
	return os;
}


// Specialization vor 'void'
template<typename GT>
class Periodic_4_hyperbolic_triangulation_vertex_base_2<GT, void> {
public:
	typedef Dummy_tds_2 									Triangulation_data_structure;
	typedef Triangulation_data_structure::Vertex_handle		Vertex_handle;
	typedef Triangulation_data_structure::Face_handle 		Face_handle;
	struct  Point {};	
	template<typename TDS2>
	struct Rebind_TDS {
    	typedef Periodic_4_hyperbolic_triangulation_vertex_base_2<GT, TDS2> Other;
  	};
};

}  // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_DS_VERTEX_BASE_2_H