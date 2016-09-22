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
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Iordan Iordanov  <Iordan.Iordanov@loria.fr>


#ifndef CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_DS_VERTEX_BASE_2_H
#define CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_DS_VERTEX_BASE_2_H

#include <CGAL/basic.h>
#include <CGAL/Dummy_tds_2.h>
#include <CGAL/Hyperbolic_octagon_word_4.h>

namespace CGAL {

template< typename GT, typename TDS = void >
class Periodic_4_hyperbolic_triangulation_ds_vertex_base_2 {
public:
	typedef TDS 							Triangulation_data_structure;
	typedef typename TDS::Vertex_handle     Vertex_handle;
	typedef typename TDS::Face_handle 		Face_handle;
	typedef typename GT::Point_2			Point;
	typedef unsigned short int 				Int;
	typedef Hyperbolic_octagon_word_4<Int, GT>		Offset;


	template <typename TDS2>
  	struct Rebind_TDS {
    	typedef Periodic_4_hyperbolic_triangulation_ds_vertex_base_2<GT, TDS2> Other;
  	};

private:
	Face_handle _face;
	int 		_idx;
	Point 		_p;
	Offset 		_o; 	// Used only during insert();
	bool 		_stored_offset;

public:
	Periodic_4_hyperbolic_triangulation_ds_vertex_base_2() :
	_face(), _stored_offset(false), _idx(-1)
	{}

	Periodic_4_hyperbolic_triangulation_ds_vertex_base_2(const Point & p) : 
	_face(), _p(p), _stored_offset(false), _idx(-1) 
	{}

  Periodic_4_hyperbolic_triangulation_ds_vertex_base_2(const Point & p, Face_handle fh) : 
  _face(fh), _p(p), _stored_offset(false), _idx(-1) 
  {}

	Periodic_4_hyperbolic_triangulation_ds_vertex_base_2(const Face_handle& fh) :
	_face(fh), _stored_offset(false), _idx(-1)
	{}

	const Face_handle& face() {
		return _face;
	} 

	void set_idx(int idx) {
		_idx = idx;
	}

	int idx() {
		return _idx;
	}

	void set_face(const Face_handle& fh) {
		_face = fh;
	}

	void set_point(const Point & p) { _p = p; }
  	
  	const Point&  point() const { return _p; }

  	// the non const version of point() is undocument
  	// but needed to make the point iterator works
  	// using Lutz projection scheme
  	Point&        point() { return _p; }


	// the following trivial is_valid allows
  	// the user of derived face base classes 
  	// to add their own purpose checking
  	bool is_valid(bool = false, int = 0) const { 
    	return face() != Face_handle();
  	}

  	// For use by the Compact_container.
  	void *   for_compact_container() const { 
  		return _face.for_compact_container(); 
  	}
  	
  	void * & for_compact_container() { 
  		return _face.for_compact_container(); 
  	}

  	void store_offset(Offset o) {
  		if (!_stored_offset) {
  			_o = o;
        //cout << " -- Stored offset " << _o << " in vertex " << _idx << endl;
  			_stored_offset = true;
  		}
  	}

  	Offset get_offset() {
        //CGAL_assertion(_stored_offset);
  		return _o;
  	}

  	void remove_offset() {
      //std::cout << " -- Clearing offset " << _o << " from vertex " << _idx << std::endl;
  		_o = Offset();
  		_stored_offset = false;
  	}

  	bool stored_offset() {
  		return _stored_offset;
  	}

};

template < class TDS >
inline 
std::istream&
operator>>(std::istream &is, Periodic_4_hyperbolic_triangulation_ds_vertex_base_2<TDS> &) {
	return is;
}

template < class TDS >
inline
std::ostream&
operator<<(std::ostream &os, Periodic_4_hyperbolic_triangulation_ds_vertex_base_2<TDS> &) {
	return os;
}


// Specialization vor 'void'
template<typename GT>
class Periodic_4_hyperbolic_triangulation_ds_vertex_base_2<GT, void> {
public:
	typedef Dummy_tds_2 									Triangulation_data_structure;
	typedef Triangulation_data_structure::Vertex_handle		Vertex_handle;
	typedef Triangulation_data_structure::Face_handle 		Face_handle;
	struct  Point {};	
	template<typename TDS2>
	struct Rebind_TDS {
    	typedef Periodic_4_hyperbolic_triangulation_ds_vertex_base_2<GT, TDS2> Other;
  	};
};

}  // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_DS_VERTEX_BASE_2_H