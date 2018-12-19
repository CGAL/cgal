// Copyright (c) 1999-2016 INRIA Sophia Antipolis, INRIA Nancy - Grand Est (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Iordan Iordanov  <Iordan.Iordanov@loria.fr>
//

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_VERTEX_BASE_2_H
#define CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_VERTEX_BASE_2_H

#include <CGAL/license/Periodic_4_hyperbolic_triangulation_2.h>

#include <CGAL/basic.h>
#include <CGAL/Triangulation_vertex_base_2.h>

namespace CGAL {

template< class GT, class Vb = CGAL::Triangulation_vertex_base_2<GT> >
class Periodic_4_hyperbolic_triangulation_vertex_base_2 : public Vb {
public:
  	typedef Vb                                          Base;
	typedef typename Vb::Triangulation_data_structure   Triangulation_data_structure;
	typedef typename Triangulation_data_structure::Face_handle 		             
        	                                            Face_handle;
	typedef typename GT::Point_2			            Point;
	typedef typename GT::Hyperbolic_translation 		Hyperbolic_translation;


	template <typename TDS2>
  	struct Rebind_TDS {
      	typedef typename Vb::template Rebind_TDS<TDS2>::Other                  Vb2;
    	typedef Periodic_4_hyperbolic_triangulation_vertex_base_2<GT, Vb2>     Other;
  	};

private:
	Hyperbolic_translation 		_tr; 	// Used only during insert();
	bool 		                _stored_translation;

public:
	Periodic_4_hyperbolic_triangulation_vertex_base_2() : 
  	Base(), _stored_translation(false)
	{}

	Periodic_4_hyperbolic_triangulation_vertex_base_2(const Point & p) : 
	Base(p), _stored_translation(false) 
	{}

  	Periodic_4_hyperbolic_triangulation_vertex_base_2(const Point & p, Face_handle fh) : 
  	Base(p, fh), _stored_translation(false) 
  	{}

	Periodic_4_hyperbolic_triangulation_vertex_base_2(const Face_handle& fh) :
	Base(), _stored_translation(false) 
	{}

	void set_translation(const Hyperbolic_translation& tr) {
		if (!_stored_translation) {
			_tr = tr;
			_stored_translation = true;
		}
	}

	Hyperbolic_translation translation() const {
		return _tr;
	}

	void clear_translation() {
		_tr = Hyperbolic_translation();
		_stored_translation = false;
	}

	bool get_translation_flag() const {
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


}  // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_VERTEX_BASE_2_H