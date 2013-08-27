// Copyright (c) 2013 Technical University Braunschweig (Germany).
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
// Author(s):  Kan Huang <huangkandiy@gmail.com>
//             

#ifndef CGAL_TRIANGULAR_EXPANSION_VISIBILITY_2_H
#define CGAL_TRIANGULAR_EXPANSION_VISIBILITY_2_H

#include <CGAL/Arrangement_2.h>
#include <stack>
#include <deque>

namespace CGAL {

template<class Arrangement_2> 
class Triangular_expansion_visibility_2 {

public:
	typedef typename Arrangement_2::Geometry_traits_2				Geometry_traits_2;
	// Currently only consider with same type for both
	typedef Arrangement_2											Input_Arrangement_2;
	typedef Arrangement_2											Output_Arrangement_2;

	typedef typename Arrangement_2::Halfedge_const_handle			Halfedge_const_handle;
	typedef typename Arrangement_2::Ccb_halfedge_const_circulator 	Ccb_halfedge_const_circulator;
	typedef typename Arrangement_2::Face_const_handle				Face_const_handle;
    typedef typename Arrangement_2::Kernel                          Kernel;
    typedef typename CGAL::Arr_linear_traits_2<Kernel>              Linear_traits_2;

	typedef typename Geometry_traits_2::Point_2						Point_2;
	typedef typename Geometry_traits_2::Ray_2						Ray_2;
	typedef typename Geometry_traits_2::Segment_2					Segment_2;
    typedef typename Geometry_traits_2::Line_2                      Line_2;
	typedef typename Geometry_traits_2::Vector_2                	Vector_2;
	typedef typename Geometry_traits_2::FT                			Number_type;



    Triangular_expansion_visibility_2() : p_arr(NULL) {}

	/*! Constructor given an arrangement and the Regularization tag. */
    Triangular_expansion_visibility_2(const Input_Arrangement_2& arr/*, Regularization_tag r_t*/): p_arr(&arr) {

    }

	bool is_attached() {
		return (p_arr != NULL);
	}

    void attach(const Input_Arrangement_2& arr) {
		p_arr = &arr;
	}

	void detach() {
		p_arr = NULL;
	}

	Input_Arrangement_2 arr() {
		return *p_arr;
	}

	void compute_visibility(const Point_2& q, 
						   const Face_const_handle face,
						   Output_Arrangement_2& out_arr
                           ) {

    }

	void compute_visibility(const Point_2& q, 
						   const Halfedge_const_handle he,
						   Output_Arrangement_2& out_arr
						   ) {

}

private:
    Input_Arrangement_2* arr;
    Line_Arrangement_2   line_arr;
    void preprocess() {

    }

    Line_2 dual_line(const Point_2& p) {
        return Line_2(p.x(), -1, -p.y());
    }

};

} // namespace CGAL

#endif
