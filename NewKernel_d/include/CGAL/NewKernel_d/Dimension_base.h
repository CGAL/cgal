// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Marc Glisse

#ifndef CGAL_KD_DIMENSION_BASE_h
#define CGAL_KD_DIMENSION_BASE_h
#include <CGAL/Dimension.h>
#include <CGAL/assertions.h>
#include <CGAL/NewKernel_d/utils.h>
namespace CGAL {
struct Store_dimension_base {
	//TODO: add some assertions
	Store_dimension_base(int dim=UNKNOWN_DIMENSION):dim_(dim){}
	int dimension()const{return dim_;}
	void set_dimension(int dim){dim_=dim;}
	private:
	int dim_;
};
template<class=Dynamic_dimension_tag>
struct Dimension_base {
	Dimension_base(int = UNKNOWN_DIMENSION){}
	int dimension() const { return UNKNOWN_DIMENSION; }
	void set_dimension(int) {}
};
template<int dim_>
struct Dimension_base<Dimension_tag<dim_> > {
	Dimension_base(){}
	Dimension_base(int CGAL_assertion_code(dim)){CGAL_assertion(dim_==dim);}
	int dimension()const{return dim_;}
	void set_dimension(int dim){CGAL_assertion(dim_==dim);}
};
}
#endif

