/*
 * Tds_output.h
 *
 *  Created on: Jul 10, 2014
 *      Author: ivovigan
 */

#ifndef TDS_OUTPUT_H_
#define TDS_OUTPUT_H_

// Copyright (c) 2014  INRIA Sophia-Antipolis (France), INRIA Lorraine LORIA.
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
// Author(s)     : Ivo Vigan

//CGAL
#include <CGAL/basic.h>

//local
#include <CGAL/Reconstruction_triangulation_2.h>


namespace CGAL {

template<class Kernel>
class Tds_output {
public:
	typedef Reconstruction_triangulation_2<Kernel> Rt_2;
	typedef typename Rt_2::Triangulation_data_structure Tds_2;

private:
	Tds_2 reconstruction_tds;

public:
	void store_marked_elements(Tds_2 tds, int nb_ignore) {
		reconstruction_tds = tds;
	}

	Tds_2 get_reconstruction_tds() {
		return reconstruction_tds;
	}

};
}



#endif /* TDS_OUTPUT_H_ */
