// Copyright (c) 2006 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Luis Peñaranda <penarand@loria.fr>

// Tests if RS is available (this is originally a librs example).

#include <cstdlib>
#include <rs_exports.h>

int main() {
	double poly[] = {1, -2, 3, -5, 6, -7};
	int deg = 5;
	double *res = (double *) malloc ((deg+1) * sizeof (double));
	int nb = 0;
	rs_init_rs ();
	rs_reset_all ();
	set_rs_precisol (23);
	rs_num_solve_u (poly,deg+1, res, &nb);
	/* we don't need to show results now */
	return 0;
}
