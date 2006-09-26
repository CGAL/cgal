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
// $URL:  $
// $Id:  $
// 
//
// Author(s)     : Luis Peñaranda <penarand@loria.fr>

// Tests if MPFI is available.

#include <iostream>
#include "gmp.h"
#include "mpfi.h"
#include "mpfi_io.h"

int main()
{
	mpfi_t a, b;	// declare two intervals
	mpfi_init_set_si (a, 100);	// initialize & set a
	mpfi_init (b);	// initialize b

	mpfi_div_si (b, a, 3);	// b=a/3

	return 0;
}
