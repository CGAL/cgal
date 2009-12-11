// Copyright (c) 2009 Inria Lorraine (France). All rights reserved.
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
// Author: Luis Peñaranda <luis.penaranda@loria.fr>

#include <CGAL/Gmpfi.h>

#ifdef NDEBUG
#  undef NDEBUG
#  include <cassert>
#  define NDEBUG 1
#endif

int main(){
        typedef CGAL::Gmpfi Gmpfi;
	Gmpfi::Precision_type new_precision=200;
	Gmpfi::Precision_type old_precision;
	Gmpfi one(1);
	old_precision=Gmpfi::set_default_precision(new_precision);
	// this does not make sense when new_precision==old_precision
	assert(old_precision!=new_precision);
	Gmpfi two(2);
	assert(one.get_precision()==old_precision);
	assert(two.get_precision()==new_precision);
	// now, test the precision in mu
        return 0;
}
