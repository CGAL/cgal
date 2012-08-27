// Copyright (c) 2010 Inria Lorraine (France). All rights reserved.
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
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#include <mpfi.h>

int main(){
        mpfi_t a;
        mpfi_init_set_si(a,-100);
        mpfi_add_ui(a,a,101);
        if(mpfi_cmp_ui(a,1)==0) {
          mpfi_clear(a);
          return 0;
        }
        else {
          return 1;
        }
}
