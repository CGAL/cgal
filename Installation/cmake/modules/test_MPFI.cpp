// Copyright (c) 2010 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
