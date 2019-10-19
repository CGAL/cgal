// Copyright (c) 2010 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#include <iostream>
#include <string>
#include <mpfi.h>

int main(){
    std::cout<<"version="<<std::string(mpfi_get_version())<<std::endl;
    return 0;
}
