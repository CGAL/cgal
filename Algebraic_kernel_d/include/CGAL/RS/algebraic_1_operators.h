// Copyright (c) 2006-2008 Inria Lorraine (France). All rights reserved.
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

#ifndef CGAL_RS_ALGEBRAIC_1_OPERATORS_H
#define CGAL_RS_ALGEBRAIC_1_OPERATORS_H

namespace CGAL{

inline
Algebraic_1 Algebraic_1::operator+()const{
        return *this;
}

inline
Algebraic_1 Algebraic_1::operator-()const{
        mpfi_t inv;
        mpfi_init2(inv,mpfi_get_prec(mpfi()));
        mpfi_neg(inv,mpfi());
        Algebraic_1 *inverse=new Algebraic_1(inv,
                                             pol().minusx(),
                                             nr(),
                                             mult(),
                                             -lefteval());
        return *inverse;
}

} // namespace CGAL

#endif  // CGAL_RS_ALGEBRAIC_1_OPERATORS_H
