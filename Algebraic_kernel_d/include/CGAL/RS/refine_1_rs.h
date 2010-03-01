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

#ifndef CGAL_RS_REFINE_1_RS_H
#define CGAL_RS_REFINE_1_RS_H

#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>
#include <CGAL/RS/algebraic_1.h>
#include <CGAL/assertions.h>

namespace RS3{

inline void refine_1(const CGAL::Algebraic_1 &a,unsigned int s=10000){
        CGAL_precondition(a.inf()<=a.sup());
        rs3_refine_u_root((mpfi_ptr)a.mpfi(),
                          a.pol().get_coefs(),
                          a.pol().get_degree(),
                          mpfi_get_prec(a.mpfi())+s,
                          0,
                          0);
        CGAL_assertion(a.inf()<=a.sup());
}

} // namespace RS3

#endif  // CGAL_RS_REFINE_1_RS_H
