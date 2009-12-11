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
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>

int main(){
        typedef CGAL::Gmpfi                     NT;
        typedef CGAL::Field_with_kth_root_tag   Tag;
        typedef CGAL::Tag_false                 Is_exact;
        CGAL::test_algebraic_structure<NT,Tag,Is_exact>();
        CGAL::test_algebraic_structure<NT,Tag,Is_exact>(NT(4),NT(6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag,Is_exact>(NT(-4),NT(6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag,Is_exact>(NT(4),NT(-6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag,Is_exact>(NT(-4),NT(-6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag,Is_exact>(NT(4),NT(6),NT(-15));
        CGAL::test_algebraic_structure<NT,Tag,Is_exact>(NT(-4),NT(6), NT(15));
        CGAL::test_algebraic_structure<NT,Tag,Is_exact>(NT(4),NT(-6),NT(-15));
        CGAL::test_algebraic_structure<NT,Tag,Is_exact>(NT(-4),NT(-6),NT(-15));
        CGAL::test_real_embeddable<NT>();
        return 0;
}
