// Copyright (c) 1997-2002  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// Author(s)     : Radu Ursu


#ifdef _MSC_VER
#  if _MSC_VER < 1310
#    error Unsupported version of VC++
#  else
#    ifdef __ICL
#      include "CGAL/icl_8_1.h"
#    else 
#      if _MSC_VER == 1310
#        include "CGAL/cl_1310.h"
#      else
#        include "CGAL/cl_1400.h"
#      endif
#    endif
#  endif

#  define CGAL_LIB_NAME CGAL
#  include "CGAL/auto_link.h"
#endif
