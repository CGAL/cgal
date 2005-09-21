// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Laurent Saboret, Pierre Alliez
//
// Purpose: Workaround a Visual C++ bug that emits C4503 "name too long" warnings


#include <CGAL/config.h>

#if defined(CGAL_CFG_LONGNAME_BUG)

#define Cartesian                       Cn
#define Simple_cartesian                SC
#define I_HalfedgeDS_iterator           IPI
#define I_HalfedgeDS_const_iterator     IPCI
#define In_place_list_iterator          IPLI
#define In_place_list_const_iterator    IPLCI

#if defined(_MSC_VER)
// Has no effect, probably bug in MSVC
#pragma warning(disable:4503)
#endif

#endif
