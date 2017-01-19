// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_SNC_SM_EXPLORER_H
#define CGAL_SNC_SM_EXPLORER_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/basic.h>
#include <CGAL/Nef_S2/SM_const_decorator.h>

namespace CGAL {

template <typename SMCDEC>
class SNC_SM_explorer : public SMCDEC {
  typedef SMCDEC                            Base;
  typedef SNC_SM_explorer<SMCDEC>           Self;
  //  typedef typename Base::Kernel    Kernel;
  //  typedef typename Kernel::Point_3 Point_3;

 public:
  SNC_SM_explorer(const Base& E) : Base(E) {}
  Self& operator=(const Self& E) {
    Base::operator=(E); 
    return *this;
  }
};

} //namespace CGAL

#endif  // CGAL_SNC_SM_EXPLORER_H
