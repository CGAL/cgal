// Copyright (c) 1998  ETH Zurich (Switzerland).
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
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Gabriele Neyer<neyer@inf.ethz.ch>
#ifndef __R_tree_internal_storage_H__
#define __R_tree_internal_storage_H__
#include<CGAL/R_tree_internal_db.h>
#include<CGAL/R_tree.h>

CGAL_BEGIN_NAMESPACE

class R_tree_internal_storage{
 public:
  typedef CGAL::R_tree_internal_db IO_tree_traits_nodes;
  typedef CGAL::R_tree_internal_db IO_tree_traits_leaves;
  const static  int IO_min_cap_nodes=2;
  const static  int IO_max_cap_nodes=4;
  const static  int IO_min_cap_leaves=2;
  const static  int IO_page_size=70;
  const static bool headerextern=false;
};

CGAL_END_NAMESPACE
#endif
