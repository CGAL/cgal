// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef   CGAL_PARTITION_OPT_CVX_DIAGONAL_LIST_H
#define   CGAL_PARTITION_OPT_CVX_DIAGONAL_LIST_H

#include <CGAL/license/Partition_2.h>


#include <utility>
#include <list>
#include <iostream>

typedef std::pair<int, int>                   Partition_opt_cvx_diagonal;
typedef std::list<Partition_opt_cvx_diagonal> Partition_opt_cvx_diagonal_list;

inline
std::ostream& operator<<(std::ostream& os,
                         const Partition_opt_cvx_diagonal_list& d)
{
   Partition_opt_cvx_diagonal_list::const_iterator it;
   for (it = d.begin(); it != d.end(); it++)
   {
      os << "(" << (*it).first << ", " << (*it).second << ") ";
   }
   return os;
}

#endif // CGAL_PARTITION_OPT_CVX_DIAGONAL_LIST_H
