// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_INDIRECT_NOT_LESS_YX_2_H
#define CGAL_INDIRECT_NOT_LESS_YX_2_H

namespace CGAL {

template <class Traits>
class Indirect_not_less_yx_2 
{
   public:
     typedef typename Traits::Less_yx_2     Less_yx_2;

     Indirect_not_less_yx_2(const Traits& traits) : 
         less_yx_2(traits.less_yx_2_object()) {}

     template <class Iterator>
     bool 
     operator()( const Iterator& p, const Iterator& q) const
     { return less_yx_2( *q, *p); }

   private:
     Less_yx_2 less_yx_2;
};

}

#endif // CGAL_INDIRECT_NOT_LESS_YX_2_H
