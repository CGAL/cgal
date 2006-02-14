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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
#ifndef CGAL_CIRCULATOR_LIST_H
#define CGAL_CIRCULATOR_LIST_H

#include <list>
#include <iostream>

namespace CGAL {

template <class Circulator>
class Circulator_list : public std::list<Circulator> 
{

public:
  Circulator_list() 
  {}

  Circulator_list(Circulator first) 
  {
      if (first == NULL) return;

      Circulator current = first;
      do 
      {
         push_back(current);
      } while (++current != first);
  }
};

template <class C>
std::ostream& operator<<(std::ostream& os, const Circulator_list<C>& c)
{
     typename Circulator_list<C>::const_iterator current;
     for (current = c.begin(); current != c.end(); current++)
     {
        os << **current << " ";
     }
     return os;
}
}

#endif // CGAL_CIRCULATOR_LIST_H
