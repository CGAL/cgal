// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#include <iostream>

namespace CGAL {


// makes *p the rightmost child of *q
template<class Traits>
void Rotation_tree_2<Traits>::set_rightmost_child(Self_iterator p, 
                                                  Self_iterator q)
{
   CGAL_assertion(q != this->end());

   if (p != this->end())
   {
      (*p).clear_right_sibling();
      if (rightmost_child(q) != this->end())
      {
         (*p).set_left_sibling(rightmost_child(q));
         (*rightmost_child(q)).set_right_sibling(p);
      }
      else
         (*p).clear_left_sibling();
      (*p).set_parent(q);
      (*q).set_rightmost_child(p);
   }
   else
   {
      (*q).clear_rightmost_child();
   }
}

// makes *p the left sibling of *q
template <class Traits>
void Rotation_tree_2<Traits>::set_left_sibling(Self_iterator p, 
                                               Self_iterator q)
{
   if (q == this->end()) return;
       
   if (p != this->end())
   {
      if (left_sibling(q) != this->end())
      {
         (*left_sibling(q)).set_right_sibling(p);
         (*p).set_left_sibling(left_sibling(q));
      }
      else
         (*p).clear_left_sibling();

      (*q).set_left_sibling(p);
      (*p).set_right_sibling(q);
      set_parent(parent(q),p);
   }
   else
   {
      if (left_sibling(q) != this->end())
         (*(*q).left_sibling()).clear_right_sibling();
      (*q).clear_left_sibling();
   }
}

// makes p the right sibling of q
template <class Traits>
void Rotation_tree_2<Traits>::set_right_sibling(Self_iterator p, 
                                                Self_iterator q)
{
   if (q == this->end()) return;
       
   if (p != this->end())
   {
      if (right_sibling(q) != this->end())
      {
         (*right_sibling(q)).set_left_sibling(p);
         (*p).set_right_sibling(right_sibling(q));
      }
      else
         (*p).clear_right_sibling();
      (*q).set_right_sibling(p);
      (*p).set_left_sibling(q);
      set_parent(parent(q),p);
   }
   else
   {
      if (right_sibling(q) != this->end())
         (*right_sibling(q)).clear_left_sibling();
      (*q).clear_right_sibling();
   }
}

// NOTE:  this function does not actually remove the node p from the
//        list; it only reorganizes the pointers so this node is not
//        in the tree structure anymore
template <class Traits>
void Rotation_tree_2<Traits>::erase(Self_iterator p)
{
   CGAL_assertion((*p).is_a_leaf());

   Self_iterator s;
   s = right_sibling(p);
   if (s != this->end())
      set_left_sibling(left_sibling(p),s);
   s = left_sibling(p);
   if (s != this->end())
      set_right_sibling(right_sibling(p),s);
   s = parent(p);
   // if p was the rightmost child of its parent, then set its left
   // sibling as the new rightmost child
   if (rightmost_child(s) == p)
      set_rightmost_child(left_sibling(p),s);
}

template <class Traits>
std::ostream& operator<<(std::ostream& os, const Rotation_tree_2<Traits>& tree)
{
    typename Rotation_tree_2<Traits>::const_iterator it;
    for (it = tree.begin(); it != tree.end(); it++)
         os << *it << " " << std::endl;
    return os;
}

}
