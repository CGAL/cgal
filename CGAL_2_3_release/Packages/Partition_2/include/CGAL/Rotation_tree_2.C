// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : include/CGAL/Rotation_tree_2.C
// package       : $CGAL_Package: Partition_2 $
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
// chapter       : Planar Polygon Partitioning
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: Rotation tree for visibility graph computation
// ============================================================================

#include <iostream>
#include <CGAL/ch_utils.h>

namespace CGAL {


// makes *p the rightmost child of *q
template<class Traits>
void Rotation_tree_2<Traits>::set_rightmost_child(Self_iterator p, 
                                                  Self_iterator q)
{
   CGAL_assertion(q != end());

   if (p != end())
   {
      (*p).clear_right_sibling();
      if (rightmost_child(q) != end())
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
   if (q == end()) return;
       
   if (p != end())
   {
      if (left_sibling(q) != end())
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
      if (left_sibling(q) != end())
         (*(*q).left_sibling()).clear_right_sibling();
      (*q).clear_left_sibling();
   }
}

// makes p the right sibling of q
template <class Traits>
void Rotation_tree_2<Traits>::set_right_sibling(Self_iterator p, 
                                                Self_iterator q)
{
   if (q == end()) return;
       
   if (p != end())
   {
      if (right_sibling(q) != end())
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
      if (right_sibling(q) != end())
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
   if (s != end())
      set_left_sibling(left_sibling(p),s);
   s = left_sibling(p);
   if (s != end())
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
