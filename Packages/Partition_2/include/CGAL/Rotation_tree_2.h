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
// file          : include/CGAL/Rotation_tree_2.h
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
// implementation: Rotation tree for vertex visibility graph computation
// ============================================================================

#ifndef  CGAL_ROTATION_TREE_H
#define  CGAL_ROTATION_TREE_H

#include <vector>
#include <CGAL/ch_utils.h>
#include <CGAL/Rotation_tree_node_2.h>

namespace CGAL {

template <class Traits_>
class Rotation_tree_2 : public std::vector< Rotation_tree_node_2<Traits_> >
{
public:
   typedef Traits_                               Traits;
   typedef Rotation_tree_node_2<Traits>          Node;
   typedef typename std::vector<Node>::iterator    Self_iterator;
   typedef typename Traits::Point_2              Point_2;


   // constructor
   template<class ForwardIterator>
   Rotation_tree_2(ForwardIterator first, ForwardIterator beyond)
   {
      typedef typename Traits::R                               R;
      typedef typename Traits::R::FT                           FT;
      typedef typename Traits::Less_xy_2                       Less_xy_2;
      typedef ch_Binary_predicate_reversor<Point_2, Less_xy_2> Greater_xy_2;
   
      for (ForwardIterator it = first; it != beyond; it++)
         push_back(*it);
   
      std::sort(begin(), end(), Greater_xy_2(Traits().less_xy_2_object()));
      std::unique(begin(), end());
   
      // b is the point with the largest x coordinate
      Node largest_x = front();
   
      // push the point p_minus_infinity
      push_back(Point_2( CGAL::to_double(largest_x.x())+1,
                         -CGAL::to_double(largest_x.y())));

      // push the point p_infinity
      push_back(Point_2(CGAL::to_double(largest_x.x())+1,
                         CGAL::to_double(largest_x.y())));
   
      _p_inf = end();  // record the iterators to these extreme points
      _p_inf--;
      _p_minus_inf = _p_inf;
      _p_minus_inf--;
   
      Self_iterator child;
      // make p_minus_inf a child of p_inf
      set_rightmost_child(_p_minus_inf, _p_inf); 
      child = begin();               // now points to p_0
      while (child != _p_minus_inf)  // make all points children of p_minus_inf
      {
         set_rightmost_child(child, _p_minus_inf);
         child++;
      }
   }


   // the point that comes first in the right-to-left ordering is first
   // in the ordering, after the auxilliary points p_minus_inf and p_inf
   Self_iterator rightmost_point_ref() 
   {
      return begin();
   }

   Self_iterator right_sibling(Self_iterator p) 
   {
      if (!(*p).has_right_sibling()) return end();
      return (*p).right_sibling();
   }

   Self_iterator left_sibling(Self_iterator p) 
   {
      if (!(*p).has_left_sibling()) return end();
      return (*p).left_sibling();
   }

   Self_iterator rightmost_child(Self_iterator p) 
   {
      if (!(*p).has_children()) return end();
      return (*p).rightmost_child();
   }

   Self_iterator parent(Self_iterator p) 
   {
      if (!(*p).has_parent()) return end();
      return (*p).parent();
   }

   bool parent_is_p_infinity(Self_iterator p) 
   {
      return parent(p) == _p_inf;
   }

   bool parent_is_p_minus_infinity(Self_iterator p) 
   {
      return parent(p) == _p_minus_inf;
   }

   // makes *p the parent of *q
   void set_parent (Self_iterator p, Self_iterator q)
   {
      CGAL_assertion(q != end());
      if (p == end())
         (*q).clear_parent();
      else
         (*q).set_parent(p);
   }

   // makes *p the rightmost child of *q
   void set_rightmost_child(Self_iterator p, Self_iterator q);

   // makes *p the left sibling of *q
   void set_left_sibling(Self_iterator p, Self_iterator q);

   // makes *p the right sibling of *q
   void set_right_sibling(Self_iterator p, Self_iterator q);

   // NOTE:  this function does not actually remove the node p from the
   //        list; it only reorganizes the pointers so this node is not
   //        in the tree structure anymore
   void erase(Self_iterator p);

private:
   Self_iterator _p_inf;
   Self_iterator _p_minus_inf;
};

}

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Rotation_tree_2.C>
#endif // CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION


#endif // CGAL_ROTATION_TREE_H
