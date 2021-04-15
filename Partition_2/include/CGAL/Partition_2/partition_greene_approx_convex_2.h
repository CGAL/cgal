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

#ifndef CGAL_GREENE_APPROX_H
#define CGAL_GREENE_APPROX_H

#include <CGAL/license/Partition_2.h>


#include<vector>
#include<algorithm>
#include<iterator>
#include<CGAL/Partition_2/Circulator_pair.h>
#include<CGAL/Partition_2/partition_y_monotone_2.h>
#include<CGAL/Partition_2/Turn_reverser.h>
#include<CGAL/IO/Tee_for_output_iterator.h>
#include<CGAL/Partition_2/partition_assertions.h>
#include<CGAL/partition_is_valid_2.h>
#include<CGAL/Partition_traits_2.h>
#include<CGAL/is_y_monotone_2.h>
#include<CGAL/Partition_2/is_degenerate_polygon_2.h>

// These things should be constant:
//   front is where you add things to a chain
//   top chain shares its back with front (top) of stack
//   bottom chain shares its back with the back (bottom) of the stack
//   bottom of the stack has lower y value than top.
// ==> chain direction is cw, make ccw polygon from front to back
// and chain direction is ccw, make ccw polygon from back to front

namespace CGAL {

template <class BidirectionalCirculator>
bool is_adjacent_to(BidirectionalCirculator new_point_ref,
                    BidirectionalCirculator old_point_ref)
{
#ifdef CGAL_GREENE_APPROX_DEBUG
   std::cout << "is_adjacent_to: new_point " << *new_point_ref
             <<             " old_point " << *old_point_ref << std::endl;
#endif
   // find the old point in original list of points
   if (*new_point_ref == *(++old_point_ref)) return true;  // check ccw
#ifdef CGAL_GREENE_APPROX_DEBUG
   std::cout << "is_adjacent_to:  next point is " << *old_point_ref
             << std::endl;
#endif
   old_point_ref--;
#ifdef CGAL_GREENE_APPROX_DEBUG
   std::cout << "is_adjacent_to:  old_point is " << *old_point_ref
             << std::endl;
#endif
   if (*new_point_ref == *(--old_point_ref)) return true;  // check cw
#ifdef CGAL_GREENE_APPROX_DEBUG
   std::cout << "is_adjacent_to:  previous points is " << *old_point_ref
             << std::endl;
#endif
   return false;
}


// erases vertices in the range [first, last) (counterclockwise)
template<class BidirectionalCirculator, class Polygon>
void erase_vertices(BidirectionalCirculator first,
                    BidirectionalCirculator last,
                    Polygon& polygon,
                    bool& update_required)
{
   typedef typename Polygon::iterator  Vertex_iterator;
   Vertex_iterator it = polygon.begin();
   Vertex_iterator temp_it;
   update_required = false;

#ifdef CGAL_GREENE_APPROX_DEBUG
   std::cout << "before erasing polygon is ";
   for (Vertex_iterator v = polygon.begin(); v != polygon.end(); v++)
   {
      std::cout << " " << *v;
   }
   std::cout << std::endl;
   std::cout << "erasing from " << *first << " to " << *last << std::endl;
#endif
   it = first.current_iterator();

   CGAL_partition_assertion (it != polygon.end());

   while ( (it != polygon.end()) && (*it != *last) )
   {
      while ( (it != polygon.end()) && (*it != *last) )
      {
#ifdef CGAL_GREENE_APPROX_DEBUG
         std::cout << "erase_vertices: erasing " << *it << std::endl;
#endif
         temp_it = it;
         // when the beginning vertex of the polygon is deleted, all the
         // circulators for that polygon become invalid and must be updated.
         if (it == polygon.begin())
            update_required = true;
         it++;
#ifdef CGAL_GREENE_APPROX_DEBUG
         if (it != polygon.end())
            std::cout << "erase_vertices: next vertex is " << *it << std::endl;
#endif
         polygon.erase(temp_it);
      }
      if (it == polygon.end())
         it = polygon.begin();
   }
   // Yes, both loops here are necessary in case the range of vertices
   // includes the last and first points.
}


template <class Polygon, class BidirectionalCirculator,
          class OutputIterator, class Traits>
void visible(Polygon& polygon,
             const BidirectionalCirculator& new_point_ref,
             Circ_pair< BidirectionalCirculator >& stack,
             Circ_pair< BidirectionalCirculator >& bottom_chain,
             Circ_pair< BidirectionalCirculator >& top_chain,
             OutputIterator&  result,
             const Traits& traits)
{

#ifdef CGAL_GREENE_APPROX_DEBUG
   std::cout << "visible: stack.back " << *stack.back()
             << " stack.before_back " << *stack.before_back() <<
                " new_point_ref " << *new_point_ref << std::endl;
#endif

    typedef typename Traits::Point_2      Point_2;
    typedef typename Traits::Left_turn_2   Left_turn_2;
    Left_turn_2    left_turn = traits.left_turn_2_object();
    Turn_reverser<Point_2, Left_turn_2> right_turn(left_turn);

    if (((bottom_chain.direction() == CLOCKWISE) &&
        right_turn(*stack.back(), *stack.before_back(), *new_point_ref)) ||
        ((bottom_chain.direction() == COUNTERCLOCKWISE) &&
        left_turn(*stack.back(), *stack.before_back(), *new_point_ref)))
    {
       typedef typename Traits::Polygon_2 new_Polygon_2;
       new_Polygon_2 new_polygon;

       bool done = false;
       bool big_angle_at_stack = false;
       bool is_visible = false;
       bool update_required;

#ifdef CGAL_GREENE_APPROX_DEBUG
       std::cout << "visible: bottom_chain.front " << *bottom_chain.front()
                 << std::endl;
#endif
       do
       {
          new_Polygon_2 new_polygon;
          stack.pop_back();
#ifdef CGAL_GREENE_APPROX_DEBUG
          std::cout << "visible: stack.back " << *stack.back() << std::endl;
#endif
          if (bottom_chain.direction() == CLOCKWISE)
          {
             std::copy(bottom_chain.front(), stack.before_back(),
                       std::back_inserter(new_polygon));
             erase_vertices(bottom_chain.before_front(), stack.back(), polygon,
                            update_required);
          }
          else
          {
             std::copy(stack.back(), bottom_chain.after_front(),
                       std::back_inserter(new_polygon));
             erase_vertices(bottom_chain.back(), bottom_chain.front(), polygon,
                            update_required);
          }
          if (!is_degenerate_polygon_2(new_polygon.vertices_begin(),
                                       new_polygon.vertices_end(), traits))
          {
             *result = new_polygon;
             result++;
          }
          bottom_chain.push_back(stack.back());
          if (stack.back() == stack.front())   // form new stack with previous
          {                                    // point and old stack top
                                               // (still on stack)
               done = true;
               typename Traits::Less_yx_2  less_yx = traits.less_yx_2_object();
               if (less_yx(*(stack.front()),*bottom_chain.front()))
               {
#ifdef CGAL_GREENE_APPROX_DEBUG
                  std::cout << "visible:  reversing stack and swapping chains "
                            << std::endl;
#endif
                  stack.push_front(bottom_chain.front());
                                               // reverse stack direction
                  stack.change_dir();
                  bottom_chain = top_chain;         // swap chains
                  top_chain.initialize(stack.front());
                  top_chain.change_dir();
#ifdef CGAL_GREENE_APPROX_DEBUG
                  std::cout << "visible:  stack is now " << *stack.back()
                            << " " << *stack.front()
                            << " dir " << int(stack.direction())
                            << std::endl;
                  std::cout << "visible:  bottom_chain is now "
                            << *bottom_chain.back() << " "
                            << *bottom_chain.front()
                            << " dir " << int(bottom_chain.direction())
                            << std::endl;
                  std::cout << "visible:  top_chain is now "
                            << *top_chain.back() << " " << *top_chain.front()
                            << " dir " << int(top_chain.direction())
                            << std::endl;
#endif
               }
               else
               {
                  stack.push_back(bottom_chain.front());
                  bottom_chain.push_back(bottom_chain.front());
               }
            }
            else // stack size must be >= 2 here
            {
#ifdef CGAL_GREENE_APPROX_DEBUG
               std::cout << "visible: stack.before_back "
                         << *stack.before_back() << std::endl;
#endif
               // angle at stack > 180
               if (bottom_chain.direction() == CLOCKWISE)
                  big_angle_at_stack = right_turn(*bottom_chain.front(),
                                                 *stack.back(),
                                                 *stack.before_back());
               else
                  big_angle_at_stack = left_turn(*bottom_chain.front(),
                                                        *stack.back(),
                                                        *stack.before_back());
               if (bottom_chain.direction() == CLOCKWISE)
                  is_visible = !right_turn(*stack.back(), *stack.before_back(),
                                          *new_point_ref);
               else
                  is_visible = !left_turn(*stack.back(), *stack.before_back(),
                                         *new_point_ref);
               // point can see stack bottom
            }
        } while (!done && !big_angle_at_stack && !is_visible);
        if (big_angle_at_stack)  // previous point is placed on bottom of stack
        {
           stack.push_back(bottom_chain.front());
           bottom_chain.push_back(bottom_chain.front());
        }
    }
}


template <class Polygon, class BidirectionalCirculator,
          class OutputIterator, class Traits>
void stack_extend(Polygon& polygon,
                  BidirectionalCirculator& point_ref,
                  Circ_pair< BidirectionalCirculator >& stack,
                  Circ_pair< BidirectionalCirculator >& top_chain,
                  OutputIterator& result, const Traits& traits)
{
#ifdef CGAL_GREENE_APPROX_DEBUG
   std::cout << "stack_extend" << std::endl;
   std::cout << "stack_extend:  stack.before_front() " << *stack.before_front()
    << " stack.front " << *stack.front() << " point_ref " << *point_ref
    << std::endl;
#endif

   typedef typename Traits::Point_2      Point_2;
   typedef typename Traits::Left_turn_2   Left_turn_2;
   Left_turn_2    left_turn = traits.left_turn_2_object();
   Turn_reverser<Point_2, Left_turn_2> right_turn(left_turn);
   if (((stack.direction() == COUNTERCLOCKWISE) &&
        right_turn(*stack.before_front(), *stack.front(), *point_ref)) ||
       ((stack.direction() == CLOCKWISE) &&
        left_turn(*stack.before_front(), *stack.front(), *point_ref)))
   {
      // new stack top becomes new first (and only) element of top chain
      stack.push_front(point_ref);
      top_chain.initialize(point_ref);
#ifdef CGAL_GREENE_APPROX_DEBUG
      std::cout << "stack_extend: top_chain.front "
                << *top_chain.front() << std::endl;
#endif
   }
   else
      change_top_chain(polygon, point_ref, stack, top_chain, result, traits);
}

template <class Polygon, class BidirectionalCirculator,
          class OutputIterator, class Traits>
void change_top_chain(Polygon& polygon,
                      BidirectionalCirculator new_point_ref,
                      Circ_pair< BidirectionalCirculator >& stack,
                      Circ_pair< BidirectionalCirculator >& top_chain,
                      OutputIterator& result, const Traits& traits)
{
   BidirectionalCirculator next_point_ref = new_point_ref;
   if (top_chain.direction() == COUNTERCLOCKWISE)
      next_point_ref++;
   else
      next_point_ref--;

#ifdef CGAL_GREENE_APPROX_DEBUG
   std::cout << "change_top_chain: top_chain.front " << *top_chain.front()
             << " new_point_ref " << *new_point_ref
             << " next_point_ref " << *next_point_ref << std::endl;
#endif
   typedef typename Traits::Point_2      Point_2;
   typedef typename Traits::Left_turn_2   Left_turn_2;
   Left_turn_2    left_turn = traits.left_turn_2_object();
   Turn_reverser<Point_2, Left_turn_2> right_turn(left_turn);
   if (((top_chain.direction() == COUNTERCLOCKWISE) &&
       !right_turn(*top_chain.front(), *new_point_ref, *next_point_ref)) ||
      ((top_chain.direction() == CLOCKWISE) &&
       !left_turn(*top_chain.front(), *new_point_ref, *next_point_ref)))
   {
      top_chain.push_front(new_point_ref);
   }
   else
   {
      typedef typename Traits::Polygon_2 new_Polygon_2;

      BidirectionalCirculator old_top_ref;
      bool done = false;
      bool big_angle_at_stack = false;
      bool small_angle_at_point = false;
      bool update_required;

      // pop off element already on top chain
      old_top_ref = stack.front();
      stack.pop_front();

      do
      {
         new_Polygon_2 new_polygon;
#ifdef CGAL_GREENE_APPROX_DEBUG
         std::cout << "change_top_chain: stack.front "
                   << *stack.front() << std::endl;
#endif
         if (top_chain.direction() == COUNTERCLOCKWISE)
         {
            std::copy(stack.front(), next_point_ref,
                      std::back_inserter(new_polygon));
            erase_vertices(old_top_ref, new_point_ref, polygon,
                           update_required);
         }
         else
         {
            std::copy(new_point_ref, stack.before_front(),
                      std::back_inserter(new_polygon));
            erase_vertices(top_chain.front(), stack.front(), polygon,
                           update_required);
            top_chain.push_front(stack.front());
         }
         if (!is_degenerate_polygon_2(new_polygon.vertices_begin(),
                                      new_polygon.vertices_end(), traits))
         {
            *result = new_polygon;
            result++;
         }
         if (stack.front() == stack.back())          // the "stack empty" case
         {
            done = true;
            stack.push_front(new_point_ref);
#ifdef CGAL_GREENE_APPROX_DEBUG
            std::cout << "change_top_chain: stack is empty. New stack top is "
                      << *stack.front() << std::endl;
            std::cout << "stack.back is " << *stack.back() << std::endl;
            std::cout << "stack.before_back is " << *stack.before_back()
                      << std::endl;
#endif
            top_chain.initialize(new_point_ref);
            // top chain should share top of stack
         }
         else    // stack size must be >= 2 here
         {
#ifdef CGAL_GREENE_APPROX_DEBUG
            std::cout << "change_top_chain: stack.front "
                      << *stack.front() << std::endl;
#endif
            if (top_chain.direction() == COUNTERCLOCKWISE)
               big_angle_at_stack = !left_turn(*stack.before_front(),
                                              *stack.front(), *new_point_ref);
            else
               big_angle_at_stack = !right_turn(*stack.before_front(),
                                             *stack.front(), *new_point_ref);
            if (top_chain.direction() == COUNTERCLOCKWISE)
               small_angle_at_point = left_turn(*stack.front(), *new_point_ref,
                                               *next_point_ref);
            else
               small_angle_at_point = right_turn(*stack.front(),
                                                 *new_point_ref,
                                                 *next_point_ref);
            if (!big_angle_at_stack && !small_angle_at_point)
            {
               old_top_ref = stack.front();
               stack.pop_front();
            }
         }
      } while (!done && !big_angle_at_stack && !small_angle_at_point);
      if (big_angle_at_stack)  // promote point to stack
      {
         stack.push_front(new_point_ref);
         top_chain.initialize(new_point_ref);
         // top chain shares top of stack
      }
      else if (small_angle_at_point)
      {
         top_chain.push_back(stack.front());
         top_chain.push_front(new_point_ref);
      }
   }
}

template <class Polygon, class BidirectionalCirculator, class OutputIterator,
          class Traits>
void change_bottom_chain(Polygon& polygon,
                         BidirectionalCirculator new_point_ref,
                         Circ_pair< BidirectionalCirculator >& stack,
                         Circ_pair< BidirectionalCirculator >& bottom_chain,
                         Circ_pair< BidirectionalCirculator >& top_chain,
                         OutputIterator& result, const Traits& traits)
{
   BidirectionalCirculator next_point_ref = new_point_ref;
   if (bottom_chain.direction() == COUNTERCLOCKWISE)
      next_point_ref++;
   else
      next_point_ref--;
#ifdef CGAL_GREENE_APPROX_DEBUG
   std::cout << "change_bottom_chain: bottom_chain.front "
             << *bottom_chain.front() << " new_point_ref " << *new_point_ref
             << " next_point_ref " << *next_point_ref << std::endl;
#endif
   typedef typename Traits::Point_2      Point_2;
   typedef typename Traits::Left_turn_2   Left_turn_2;
   Left_turn_2    left_turn = traits.left_turn_2_object();
   Turn_reverser<Point_2, Left_turn_2> right_turn(left_turn);
   if (((bottom_chain.direction() == CLOCKWISE) &&
        !left_turn(*bottom_chain.front(), *new_point_ref,
                          *next_point_ref)) ||
       ((bottom_chain.direction() == COUNTERCLOCKWISE) &&
        !right_turn(*bottom_chain.front(), *new_point_ref, *next_point_ref)))
   {
      bottom_chain.push_front(new_point_ref);
   }
   else
   {
      typedef typename Traits::Polygon_2 new_Polygon_2;

      bool done = false;
      bool small_angle_at_point = false;
      bool update_required;

      do
      {
         new_Polygon_2 new_polygon;
         stack.pop_back();  // remove old bottom of stack
#ifdef CGAL_GREENE_APPROX_DEBUG
         std::cout << "change_bottom_chain: stack.front " << *stack.front()
                   << " stack.back " << *stack.back() << std::endl;
#endif
         if (bottom_chain.direction() == CLOCKWISE)
         {
            std::copy(new_point_ref, stack.before_back(),
                      std::back_inserter(new_polygon));
            erase_vertices(bottom_chain.front(), stack.back(), polygon,
                           update_required);
         }
         else
         {
            std::copy(stack.back(), next_point_ref,
                      std::back_inserter(new_polygon));
            erase_vertices(bottom_chain.back(), new_point_ref, polygon,
                           update_required);
         }
         if (!is_degenerate_polygon_2(new_polygon.vertices_begin(),
                                      new_polygon.vertices_end(), traits))
         {
            *result = new_polygon;
            result++;
         }
         bottom_chain.initialize(stack.back());
         if (stack.back() == stack.front())   // form new stack with new point
         {                                // and old stack top (still on stack)
            done = true;
            typename Traits::Less_yx_2  less_yx = traits.less_yx_2_object();
            if (less_yx(*(stack.front()),*new_point_ref))
            {
#ifdef CGAL_GREENE_APPROX_DEBUG
               std::cout << "change_bottom_chain: reversing stack and "
                         << "swapping chains" << std::endl;
#endif
               stack.push_front(new_point_ref);  // reverse stack direction
               stack.change_dir();
               bottom_chain = top_chain;     // swap the chains
               top_chain.initialize(stack.front());
               top_chain.change_dir();
#ifdef CGAL_GREENE_APPROX_DEBUG
               std::cout << "change_bottom_chain:  stack is now "
                    << *stack.back() << " " << *stack.front()
                    << " dir " << int(stack.direction()) << std::endl;
               std::cout << "change_bottom_chain:  bottom_chain is now "
                    << *bottom_chain.back() << " " << *bottom_chain.front()
                    << " dir " << int(bottom_chain.direction()) << std::endl;
               std::cout << "change_bottom_chain:  top_chain is now "
                    << *top_chain.back() << " " << *top_chain.front()
                    << " dir " << int(top_chain.direction()) << std::endl;
#endif
            }
            else
               stack.push_back(new_point_ref);
         }
         else    // stack size must be >= 2 here
         {
            // angle at new point is < 180
            if (bottom_chain.direction() == CLOCKWISE)
               small_angle_at_point = right_turn(*stack.back(), *new_point_ref,
                                                *next_point_ref);
            else
               small_angle_at_point = left_turn(*stack.back(), *new_point_ref,
                                               *next_point_ref);
         }
      } while (!done && !small_angle_at_point);

      if (small_angle_at_point)
      {
         bottom_chain.push_back(stack.back());
         bottom_chain.push_front(new_point_ref);
      }
   }
}

template <class Polygon, class BidirectionalCirculator,
          class OutputIterator, class Traits>
void make_polygons_from_stack(Polygon& polygon,
                            const BidirectionalCirculator& high_point_ref,
                            Circ_pair< BidirectionalCirculator >& stack,
                            Circ_pair< BidirectionalCirculator >& bottom_chain,
                            OutputIterator& result, const Traits& traits)
{
   bool update_required;
   // make polygons by connecting the high point to every point on the stack
   // except the first and the last.
#ifdef CGAL_GREENE_APPROX_DEBUG
   std::cout << "make_polygons_from_stack: high_point_ref " << *high_point_ref
        << " stack.back " << *stack.back()
        << " stack.before_back " << *stack.before_back() << std::endl;
#endif
   typedef typename Traits::Polygon_2 new_Polygon_2;
   new_Polygon_2 new_polygon;
   BidirectionalCirculator next_point_ref = high_point_ref;
   if (bottom_chain.direction() == COUNTERCLOCKWISE)
      next_point_ref++;

   stack.pop_back();
   while (stack.front() != stack.back())
   {
       new_Polygon_2 new_polygon;
#ifdef CGAL_GREENE_APPROX_DEBUG
       std::cout << "make_polygons_from_stack: stack.back " << *stack.back()
                 << std::endl;
#endif
       if (bottom_chain.direction() == CLOCKWISE)
       {
         std::copy(high_point_ref, stack.before_back(),
                   std::back_inserter(new_polygon));
         erase_vertices(bottom_chain.front(), stack.back(), polygon,
                        update_required);
         bottom_chain.initialize(stack.back());
       }
       else
       {
         std::copy(stack.back(), next_point_ref,
                   std::back_inserter(new_polygon));
         erase_vertices(bottom_chain.back(), high_point_ref, polygon,
                        update_required);
         bottom_chain.push_back(stack.back());
       }
       if (!is_degenerate_polygon_2(new_polygon.vertices_begin(),
                                    new_polygon.vertices_end(), traits))
       {
          *result = new_polygon;
          result++;
       }
       stack.pop_back();
   }
   // add remaining points from the top chain if there is more than one
   std::copy(polygon.begin(), polygon.end(), std::back_inserter(new_polygon));
   if (!is_degenerate_polygon_2(new_polygon.vertices_begin(),
                                new_polygon.vertices_end(), traits))
   {
      *result = new_polygon;
      result++;
   }
}

template<class BidirectionalCirculator, class Traits>
void find_smallest_yx(BidirectionalCirculator& first, const Traits& traits)
{
   BidirectionalCirculator current = first;
   current++;
   // find out which direction to go
   typename Traits::Less_yx_2     less_yx = traits.less_yx_2_object();
   if (less_yx(*current, *first))   // go foward
   {
      do
      {
         first++;
         current++;
      } while (less_yx(*current, *first));
   }
   else
   {
      current = first;
      current--;
      if (less_yx(*current, *first))    // go backward
      {
        do
        {
          first--;
          current--;
        } while (less_yx(*current, *first));
      }
   } // otherwise both previous and next are larger, so no need to move
}

template <class BidirectionalCirculator, class Traits>
bool second_point_is_next(BidirectionalCirculator first, const Traits& traits)
{
   BidirectionalCirculator next = first;
   next++;
   BidirectionalCirculator previous = first;
   previous--;
   typename Traits::Less_yx_2        less_yx = traits.less_yx_2_object();
   return less_yx(*next, *previous);
}

template <class BidirectionalCirculator, class Traits>
BidirectionalCirculator next_vertex(BidirectionalCirculator& ccw_current,
                                    BidirectionalCirculator& cw_current,
                                    const Traits& traits)
{
  BidirectionalCirculator ccw_next = ccw_current;
  ccw_next++;
  BidirectionalCirculator cw_next = cw_current;
  cw_next--;
#ifdef CGAL_GREENE_APPROX_DEBUG
  std::cout << "next_vertex: ccw_next " << *ccw_next << " cw_next " << *cw_next
            << std::endl;
#endif

  if (ccw_next == cw_next)
  {
     cw_current = cw_next;
     ccw_current = ccw_next;
     return cw_next;
  }
  else
  {
     typename Traits::Less_yx_2     less_yx = traits.less_yx_2_object();
     if (less_yx(*ccw_next,*cw_next))
     {
        ccw_current = ccw_next;
        return ccw_next;
     }
     else
     {
        cw_current = cw_next;
        return cw_next;
     }
  }
}


template <class ForwardIterator, class OutputIterator, class Traits>
void ga_convex_decomposition(ForwardIterator first, ForwardIterator beyond,
                             OutputIterator& result, const Traits& traits)
{
   if (first == beyond) return;

   typedef typename Traits::Point_2                        Point_2;
   typedef std::list<Point_2>                              Vertex_list;
   typedef Circulator_from_container<Vertex_list>         Vertex_circulator;

   Vertex_list  polygon(first, beyond);

   CGAL_partition_precondition(
    orientation_2(polygon.begin(), polygon.end(), traits) == COUNTERCLOCKWISE);
   CGAL_partition_precondition(
    is_y_monotone_2(polygon.begin(), polygon.end(), traits));

   Vertex_circulator point_ref(&polygon);
   Vertex_circulator circ = point_ref;

#ifdef CGAL_GREENE_APPROX_DEBUG
   std::cout << "before find_smallest_yx " << std::endl;
   do
   {
      std::cout << *circ << std::endl;
   } while (++circ != point_ref);
#endif

   find_smallest_yx(point_ref, traits);

   circ = point_ref;
#ifdef CGAL_GREENE_APPROX_DEBUG
   std::cout << "after find_smallest_yx " << std::endl;
   do
   {
      std::cout << *circ << std::endl;
   } while (++circ != point_ref);
#endif

   Vertex_circulator ccw_chain_ref = point_ref;
   Vertex_circulator cw_chain_ref = point_ref;

   Circ_pair< Vertex_circulator > stack(point_ref, COUNTERCLOCKWISE);
   Circ_pair< Vertex_circulator > bottom_chain(point_ref, CLOCKWISE);
   Circ_pair< Vertex_circulator > top_chain(point_ref, COUNTERCLOCKWISE);

   // discover which is the direction to the next point
   if (second_point_is_next(point_ref, traits))
   {
     ccw_chain_ref++;
     stack.push_front(ccw_chain_ref);
     stack.set_direction(COUNTERCLOCKWISE);
     top_chain.initialize(ccw_chain_ref);
     top_chain.set_direction(COUNTERCLOCKWISE);
     bottom_chain.set_direction(CLOCKWISE);
   }
   else
   {
     cw_chain_ref--;
     stack.push_front(cw_chain_ref);
     stack.set_direction(CLOCKWISE);
     top_chain.initialize(cw_chain_ref);
     top_chain.set_direction(CLOCKWISE);
     bottom_chain.set_direction(COUNTERCLOCKWISE);
   }
#ifdef CGAL_GREENE_APPROX_DEBUG
   std::cout << "after inserting first two points: ccw_chain_ref "
             << *ccw_chain_ref << " cw_chain_ref " << *cw_chain_ref
             << std::endl;
#endif

   while (ccw_chain_ref != cw_chain_ref)
   {
      point_ref = next_vertex(ccw_chain_ref, cw_chain_ref, traits);

#ifdef CGAL_GREENE_APPROX_DEBUG
      std::cout << "after next_vertex: ccw_chain_ref " << *ccw_chain_ref
                << " cw_chain_ref " << *cw_chain_ref << std::endl;
      std::cout << "current point: " << *point_ref << std::endl;
#endif
      if (is_adjacent_to(point_ref, bottom_chain.front()))
         visible(polygon, point_ref, stack, bottom_chain, top_chain,
                 result, traits);
      if (ccw_chain_ref == cw_chain_ref)
      {
         make_polygons_from_stack(polygon, point_ref, stack,
                                  bottom_chain, result, traits);
         return;
      }
      if (is_adjacent_to(point_ref, stack.front()))
         stack_extend(polygon, point_ref, stack, top_chain,
                      result, traits);
      else if (is_adjacent_to(point_ref, top_chain.front()))
         change_top_chain(polygon, point_ref, stack, top_chain,
                          result, traits);
      else
         change_bottom_chain(polygon, point_ref, stack, bottom_chain,
                             top_chain, result, traits);
   }
}



// should change this so it is like the others in that it adds diagonals
// to a partition polygon and then calls partition
//   have to be able to go over the vertices again in order to check the
//   result, so they have to be stored somewhere else.
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator partition_greene_approx_convex_2(InputIterator first,
                                                InputIterator beyond,
                                                OutputIterator result,
                                                const Traits& traits)
{
   if (first == beyond) return result;

   typedef typename Traits::Polygon_2                        Polygon_2;

#if defined(CGAL_PARTITION_NO_POSTCONDITIONS) || \
    defined(CGAL_NO_POSTCONDITIONS) || defined(NDEBUG)
   OutputIterator res(result);
#else

   Tee_for_output_iterator<OutputIterator, Polygon_2>    res(result);
#endif // no postconditions

   Polygon_2  polygon(first, beyond);
   CGAL_partition_precondition(
      orientation_2(polygon.vertices_begin(), polygon.vertices_end(),
                                               traits) == COUNTERCLOCKWISE);

   std::list<Polygon_2> MP_list;
   typename std::list<Polygon_2>::iterator MP_it;

   y_monotone_partition_2(polygon.vertices_begin(), polygon.vertices_end(),
                          std::back_inserter(MP_list), traits);

   for(MP_it = MP_list.begin(); MP_it != MP_list.end(); MP_it++)
   {
      ga_convex_decomposition((*MP_it).vertices_begin(),
                              (*MP_it).vertices_end(), res, traits);
   }

   CGAL_partition_postcondition(
       convex_partition_is_valid_2(polygon.vertices_begin(),
                                   polygon.vertices_end(),
                                   res.output_so_far_begin(),
                                   res.output_so_far_end(), traits));

#if defined(CGAL_PARTITION_NO_POSTCONDITIONS) || \
    defined(CGAL_NO_POSTCONDITIONS) || defined(NDEBUG)
   return res;
#else
   return res.to_output_iterator();
#endif // no postconditions
}

template <class InputIterator, class OutputIterator>
inline
OutputIterator partition_greene_approx_convex_2(InputIterator first,
                                                InputIterator beyond,
                                                OutputIterator result)
{
   typedef typename std::iterator_traits<InputIterator>::value_type Point_2;
   typedef typename Kernel_traits<Point_2>::Kernel   K;
   return partition_greene_approx_convex_2(first, beyond, result,
                                           Partition_traits_2<K>());
}

}

#endif // CGAL_GREENE_APPROX_H
