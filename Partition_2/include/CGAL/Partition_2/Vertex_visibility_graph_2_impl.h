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

namespace CGAL {


/*
// ??? need to finish this ???
template <class Traits>
template <class ForwardIterator>
bool 
Vertex_visibility_graph_2<Traits>::is_valid(ForwardIterator first, 
                                     ForwardIterator beyond)
{
   std::vector<Point_2> vertices(first, beyond);
   bool edge_there[vertices.size()];

   // for each edge in the graph determine if it is either an edge of the
   // polygon or, if not, if it intersects the polygon in the interior of the
   // edge.
   for (iterator e_it = edges.begin(); e_it != edges.end(); e_it++)
   {
      Segment_2 s = construct_segment_2((*e_it).first, (*e_it).second);
      if (is_an_edge(*e_it))
         edge_there[edge_num] = true;
      else if (do_intersect_in_interior(s, first, beyond))
      return false;
   }
   // check if all the edges of the polygon are present
   //
   // ??? how do you check if there are missing edges ???
}

*/

// want to determine, for each vertex p of the polygon, the line segment
// immediately below it.  For vertical edges, the segment below is not the
// one that begins at the other endpoint of the edge.
template <class Traits>
void 
Vertex_visibility_graph_2<Traits>::initialize_vertex_map(
                              const Polygon& polygon, Vertex_map& vertex_map)
{
   typedef typename Vertex_map::value_type           Map_pair;

   // Create an event list that is a list of circulators for the polygon
   Iterator_list<Polygon_const_iterator>     
                           iterator_list(polygon.begin(), polygon.end());

   // Sort the event list (iterators to points) from left to right 
   // (using less_xy)
   iterator_list.sort(Indirect_less_xy_2<Traits>());
   // Create an ordered list of edge endpoints (iterators), initially empty
   typedef std::set< Point_pair, Segment_less_yx_2 > Ordered_edge_set;
   typedef typename Ordered_edge_set::iterator       Ordered_edge_set_iterator;

   Ordered_edge_set              ordered_edges;
   Ordered_edge_set_iterator     edge_it;
   Vertex_map_iterator   vm_it;
   Vertex_map_iterator   vis_it;

   Polygon_const_iterator event_it;
   Polygon_const_iterator next_endpt;
   Polygon_const_iterator prev_endpt;

   // initialize the map by associating iterators and points and indicating 
   // that no points can see anything.
   for (Polygon_const_iterator it = polygon.begin();it != polygon.end();it++)
   {
      vertex_map.insert(Map_pair(*it, Iterator_pair(it, polygon.end())));
   }
      
   // now go through the events in sorted order.  
   while (!iterator_list.empty())
   {
      event_it = iterator_list.front();
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
      std::cout << "event = " << *event_it << std::endl;     
#endif
      next_endpt = event_it; next_endpt++;
      if (next_endpt == polygon.end()) next_endpt = polygon.begin();
      iterator_list.pop_front();

      // the first edge that is not less than (below) this edge, so ...
      edge_it = ordered_edges.lower_bound(Point_pair(*event_it,*next_endpt));

      // ...if there is no edge below this one then nothing is visible,
      // otherwise....
      if (edge_it != ordered_edges.begin())
      {
         edge_it--; // ...the first visible edge is the previous edge

         // find the event point in the vertex map
         vm_it = vertex_map.find(*event_it);

         // Find the entry for the edge's first endpoint in the vertex map.
         vis_it = vertex_map.find((*edge_it).first);
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         std::cout << "the potential visibility point is " << (*vis_it).first 
                   << endl;
#endif
         // an edge that ends at this event point cannot be below this 
         // endpoint
         if (!is_next_to(polygon, (*vis_it).second.first, event_it))
         {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
            cout << "the edge beginning at  " << *(*vis_it).second.first 
                 << " is visible" << endl;
#endif
            // set the visibility iterator for this point to the iterator
            // corresponding to the edge endpoint that is to the left of
            // the vertical line
            if (less_xy_2((*vis_it).first,  (*vm_it).first))
            {
               Polygon_const_iterator next_vtx = (*vis_it).second.first;
               next_vtx++; 
               if (next_vtx == polygon.end()) next_vtx = polygon.begin();
               (*vm_it).second.second = next_vtx;
            }
            else
               (*vm_it).second.second = (*vis_it).second.first;
         }
         // skip over the edge that ends at this event point. If there 
         // is another edge above this event's edge then it is visible. 
         // since it can't also end at the event point.
         else if (edge_it != ordered_edges.begin() &&
                  --edge_it != ordered_edges.begin())
         {
            vis_it = vertex_map.find((*edge_it).first);
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
            std::cout << "the edge beginning at  " << *(*vis_it).second.first 
                      << " is visible" << endl;
#endif
            // set the visibility iterator for this point to the iterator
            // corresponding to the edge endpoint that is to the left of
            // the vertical line
            if (less_xy_2((*vis_it).first,  (*vm_it).first))
            {
                Polygon_const_iterator next_vtx = (*vis_it).second.first;
                next_vtx++; 
                if (next_vtx == polygon.end()) next_vtx = polygon.begin();
                (*vm_it).second.second = next_vtx;
            }
            else
                (*vm_it).second.second = (*vis_it).second.first;
         }
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         else
            std::cout << "nothing is visible " << endl;
#endif
      }
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
      else
         cout << "nothing is visible " << endl;
#endif
      prev_endpt = event_it; 
      if (prev_endpt == polygon.begin()) 
         prev_endpt = polygon.end();
      prev_endpt--;
      // if the other endpoint of the next edge is to the right of the
      // sweep line, then insert this edge
      if (less_xy_2(*event_it, *next_endpt))
      {
         ordered_edges.insert(Point_pair(*event_it,*next_endpt));
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
             cout << "inserting edge from "
                  << *event_it << " to " << *next_endpt << endl;
#endif
      }
      else // other endpoint not to the right, so erase it 
      {
         ordered_edges.erase(Point_pair(*event_it,*next_endpt));
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         std::cout << "erasing edge from "
                   << *event_it << " to " << *next_endpt << endl;
#endif
      }
      // if the other endpoint of the previous edge is to the right of the
      // sweep line, insert it
      if (less_xy_2(*event_it, *prev_endpt))
      {
         ordered_edges.insert(Point_pair(*prev_endpt,*event_it));
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         cout << "inserting edge from "
              << *prev_endpt << " to " << *event_it << endl;
#endif
       }
       else // other endpoint is not to the right, so erase it
       {
          ordered_edges.erase(Point_pair(*prev_endpt,*event_it));
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
          std::cout << "erasing edge from "
                    << *prev_endpt << " to " << *event_it << endl;
#endif
       }
   }
}

// determines if one makes a left turn going from p to q to q's parent.
// if q's parent is p_infinity, then a left turn is made when p's x value
// is less than q's x value or the x values are the same and p's y value is
// less than q's.
// if p, q, and q's parent are collinear, then one makes a "left turn"
// if q is between p and q's parent (since this means that p can't see 
// q's parent and thus should not become a child of that node)
template <class Traits>
bool 
Vertex_visibility_graph_2<Traits>::left_turn_to_parent(
                                   Tree_iterator p, 
                                   Tree_iterator q, 
                                   Tree& tree)
{
   typedef typename Traits::Point_2 Point;
   if (tree.parent_is_p_infinity(q)) 
   {
      return (less_xy_2(Point(*p), Point(*q)));
   }
   else if (orientation_2(Point(*p), Point(*q), Point(*q->parent())) == COLLINEAR &&
            (collinear_ordered_2(Point(*p), Point(*q), Point(*q->parent())) ||
             collinear_ordered_2(Point(*p), Point(*q), Point(*q->parent()))))
      
   {
      return true;
   }
   else
   {
      return left_turn_2(Point(*p), Point(*q), Point(*q->parent()));
   }
}

// returns true if the diagonal from p to q cuts the interior angle at p

template <class Traits>
bool 
Vertex_visibility_graph_2<Traits>::diagonal_in_interior(
                             const Polygon& polygon, 
                             Polygon_const_iterator p,
                             Polygon_const_iterator q)
{
   Turn_reverser<Point_2, Left_turn_2> right_turn(left_turn_2);
   Polygon_const_iterator before_p;
   if (p == polygon.begin())
      before_p = polygon.end();
   else 
      before_p = p;
   before_p--;
   Polygon_const_iterator after_p = p; after_p++;
   if (after_p == polygon.end()) after_p = polygon.begin();

   if (right_turn(*before_p, *p, *after_p))
   {
      if (right_turn(*before_p, *p, *q) && right_turn(*q, *p, *after_p))
         return false;
   }
   else // left turn or straight at vertex
   {
/*
      // p should not be able to see q through its own edge
      if (are_strictly_ordered_along_line(*p, *after_p, *q))
         return false;
*/
      if (right_turn(*before_p, *p, *q) || right_turn(*q, *p, *after_p))
         return false;
   }
   return true;
}
 

// returns true if the looker can see the point_to_see 
template <class Traits>
bool Vertex_visibility_graph_2<Traits>::point_is_visible(
                                           const Polygon& polygon, 
                                           Polygon_const_iterator point_to_see,
                                           Vertex_map_iterator looker)
{
   // Collect pointers to the current visibility segments for the looker
   // (the current visibility point and the two vertices flanking this vertex)
   Polygon_const_iterator vis_endpt = (*looker).second.second;
   Polygon_const_iterator next_vis_endpt = vis_endpt; next_vis_endpt++;
   if (next_vis_endpt == polygon.end()) next_vis_endpt = polygon.begin();
   Polygon_const_iterator prev_vis_endpt;
   if (vis_endpt == polygon.begin())
      prev_vis_endpt = polygon.end();
   else
      prev_vis_endpt = vis_endpt;
   prev_vis_endpt--;

#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
     cout << "looker is " << (*looker).first << " point to see is " 
          << *point_to_see;
     cout << " visibility points are prev: " << *prev_vis_endpt
          << " vis: " << *vis_endpt << " next: " << *next_vis_endpt << endl;
#endif

    // if the point to see is the current visibility point or if the looker 
    // and the point to see flank the old visibility point, they are visible 
    // to each other since it is known at this point that the edge from 
    // the looker to the point to see goes through the interior of the polygon
    if ((*looker).second.second == point_to_see)
        
    {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
       std::cout << "looker sees point" << std::endl;
#endif
       return true;
    }
    else if (((*looker).second.first == prev_vis_endpt && 
              point_to_see == next_vis_endpt) ||
             ((*looker).second.first == next_vis_endpt && 
              point_to_see == prev_vis_endpt))
    {
       if (orientation_2(*prev_vis_endpt, *vis_endpt, *next_vis_endpt) == 
           COLLINEAR &&
           (collinear_ordered_2((*looker).first, *vis_endpt, *point_to_see) ||
            collinear_ordered_2(*point_to_see, *vis_endpt, (*looker).first)))
       {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
          cout << "looker does NOT see point" << endl;
#endif
          return false;
       }
       else                      
       {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
          cout << "looker sees point" << endl;
#endif
          return true;
       }
    }
    else if ((*looker).second.first == prev_vis_endpt ||
             point_to_see == prev_vis_endpt)
    // point to see or looker is not adjacent to old visibility, so check 
    // intersection with next visibility segment
    {
       if (orientation_2(*vis_endpt, *next_vis_endpt, (*looker).first) !=
           orientation_2(*vis_endpt, *next_vis_endpt, *point_to_see) &&
           orientation_2((*looker).first, *point_to_see, *vis_endpt) !=
           orientation_2((*looker).first, *point_to_see, *next_vis_endpt))
       {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
          cout << "looker does NOT see point" << endl;
#endif
          return false;
       }
       else
       {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
          cout << "looker sees point" << endl;
#endif
         return true;
       }
    }
    else if ((*looker).second.first == next_vis_endpt ||
             point_to_see == next_vis_endpt)
    // point to see or looker is not adjacent to old visibility, so check 
    // intersection with previous visibility segment
    {
       if (orientation_2(*vis_endpt, *prev_vis_endpt, (*looker).first) !=
           orientation_2(*vis_endpt, *prev_vis_endpt, *point_to_see) &&
           orientation_2((*looker).first, *point_to_see, *vis_endpt) !=
           orientation_2((*looker).first, *point_to_see, *prev_vis_endpt))
       {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
          cout << "looker does NOT see point" << endl;
#endif
          return false;
       }
       else
       {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         cout << "looker sees point" << endl;
#endif
         return true;
       }
    }
    else 
    // neither is adjacent to the old visibility point so check intersection
    // with both visibility segments
    {
       if (orientation_2(*vis_endpt, *next_vis_endpt, (*looker).first) !=
           orientation_2(*vis_endpt, *next_vis_endpt, *point_to_see) &&
           orientation_2((*looker).first, *point_to_see, *vis_endpt) !=
           orientation_2((*looker).first, *point_to_see, *next_vis_endpt))
       {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         cout << "looker does NOT see point" << endl;
#endif
         return false;
       }
       else if (orientation_2(*vis_endpt, *prev_vis_endpt, (*looker).first) !=
                orientation_2(*vis_endpt, *prev_vis_endpt, *point_to_see) &&
                orientation_2((*looker).first, *point_to_see, *vis_endpt) !=
                orientation_2((*looker).first, *point_to_see, *prev_vis_endpt))
       {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         cout << "looker does NOT see point" << endl;
#endif
         return false;
       }
       else // no intersection
       {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         cout << "looker sees point" << endl;
#endif
         return true;
       }
   }
}
   
template <class Traits>
void Vertex_visibility_graph_2<Traits>::update_visibility(
                                                      Vertex_map_iterator p_it,
                                                      Vertex_map_iterator q_it,
                                                      const Polygon& polygon,
                                                      int are_adjacent)
{
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
      std::cout << "updating visibility:  " << std::endl;
#endif
   Polygon_const_iterator prev_q; 
   Polygon_const_iterator turn_q; 
   if ((*q_it).second.first == polygon.begin()) 
      prev_q = polygon.end();
   else
      prev_q = (*q_it).second.first; 
   prev_q--;

   // determine if the vertex before or after q is the one that will
   // be encountered next when moving in the direction from p to q.
   if (prev_q == (*p_it).second.first)
   {
      turn_q = (*q_it).second.first; 
      turn_q++;
      if (turn_q == polygon.end()) turn_q = polygon.begin();
   }
   else
      turn_q = prev_q;

#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
   std::cout << "prev_q = " << *prev_q  << " turn_q = " << *turn_q 
             << std::endl;
#endif

   if (are_adjacent)
   {
      if (orientation_2((*p_it).first, (*q_it).first, *turn_q) == RIGHT_TURN)
      {
         (*p_it).second.second = (*q_it).second.second; // p sees what q sees
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         std::cout << "adjacent with right turn; p now sees what q sees" 
                   << std::endl;
#endif
      }
      else // turn left or go straight
      {
         (*p_it).second.second = (*q_it).second.first;  // p sees q
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         std::cout << "adjacent and NOT right turn; p now sees q " 
                   << std::endl;
#endif
      }
   }
   // if not adjacent, the edge was an interior one and so the "next" vertex
   // (the turn vertex) has to be the one that follows q.
   //
   // if Segment(q) == vis(p), i.e., p already sees q's segment
   else if ((*q_it).second.first == (*p_it).second.second ||
            prev_q == (*p_it).second.second)
   {
      turn_q = (*q_it).second.first; turn_q++;
      if (turn_q == polygon.end()) turn_q = polygon.begin();

#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
      std::cout << "prev_q = " << *prev_q  << " turn_q = " << *turn_q 
                << std::endl;
#endif
      // q sees nothing or there is not a right turn to the point after q
      if ((*q_it).second.second == polygon.end() || 
          orientation_2((*p_it).first, (*q_it).first, *turn_q) != RIGHT_TURN)
      {
         (*p_it).second.second = (*q_it).second.first; // p sees q
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         std::cout << "p sees q's segment, q sees nothing and not right to "
                   << " next point; p sees q " << std::endl;
#endif
      }
      else
      {
         (*p_it).second.second = (*q_it).second.second; // p sees what q sees
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         std::cout << "p sees q's segment, q sees something;"
                   << " p sees what q sees" << std::endl;
#endif 
      }
   }
   // Before(p,q,vis(p)) == true if q lies nearer to p than segment vis(p)
   // NOTE:  it is known that p is always looking to the right.
   else  if ((*p_it).second.second != polygon.end())  // if p sees something
   {
      Polygon_const_iterator next_v_p = (*p_it).second.second; next_v_p++;
      if (next_v_p == polygon.end()) next_v_p = polygon.begin();

      // don't need to do this for the previous visibility point since
      // if it were closer to p than q when looking from p to q, q would
      // not be visible.
      Segment_2 next_seg = construct_segment_2(*(*p_it).second.second, 
                                               *next_v_p);
      Ray_2 ray = construct_ray_2((*p_it).first, (*q_it).first);
      Segment_2 i_seg;
      Point_2 i_point;
         
      Object_2 next_result = intersect_2(next_seg, ray);
         
      if (assign_2(i_point, next_result)) 
      {
         if (collinear_ordered_2((*p_it).first, (*q_it).first, i_point))
         {
            (*p_it).second.second = (*q_it).second.first;
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
            std::cout << "p sees something in direction of q, but q is closer;"
                      << " p sees q" << std::endl;
#endif 
         }
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         else 
         {
            std::cout << "p sees something in direction of q that's closer "
                      << "than q; p doesn't see  q" << std::endl;
         }
#endif 
      }
      else if (assign_2(i_seg, next_result))
      {
         if (collinear_ordered_2((*p_it).first,(*q_it).first,i_seg.source()) &&
             collinear_ordered_2((*p_it).first,(*q_it).first,i_seg.target()))
         {
            (*p_it).second.second = (*q_it).second.first;
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
            std::cout << "p sees something in direction of q, but q is closer;"
                      << " p sees q" << std::endl;
#endif 
         }
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         else 
         {
            std::cout << "p sees something in direction of q that's closer "
                      << " than q; p doesn't see  q" << std::endl;
         }
#endif 
      }
      else
      {
         (*p_it).second.second = (*q_it).second.first;
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         std::cout << "p doesn't see something in direction of q; p sees q" 
                   << std::endl;
#endif 
      }
   }
   else // p sees what q sees
   {
      (*p_it).second.second = (*q_it).second.first;
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
      std::cout << "p sees nothing; p sees what q sees" << std::endl;
#endif 
   }
}

template <class Traits>
void Vertex_visibility_graph_2<Traits>::update_collinear_visibility(
                                                    Vertex_map_iterator p_it,
                                                    Vertex_map_iterator q_it,
                                                    const Polygon& polygon)
{
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
   std::cout << "updating collinear visibility" << std::endl;
#endif
   Polygon_const_iterator prev_q; 
   if ((*q_it).second.first == polygon.begin()) 
      prev_q = polygon.end();
   else
      prev_q = (*q_it).second.first;
   prev_q--;

   Polygon_const_iterator next_q = (*q_it).second.first; 
   next_q++;
   if (next_q == polygon.end()) next_q = polygon.begin();

#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
   std::cout << "q's neighbors are: prev " << *prev_q 
             << " q " << (*q_it).first 
             << " next " << *next_q << std::endl;
#endif

   // if the point before q is above the line containing p and q, make
   // this p's visibility point
   if (left_turn_2((*p_it).first, (*q_it).first, *prev_q))
   {
      if (point_is_visible(polygon, prev_q, p_it))
         (*p_it).second.second = prev_q;
   }
   // check the same thing for the point after q and, if it is still visible
   // (even after possibly updating the visibility in the above if) the
   // update again.
   if (left_turn_2((*p_it).first, (*q_it).first, *next_q))
   {
      if (point_is_visible(polygon, next_q, p_it))
         (*p_it).second.second = next_q;
   }
}

   // The segment between points p and q is a potential visibility edge
   // This function determines if the edge should be added or not (based
   // on p's current visibility point) and updates p's visibility point
   // where appropriate
template <class Traits>
void Vertex_visibility_graph_2<Traits>::handle(Tree_iterator p, 
                                        Tree_iterator q, 
                                        const Polygon& polygon,
                                        Vertex_map& vertex_map)
{
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
      std::cout << "Handling edge from " << (*p).x() << " " << (*p).y() 
                << " to " << (*q).x() << " " << (*q).y() << std::endl;
#endif
   Vertex_map_iterator p_it = vertex_map.find(*p);
   Vertex_map_iterator q_it = vertex_map.find(*q);
   CGAL_assertion (p_it != vertex_map.end());
   CGAL_assertion (q_it != vertex_map.end());
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
   std::cout << "p currently sees : ";
   if ((*p_it).second.second != polygon.end())
      std::cout << *((*p_it).second.second) << endl;
   else
      std::cout << " NADA" << endl;
#endif

   // if p and q are adjacent
   if (are_adjacent(polygon, (*p_it).second.first, (*q_it).second.first))
   {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
      cout << "are adjacent" << endl;
#endif
      insert_edge(Point_pair(*p,*q));
      update_visibility(p_it, q_it, polygon, 1);
   }
   else 
   {
      bool interior_at_p = diagonal_in_interior(polygon, (*p_it).second.first,
                                                (*q_it).second.first);
      bool interior_at_q = diagonal_in_interior(polygon, (*q_it).second.first,
                                                (*p_it).second.first);
      // line of site is through the interior of the polygon
      if (interior_at_p && interior_at_q)
      {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         cout << "both interior" << endl;
#endif
         // if p sees something and q is visible only through collinear
         // points then update p's visibility if one of the points adjacent
         // to q is above the line unless p's current visibility point 
         // obscures the view.
         if ((*p_it).second.second != polygon.end() &&
             are_strictly_ordered_along_line_2((*p_it).first,
                                             *(*p_it).second.second,
                                                (*q_it).first))
         {
            update_collinear_visibility(p_it, q_it, polygon);
         }
         // p current sees nothing or q is visible to p
         else if ((*p_it).second.second == polygon.end() ||
                  point_is_visible(polygon, (*q_it).second.first, p_it))
         {
            insert_edge(Point_pair(*p,*q));
            update_visibility(p_it, q_it, polygon, 0);
         }
      }
      else if (!interior_at_p && !interior_at_q) // both points exterior
      {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
            cout << "both exterior" << endl;
#endif
         // p currently sees nothing or q is visible to p
         if ((*p_it).second.second == polygon.end() ||
             point_is_visible(polygon, (*q_it).second.first, p_it))
         {
            (*p_it).second.second = (*q_it).second.first;
         }
      }
   }
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
   std::cout << "p now sees : ";
   if ((*p_it).second.second != polygon.end())
      std::cout << *((*p_it).second.second) << endl;
   else
      std::cout << " NADA" << endl;
#endif
}

}
