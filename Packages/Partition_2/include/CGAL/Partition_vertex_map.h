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
// file          : include/CGAL/Partition_vertex_map.h
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
// implementation: Map used to test for validity of polygon partition
// ============================================================================

#ifndef CGAL_PARTITION_VERTEX_MAP_H
#define CGAL_PARTITION_VERTEX_MAP_H

#include <map>
#include <iostream>
#include <CGAL/circulator.h>
#include <CGAL/Indirect_less_xy_2.h>
#include <cassert>

namespace CGAL {

const int PARTITION_VMAP_UNSHARED_EDGE = -1;

template <class Traits>
class Partition_vertex_map;

template <class Iterator>
class Edge_info;

template <class Iterator, class Traits>
class CW_indirect_edge_info_compare
{
public:
   typedef typename Traits::Leftturn_2              Leftturn_2;
   typedef typename Traits::Less_xy_2               Less_xy_2;
   typedef typename Traits::Point_2                 Point_2;
   typedef Edge_info<Iterator>                      Edge_info;

   CW_indirect_edge_info_compare (Iterator v_it) : vertex_it(v_it),
      left_turn(Traits().leftturn_2_object()),
      less_xy(Traits().less_xy_2_object())
   {}

   bool operator()(Edge_info e1, Edge_info e2)
   {
      bool e1_less = less_xy((*e1.endpoint()), *vertex_it);
      bool e2_less = less_xy((*e2.endpoint()), *vertex_it);
      bool e1_to_e2_left_turn = left_turn((*e1.endpoint()), *vertex_it, 
                                (*e2.endpoint()));

      // if both edges are on the same side of the vertical line through 
      // _vertex then e1 comes before e2 (in CW order from the vertical line)
      // if one  makes a left turn going from e1 to e2
      if (e1_less == e2_less)
          return e1_to_e2_left_turn;
      else // e1 comes first if it is to the right of the vertical line
          return !e1_less;
   }

private:
   Iterator          vertex_it;  
   Leftturn_2        left_turn;
   Less_xy_2         less_xy;
};


template <class Iterator>
class Edge_info
{
public:
   Edge_info() {}
   Edge_info(Iterator e_ref, int p_num1, int p_num2) : _endpoint_ref(e_ref),
      _poly_num1(p_num1), _poly_num2(p_num2)
   { }

   Edge_info(Iterator e_ref, int p_num1) : _endpoint_ref(e_ref),
      _poly_num1(p_num1), _poly_num2(PARTITION_VMAP_UNSHARED_EDGE)
   { }
 
   void set_poly_num1(int p_num)
   {
      _poly_num1 = p_num;
   }

   void set_poly_num2(int p_num)
   {
      _poly_num2 = p_num;
   }

   void set_endpoint(Iterator e_ref)
   {
      _endpoint_ref = e_ref;
   }

   bool same_edge(Iterator e_ref)
   {
      return e_ref == endpoint();
   }

   Iterator endpoint() const { return _endpoint_ref; }
   int poly_num1() const { return _poly_num1; }
   int poly_num2() const { return _poly_num2; }


private:
   Iterator _endpoint_ref;
   int _poly_num1;
   int _poly_num2;
};


template <class Traits>
class Edge_list : public std::list< 
                  Edge_info<typename Traits::Polygon_2::Vertex_iterator> >
{
public:
   typedef typename Traits::Point_2                      Point_2;
   typedef typename Traits::Orientation_2                Orientation_pred;
   typedef typename Traits::Polygon_2::Vertex_iterator   Vertex_iterator;
   typedef Edge_info<Vertex_iterator>                    Edge;
   typedef typename std::list<Edge>::iterator            Self_iterator;
   typedef typename std::list<Edge>::const_iterator      Self_const_iterator;
   typedef Circulator_from_iterator<Self_const_iterator> Self_const_circulator;

   void insert_next(Vertex_iterator endpoint_ref, int num)
   {
       Self_iterator e_it;

       for (e_it = begin(); 
            e_it != end() && (*e_it).endpoint() != endpoint_ref;
            e_it++) 
       {}

       if (e_it != end())
            (*e_it).set_poly_num2(num);
       else
          push_back(Edge(endpoint_ref, num));
   }

   void insert_prev(Vertex_iterator endpoint_ref, int num)
   {
       Self_iterator e_it;

       for (e_it = begin(); 
            e_it != end() && (*e_it).endpoint() != endpoint_ref;
            e_it++) 
       {}

       if (e_it != end())
            (*e_it).set_poly_num2(num);
       else
          push_front(Edge(endpoint_ref, num));
   }

   // PRE: polygons must be simple
   bool edges_overlap(Vertex_iterator vertex_it) 
   {

#ifdef CGAL_PARTITION_CHECK_DEBUG
       std::cout << "before sort: edges for " << *vertex_it << std::endl;
       std::cout << *this << std::endl;
#endif

       int num_unshared = 0;

       // Don't want to sort the edges for vertices of degree 2 because they
       // are already in CCW order (since the partition polygons were in CCW
       // order), and this is what you need when you construct the union 
       // polygon.
       if (size() > 2)
        sort(CW_indirect_edge_info_compare<Vertex_iterator,Traits>(vertex_it));

#ifdef CGAL_PARTITION_CHECK_DEBUG
       std::cout << "after sort: edges for " << *vertex_it  << std::endl;
       std::cout << *this << std::endl;
#endif

       Self_const_iterator prev_e_it = begin();
       Self_const_iterator e_it;

       for (e_it = begin(); e_it != end(); e_it++)
       {
          if ((*e_it).poly_num1() == PARTITION_VMAP_UNSHARED_EDGE) 
             num_unshared++;
          if ((*e_it).poly_num2() == PARTITION_VMAP_UNSHARED_EDGE) 
             num_unshared++;

          if ((*prev_e_it).poly_num1() != (*e_it).poly_num1() &&
              (*prev_e_it).poly_num1() != (*e_it).poly_num2() &&
              (*prev_e_it).poly_num2() != (*e_it).poly_num1() &&
              (*prev_e_it).poly_num2() != (*e_it).poly_num2())
          {
             return true;
          }
          prev_e_it = e_it;
       }
       if ((*prev_e_it).poly_num1() != (*begin()).poly_num1() &&
           (*prev_e_it).poly_num1() != (*begin()).poly_num2() &&
           (*prev_e_it).poly_num2() != (*begin()).poly_num1() &&
           (*prev_e_it).poly_num2() != (*begin()).poly_num2())
       {
          return true;
       }

       return (num_unshared > 2);
   }


   // NOTE:  the edges here are sorted in CW order so the next CCW edge
   //        comes BEFORE the edge with endpoint v_it in the sorted list
   Edge next_ccw_edge_info(Vertex_iterator v_it) const
   {
      Self_const_circulator first_e(begin(), end(), begin());
      Self_const_circulator e_circ = first_e;

      do
      {
         if ((*e_circ).endpoint() == v_it)
         {
            e_circ--; // go to the previous endpoint 
            return *e_circ;
         }
      }
      while (++e_circ != first_e);
      return *first_e; // shouldn't get here unless v_it is not in list
   }
};

template <class Traits>
std::ostream& operator<<(std::ostream& os, const Edge_list<Traits>& edges) 
{
   typename Edge_list<Traits>::const_iterator  e_it;

   for (e_it = edges.begin(); e_it != edges.end(); e_it++)
   {
          os << "   " << (*(*e_it).endpoint())
             << " from poly #" << (*e_it).poly_num1() 
             << " and poly #" << (*e_it).poly_num2() 
             << std::endl;
   }
   return os;
}

template <class Traits>
class Partition_vertex_map : 
                public std::map<typename Traits::Polygon_2::Vertex_iterator,
                                Edge_list<Traits>, 
                                Indirect_less_xy_2<Traits> >
{
public:

   typedef typename std::map<typename Traits::Polygon_2::Vertex_iterator,
                             Edge_list<Traits>,
                             Indirect_less_xy_2<Traits> >::iterator
                                                       Self_iterator;
   typedef typename Traits::Point_2                    Point_2;
   typedef typename Traits::Polygon_2::Vertex_iterator Vertex_iterator;

   Partition_vertex_map() {}

   template <class InputIterator>
   Partition_vertex_map(InputIterator first_poly, InputIterator last_poly)
   {  build(first_poly, last_poly); }
   
   template <class InputIterator>
   void build(InputIterator poly_first, InputIterator poly_last)
   {
      typedef std::pair<Self_iterator, bool>          Location_pair;
      typedef Edge_list<Traits>                       Edge_list;
      typedef std::pair<Vertex_iterator, Edge_list>   P_Vertex;
   
      Location_pair v_loc_pair;
      Location_pair begin_v_loc_pair;
      Location_pair prev_v_loc_pair;
   
      Vertex_iterator vtx_begin;
      Vertex_iterator vtx_end;
      Vertex_iterator v_it;
   
      int poly_num = 0;
      for (; poly_first != poly_last; poly_first++, poly_num++)
      {
        vtx_begin = (*poly_first).vertices_begin();
        vtx_end = (*poly_first).vertices_end();
        begin_v_loc_pair = insert(P_Vertex(vtx_begin, Edge_list()));
        prev_v_loc_pair = begin_v_loc_pair;
        v_it = vtx_begin;
        for (v_it++; v_it != vtx_end; v_it++)
        {
           v_loc_pair = insert(P_Vertex(v_it, Edge_list()));
           insert_next_edge(prev_v_loc_pair.first, v_loc_pair.first, poly_num);
           insert_prev_edge(v_loc_pair.first, prev_v_loc_pair.first, poly_num);
           prev_v_loc_pair = v_loc_pair;
        }
        insert_next_edge(prev_v_loc_pair.first, begin_v_loc_pair.first, 
                         poly_num);
        insert_prev_edge(begin_v_loc_pair.first, prev_v_loc_pair.first, 
                         poly_num);
      }
   }


   void insert_next_edge(Self_iterator& v1_ref, Self_iterator& v2_ref, 
                         int num)
   {
      (*v1_ref).second.insert_next((*v2_ref).first, num);
   }

   void insert_prev_edge(Self_iterator& v1_ref, Self_iterator& v2_ref, 
                         int num)
   {
      (*v1_ref).second.insert_prev((*v2_ref).first, num);
   }

   bool polygons_overlap()
   {
      Self_iterator v_it;
      for (v_it = begin(); v_it != end(); v_it++)
      {
         if ((*v_it).second.edges_overlap((*v_it).first)) return true;
      }
      return false;
   }

   template <class OutputIterator>
   OutputIterator union_vertices(OutputIterator result)
   {
       if (empty()) return result;
   
       Self_iterator m_it = begin();
       Vertex_iterator v_it;
       Vertex_iterator first_v_it;
       Vertex_iterator prev_v_it;
       Vertex_iterator next_v_it;
   
       // find a vertex with degree 2 (there must be at least one)
       while (m_it != end() && (*m_it).second.size() != 2)
          m_it++;
       CGAL_assertion (m_it != end());

       // insert this vertex and the two around it
       first_v_it = prev_v_it = (*(*m_it).second.begin()).endpoint();
#ifdef CGAL_PARTITION_CHECK_DEBUG
       std::cout << "union_vertices: inserting " << (*prev_v_it) << std::endl;
#endif
       *result = *prev_v_it;
       result++;
#ifdef CGAL_PARTITION_CHECK_DEBUG
       std::cout << "union_vertices: inserting "<< *(*m_it).first << std::endl;
#endif
       *result = *(*m_it).first;
       result++;
       next_v_it = (*m_it).second.back().endpoint();
#ifdef CGAL_PARTITION_CHECK_DEBUG
       std::cout << "union_vertices: inserting " << *next_v_it << std::endl;
#endif
       *result = *next_v_it;
       result++;

       // find the map iterator corresponding to the next vertex 
       prev_v_it  = (*m_it).first;
       v_it = next_v_it;
       m_it = find(v_it);
       
       while (v_it != first_v_it && m_it != end())
       {
#ifdef CGAL_PARTITION_CHECK_DEBUG
          std::cout << "union_vertices: prev_v_it " << (*prev_v_it)
                    << " v_it " << (*v_it) << " next_v_it "
                    << (*next_v_it) << std::endl;
#endif
          // Don't want to sort the edges for vertices of degree 2 because they
          // are already in CCW order (since the partition polygons were in CCW
          // order), and this is what you need to begin the construction 
          // of the union polygon.
          if ((*m_it).second.size() > 2)
           (*m_it).second.sort(
             CW_indirect_edge_info_compare<Vertex_iterator,Traits>(
                                                                (*m_it).first));
          // find the previous vertex in this vertex's list
          next_v_it=(*m_it).second.next_ccw_edge_info(prev_v_it).endpoint();
          if (next_v_it != first_v_it)
          {
#ifdef CGAL_PARTITION_CHECK_DEBUG
             std::cout << "union_vertices: inserting "
                       << *next_v_it << std::endl;
#endif
             *result = *next_v_it;
             result++;
          }
          prev_v_it  = v_it;
          v_it = next_v_it;
          m_it = find(v_it);
          CGAL_assertion (m_it == end() || (*m_it).first == v_it);
       }
#ifdef CGAL_PARTITION_CHECK_DEBUG
       if (v_it == first_v_it)
          std::cout << "union_vertices: stopped because first was reached "
                    << std::endl;
       else
          std::cout << "union_vertices: stopped because end was reached "
                    << std::endl;
#endif
       return result;
   }
};

}


#endif // CGAL_PARTITION_VERTEX_MAP_H
