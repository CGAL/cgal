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
// 
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_PARTITION_VERTEX_MAP_H
#define CGAL_PARTITION_VERTEX_MAP_H

#include <map>
#include <iostream>
#include <CGAL/circulator.h>
#include <CGAL/assertions.h>

#include <sstream>

namespace CGAL {

const int PARTITION_VMAP_UNSHARED_EDGE = -1;

template<class Traits_>
class Vertex_info 
{
public:

  typedef Traits_                                       Traits;
  typedef typename Traits::Polygon_2                    Polygon_2 ; 
  typedef typename Traits::Polygon_2::Vertex_const_iterator   Vertex_const_iterator;
  typedef typename Traits::Less_xy_2                    Less_xy_2;

  Vertex_info ( Vertex_const_iterator const& vx_it, Polygon_2 const* poly_ptr )
   :
    m_vx_it(vx_it)
   ,m_poly_ptr(poly_ptr)
  {}

  Vertex_const_iterator  vertex_it() const { return m_vx_it ; }
  Polygon_2 const* poly_ptr () const { return m_poly_ptr ; }

  friend bool operator == ( Vertex_info const& a, Vertex_info const& b )
  {
    return a.poly_ptr() == b.poly_ptr() && a.vertex_it() == b.vertex_it() ;   
  }

  friend bool operator != ( Vertex_info const& a, Vertex_info const& b ) { return !(a==b); }

  friend bool operator < ( Vertex_info const& a, Vertex_info const& b )
  {
    return Traits().less_xy_2_object()(*a.vertex_it(), *b.vertex_it());
  }

private:

  Vertex_const_iterator  m_vx_it   ;
  Polygon_2 const* m_poly_ptr ;
} ;

template <class Traits_>
class Edge_info
{
public:

 typedef Traits_                    Traits;
 typedef CGAL::Vertex_info<Traits>  Vertex_info;

public:

   Edge_info(Vertex_info e_ref, int p_num1) : _endpoint_ref(e_ref),
      _poly_num1(p_num1), _poly_num2(PARTITION_VMAP_UNSHARED_EDGE)
   { }

   void set_poly_num2(int p_num)
   {
      _poly_num2 = p_num;
   }

   Vertex_info endpoint() const { return _endpoint_ref; }
   int poly_num1() const { return _poly_num1; }
   int poly_num2() const { return _poly_num2; }


private:
   Vertex_info _endpoint_ref;
   int _poly_num1;
   int _poly_num2;
};

template <class Traits>
class CW_indirect_edge_info_compare
{
public:

   typedef CGAL::Vertex_info<Traits> Vertex_info ;
   typedef CGAL::Edge_info<Traits>   Edge_info;

   typedef typename Vertex_info::Vertex_const_iterator Vertex_const_iterator ;

   typedef typename Traits::Left_turn_2 Left_turn_2;
   typedef typename Traits::Less_xy_2   Less_xy_2;
   typedef typename Traits::Point_2     Point_2;

   CW_indirect_edge_info_compare (Vertex_const_iterator v_info) : vertex_it(v_info),
      left_turn(Traits().left_turn_2_object()),
      less_xy(Traits().less_xy_2_object())
   {}

   bool operator()(Edge_info e1, Edge_info e2)
   {
      bool e1_less = less_xy((*e1.endpoint().vertex_it()), *vertex_it);
      bool e2_less = less_xy((*e2.endpoint().vertex_it()), *vertex_it);
      bool e1_to_e2_left_turn = left_turn((*e1.endpoint().vertex_it()), *vertex_it, 
                                (*e2.endpoint().vertex_it()));

      // if both edges are on the same side of the vertical line through 
      // _vertex then e1 comes before e2 (in CW order from the vertical line)
      // if one  makes a left turn going from e1 to e2
      if (e1_less == e2_less)
          return e1_to_e2_left_turn;
      else // e1 comes first if it is to the right of the vertical line
          return !e1_less;
   }

private:
   Vertex_const_iterator vertex_it;  
   Left_turn_2     left_turn;
   Less_xy_2       less_xy;
};




namespace Partition_2 {

template <class Traits_>
class Edge_list 
{
public:

  typedef Traits_                                       Traits;
  typedef Edge_list<Traits>                             Self;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::Orientation_2                Orientation_pred;
  typedef typename Traits::Polygon_2                    Polygon_2 ; 
  typedef typename Traits::Polygon_2::Vertex_const_iterator   Vertex_const_iterator;

  typedef CGAL::Vertex_info<Traits> Vertex_info;
  typedef CGAL::Edge_info<Traits>   Edge_info;

  typedef std::list<Edge_info> List ;

  typedef typename List::iterator       Self_iterator;
  typedef typename List::const_iterator Self_const_iterator;
  typedef typename List::size_type      size_type ;

  typedef Circulator_from_iterator<Self_const_iterator> Self_const_circulator;

  Self_const_iterator begin() const { return m_list.begin() ; }
  Self_iterator       begin()       { return m_list.begin() ; }
  Self_const_iterator end  () const { return m_list.end  () ; }
  Self_iterator       end  ()       { return m_list.end  () ; }

  size_type size() const { return m_list.size() ; }

  Edge_info const& front() const { return m_list.front() ; }
  Edge_info      & front()       { return m_list.front() ; }

  Edge_info const& back() const { return m_list.back() ; }
  Edge_info      & back()       { return m_list.back() ; }

  template<class Compare> void sort ( Compare c ) { m_list.sort(c); }

  void insert_next(Vertex_info endpoint_ref, int num)
  {

    Self_iterator e_it;

    for (e_it = m_list.begin();  e_it != m_list.end() && e_it->endpoint() != endpoint_ref ; e_it++) 
    {
    }

    if (e_it != m_list.end())
    {
      (*e_it).set_poly_num2(num);
    }
    else
    {
      m_list.push_back(Edge_info(endpoint_ref, num));
    }
  }

  void insert_prev(Vertex_info endpoint_ref, int num)
  {

    Self_iterator e_it;

    for (e_it = m_list.begin(); e_it != m_list.end() &&  e_it->endpoint() != endpoint_ref ; e_it++) 
    {
    }

    if (e_it != m_list.end())
    {
      (*e_it).set_poly_num2(num);
    }
    else
    {
       m_list.push_front(Edge_info(endpoint_ref, num));
    }
  }

  // PRE: polygons must be simple
  bool edges_overlap(Vertex_const_iterator vertex_it) 
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
    if (m_list.size() > 2)
    {
      m_list.sort(CW_indirect_edge_info_compare<Traits>(vertex_it));
    }

#ifdef CGAL_PARTITION_CHECK_DEBUG
    std::cout << "after sort: edges for " << *vertex_it  << std::endl;
    std::cout << *this << std::endl;
#endif

    Self_const_iterator prev_e_it = m_list.begin();
    Self_const_iterator e_it;

    for (e_it = m_list.begin(); e_it != m_list.end(); e_it++)
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
    if ((*prev_e_it).poly_num1() != (*m_list.begin()).poly_num1() &&
        (*prev_e_it).poly_num1() != (*m_list.begin()).poly_num2() &&
        (*prev_e_it).poly_num2() != (*m_list.begin()).poly_num1() &&
        (*prev_e_it).poly_num2() != (*m_list.begin()).poly_num2())
    {
       return true;
    }

    return (num_unshared > 2);
  }


  // NOTE:  the edges here are sorted in CW order so the next CCW edge
  //        comes BEFORE the edge with endpoint v_info in the sorted list
  Edge_info next_ccw_edge_info(Vertex_info v_info) const
  {
    Self_const_circulator first_e(m_list.begin(), m_list.end(), m_list.begin());
    Self_const_circulator e_circ = first_e;

    do
    {
      if ((*e_circ).endpoint() == v_info)
      {
        e_circ--; // go to the previous endpoint 
        return *e_circ;
      }
    }
    while (++e_circ != first_e);
    return *first_e; // shouldn't get here unless v_info is not in list
  }

private :

  List m_list ;
};


template <class Traits>
std::ostream& operator<<(std::ostream& os, const Edge_list<Traits>& edges) 
{
   typename Edge_list<Traits>::const_iterator  e_it;

   for (e_it = edges.begin(); e_it != edges.end(); e_it++)
   {
    os << "edge with endpoint (" << (*(*e_it).endpoint().vertex_it())
             << ") from poly #" << (*e_it).poly_num1() 
             << " and poly #" << (*e_it).poly_num2() 
             << std::endl;
   }
   return os;
}

} // namesapce Partition_2

template <class Traits_>
class Partition_vertex_map  
{
public:

   typedef Traits_                         Traits;
   typedef CGAL::Vertex_info<Traits>       Vertex_info;
   typedef CGAL::Edge_info<Traits>         Edge_info;
   typedef Partition_2::Edge_list<Traits>  Edge_list;  

   typedef Partition_vertex_map<Traits> Self;

   typedef std::map<Vertex_info, Edge_list> Map ;

   typedef typename Map::const_iterator Self_const_iterator;
   typedef typename Map::iterator       Self_iterator;

   typedef typename Traits::Point_2   Point_2;
   typedef typename Traits::Polygon_2 Polygon_2 ;

   typedef typename Polygon_2::Vertex_const_iterator Vertex_const_iterator;

   Partition_vertex_map() {}

   template <class InputIterator>
   Partition_vertex_map(InputIterator first_poly, InputIterator last_poly)
   {  _build(first_poly, last_poly); }
  
   Self_const_iterator begin() const { return m_map.begin() ; }
   Self_iterator       begin()       { return m_map.begin() ; }
   Self_const_iterator end  () const { return m_map.end  () ; }
   Self_iterator       end  ()       { return m_map.end  () ; }

   bool polygons_overlap()
   {
      Self_iterator v_info;
      for (v_info = m_map.begin(); v_info != m_map.end(); v_info++)
      {
         if ((*v_info).second.edges_overlap((*v_info).first.vertex_it())) 
            return true;
      }
      return false;
   }

   template <class OutputIterator>
   OutputIterator union_vertices(OutputIterator result)
   {
       if (m_map.empty()) 
         return result;
   
       Self_iterator m_it = m_map.begin();
   
       // find a vertex with degree 2 (there must be at least one)
       while (m_it != m_map.end() && (*m_it).second.size() != 2)
          m_it++;
       CGAL_assertion (m_it != m_map.end());

       // insert this vertex and the two around it
       Vertex_info first_v_info = (*m_it).second.front().endpoint();
       Vertex_info prev_v_info = first_v_info ;
#ifdef CGAL_PARTITION_CHECK_DEBUG
       std::cout << "union_vertices: inserting " << (*prev_v_info.vertex_it()) << std::endl;
#endif
       *result = *prev_v_info.vertex_it();
       result++;
#ifdef CGAL_PARTITION_CHECK_DEBUG
       std::cout << "union_vertices: inserting "<< *(*m_it).first.vertex_it() << std::endl;
#endif
       *result = *(*m_it).first.vertex_it();
       result++;
       Vertex_info next_v_info = (*m_it).second.back().endpoint();
#ifdef CGAL_PARTITION_CHECK_DEBUG
       std::cout << "union_vertices: inserting " << *next_v_info.vertex_it() << std::endl;
#endif
       *result = *next_v_info.vertex_it();
       result++;

       // find the map iterator corresponding to the next vertex 
       prev_v_info  = (*m_it).first;
       Vertex_info v_info = next_v_info;
       m_it = m_map.find(v_info);
       
       while (v_info != first_v_info && m_it != m_map.end())
       {
#ifdef CGAL_PARTITION_CHECK_DEBUG
          std::cout << "union_vertices: prev_v_info " << (*prev_v_info.vertex_it())
                    << " v_info " << (*v_info.vertex_it()) << " next_v_info "
                    << (*next_v_info.vertex_it()) << std::endl;
#endif
    // Don't want to sort the edges for vertices of degree 2 because they
    // are already in CCW order (since the partition polygons were in CCW
    // order), and this is what you need to begin the construction 
    // of the union polygon.
          if ((*m_it).second.size() > 2)
          {
            (*m_it).second.sort(
              CW_indirect_edge_info_compare<Traits>((*m_it).first.vertex_it()));
       	  }

          // find the previous vertex in this vertex's list
          next_v_info=(*m_it).second.next_ccw_edge_info(prev_v_info).endpoint();
          if (next_v_info != first_v_info)
          {
#ifdef CGAL_PARTITION_CHECK_DEBUG
             std::cout << "union_vertices: inserting "
                       << *next_v_info.vertex_it() << std::endl;
#endif
             *result = *next_v_info.vertex_it();
             result++;
          }
          prev_v_info  = v_info;
          v_info = next_v_info;
          m_it = m_map.find(v_info);
          CGAL_assertion (m_it == m_map.end() || (*m_it).first == v_info);
       }
#ifdef CGAL_PARTITION_CHECK_DEBUG
       if (v_info == first_v_info)
          std::cout << "union_vertices: stopped because first was reached "
                    << std::endl;
       else
          std::cout << "union_vertices: stopped because end was reached "
                    << std::endl;
#endif
       return result;
   }

private :

   template <class InputIterator>
   void _build(InputIterator poly_first, InputIterator poly_last)
   {

      typedef std::pair<Self_iterator, bool>          Location_pair;
      typedef std::pair<Vertex_info, Edge_list>   P_Vertex;
   
      Location_pair v_loc_pair;
      Location_pair begin_v_loc_pair;
      Location_pair prev_v_loc_pair;
   
      Vertex_const_iterator vtx_begin;
      Vertex_const_iterator vtx_end;
      Vertex_const_iterator v_it;
   
      int poly_num = 0;
      for (; poly_first != poly_last; poly_first++, poly_num++)
      {

        Polygon_2 const* poly_ptr = &(*poly_first);

        vtx_begin = (*poly_first).vertices_begin();
        vtx_end   = (*poly_first).vertices_end();
        begin_v_loc_pair = m_map.insert(P_Vertex( Vertex_info(vtx_begin,poly_ptr), Edge_list()));
  
        prev_v_loc_pair = begin_v_loc_pair;
        v_it = vtx_begin;

        for (v_it++; v_it != vtx_end; v_it++)
        {

           v_loc_pair = m_map.insert(P_Vertex( Vertex_info(v_it,poly_ptr), Edge_list()));

           insert_next_edge(prev_v_loc_pair.first,  v_loc_pair.first, poly_num);

           insert_prev_edge(v_loc_pair.first, prev_v_loc_pair.first, poly_num);

           prev_v_loc_pair = v_loc_pair;
        }
        insert_next_edge(prev_v_loc_pair.first, begin_v_loc_pair.first, 
                         poly_num);
        insert_prev_edge(begin_v_loc_pair.first, prev_v_loc_pair.first, 
                         poly_num);
      }
   }


   void insert_next_edge(Self_iterator& v1_ref, Self_iterator& v2_ref, int num)
   {
      (*v1_ref).second.insert_next((*v2_ref).first, num);
   }

   void insert_prev_edge(Self_iterator& v1_ref, Self_iterator& v2_ref, int num)
   {
      (*v1_ref).second.insert_prev((*v2_ref).first, num);
   }

private :

  Map m_map ;
};

}

#endif // CGAL_PARTITION_VERTEX_MAP_H
