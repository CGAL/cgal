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
// file          : include/CGAL/Partition_vertex_map.C
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

#include <CGAL/Partition_vertex_map.h>

namespace CGAL {

template <class Traits>
template <class InputIterator>
void 
Partition_vertex_map<Traits>::build(InputIterator poly_first, 
                                    InputIterator poly_last)
{
   typedef typename Traits::Polygon_2::Vertex_const_iterator
                                                   Poly_vtx_const_iterator;
   typedef std::pair<Self_iterator, bool>          Location_pair;
   typedef Edge_list<Traits>                       Edge_list;
   typedef typename Traits::Point_2                Point_2;
   typedef std::pair<Point_2, Edge_list>           P_Vertex;


   Location_pair v_loc_pair;
   Location_pair begin_v_loc_pair; 
   Location_pair prev_v_loc_pair;

   Poly_vtx_const_iterator begin;  
   Poly_vtx_const_iterator end;  
   Poly_vtx_const_iterator v_it;  

   int poly_num = 0;
   for (; poly_first != poly_last; poly_first++, poly_num++)
   {
      begin = (*poly_first).vertices_begin();
      end = (*poly_first).vertices_end();
      begin_v_loc_pair= insert(P_Vertex(*begin, Edge_list()));
      prev_v_loc_pair = begin_v_loc_pair;
      v_it = begin;
      for (v_it++; v_it != end; v_it++)
      {
         v_loc_pair = insert(P_Vertex(*v_it, Edge_list()));
         insert_next_edge(prev_v_loc_pair.first, v_loc_pair.first, poly_num);
         insert_prev_edge(v_loc_pair.first, prev_v_loc_pair.first, poly_num);
         prev_v_loc_pair = v_loc_pair;
      }
      insert_next_edge(prev_v_loc_pair.first, begin_v_loc_pair.first, poly_num);
      insert_prev_edge(begin_v_loc_pair.first, prev_v_loc_pair.first, poly_num);
   }
}


template <class Traits>
template <class OutputIterator>
OutputIterator
Partition_vertex_map<Traits>::union_vertices(OutputIterator result) 
{
    if (empty()) return result;

    Self_iterator first = begin();
    Self_iterator v_it = first;
    Self_iterator prev_v_it;
    bool inserting = false;
    Self_iterator next_v_it;

    do 
    {
       // Don't want to sort the edges for vertices of degree 2 because they
       // are already in CCW order (since the partition polygons were in CCW
       // order), and this is what you need when to begin the construction of
       // the union polygon.
       if ((*v_it).second.size() > 2)
       {
        (*v_it).second.sort(
         CW_indirect_edge_info_compare<Self_iterator,Traits>((*v_it).first));
       }
       if (!inserting)
       {
           if ((*v_it).second.size() == 2)
           {
              inserting = true;
              // insert this vertex and the two around it  
              first = prev_v_it = (*(*v_it).second.begin()).endpoint();
#ifdef CGAL_PARTITION_CHECK_DEBUG
              std::cout << "union_vertices: inserting " 
                        << (*prev_v_it).first << std::endl;
#endif
              *result = (*prev_v_it).first;
              result++;
#ifdef CGAL_PARTITION_CHECK_DEBUG
              std::cout << "union_vertices: inserting " 
                        << (*v_it).first << std::endl;
#endif
              *result = (*v_it).first;
              result++;
              next_v_it = (*v_it).second.last_edge_info().endpoint();
#ifdef CGAL_PARTITION_CHECK_DEBUG
              std::cout << "union_vertices: inserting " 
                        << (*next_v_it).first << std::endl;
#endif
              *result = (*next_v_it).first;
              result++;
           }
           else 
           {
              next_v_it = v_it;
              next_v_it++;
           }
       }
       else 
       {
          // find the previous vertex in this vertex's list
          next_v_it =(*v_it).second.next_ccw_edge_info(prev_v_it).endpoint();
          if (next_v_it != first) 
          {
#ifdef CGAL_PARTITION_CHECK_DEBUG
             std::cout << "union_vertices: inserting " 
                       << (*next_v_it).first << std::endl;
#endif
             *result = (*next_v_it).first;
             result++;
          }
       }
       prev_v_it  = v_it;
       v_it = next_v_it;
#ifdef CGAL_PARTITION_CHECK_DEBUG
       std::cout << "union_vertices: prev_v_it " << (*prev_v_it).first 
                 << " v_it " << (*v_it).first << " next_v_it " 
                 << (*next_v_it).first << std::endl;
#endif
    }
    while (v_it != first && v_it != end());
#ifdef CGAL_PARTITION_CHECK_DEBUG
    if (v_it == first)
       std::cout << "union_vertices: stopped because first was reached " 
                 << std::endl;
    else
       std::cout << "union_vertices: stopped because end was reached " 
                 << std::endl;
#endif
    return result;
}

}
