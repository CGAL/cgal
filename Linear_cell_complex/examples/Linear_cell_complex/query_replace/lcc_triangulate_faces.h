// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of 3d-query-replace.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef LCC_TRIANGULATE_FACES_H
#define LCC_TRIANGULATE_FACES_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Linear_cell_complex_operations.h>
///////////////////////////////////////////////////////////////////////////////
template<typename Face_handle>
bool is_external(Face_handle fh)
{
  return fh->info().is_external;
}

template<typename Face_handle>
int number_of_existing_edge(Face_handle fh)
{
  unsigned res=0;
  for(int i=0; i<3; ++i)
  { if(fh->info().exist_edge[i]) ++res; }
  return res;
}

template<typename Face_handle>
int get_free_edge(Face_handle fh)
{
  CGAL_assertion( number_of_existing_edge(fh)==2 );
  for(int i=0; i<3; ++i)
  { if(!fh->info().exist_edge[i]) return i; }

  CGAL_assertion(false);
  return -1;
}
///////////////////////////////////////////////////////////////////////////////
/// Triangulate the face containing dart d1.
template<typename LCC>
void constrained_delaunay_triangulation(LCC &lcc, typename LCC::Dart_handle d1)
{
  struct Vertex_info
  {
    typename LCC::Dart_handle dh;
    typename LCC::Vector v;
  };

  struct Face_info {
    bool exist_edge[3];
    bool is_external;
    bool is_process;
  };

  typedef CGAL::Projection_traits_3<CGAL::Exact_predicates_inexact_constructions_kernel> P_traits;
  typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, P_traits> Vb;

  typedef CGAL::Triangulation_face_base_with_info_2<Face_info,P_traits> Fb1;

  typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>    Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                   TDS;
  typedef CGAL::Exact_predicates_tag                                    Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS,
                                                     Itag>              CDT;

  if(lcc.template beta<1,1,1>(d1)==d1) { return; } // The face is already triangulated

  typename LCC::Vector normal=CGAL::compute_normal_of_cell_2(lcc,d1);
  P_traits cdt_traits(normal);
  CDT cdt(cdt_traits);

  //inserting the constraints edge by edge
  typename LCC::template Dart_of_orbit_range<1>::iterator
    it(lcc.template darts_of_orbit<1>(d1).begin());

  typename CDT::Vertex_handle previous=LCC::null_handle, first=LCC::null_handle,
    vh=LCC::null_handle;

   for(typename LCC::template Dart_of_orbit_range<1>::iterator
         itend(lcc.template darts_of_orbit<1>(d1).end()); it!=itend; ++it)
   {
     vh=cdt.insert(lcc.point(it));
     vh->info().dh=it;
     if(first==nullptr)
     { first=vh; }
     if(previous!=nullptr)
     {
       if(previous!=vh)
       { cdt.insert_constraint(previous,vh); }
     }

     previous=vh;
   }
   cdt.insert_constraint(previous,first);
   CGAL_assertion(cdt.is_valid());

   // sets mark is_external
   for(typename CDT::All_faces_iterator fit=cdt.all_faces_begin(),
         fitend=cdt.all_faces_end(); fit!=fitend; ++fit)
   {
     fit->info().is_external  = true;
     fit->info().is_process   = false;
     fit->info().exist_edge[0]=false;
     fit->info().exist_edge[1]=false;
     fit->info().exist_edge[2]=false;
   }

   std::queue<typename CDT::Face_handle> face_queue;
   typename CDT::Face_handle face_internal = nullptr;

   face_queue.push(cdt.infinite_vertex()->face());
   while(!face_queue.empty())
   {
     typename CDT::Face_handle fh=face_queue.front();
     face_queue.pop();
     if(!fh->info().is_process)
     {
       fh->info().is_process = true;
       for(int i = 0; i<3; ++i)
       {
         if(!cdt.is_constrained(std::make_pair(fh, i)))
         {
           face_queue.push(fh->neighbor(i));
         }
         else if(face_internal==nullptr)
         {
           face_internal=fh->neighbor(i);
         }
       }
     }
   }
   if(face_internal!=nullptr)
   { face_queue.push(face_internal); }
   
   while(!face_queue.empty())
   {
     typename CDT::Face_handle fh=face_queue.front();
     face_queue.pop();
     if(!fh->info().is_process)
     {
       fh->info().is_process =true;
       fh->info().is_external=false;
       for(int i = 0; i <3; ++i)
       {
         if(!cdt.is_constrained(std::make_pair(fh, i)))
         {
           face_queue.push(fh->neighbor(i));
         }
       }
     }
   }

   for(typename CDT::Finite_edges_iterator eit = cdt.finite_edges_begin(),
         eitend=cdt.finite_edges_end(); eit!=eitend; ++eit)
   {
     typename CDT::Face_handle fh=eit->first;
     int index=eit->second;
     typename CDT::Face_handle opposite_fh=fh->neighbor(index);
     if(cdt.is_constrained(std::make_pair(fh, index)))
     {
       fh->info().exist_edge[index]=true;
       opposite_fh->info().exist_edge[cdt.mirror_index(fh,index)]=true;

       if(!fh->info().is_external && number_of_existing_edge(fh)==2)
         face_queue.push(fh);
       if(!opposite_fh->info().is_external &&
          number_of_existing_edge(opposite_fh)==2)
         face_queue.push(opposite_fh);
     }
   }

   while(!face_queue.empty())
   {
     typename CDT::Face_handle fh=face_queue.front();
     face_queue.pop();
     CGAL_assertion(number_of_existing_edge(fh)>=2); // i.e. ==2 or ==3
     CGAL_assertion(!fh->info().is_external);

     if(number_of_existing_edge(fh)==2)
     {
       int index=get_free_edge(fh);
       typename CDT::Face_handle opposite_fh=fh->neighbor(index);

       CGAL_assertion( !fh->info().exist_edge[index] );
       CGAL_assertion( !opposite_fh->info().
                       exist_edge[cdt.mirror_index(fh,index)] );
       // triangle is (vc, vb, va)
       const typename CDT::Vertex_handle va = fh->vertex(cdt. cw(index));
       const typename CDT::Vertex_handle vb = fh->vertex(cdt.ccw(index));
       const typename CDT::Vertex_handle vc = fh->vertex(index);

       typename LCC::Dart_handle dd1 = nullptr;
       for(typename LCC::template Dart_of_cell_range<0, 2>::iterator
             iti=lcc.template darts_of_cell<0, 2>(va->info().dh).begin();
           dd1==nullptr && iti.cont(); ++iti)
       {
         if(lcc.point(lcc.template beta<1>(iti))==vc->point())
         { dd1=iti; }
       }

       typename LCC::Dart_handle dd2 = nullptr;
       for(typename LCC::template Dart_of_cell_range<0, 2>::iterator
             iti=lcc.template darts_of_cell<0, 2>(vb->info().dh).begin();
           dd2==nullptr && iti.cont(); ++iti)
       {
         if(lcc.point(lcc.template beta<0>(iti))==vc->point())
         { dd2=iti; }
       }

       //       assert(((lcc.beta<0,0>(dd1)==dd2) || lcc.beta<1,1>(dd1)==dd2));

       typename LCC::Dart_handle ndart=lcc.insert_cell_1_in_cell_2(dd1, dd2);
       va->info().dh=lcc.template beta<2>(ndart);

       fh->info().exist_edge[index]=true;
       opposite_fh->info().exist_edge[cdt.mirror_index(fh,index)]=true;

       if(!opposite_fh->info().is_external &&
          number_of_existing_edge(opposite_fh)==2)
       { face_queue.push(opposite_fh); }
     }
   }
}
///////////////////////////////////////////////////////////////////////////////
/// Triangulate all marked faces (each face having at least one marked dart)
template<typename LCC>
void triangulate_marked_faces(LCC &lcc, typename LCC::size_type amark)
{
  // We are going to call constrained_delaunay_triangulation several time for
  // a same face when it has all its darts marked. But only the first call will
  // triangulate the face, the other ones will do nothing since the faces are
  // already triangulated. It is maybe faster to mark darts of the face in order
  // to avoid these successive calls, but not sure since we must unmark the
  // marked darts in another loop.
  for(auto it=lcc.darts().begin(); it!=lcc.darts().end(); ++it)
  {
    if(lcc.is_marked(it, amark))
    { constrained_delaunay_triangulation(lcc, it); }
  }
}
///////////////////////////////////////////////////////////////////////////////
/// Triangulate all faces of the lcc.
template<typename LCC>
void triangulate_all_faces(LCC &lcc)
{
  for(auto it=lcc.darts().begin(); it!=lcc.darts().end(); ++it)
  { constrained_delaunay_triangulation(lcc, it); }
}
///////////////////////////////////////////////////////////////////////////////
#endif // LCC_TRIANGULATE_FACES_H
