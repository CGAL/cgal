// Copyright (c) 2017 
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     :

#ifndef CGAL_BUFFER_FOR_VAO_H
#define CGAL_BUFFER_FOR_VAO_H

#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/IO/Color.h>

#include <vector>
#include <cstdlib>
#include <queue>
#include <boost/unordered_map.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
typedef Local_kernel::Point_3  Local_point;
typedef Local_kernel::Vector_3 Local_vector;

//------------------------------------------------------------------------------
namespace internal {
  template <class Point, class Vector>
  void newell_single_step_3(const Point& p, const Point& q, Vector& n)
  {
    // Compute normal of the face by using Newell's method: for each edge PQ
    // Nx += (Py - Qy) * (Pz + Qz);
    // Ny += (Pz - Qz) * (Px + Qx);
    // Nz += (Px - Qx) * (Py + Qy);
    n = Vector(n.x()+((p.y()-q.y())*(p.z()+q.z())),
               n.y()+((p.z()-q.z())*(p.x()+q.x())),
               n.z()+((p.x()-q.x())*(p.y()+q.y())));
    }

  Local_vector compute_normal_of_face(const std::vector<Local_point>& points)
  {
    Local_vector normal(CGAL::NULL_VECTOR);
    unsigned int nb = 0;
    for (std::size_t i=0; i<points.size(); ++i)
    {
      newell_single_step_3(points[i], points[(i+1)%points.size()], normal);
      ++nb;
    }
    
    assert(nb>0);
    return (typename Local_kernel::Construct_scaled_vector_3()(normal, 1.0/nb));
  }

  template<int dim>
  struct Geom_utils;

  template<>
  struct Geom_utils<3>
  {
    template<typename KPoint>
    static Local_point get_local_point(const KPoint& p)
    {
      CGAL::Cartesian_converter<typename CGAL::Kernel_traits<KPoint>::Kernel, Local_kernel> converter;
      return converter(p);
    }

    template<typename KVector>
    static Local_vector get_local_vector(const KVector& v)
    {
      CGAL::Cartesian_converter<typename CGAL::Kernel_traits<KVector>::Kernel, Local_kernel> converter;
      return converter(v);
    }
  };

  template<>
  struct Geom_utils<2>
  {
    template<typename KPoint>
    static Local_point get_local_point(const KPoint& p)
    {
      CGAL::Cartesian_converter<typename CGAL::Kernel_traits<KPoint>::Kernel, Local_kernel> converter;
      return Local_point(converter(p.x()),0,converter(p.y()));
    }

    template<typename KVector>
    static Local_vector get_local_vector(const KVector& v)
    {
      CGAL::Cartesian_converter<typename CGAL::Kernel_traits<KVector>::Kernel, Local_kernel> converter;
      return Local_vector(converter(v.x()),0,converter(v.y()));
    }
  };

  template<typename KPoint>
  Local_point get_local_point(const KPoint& p)
  {
    return Geom_utils<CGAL::Ambient_dimension<KPoint>::value>::get_local_point(p);
  }

  template<typename KVector>
  Local_vector get_local_vector(const KVector& v)
  {
    return Geom_utils<CGAL::Ambient_dimension<KVector>::value>::get_local_vector(v);
  }
} // End namespace internal
//------------------------------------------------------------------------------
template<typename T>
class Buffer_for_vao
{
public:
  Buffer_for_vao(std::vector<T>* pos,
                 CGAL::Bbox_3* bbox=NULL,
                 std::vector<T>* color=NULL,
                 std::vector<T>* flat_normal=NULL,
                 std::vector<T>* gourod_normal=NULL,
                 std::vector<std::size_t>& indices=NULL) :
    m_pos_buffer(pos),
    m_bb(bbox),
    m_color_buffer(color),
    m_flat_normal_buffer(flat_normal),
    m_gourod_normal_buffer(gourod_normal),
    m_index_buffer(indices),
    m_face_started(false)
  {}

  void clear()
  {
    m_pos_buffer->clear();
    if (m_color_buffer!=NULL)         { m_color_buffer->clear(); }
    if (m_flat_normal_buffer!=NULL)   { m_flat_normal_buffer->clear(); }
    if (m_gourod_normal_buffer!=NULL) { m_gourod_normal_buffer->clear(); }
    if (m_index_buffer!=NULL)         { m_index_buffer->clear(); }
  }

  bool is_empty() const
  {
    return m_pos_buffer->empty() &&
      (m_color_buffer!=NULL || m_color_buffer->empty()) &&
      (m_flat_normal_buffer!=NULL || m_flat_normal_buffer->empty()) &&
      (m_gourod_normal_buffer!=NULL || m_gourod_normal_buffer->empty()) &&
      (m_index_buffer!=NULL || m_index_buffer->empty());
  }

  bool has_color() const
  { return m_color_buffer!=NULL; }

  bool has_flat_normal() const
  { return m_flat_normal_buffer!=NULL; }
  
  bool has_gourod_normal() const
  { return m_gourod_normal_buffer!=NULL; }

  bool has_indices() const
  { return m_index_buffer!=NULL; }

  // 1.1) Add a point, without color.
  template<typename KPoint>
  void add_point(const KPoint& kp)
  {
    CGAL_assertion(m_pos_buffer!=NULL);
    Local_point p=internal::get_local_point(kp);
    m_pos_buffer->push_back(p.x());
    m_pos_buffer->push_back(p.y());
    m_pos_buffer->push_back(p.z());

    if (m_bb!=NULL)
    { (*m_bb)=(*m_bb)+p.bbox(); }
  }

  // 1.2) Add a point, with color.
  template<typename KPoint>
  void add_point(const KPoint& kp, const CGAL::Color& c)
  {
    CGAL_assertion(!is_data_indexed());
    add_point(kp);
    add_color(c);
  }

  // 1.3) Add an indexed point, without color.
  template<typename KPoint, typename KVector>
  void add_indexed_point(const KPoint& kp,
                         unsigned int index)
  {
    CGAL_assertion(is_data_indexed());
    Local_point p = internal::get_local_point(kp);
    m_points_vector.push_back(p);
    m_point_index_map[p]=index;
  }

  // 1.3) Add an indexed point, with color.
  template<typename KPoint, typename KVector>
  void add_indexed_point(const KPoint& kp,
                         unsigned int index,
                         const CGAL::Color& color)
  {
    CGAL_assertion(is_data_indexed());
    add_indexed_point(kp, index);
    add_color(color);
  }

  // 2.1) Add a segment, without color.
  template<typename KPoint>
  void add_segment(const KPoint& kp1, const KPoint& kp2)
  {
    if (is_data_indexed())
    {
      add_point(kp1);
      add_point(kp2);
    }
    else
    {
      m_point_index->push_back(m_point_index_map[p1]);
      m_point_index->push_back(m_point_index_map[p2]);
    }
  }
  
  // 2.2) Add a segment, with color.
  template<typename KPoint>
  void add_segment(const KPoint& kp1, const KPoint& kp2, const CGAL::Color& c)
  {
    add_segment(kp1, kp2);
    add_color(c);
    add_color(c);
  }

  // 3.1) Add a face, without color, without normal.
  bool is_a_face_started() const
  { return m_face_started; }
  
  void face_begin()
  { face_begin_internal(false, false); }

  // 3.2) Add a face, with a color, without normal.
  void face_begin(const CGAL::Color& c)
  {
    m_color_of_face=c;
    face_begin_internal(true, false);
  }

  // 3.3) Add a face, without a color, with a normal.
  template<typename KNormal>
  void face_begin(const KNormal& kv)
  {
    m_normal_of_face=internal::get_local_vector(kv);
    face_begin_internal(false, true);
  }

  // 3.3) Add a face, with a color and with a normal.
  template<typename KNormal>
  void face_begin(const CGAL::Color& c, const KNormal& kv)
  {
    m_color_of_face=c;
    m_normal_of_face=internal::get_local_vector(kv);
    face_begin_internal(true, true);
  }

  /// Add a point at the end of the current face
  /// @param p the point to add
  /// @p_normal the vertex normal at this point (for Gourod shading)
  template<typename KPoint, typename KVector>
  bool add_point_in_face(const KPoint& p, const KVector& p_normal)
  {
    if (!is_a_face_started()) return false;
    
    if (add_point_in_face(p))
    {
      m_vertex_normals_for_face.push_back(internal::get_local_vector(p_normal));
      return true;
    }
    return false;
  }

  /// Add a point at the end of the current face, without giving the vertex normal.
  /// When this method is used, it is not possible to use the Gourod shading.
  /// @param p the point to add
  template<typename KPoint>
  bool add_point_in_face(const KPoint& kp)
  {
    if (!is_a_face_started()) return false;

    Local_point p=internal::get_local_point(kp);
    if (m_points_of_face.empty() || m_points_of_face.back()!=p) // TODO test if the distance between prev point and kp is smaller than an epsilon (?? not sure ??)
    {
      m_points_of_face.push_back(p);
      return true;
    }
    return false;
  }
  
  /// End the face: compute the triangulation.
  void face_end()
  {
    if (!is_a_face_started()) return;

    if (m_points_of_face.size()<3)
    {
      std::cout<<"PB: you try to triangulate a face with "<<m_points_of_face.size()<<" vertices."
               <<std::endl;
      
      m_face_started=false;
      m_points_of_face.clear();
      m_vertex_normals_for_face.clear();

      return;
    }
    
    Local_vector normal=(m_started_face_has_normal?m_normal_of_face:
                         internal::compute_normal_of_face(m_points_of_face));

    if (m_points_of_face.size()==3) // Triangle: no need to triangulate
    { triangular_face_end_internal(normal); }
    else if (is_current_face_convex())
    {
      if (m_points_of_face.size()==4)
      { convex_quadrangular_face_end_internal(normal); }
      else
      { convex_face_end_internal(normal); }
    }
    else
    { // Non convex and more than 3 points: we triangulate
      nonconvex_face_end_internal(normal);
    }

    m_face_started=false;
    m_points_of_face.clear();
    m_vertex_normals_for_face.clear();
  }

protected:
  void face_begin_internal(bool has_color, bool has_normal)
  {
    if (is_a_face_started())
    {
      std::cerr<<"You cannot start a new face before to finish the previous one."<<std::endl;
      return;
    }
    
    m_face_started=true;
    m_started_face_is_colored=has_color;
    m_started_face_has_normal=has_normal;
  }

  void triangular_face_end_internal(const Local_vector& normal)
  {
    for (int i=0; i<3; ++i)
    {
        // Add the position of the point, possibly with its color
      if (m_started_face_is_colored)
      { add_point(m_points_of_face[i], m_color_of_face); }
      else
      { add_point(m_points_of_face[i]); }
      
      // Add the flat normal
      add_flat_normal(normal);
      
      // Its smooth normal (if given by the user)
      if (m_vertex_normals_for_face.size()==3)
      { // Here we have 3 vertex normals; we can use Gourod
        add_gourod_normal(m_vertex_normals_for_face[i]);
      }
      else
      { // Here user does not provide all vertex normals: we use face normal istead
        // and thus we will not be able to use Gourod
        add_gourod_normal(normal);
      }
    }
  }
  
  void convex_quadrangular_face_end_internal(const Local_vector& normal)
  {
    if(faces_vector)
    {
      add_point(points_of_face[0], *faces_vector);
      add_point(points_of_face[1], *faces_vector);
      add_point(points_of_face[2], *faces_vector);
      
      add_point(points_of_face[0], *faces_vector);
      add_point(points_of_face[2], *faces_vector);
      add_point(points_of_face[3], *faces_vector);
    }
    else if(idx_faces_vector)
    {
      idx_faces_vector->push_back(point_index_map[ points_of_face[0] ]);
      idx_faces_vector->push_back(point_index_map[ points_of_face[1] ]);
      idx_faces_vector->push_back(point_index_map[ points_of_face[2] ]);
      
      idx_faces_vector->push_back(point_index_map[ points_of_face[0] ]);
      idx_faces_vector->push_back(point_index_map[ points_of_face[2] ]);
      idx_faces_vector->push_back(point_index_map[ points_of_face[3] ]);
    }
    if(!m_data_is_indexed)
      for(int i=0; i<6; ++i)
      {
        if (face_color_vector && m_started_face_is_colored)
          add_color(color_of_face, *face_color_vector);
        if(flat_normal_vector)
          add_normal(normal_for_face, *flat_normal_vector);
        if (smooth_normal_vector
            && vertex_normals_for_face.size()==4)
          add_normal(vertex_normals_for_face[0], *smooth_normal_vector);
        else if (smooth_normal_vector)
          add_normal(normal_for_face, *smooth_normal_vector);
      }
  }
  
  void convex_face_end_internal(const Local_vector& normal)
  {
    Local_point p0,p1,p2;
    std::size_t id;
    for(id = 1; id<points_of_face.size()-1; ++id)
    {
      
      p0 = points_of_face[0];
      p1 = points_of_face[id];
      p2 = points_of_face[id+1];
      if(faces_vector)
      {
        add_point(p0, *faces_vector);
        add_point(p1, *faces_vector);
        add_point(p2, *faces_vector);
      }
      else if(idx_faces_vector)
      {
        idx_faces_vector->push_back(point_index_map[p0]);
        idx_faces_vector->push_back(point_index_map[p1]);
        idx_faces_vector->push_back(point_index_map[p2]);
      }
      
      if(!m_data_is_indexed)
      {
        if (face_color_vector && m_started_face_is_colored)
          for (int i = 0; i<6; ++i)
          {
            add_color(color_of_face, *face_color_vector);
          }
        if(flat_normal_vector)
          add_normal(normal_for_face, *flat_normal_vector);
        if (smooth_normal_vector
            && vertex_normals_for_face.size()==4)
          add_normal(vertex_normals_for_face[0], *smooth_normal_vector);
        else if (smooth_normal_vector)
          add_normal(normal_for_face, *smooth_normal_vector);
      }
    }
  }
  
  void nonconvex_face_end_internal(const Local_vector& normal)
  {
    try
    {
      P_traits cdt_traits(normal);
      CDT cdt(cdt_traits);
      
      bool with_vertex_normal=(m_vertex_normals_for_face.size()==m_points_of_face.size());
      
      // (1) We insert all the edges as contraint in the CDT.
      typename CDT::Vertex_handle previous=NULL, first=NULL;
      for (int i=0; i<m_points_of_face.size(); ++i)
      {
        typename CDT::Vertex_handle vh = cdt.insert(m_points_of_face[i]);
        if(first==NULL)
        { first=vh; }
        
        if (with_vertex_normal)
        { vh->info().v=m_vertex_normals_for_face[i]; }
        else
        { vh->info().v=normal; }
        
        if(previous!=NULL && previous!=vh)
        { cdt.insert_constraint(previous, vh); }
        previous=vh;
      }
      
      if (previous!=NULL && previous!=first)
        cdt.insert_constraint(previous, first);
      
      // (2) We mark all external triangles
      // (2.1) We initialize is_external and is_process values 
      for(typename CDT::All_faces_iterator fit = cdt.all_faces_begin(),
            fitend = cdt.all_faces_end(); fit!=fitend; ++fit)
      {
        fit->info().is_external = true;
        fit->info().is_process = false;
      }
      // (2.2) We check if the facet is external or internal
      std::queue<typename CDT::Face_handle> face_queue;
      typename CDT::Face_handle face_internal = NULL;
      if (cdt.infinite_vertex()->face()!=NULL)
        face_queue.push(cdt.infinite_vertex()->face());
      while(! face_queue.empty() )
      {
        typename CDT::Face_handle fh = face_queue.front();
        face_queue.pop();
        if(!fh->info().is_process)
        {
          fh->info().is_process = true;
          for(int i=0; i<3; ++i)
          {
            if(!cdt.is_constrained(std::make_pair(fh, i)))
            {
              if (fh->neighbor(i)!=NULL)
                face_queue.push(fh->neighbor(i));
            }
            else if (face_internal==NULL)
            {
              face_internal = fh->neighbor(i);
            }
          }
        }
      }
      
      if ( face_internal!=NULL )
        face_queue.push(face_internal);
      
      while(! face_queue.empty() )
      {
        typename CDT::Face_handle fh = face_queue.front();
        face_queue.pop();
        if(!fh->info().is_process)
        {
          fh->info().is_process = true;
          fh->info().is_external = false;
          for(int i=0; i<3; ++i)
          {
            if(!cdt.is_constrained(std::make_pair(fh, i)))
            {
              if (fh->neighbor(i)!=NULL)
                face_queue.push(fh->neighbor(i));
            }
          }
        }
      }
      
      // (3) Now we iterates on the internal faces to add the vertices to the
      //     positions and the normals to the appropriate vectors
      for(typename CDT::Finite_faces_iterator ffit=cdt.finite_faces_begin(),
            ffitend = cdt.finite_faces_end(); ffit!=ffitend; ++ffit)
      {
        if(!ffit->info().is_external)
        {
          for(int i=0; i<3; ++i)
          {
            // Add the position of the point, possibly with its color
            if (m_started_face_is_colored)
            { add_point(ffit->vertex(i)->point(), m_color_of_face); }
            else
            { add_point(ffit->vertex(i)->point()); }
            
            // Add the flat normal
            add_flat_normal(normal);
            
            // Its smoth normal (if given by the user)
            add_gourod_normal(ffit->vertex(i)->info().v);
          }
        }
      }
    }
    catch(...)
    { // Triangulation crash: the face is not filled
      std::cout<<"Catch: face not filled."<<std::endl;
    }
  }
  
    bool is_current_face_convex()(const Local_vector& N)
   {

     typedef Local_kernel::Orientation Orientation;
     Orientation orientation;
     Local_vector normal = N;
     bool normal_is_ok;
     std::size_t id = 0;
     do{
       normal_is_ok = true;
       //Initializes the facet orientation

       Local_point S,T;
       S = facet[id];
       T = facet[id+1];
       Local_vector V1 = Local_vector((T-S).x(), (T-S).y(), (T-S).z());
       S = T;
       T = facet[id+2];
       Local_vector V2 = Local_vector((T-S).x(), (T-S).y(), (T-S).z());

       if(normal == Local_vector(0,0,0))
         normal_is_ok = false;
       if(normal_is_ok)
       {
         orientation = Local_kernel::Orientation_3()(V1, V2, normal);
         if( orientation == CGAL::COPLANAR )
           normal_is_ok = false;
       }
       //Checks if the normal is good : if the normal is null
       // or if it is coplanar to the facet, we need another one.
       if(!normal_is_ok)
       {
         normal = CGAL::cross_product(V1,V2);
       }

     }while( ++id != facet.size() && !normal_is_ok);
     //if no good normal can be found, stop here.
     if (!normal_is_ok)
       return false;

     //computes convexness

     //re-initializes he_end;

     for(id=0; id<facet.size(); ++id)
     {
       Local_point S,T;
       S = facet[id%facet.size()];
       T = facet[(id+1)%facet.size()];
       Local_vector V1 = Local_vector((T-S).x(), (T-S).y(), (T-S).z());
       S = T;
       T = facet[(id+2)%facet.size()];
       Local_vector V2 = Local_vector((T-S).x(), (T-S).y(), (T-S).z());
       Orientation res = Local_kernel::Orientation_3()(V1, V2, normal) ;

       if(res!= orientation && res != CGAL::ZERO)
         return false;
     }
     return true;
   }

  void add_color(const CGAL::Color& acolor)
  {
    if (m_color_buffer!=NULL)
    {
      m_color_buffer->push_back((float)acolor.red()/(float)255);
      m_color_buffer->push_back((float)acolor.green()/(float)255);
      m_color_buffer->push_back((float)acolor.blue()/(float)255);
    }
  }
  
  template<typename KVector>
  void add_normal_in_buffer(const KVector& kv, std::vector<float>* buffer)
  {
    if (buffer!=NULL)
    {
      Local_vector n=internal::get_local_vector(kv);
      buffer->push_back(n.x());
      buffer->push_back(n.y());
      buffer->push_back(n.z());
    }
  }

  template<typename KVector>
  void add_flat_normal(const KVector& kv)
  { add_normal_in_buffer(kv, m_flat_normal_buffer); }

  template<typename KVector>
  void add_gourod_normal(const KVector& kv)
  { add_normal_in_buffer(kv, m_gourod_normal_buffer); }

protected:
  // Types usefull for triangulation
  struct Vertex_info
  {
    Local_vector v;
  };

  struct Face_info
  {
    bool exist_edge[3];
    bool is_external;
    bool is_process;
  };

  typedef CGAL::Triangulation_2_projection_traits_3<CGAL::Exact_predicates_inexact_constructions_kernel> P_traits;
  typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, P_traits> Vb;
  typedef CGAL::Triangulation_face_base_with_info_2<Face_info, P_traits>     Fb1;
  typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>         Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                        TDS;
  typedef CGAL::Exact_predicates_tag                                         Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS, Itag>    CDT;

protected:
  std::vector<float>* m_pos_buffer;
  std::vector<float>* m_color_buffer;
  std::vector<float>* m_flat_normal_buffer;
  std::vector<float>* m_gourod_normal_buffer;
  std::vector<unsigned int>* m_index_buffer;

  CGAL::Bbox_3* m_bb;
  
  // Local variables, used when we started a new face.
  bool m_face_started;
  bool m_started_face_is_colored;
  bool m_started_face_has_normal;
  std::vector<Local_point> m_points_of_face;
  std::vector<Local_vector> m_vertex_normals_for_face;
  CGAL::Color m_color_of_face;
  Local_vector m_normal_of_face;
};

#endif // CGAL_BUFFER_FOR_VAO_H
