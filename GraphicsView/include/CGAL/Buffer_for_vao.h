// Copyright (c) 2018  GeometryFactory Sarl (France).
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
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_VBO_BUFFER_FILLER_H
#define CGAL_VBO_BUFFER_FILLER_H

#include <CGAL/license/GraphicsView.h>

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

namespace CGAL
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
  typedef Local_kernel::Point_3  Local_point;
  typedef Local_kernel::Vector_3 Local_vector;

//------------------------------------------------------------------------------
namespace internal
{
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

  inline
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
    return (Local_kernel::Construct_scaled_vector_3()(normal, 1.0/nb));
  }

  ////////////////////////////////////////////////////////////////
  // Structs to transform any CGAL point/vector into a Local_point/Local_vector
  template<typename K>
  struct Geom_utils
  {
    static Local_point get_local_point(const typename K::Point_2& p)
    {
      CGAL::Cartesian_converter<K, Local_kernel> converter;
      return Local_point(converter(p.x()), 0, converter(p.y()));
    }
    static Local_point get_local_point(const typename K::Weighted_point_2& p)
    {
      typename K::Point_2 lp(p);
      return Geom_utils<K>::get_local_point(lp);
    }
    static Local_point get_local_point(const typename K::Point_3& p)
    {
      CGAL::Cartesian_converter<K, Local_kernel> converter;
      return converter(p);
    }
    static Local_point get_local_point(const typename K::Weighted_point_3& p)
    {
      typename K::Point_3 lp(p);
      return Geom_utils<K>::get_local_point(lp);
    }
    static Local_vector get_local_vector(const typename K::Vector_2& v)
    {
      CGAL::Cartesian_converter<K, Local_kernel> converter;
      return Local_vector(converter(v.x()), 0, converter(v.y()));
    }
    static Local_vector get_local_vector(const typename K::Vector_3& v)
    {
      CGAL::Cartesian_converter<K, Local_kernel> converter;
      return converter(v);
    }
  };

  // Specialization for Local_kernel, because there is no need of convertion here.
  template<>
  struct Geom_utils<Local_kernel>
  {
    static Local_point get_local_point(const Local_kernel::Point_2& p)
    { return Local_point(p.x(), 0, p.y()); }
    static Local_point get_local_point(const Local_kernel::Weighted_point_2& p)
    { return Local_point(p.point().x(), 0, p.point().y());}
    static const Local_point & get_local_point(const Local_kernel::Point_3& p)
    { return p; }
    static Local_point get_local_point(const Local_kernel::Weighted_point_3& p)
    { return Local_point(p);}
    static Local_vector get_local_vector(const Local_kernel::Vector_2& v)
    { return Local_vector(v.x(), 0, v.y()); }
    static const Local_vector& get_local_vector(const Local_kernel::Vector_3& v)
    { return v; }
  };    

  ////////////////////////////////////////////////////////////////
  // Global function to simplify function calls.
  template<typename KPoint>
  Local_point get_local_point(const KPoint& p)
  {
    return Geom_utils<typename CGAL::Kernel_traits<KPoint>::Kernel>::
      get_local_point(p);
  }
  template<typename KVector>
  Local_vector get_local_vector(const KVector& v)
  {
    return Geom_utils<typename CGAL::Kernel_traits<KVector>::Kernel>::
      get_local_vector(v);
  }
} // End namespace internal

//------------------------------------------------------------------------------
template<typename BufferType=float, typename IndexType=std::size_t>
class Buffer_for_vao
{
public:
  Buffer_for_vao(std::vector<BufferType>* pos=NULL,
                 std::vector<IndexType>* indices=NULL,
                 CGAL::Bbox_3* bbox=NULL,
                 std::vector<BufferType>* color=NULL,
                 std::vector<BufferType>* flat_normal=NULL,
                 std::vector<BufferType>* gouraud_normal=NULL) :
    m_pos_buffer(pos),
    m_index_buffer(indices),
    m_color_buffer(color),
    m_flat_normal_buffer(flat_normal),
    m_gouraud_normal_buffer(gouraud_normal),
    m_bb(bbox),
    m_face_started(false)
  {}

  void clear()
  {
    if (m_pos_buffer!=NULL)            { m_pos_buffer->clear(); }
    if (m_color_buffer!=NULL)          { m_color_buffer->clear(); }
    if (m_index_buffer!=NULL)          { m_index_buffer->clear(); }
    if (m_flat_normal_buffer!=NULL)    { m_flat_normal_buffer->clear(); }
    if (m_gouraud_normal_buffer!=NULL) { m_gouraud_normal_buffer->clear(); }
  }

  bool is_empty() const
  {
    return
      (m_pos_buffer!=NULL && m_pos_buffer->empty()) &&
      (m_color_buffer!=NULL || m_color_buffer->empty()) &&
      (m_flat_normal_buffer!=NULL || m_flat_normal_buffer->empty()) &&
      (m_gouraud_normal_buffer!=NULL || m_gouraud_normal_buffer->empty()) &&
      (m_index_buffer!=NULL || m_index_buffer->empty());
  }

  bool has_position() const
  { return m_pos_buffer!=NULL; }

  bool has_indices() const
  { return m_index_buffer!=NULL; }

  bool has_color() const
  { return m_color_buffer!=NULL; }

  bool has_flat_normal() const
  { return m_flat_normal_buffer!=NULL; }
  
  bool has_gouraud_normal() const
  { return m_gouraud_normal_buffer!=NULL; }

  // 1.1) Add a point, without color. Return the index of the added point.
  template<typename KPoint>
  std::size_t add_point(const KPoint& kp)
  {
    if (!has_position()) return (std::size_t)-1;
    
    Local_point p=internal::get_local_point(kp);
    add_point_in_buffer(p, *m_pos_buffer);

    if (m_bb!=NULL)
    { (*m_bb)=(*m_bb)+p.bbox(); }

    return m_pos_buffer->size()-3;
  }

  // 1.2) Add a point, with color.
  template<typename KPoint>
  void add_point(const KPoint& kp, const CGAL::Color& c)
  {
    add_point(kp);
    add_color(c);
  }

  // 1.3) Add an indexed point, without color.
  template<typename T>
  void add_indexed_point(T index)
  {
    if (!has_indices()) return;
    m_index_buffer->push_back((IndexType)index);
  }

  // 2.1) Add a segment, without color.
  template<typename KPoint>
  void add_segment(const KPoint& kp1, const KPoint& kp2)
  {
    add_point(kp1);
    add_point(kp2);
  }
  
  // 2.2) Add a segment, with color.
  template<typename KPoint>
  void add_segment(const KPoint& kp1, const KPoint& kp2, const CGAL::Color& c)
  {
    add_segment(kp1, kp2);
    add_color(c);
    add_color(c);
  }

  // 2.3) Add an indexed segment, without color.
  template<typename T>
  void add_indexed_segment(T index1, T index2)
  {
    add_indexed_point(index1);
    add_indexed_point(index2);
  }

  /// @return true iff a face has begun.
  bool is_a_face_started() const
  { return m_face_started; }
  
  // 3.1) Add a face, without color, without normal.
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

  /// Add a point at the end of the current face, without giving the vertex normal.
  /// When this method is used, it is not possible to use the Gouraud shading.
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
  
  /// Add a point at the end of the current face
  /// @param p the point to add
  /// @p_normal the vertex normal at this point (for Gouraud shading)
  template<typename KPoint, typename KVector>
  bool add_point_in_face(const KPoint& kp, const KVector& p_normal)
  {
    if (add_point_in_face(kp))
    {
      m_vertex_normals_for_face.push_back(internal::get_local_vector(p_normal));
      return true;
    }
    return false;
  }

  /// Add an indexed point at the end of the current face, without giving the vertex normal.
  /// When Indexation is used, it is not possible to use flat shading or multiple colors 
  /// for face sor edges.
  /// Note that we still need the point itself, in order to triangulate the face when necessary.
  template<typename T, typename KPoint>
  bool add_indexed_point_in_face(T index, const KPoint& kp)
  {
    if (add_point_in_face(kp))
    {
      m_indices_of_points_of_face.push_back(index);
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
      std::cerr<<"PB: you try to triangulate a face with "<<m_points_of_face.size()<<" vertices."
               <<std::endl;
      
      m_face_started=false;
      m_points_of_face.clear();
      m_vertex_normals_for_face.clear();

      return;
    }

    if (m_indices_of_points_of_face.size()>0 &&
        m_indices_of_points_of_face.size()!=m_points_of_face.size())
    {
      std::cerr<<"PB: you mixed some add_point_in_face(...) and some add_indexed_point_in_face(...)"
               <<" for a same face. Indices for this face are ignored."<<std::endl;
      m_indices_of_points_of_face.clear();
    }

    if (m_vertex_normals_for_face.size()>0 &&
        m_vertex_normals_for_face.size()!=m_points_of_face.size())
    {
      std::cerr<<"PB: you only gave some vertex normals (and not all) for a same face. "
               <<"All vertex normal are ignored and thus it is not possible to use Gouraud "
               <<"shading for this face."
               <<std::endl;
      m_vertex_normals_for_face.clear();
    }
    
    Local_vector normal=(m_started_face_has_normal?m_normal_of_face:
                         internal::compute_normal_of_face(m_points_of_face));

    if (m_points_of_face.size()==3)
    { triangular_face_end_internal(normal); } // Triangle: no need to triangulate
    else if (is_current_face_convex(normal))
    {
      if (m_points_of_face.size()==4)
      { convex_quadrangular_face_end_internal(normal); } // Convex quad
      else 
      { convex_face_end_internal(normal); } // Convex face with > 4 vertices
    }
    else
    { // Non convex and more than 3 points: we triangulate
      nonconvex_face_end_internal(normal);
    }

    m_face_started=false;
  }

  /// adds `kp` coordinates to `buffer`
  template<typename KPoint>
  static void add_point_in_buffer(const KPoint& kp, std::vector<float>& buffer)
  {
    Local_point p=internal::get_local_point(kp);
    buffer.push_back(p.x());
    buffer.push_back(p.y());
    buffer.push_back(p.z());
  }

  /// adds `kv` coordinates to `buffer`
  template<typename KVector>
  static void add_normal_in_buffer(const KVector& kv, std::vector<float>& buffer)
  {
    Local_vector n=internal::get_local_vector(kv);
    buffer.push_back(n.x());
    buffer.push_back(n.y());
    buffer.push_back(n.z());
  }

  ///adds `acolor` RGB components to `buffer`
  static void add_color_in_buffer(const CGAL::Color& acolor, std::vector<float>& buffer)
  {
    buffer.push_back((float)acolor.red()/(float)255);
    buffer.push_back((float)acolor.green()/(float)255);
    buffer.push_back((float)acolor.blue()/(float)255);
  }

  /// @return true iff the points of 'facet' form a convex face
  static bool is_facet_convex(const std::vector<Local_point>& facet,
                              const Local_vector& normal)
  {
    Local_kernel::Orientation orientation, local_orientation;
    std::size_t id=0;
    do
    {
      const Local_point& S=facet[id];
      const Local_point& T=facet[(id+1==facet.size())?0:id+1];
      Local_vector V1=Local_vector((T-S).x(), (T-S).y(), (T-S).z());

      const Local_point& U=facet[(id+2==facet.size())?0:id+2];
      Local_vector V2=Local_vector((U-T).x(), (U-T).y(), (U-T).z());

      orientation = Local_kernel::Orientation_3()(V1, V2, normal);
      // Is it possible that orientation==COPLANAR ? Maybe if V1 or V2 is very small ?
    }
    while(++id!=facet.size() &&
          (orientation==CGAL::COPLANAR || orientation==CGAL::ZERO));

    //Here, all orientations were COPLANAR. Not sure this case is possible,
    // but we stop here.
    if (orientation==CGAL::COPLANAR || orientation==CGAL::ZERO)
    { return false; }

    // Now we compute convexness
    for(id=0; id<facet.size(); ++id)
    {
      const Local_point& S=facet[id];
      const Local_point& T=facet[(id+1==facet.size())?0:id+1];
      Local_vector V1=Local_vector((T-S).x(), (T-S).y(), (T-S).z());
      
      const Local_point& U=facet[(id+2==facet.size())?0:id+2];
      Local_vector V2=Local_vector((U-T).x(), (U-T).y(), (U-T).z());
      
      local_orientation=Local_kernel::Orientation_3()(V1, V2, normal) ;

      if(local_orientation!=CGAL::ZERO && local_orientation!=orientation)
      { return false; }
    }
    return true;
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

    m_points_of_face.clear();
    m_vertex_normals_for_face.clear();
    m_indices_of_points_of_face.clear();
  }

  void triangular_face_end_internal(const Local_vector& normal)
  {
    for (int i=0; i<3; ++i)
    {
      // If user gave vertex indices
      if (m_indices_of_points_of_face.size()>0)
      { 
	add_indexed_point(m_indices_of_points_of_face[i]); 
	}
      else
      {
        add_point(m_points_of_face[i]); // Add the position of the point
      if (m_started_face_is_colored)
      { add_color(m_color_of_face); } // Add the color      
      add_flat_normal(normal); // Add the flat normal
      // Its smooth normal (if given by the user)
      if (m_vertex_normals_for_face.size()>0)
      { // Here we have 3 vertex normals; we can use Gouraud
        add_gouraud_normal(m_vertex_normals_for_face[i]);
      }
      else
      { // Here user does not provide all vertex normals: we use face normal istead
        // and thus we will not be able to use Gouraud
        add_gouraud_normal(normal);
      }
      }
    }
  }
  
  void convex_quadrangular_face_end_internal(const Local_vector& normal)
  {
    // Add indices when they exist
    if (m_indices_of_points_of_face.size()>0)
    {
      add_indexed_point(m_indices_of_points_of_face[0]);
      add_indexed_point(m_indices_of_points_of_face[1]);
      add_indexed_point(m_indices_of_points_of_face[2]);
      
      add_indexed_point(m_indices_of_points_of_face[0]);
      add_indexed_point(m_indices_of_points_of_face[2]);
      add_indexed_point(m_indices_of_points_of_face[3]);
    }
    else
    {
      // (1) add points
      add_point(m_points_of_face[0]);
      add_point(m_points_of_face[1]);
      add_point(m_points_of_face[2]);

      add_point(m_points_of_face[0]);
      add_point(m_points_of_face[2]);
      add_point(m_points_of_face[3]);

      // (2) Add flat and smooth normals and color
      for(unsigned int i=0; i<6; ++i)
      {
        if (m_started_face_is_colored)
        { add_color(m_color_of_face); }

        add_flat_normal(normal);

        if (m_vertex_normals_for_face.size()==0)
        { add_gouraud_normal(normal); }   
      }

      if (m_vertex_normals_for_face.size()>0)
      {
        add_gouraud_normal(m_vertex_normals_for_face[0]);
        add_gouraud_normal(m_vertex_normals_for_face[1]);
        add_gouraud_normal(m_vertex_normals_for_face[2]);
        
        add_gouraud_normal(m_vertex_normals_for_face[0]);
        add_gouraud_normal(m_vertex_normals_for_face[2]);
        add_gouraud_normal(m_vertex_normals_for_face[3]);
      }
    }
  }
  
  void convex_face_end_internal(const Local_vector& normal)
  {
    for(std::size_t i=1; i<m_points_of_face.size()-1; ++i)
    {
      // Add indices when they exist
      if (m_indices_of_points_of_face.size()>0)
      {
        add_indexed_point(m_indices_of_points_of_face[0]);
        add_indexed_point(m_indices_of_points_of_face[i]);
        add_indexed_point(m_indices_of_points_of_face[i+1]);
      }
      else
      {
        Local_point& p0 = m_points_of_face[0];
        Local_point& p1 = m_points_of_face[i];
        Local_point& p2 = m_points_of_face[i+1];

        // (1) add points
        add_point(p0);
        add_point(p1);
        add_point(p2);

        // (2) Add flat normal and color
        for(unsigned int j=0; j<3; ++j)
        {
          if (m_started_face_is_colored)
          { add_color(m_color_of_face); }

          add_flat_normal(normal);

          if (m_vertex_normals_for_face.size()==0)
          { add_gouraud_normal(normal); } // No smooth normal, we use the flat one instead
        }

        // (3) Add smooth normals if they exist
        if (m_vertex_normals_for_face.size()>0)
        {
          add_gouraud_normal(m_vertex_normals_for_face[0]);
          add_gouraud_normal(m_vertex_normals_for_face[i]);
          add_gouraud_normal(m_vertex_normals_for_face[i+1]);
        }
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
      for (unsigned int i=0; i<m_points_of_face.size(); ++i)
      {
        typename CDT::Vertex_handle vh = cdt.insert(m_points_of_face[i]);
        if(first==NULL)
        { first=vh; }
        
        if (with_vertex_normal)
        { vh->info().v=m_vertex_normals_for_face[i]; }
        else
        { vh->info().v=normal; }

        if (m_indices_of_points_of_face.size()>0)
        { vh->info().index=m_indices_of_points_of_face[i]; }
        
        if(previous!=NULL && previous!=vh)
        { cdt.insert_constraint(previous, vh); }
        previous=vh;
      }
      
      if (previous!=NULL && previous!=first)
      { cdt.insert_constraint(previous, first); }
      
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
      { face_queue.push(cdt.infinite_vertex()->face()); }
      while(!face_queue.empty())
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
              { face_queue.push(fh->neighbor(i)); }
            }
            else if (face_internal==NULL)
            {
              face_internal = fh->neighbor(i);
            }
          }
        }
      }
      
      if ( face_internal!=NULL )
      { face_queue.push(face_internal); }
      
      while(!face_queue.empty())
      {
        typename CDT::Face_handle fh = face_queue.front();
        face_queue.pop();
        if(!fh->info().is_process)
        {
          fh->info().is_process = true;
          fh->info().is_external = false;
          for(unsigned int i=0; i<3; ++i)
          {
            if(!cdt.is_constrained(std::make_pair(fh, i)))
            {
              if (fh->neighbor(i)!=NULL)
              { face_queue.push(fh->neighbor(i)); }
            }
          }
        }
      }
      
      // (3) Now we iterates on the internal faces to add the vertices 
      //     and the normals to the appropriate vectors
      for(typename CDT::Finite_faces_iterator ffit=cdt.finite_faces_begin(),
            ffitend = cdt.finite_faces_end(); ffit!=ffitend; ++ffit)
      {
        if(!ffit->info().is_external)
        {
          for(unsigned int i=0; i<3; ++i)
          {            
            // Add indices when they exist
            if (m_indices_of_points_of_face.size()>0)
            { add_indexed_point(ffit->vertex(i)->info().index); }
            else
            {
              // (1) add point
              add_point(ffit->vertex(i)->point());
              // (2) Add face color
              if (m_started_face_is_colored)
              { add_color(m_color_of_face); }

              // (3) Add flat normal
              add_flat_normal(normal);

              // (4) Add smooth normals (or flat if smooth normals do not exist)
              add_gouraud_normal(ffit->vertex(i)->info().v);
            }
          }
        }
      }
    }
    catch(...)
    { // Triangulation crash: the face is not filled
      std::cerr<<"Catch: face not filled."<<std::endl;
    }
  }

  /// @return true iff the current face is convex
  bool is_current_face_convex(const Local_vector& N) const
  {
    return is_facet_convex(m_points_of_face, N);
  }

  void add_color(const CGAL::Color& acolor)
  {
    if (m_color_buffer!=NULL)
    { add_color_in_buffer(acolor, *m_color_buffer); }
  }
  
  template<typename KVector>
  void add_flat_normal(const KVector& kv)
  {
    if(m_flat_normal_buffer != NULL)
    { add_normal_in_buffer(kv, *m_flat_normal_buffer); }
  }

  template<typename KVector>
  void add_gouraud_normal(const KVector& kv)
  {
    if(m_gouraud_normal_buffer != NULL)
    { add_normal_in_buffer(kv, *m_gouraud_normal_buffer); }
  }

protected:
  // Types usefull for triangulation
  struct Vertex_info
  {
    Local_vector v;
    std::size_t index;
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
  std::vector<BufferType>* m_pos_buffer;
  std::vector<IndexType>* m_index_buffer;
  std::vector<BufferType>* m_color_buffer;
  std::vector<BufferType>* m_flat_normal_buffer;
  std::vector<BufferType>* m_gouraud_normal_buffer;

  CGAL::Bbox_3* m_bb;
  
  // Local variables, used when we started a new face.
  bool m_face_started;
  bool m_started_face_is_colored;
  bool m_started_face_has_normal;
  std::vector<Local_point> m_points_of_face;
  std::vector<Local_vector> m_vertex_normals_for_face;
  std::vector<std::size_t> m_indices_of_points_of_face;
  CGAL::Color m_color_of_face;
  Local_vector m_normal_of_face;
};

} // End namespace CGAL

#endif // CGAL_VBO_BUFFER_FILLER_H
