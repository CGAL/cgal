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
// Author(s)     : Simon Giraudot

#ifndef CGAL_SURFACE_MESH_IO_PLY
#define CGAL_SURFACE_MESH_IO_PLY

#include <CGAL/IO/PLY.h>

namespace CGAL {

namespace internal {

#if !defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE) && !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES)
namespace PLY {

template <typename Point>
class Surface_mesh_filler
{
public:
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  typedef typename Kernel::Vector_3 Vector;
  typedef Surface_mesh<Point> Surface_mesh;
  typedef typename Surface_mesh::Vertex_index Vertex_index;
  typedef typename Surface_mesh::Face_index Face_index;

private:

  struct Abstract_ply_property_to_surface_mesh_property
  {
    virtual ~Abstract_ply_property_to_surface_mesh_property() { }
    virtual void assign (PLY_element& element, std::size_t index) = 0;
  };

  template <typename Simplex, typename Type>
  class PLY_property_to_surface_mesh_property : public Abstract_ply_property_to_surface_mesh_property
  {
    typedef typename Surface_mesh::template Property_map<Simplex, Type> Map;
    Map m_map;
    std::string m_name;
  public:
    PLY_property_to_surface_mesh_property (Surface_mesh& sm, const std::string& name)
      : m_name (name)
    {
      m_map = sm.template add_property_map<Simplex, Type>(prefix(Simplex()) + name).first;
    }
    
    virtual void assign (PLY_element& element, std::size_t index)
    {
      Type t{};
      element.assign (t, m_name.c_str());
      put(m_map, Simplex(index), t);
    }

    std::string prefix(Vertex_index) const { return "v:"; }
    std::string prefix(Face_index) const { return "f:"; }
  };
  
  Surface_mesh& m_mesh;
  bool m_use_floats;
  int m_normals;
  typename Surface_mesh::template Property_map<Vertex_index, Vector> m_normal_map;
  int m_vcolors;
  typename Surface_mesh::template Property_map<Vertex_index, CGAL::Color> m_vcolor_map;
  int m_fcolors;
  typename Surface_mesh::template Property_map<Face_index, CGAL::Color> m_fcolor_map;
  bool m_use_int32_t;
  std::string m_index_tag;
  std::vector<Abstract_ply_property_to_surface_mesh_property*> m_vertex_properties;
  std::vector<Abstract_ply_property_to_surface_mesh_property*> m_face_properties;

public:
  
  Surface_mesh_filler (Surface_mesh& mesh)
    : m_mesh (mesh), m_use_floats (false), m_normals(0), m_vcolors(0), m_fcolors(0)
  { }

  ~Surface_mesh_filler()
  {
    for (std::size_t i = 0; i < m_vertex_properties.size(); ++ i)
      delete m_vertex_properties[i];
    for (std::size_t i = 0; i < m_face_properties.size(); ++ i)
      delete m_face_properties[i];
  }

  bool has_simplex_specific_property (internal::PLY::PLY_read_number* property, Vertex_index)
  {
    const std::string& name = property->name();
    if (name == "x" ||
        name == "y" ||
        name == "z")
    {
      if (dynamic_cast<PLY_read_typed_number<float>*>(property))
        m_use_floats = true;
      return true;
    }
    if (name == "nx" ||
        name == "ny" ||
        name == "nz")
    {
      ++ m_normals;
      if (m_normals == 3)
        m_normal_map = m_mesh.template add_property_map<Vertex_index, Vector>("v:normal").first;
      return true;
    }
    if (name == "red" ||
        name == "green" ||
        name == "blue")
    {
      ++ m_vcolors;
      if (m_vcolors == 3)
        m_vcolor_map = m_mesh.template add_property_map<Vertex_index, CGAL::Color>("v:color").first;
      return true;
    }
    return false;
  }

  bool has_simplex_specific_property (internal::PLY::PLY_read_number* property, Face_index)
  {
    const std::string& name = property->name();
    if (name == "vertex_indices" || name == "vertex_index")
    {
      m_index_tag  = name;
      m_use_int32_t = dynamic_cast<PLY_read_typed_list<boost::int32_t>*>(property);
      CGAL_assertion (dynamic_cast<PLY_read_typed_list<boost::uint32_t>*>(property));
      return true;
    }
    if (name == "red" ||
        name == "green" ||
        name == "blue")
    {
      ++ m_fcolors;
      if (m_fcolors == 3)
        m_fcolor_map = m_mesh.template add_property_map<Face_index, CGAL::Color>("f:color").first;
      return true;
    }
    
    return false;
  }

  void instantiate_vertex_properties (PLY_element& element)
  {
    instantiate_properties<Vertex_index> (element, m_vertex_properties);
  }
  
  void instantiate_face_properties (PLY_element& element)
  {
    instantiate_properties<Face_index> (element, m_face_properties);
  }
  
  template <typename Simplex>
  void instantiate_properties  (PLY_element& element,
                                std::vector<Abstract_ply_property_to_surface_mesh_property*>& properties)
  {
    for (std::size_t j = 0; j < element.number_of_properties(); ++ j)
    {
      internal::PLY::PLY_read_number* property = element.property(j);

      if (has_simplex_specific_property (property, Simplex()))
        continue;
      
      const std::string& name = property->name();
      
      if (dynamic_cast<PLY_read_typed_number<boost::int8_t>*>(property))
      {
        properties.push_back
          (new PLY_property_to_surface_mesh_property<Simplex, boost::int8_t>(m_mesh,
                                                                             name));
      }
      else if (dynamic_cast<PLY_read_typed_number<boost::uint8_t>*>(property))
      {
        properties.push_back
          (new PLY_property_to_surface_mesh_property<Simplex, boost::uint8_t>(m_mesh,
                                                                              name));
      }
      else if (dynamic_cast<PLY_read_typed_number<boost::int16_t>*>(property))
      {
        properties.push_back
          (new PLY_property_to_surface_mesh_property<Simplex, boost::int16_t>(m_mesh,
                                                                              name));
      }
      else if (dynamic_cast<PLY_read_typed_number<boost::uint16_t>*>(property))
      {
        properties.push_back
          (new PLY_property_to_surface_mesh_property<Simplex, boost::uint16_t>(m_mesh,
                                                                               name));
      }
      else if (dynamic_cast<PLY_read_typed_number<boost::int32_t>*>(property))
      {
        properties.push_back
          (new PLY_property_to_surface_mesh_property<Simplex, boost::int32_t>(m_mesh,
                                                                              name));
      }
      else if (dynamic_cast<PLY_read_typed_number<boost::uint32_t>*>(property))
      {
        properties.push_back
          (new PLY_property_to_surface_mesh_property<Simplex, boost::uint32_t>(m_mesh,
                                                                               name));
      }
      else if (dynamic_cast<PLY_read_typed_number<float>*>(property))
      {
        properties.push_back
          (new PLY_property_to_surface_mesh_property<Simplex, float>(m_mesh,
                                                                     name));
      }
      else if (dynamic_cast<PLY_read_typed_number<double>*>(property))
      {
        properties.push_back
          (new PLY_property_to_surface_mesh_property<Simplex, double>(m_mesh,
                                                                      name));
      }
    }
  }

  void process_vertex_line (PLY_element& element)
  {
    Vertex_index vi;
    
    if (m_use_floats)
      process_line<float>(element, vi);
    else
      process_line<double>(element, vi);

    for (std::size_t i = 0; i < m_vertex_properties.size(); ++ i)
      m_vertex_properties[i]->assign (element, vi);
  }

  template <typename FT>
  void process_line (PLY_element& element, Vertex_index& vi)
  {
    FT x = (FT)0.,y = (FT)0., z = (FT)0.,
      nx = (FT)0., ny = (FT)0., nz = (FT)0.;
    element.assign (x, "x");
    element.assign (y, "y");
    element.assign (z, "z");
    Point point (x, y, z);
    vi = m_mesh.add_vertex(point);
    
    if (m_normals == 3)
    {
      element.assign (nx, "nx");
      element.assign (ny, "ny");
      element.assign (nz, "nz");
      Vector normal (nx, ny, nz);
      m_normal_map[vi] = normal;
    }

    if (m_vcolors == 3)
    {
      unsigned char r, g, b;
      element.assign (r, "red");
      element.assign (g, "green");
      element.assign (b, "blue");
      m_vcolor_map[vi] = CGAL::Color (r, g, b);
    }
  }

  bool process_face_line (PLY_element& element)
  {
    Face_index fi;

    if (m_use_int32_t)
      process_line<boost::int32_t>(element, fi);
    else
      process_line<boost::uint32_t>(element, fi);

    if (fi == Surface_mesh::null_face())
      return false;

    for (std::size_t i = 0; i < m_face_properties.size(); ++ i)
      m_face_properties[i]->assign (element, fi);

    return true;
  }

  template <typename IntType>
  void process_line (PLY_element& element, Face_index& fi)
  {
    std::vector<IntType> indices;
    element.assign (indices, m_index_tag.c_str());
    std::vector<Vertex_index> vertices;
    vertices.reserve(indices.size());
    for (std::size_t i = 0; i < indices.size(); ++ i)
      vertices.push_back (Vertex_index (indices[i]));

    fi = m_mesh.add_face(vertices);

    if (m_fcolors == 3)
    {
      unsigned char r, g, b;
      element.assign (r, "red");
      element.assign (g, "green");
      element.assign (b, "blue");
      m_fcolor_map[fi] = CGAL::Color (r, g, b);
    }
  }

};


} // namespace PLY
#endif

} // namespace internal

} // namespace CGAL


#endif // CGAL_SURFACE_MESH_IO_PLY
