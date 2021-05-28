// Copyright (c) 2018  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_SURFACE_MESH_IO_PLY_H
#define CGAL_SURFACE_MESH_IO_PLY_H

#include <CGAL/license/Surface_mesh.h>

#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#if !defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE) && !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES)

#include <CGAL/IO/PLY.h>

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#define CGAL_DEPRECATED
#endif

namespace CGAL {
namespace IO {
namespace internal {

template <typename Point>
class Surface_mesh_filler
{
public:
  typedef typename Kernel_traits<Point>::Kernel                             Kernel;
  typedef typename Kernel::Vector_3                                         Vector;
  typedef CGAL::Surface_mesh<Point>                                         Surface_mesh;
  typedef typename Surface_mesh::size_type                                  size_type;
  typedef typename Surface_mesh::Vertex_index                               Vertex_index;
  typedef typename Surface_mesh::Face_index                                 Face_index;
  typedef typename Surface_mesh::Edge_index                                 Edge_index;
  typedef typename Surface_mesh::Halfedge_index                             Halfedge_index;

private:
  struct Abstract_ply_property_to_surface_mesh_property
  {
    virtual ~Abstract_ply_property_to_surface_mesh_property() { }
    virtual void assign(PLY_element& element, size_type index) = 0;
  };

  template <typename Simplex, typename Type>
  class PLY_property_to_surface_mesh_property
      : public Abstract_ply_property_to_surface_mesh_property
  {
    typedef typename Surface_mesh::template Property_map<Simplex, Type> Map;
    Map m_map;
    std::string m_name;

  public:
    PLY_property_to_surface_mesh_property(Surface_mesh& sm, const std::string& name)
      : m_name(name)
    {
      m_map = sm.template add_property_map<Simplex, Type>(prefix(Simplex()) + name).first;
    }

    virtual void assign(PLY_element& element, size_type index)
    {
      Type t{};
      element.assign(t, m_name.c_str());
      put(m_map, Simplex(index), t);
    }

    std::string prefix(Vertex_index) const { return "v:"; }
    std::string prefix(Face_index) const { return "f:"; }
    std::string prefix(Edge_index) const { return "e:"; }
    std::string prefix(Halfedge_index) const { return "h:"; }
  };

  Surface_mesh& m_mesh;
  std::vector<Vertex_index> m_map_v2v;
  bool m_use_floats;
  int m_normals;
  typename Surface_mesh::template Property_map<Vertex_index, Vector> m_normal_map;
  int m_vcolors;
  typename Surface_mesh::template Property_map<Vertex_index, CGAL::IO::Color> m_vcolor_map;
  int m_fcolors;
  typename Surface_mesh::template Property_map<Face_index, CGAL::IO::Color> m_fcolor_map;
  bool m_use_int32_t;
  std::string m_index_tag;
  std::vector<Abstract_ply_property_to_surface_mesh_property*> m_vertex_properties;
  std::vector<Abstract_ply_property_to_surface_mesh_property*> m_face_properties;
  std::vector<Abstract_ply_property_to_surface_mesh_property*> m_edge_properties;
  std::vector<Abstract_ply_property_to_surface_mesh_property*> m_halfedge_properties;

public:
  Surface_mesh_filler(Surface_mesh& mesh)
    : m_mesh(mesh), m_use_floats(false), m_normals(0), m_vcolors(0), m_fcolors(0)
  { }

  ~Surface_mesh_filler()
  {
    for(std::size_t i = 0; i < m_vertex_properties.size(); ++i)
      delete m_vertex_properties[i];
    for(std::size_t i = 0; i < m_face_properties.size(); ++i)
      delete m_face_properties[i];
    for(std::size_t i = 0; i < m_edge_properties.size(); ++i)
      delete m_edge_properties[i];
    for(std::size_t i = 0; i < m_halfedge_properties.size(); ++i)
      delete m_halfedge_properties[i];
  }

  bool has_simplex_specific_property(internal::PLY_read_number* property, Vertex_index)
  {
    const std::string& name = property->name();
    if(name == "x" ||
       name == "y" ||
       name == "z")
    {
      if(dynamic_cast<PLY_read_typed_number<float>*>(property))
        m_use_floats = true;
      return true;
    }
    if(name == "nx" ||
       name == "ny" ||
       name == "nz")
    {
      ++ m_normals;
      if(m_normals == 3)
        m_normal_map = m_mesh.template add_property_map<Vertex_index, Vector>("v:normal").first;
      return true;
    }
    if(name == "red" ||
       name == "green" ||
       name == "blue")
    {
      ++ m_vcolors;
      if(m_vcolors == 3)
        m_vcolor_map = m_mesh.template add_property_map<Vertex_index, CGAL::IO::Color>("v:color").first;
      return true;
    }
    return false;
  }

  bool has_simplex_specific_property(internal::PLY_read_number* property, Face_index)
  {
    const std::string& name = property->name();
    if(name == "vertex_indices" || name == "vertex_index")
    {
      CGAL_assertion(dynamic_cast<PLY_read_typed_list<boost::int32_t>*>(property)
                     || dynamic_cast<PLY_read_typed_list<boost::uint32_t>*>(property));
      m_index_tag  = name;
      m_use_int32_t = dynamic_cast<PLY_read_typed_list<boost::int32_t>*>(property);
      return true;
    }
    if(name == "red" ||
       name == "green" ||
       name == "blue")
    {
      ++ m_fcolors;
      if(m_fcolors == 3)
        m_fcolor_map = m_mesh.template add_property_map<Face_index, CGAL::IO::Color>("f:color").first;
      return true;
    }

    return false;
  }

  bool has_simplex_specific_property(internal::PLY_read_number* property, Edge_index)
  {
    const std::string& name = property->name();
    if(name == "v0" || name == "v1")
      return true;
    return false;
  }

  bool has_simplex_specific_property(internal::PLY_read_number* property, Halfedge_index)
  {
    const std::string& name = property->name();
    if(name == "source" || name == "target")
      return true;
    return false;
  }

  void instantiate_vertex_properties(PLY_element& element)
  {
    m_map_v2v.reserve(element.number_of_items());
    instantiate_properties<Vertex_index>(element, m_vertex_properties);
  }

  void instantiate_face_properties(PLY_element& element)
  {
    instantiate_properties<Face_index>(element, m_face_properties);
  }

  void instantiate_edge_properties(PLY_element& element)
  {
    instantiate_properties<Edge_index>(element, m_edge_properties);
  }

  void instantiate_halfedge_properties(PLY_element& element)
  {
    instantiate_properties<Halfedge_index>(element, m_halfedge_properties);
  }

  template <typename Simplex>
  void instantiate_properties(PLY_element& element,
                              std::vector<Abstract_ply_property_to_surface_mesh_property*>& properties)
  {
    for(std::size_t j = 0; j < element.number_of_properties(); ++ j)
    {
      internal::PLY_read_number* property = element.property(j);

      if(has_simplex_specific_property(property, Simplex()))
        continue;

      const std::string& name = property->name();

      if(dynamic_cast<PLY_read_typed_number<boost::int8_t>*>(property))
      {
        properties.push_back(new PLY_property_to_surface_mesh_property<Simplex, boost::int8_t>(m_mesh, name));
      }
      else if(dynamic_cast<PLY_read_typed_number<boost::uint8_t>*>(property))
      {
        properties.push_back(new PLY_property_to_surface_mesh_property<Simplex, boost::uint8_t>(m_mesh, name));
      }
      else if(dynamic_cast<PLY_read_typed_number<boost::int16_t>*>(property))
      {
        properties.push_back(new PLY_property_to_surface_mesh_property<Simplex, boost::int16_t>(m_mesh, name));
      }
      else if(dynamic_cast<PLY_read_typed_number<boost::uint16_t>*>(property))
      {
        properties.push_back(new PLY_property_to_surface_mesh_property<Simplex, boost::uint16_t>(m_mesh, name));
      }
      else if(dynamic_cast<PLY_read_typed_number<boost::int32_t>*>(property))
      {
        properties.push_back(new PLY_property_to_surface_mesh_property<Simplex, boost::int32_t>(m_mesh, name));
      }
      else if(dynamic_cast<PLY_read_typed_number<boost::uint32_t>*>(property))
      {
        properties.push_back(new PLY_property_to_surface_mesh_property<Simplex, boost::uint32_t>(m_mesh, name));
      }
      else if(dynamic_cast<PLY_read_typed_number<float>*>(property))
      {
        properties.push_back(new PLY_property_to_surface_mesh_property<Simplex, float>(m_mesh, name));
      }
      else if(dynamic_cast<PLY_read_typed_number<double>*>(property))
      {
        properties.push_back(new PLY_property_to_surface_mesh_property<Simplex, double>(m_mesh, name));
      }
    }
  }

  void process_vertex_line(PLY_element& element)
  {
    Vertex_index vi;

    if(m_use_floats)
      process_line<float>(element, vi);
    else
      process_line<double>(element, vi);

    for(std::size_t i = 0; i < m_vertex_properties.size(); ++i)
      m_vertex_properties[i]->assign(element, vi);
  }

  template <typename FT>
  void process_line(PLY_element& element, Vertex_index& vi)
  {
    FT x = 0, y = 0, z = 0,
        nx = 0, ny = 0, nz = 0;
    element.assign(x, "x");
    element.assign(y, "y");
    element.assign(z, "z");
    Point point(x, y, z);
    vi = m_mesh.add_vertex(point);
    m_map_v2v.push_back(vi);

    if(m_normals == 3)
    {
      element.assign(nx, "nx");
      element.assign(ny, "ny");
      element.assign(nz, "nz");
      Vector normal(nx, ny, nz);
      m_normal_map[vi] = normal;
    }

    if(m_vcolors == 3)
    {
      unsigned char r, g, b;
      element.assign(r, "red");
      element.assign(g, "green");
      element.assign(b, "blue");
      m_vcolor_map[vi] = CGAL::IO::Color(r, g, b);
    }
  }

  bool process_face_line(PLY_element& element)
  {
    Face_index fi = m_mesh.null_face();

    if(m_use_int32_t)
      process_line<boost::int32_t>(element, fi);
    else
      process_line<boost::uint32_t>(element, fi);

    if(fi == Surface_mesh::null_face())
      return false;

    for(std::size_t i = 0; i < m_face_properties.size(); ++i)
      m_face_properties[i]->assign(element, fi);

    return true;
  }

  template <typename IntType>
  void process_line(PLY_element& element, Face_index& fi)
  {
    std::vector<IntType> indices;
    element.assign(indices, m_index_tag.c_str());
    std::vector<Vertex_index> vertices;
    vertices.reserve(indices.size());
    for(std::size_t i = 0; i < indices.size(); ++i)
      vertices.push_back(m_map_v2v[std::size_t(indices[i])]);

    fi = m_mesh.add_face(vertices);
    if(fi == m_mesh.null_face())
      return;

    if(m_fcolors == 3)
    {
      unsigned char r, g, b;
      element.assign(r, "red");
      element.assign(g, "green");
      element.assign(b, "blue");
      m_fcolor_map[fi] = CGAL::IO::Color(r, g, b);
    }
  }

  bool process_edge_line(PLY_element& element)
  {
    Edge_index ei = m_mesh.null_edge();

    if(m_use_int32_t)
      process_line<boost::int32_t>(element, ei);
    else
      process_line<boost::uint32_t>(element, ei);

    if(ei == Surface_mesh::null_edge())
      return false;

    for(std::size_t i = 0; i < m_edge_properties.size(); ++i)
      m_edge_properties[i]->assign(element, ei);

    return true;
  }

  template <typename IntType>
  void process_line(PLY_element& element, Edge_index& ei)
  {
    IntType v0, v1;
    element.assign(v0, "v0");
    element.assign(v1, "v1");

    Halfedge_index hi = m_mesh.halfedge(m_map_v2v[std::size_t(v0)],
        m_map_v2v[std::size_t(v1)]);
    if(hi == m_mesh.null_halfedge())
      return;

    ei = m_mesh.edge(hi);
  }

  bool process_halfedge_line(PLY_element& element)
  {
    Halfedge_index hi = m_mesh.null_halfedge();

    if(m_use_int32_t)
      process_line<boost::int32_t>(element, hi);
    else
      process_line<boost::uint32_t>(element, hi);

    if(hi == Surface_mesh::null_halfedge())
      return false;

    for(std::size_t i = 0; i < m_halfedge_properties.size(); ++i)
      m_halfedge_properties[i]->assign(element, hi);

    return true;
  }

  template <typename IntType>
  void process_line(PLY_element& element, Halfedge_index& hi)
  {
    IntType source, target;
    element.assign(source, "source");
    element.assign(target, "target");
    hi = m_mesh.halfedge(m_map_v2v[std::size_t(source)], m_map_v2v[std::size_t(target)]);
  }
};

template <typename Point>
bool fill_simplex_specific_header(std::ostream& os,
                                  const Surface_mesh<Point>& sm,
                                  std::vector<Abstract_property_printer<
                                    typename Surface_mesh<Point>::Vertex_index>*>& printers,
                                  const std::string& prop)
{
  typedef Surface_mesh<Point>                                     SMesh;
  typedef typename SMesh::Vertex_index                            VIndex;
  typedef typename Kernel_traits<Point>::Kernel                   Kernel;
  typedef typename Kernel::FT                                     FT;
  typedef typename Kernel::Vector_3                               Vector;
  typedef typename SMesh::template Property_map<VIndex, Point>    Point_map;
  typedef typename SMesh::template Property_map<VIndex, Vector>   Vector_map;
  typedef typename SMesh::template Property_map<VIndex, Color>    Vcolor_map;

  if(prop == "v:connectivity" || prop == "v:removed")
    return true;

  if(prop == "v:point")
  {
    if(boost::is_same<FT, float>::value)
    {
      os << "property float x" << std::endl
         << "property float y" << std::endl
         << "property float z" << std::endl;
    }
    else
    {
      os << "property double x" << std::endl
         << "property double y" << std::endl
         << "property double z" << std::endl;
    }
    printers.push_back(new Property_printer<VIndex, Point_map>(sm.points()));
    return true;
  }

  bool okay = false;
  if(prop == "v:normal")
  {
    Vector_map pmap;
    boost::tie(pmap, okay) = sm.template property_map<VIndex, Vector>(prop);
    if(okay)
    {
      if(boost::is_same<FT, float>::value)
      {
        os << "property float nx" << std::endl
           << "property float ny" << std::endl
           << "property float nz" << std::endl;
      }
      else
      {
        os << "property double nx" << std::endl
           << "property double ny" << std::endl
           << "property double nz" << std::endl;
      }
      printers.push_back(new Property_printer<VIndex, Vector_map>(pmap));
      return true;
    }
  }

  if(prop == "v:color")
  {
    Vcolor_map pmap;
    boost::tie(pmap, okay) = sm.template property_map<VIndex, Color>(prop);
    if(okay)
    {
      os << "property uchar red" << std::endl
         << "property uchar green" << std::endl
         << "property uchar blue" << std::endl
         << "property uchar alpha" << std::endl;

      printers.push_back(new Property_printer<VIndex, Vcolor_map>(pmap));
      return true;
    }
  }

  return false;
}

template <typename Point>
bool fill_simplex_specific_header(std::ostream& os,
                                  const Surface_mesh<Point>& sm,
                                  std::vector<Abstract_property_printer<
                                  typename Surface_mesh<Point>::Face_index>*>& printers,
                                  const std::string& prop)
{
  typedef Surface_mesh<Point>                                   SMesh;
  typedef typename SMesh::Face_index                            FIndex;
  typedef typename SMesh::template Property_map<FIndex, Color>  Fcolor_map;

  if(prop == "f:connectivity" || prop == "f:removed")
    return true;

  bool okay = false;
  if(prop == "f:color")
  {
    Fcolor_map pmap;
    boost::tie(pmap, okay) = sm.template property_map<FIndex, Color>(prop);
    if(okay)
    {
      os << "property uchar red" << std::endl
         << "property uchar green" << std::endl
         << "property uchar blue" << std::endl
         << "property uchar alpha" << std::endl;

      printers.push_back(new Property_printer<FIndex, Fcolor_map>(pmap));
      return true;
    }
  }
  return false;
}

template <typename Point>
bool fill_simplex_specific_header(std::ostream&,
                                  const Surface_mesh<Point>& ,
                                  std::vector<Abstract_property_printer<
                                  typename Surface_mesh<Point>::Edge_index>*>& ,
                                  const std::string& prop)
{
  if(prop == "e:removed")
    return true;

  return false;
}

template <typename Point>
bool fill_simplex_specific_header(std::ostream&,
                                  const Surface_mesh<Point>& ,
                                  std::vector<Abstract_property_printer<
                                  typename Surface_mesh<Point>::Halfedge_index>*>& ,
                                  const std::string& prop)
{
  if(prop == "h:connectivity")
    return true;

  return false;
}

template <typename Point>
std::string get_property_raw_name(const std::string& prop, typename Surface_mesh<Point>::Vertex_index)
{
  std::string name = prop;
  if(name.rfind("v:",0) == 0)
    name = std::string(prop.begin() + 2, prop.end());
  return name;
}

template <typename Point>
std::string get_property_raw_name(const std::string& prop, typename Surface_mesh<Point>::Face_index)
{
  std::string name = prop;
  if(name.rfind("f:",0) == 0)
    name = std::string(prop.begin() + 2, prop.end());
  return name;
}

template <typename Point>
std::string get_property_raw_name(const std::string& prop, typename Surface_mesh<Point>::Edge_index)
{
  std::string name = prop;
  if(name.rfind("e:",0) == 0)
    name = std::string(prop.begin() + 2, prop.end());
  return name;
}

template <typename Point>
std::string get_property_raw_name(const std::string& prop, typename Surface_mesh<Point>::Halfedge_index)
{
  std::string name = prop;
  if(name.rfind("h:", 0) == 0)
    name = std::string(prop.begin() + 2, prop.end());
  return name;
}

template <typename Point, typename Simplex>
void fill_header(std::ostream& os, const Surface_mesh<Point>& sm,
                 std::vector<Abstract_property_printer<Simplex>*>& printers)
{
  typedef Surface_mesh<Point>                                             SMesh;
  typedef typename SMesh::template Property_map<Simplex, boost::int8_t>   Int8_map;
  typedef typename SMesh::template Property_map<Simplex, boost::uint8_t>  Uint8_map;
  typedef typename SMesh::template Property_map<Simplex, boost::int16_t>  Int16_map;
  typedef typename SMesh::template Property_map<Simplex, boost::uint16_t> Uint16_map;
  typedef typename SMesh::template Property_map<Simplex, boost::int32_t>  Int32_map;
  typedef typename SMesh::template Property_map<Simplex, boost::uint32_t> Uint32_map;
  typedef typename SMesh::template Property_map<Simplex, boost::int64_t>  Int64_map;
  typedef typename SMesh::template Property_map<Simplex, boost::uint64_t> Uint64_map;
  typedef typename SMesh::template Property_map<Simplex, float>           Float_map;
  typedef typename SMesh::template Property_map<Simplex, double>          Double_map;

  std::vector<std::string> prop = sm.template properties<Simplex>();

  for(std::size_t i = 0; i < prop.size(); ++i)
  {
    if(fill_simplex_specific_header(os, sm, printers, prop[i]))
      continue;

    // Cut the "v:" prefix
    std::string name = get_property_raw_name<Point>(prop[i], Simplex());

    bool okay = false;
    {
      Int8_map pmap;
      boost::tie(pmap, okay) = sm.template property_map<Simplex,boost::int8_t>(prop[i]);
      if(okay)
      {
        os << "property char " << name << std::endl;
        printers.push_back(new internal::Char_property_printer<Simplex,Int8_map>(pmap));
        continue;
      }
    }
    {
      Uint8_map pmap;
      boost::tie(pmap, okay) = sm.template property_map<Simplex,boost::uint8_t>(prop[i]);
      if(okay)
      {
        os << "property uchar " << name << std::endl;
        printers.push_back(new internal::Char_property_printer<Simplex,Uint8_map>(pmap));
        continue;
      }
    }
    {
      Int16_map pmap;
      boost::tie(pmap, okay) = sm.template property_map<Simplex,boost::int16_t>(prop[i]);
      if(okay)
      {
        os << "property short " << name << std::endl;
        printers.push_back(new internal::Simple_property_printer<Simplex,Int16_map>(pmap));
        continue;
      }
    }
    {
      Uint16_map pmap;
      boost::tie(pmap, okay) = sm.template property_map<Simplex,boost::uint16_t>(prop[i]);
      if(okay)
      {
        os << "property ushort " << name << std::endl;
        printers.push_back(new internal::Simple_property_printer<Simplex,Uint16_map>(pmap));
        continue;
      }
    }
    {
      Int32_map pmap;
      boost::tie(pmap, okay) = sm.template property_map<Simplex,boost::int32_t>(prop[i]);
      if(okay)
      {
        os << "property int " << name << std::endl;
        printers.push_back(new internal::Simple_property_printer<Simplex,Int32_map>(pmap));
        continue;
      }
    }
    {
      Uint32_map pmap;
      boost::tie(pmap, okay) = sm.template property_map<Simplex,boost::uint32_t>(prop[i]);
      if(okay)
      {
        os << "property uint " << name << std::endl;
        printers.push_back(new internal::Simple_property_printer<Simplex,Uint32_map>(pmap));
        continue;
      }
    }
    {
      Int64_map pmap;
      boost::tie(pmap, okay) = sm.template property_map<Simplex,boost::int64_t>(prop[i]);
      if(okay)
      {
        os << "property int " << name << std::endl;
        printers.push_back(new internal::Simple_property_printer<Simplex,Int64_map,boost::int32_t>(pmap));
        continue;
      }
    }
    {
      Uint64_map pmap;
      boost::tie(pmap, okay) = sm.template property_map<Simplex,boost::uint64_t>(prop[i]);
      if(okay)
      {
        os << "property uint " << name << std::endl;
        printers.push_back(new internal::Simple_property_printer<Simplex,Uint64_map,boost::uint32_t>(pmap));
        continue;
      }
    }
    {
      Float_map pmap;
      boost::tie(pmap, okay) = sm.template property_map<Simplex,float>(prop[i]);
      if(okay)
      {
        os << "property float " << name << std::endl;
        printers.push_back(new internal::Simple_property_printer<Simplex,Float_map>(pmap));
        continue;
      }
    }
    {
      Double_map pmap;
      boost::tie(pmap, okay) = sm.template property_map<Simplex,double>(prop[i]);
      if(okay)
      {
        os << "property double " << name << std::endl;
        printers.push_back(new internal::Simple_property_printer<Simplex,Double_map>(pmap));
        continue;
      }
    }
  }
}

} // namespace internal

/// \ingroup PkgSurfaceMeshIOFuncPLY
///
/// \attention When reading a binary file, the flag `std::ios::binary` flag must be set during the creation of the `ifstream`.
///
/// \brief extracts the surface mesh from an input stream in the \ref IOStreamPLY
///        and appends it to the surface mesh `sm`.
///
/// - the operator reads the vertex `point` property and the face
///   `vertex_index` (or `vertex_indices`) property;
/// - if three PLY properties `nx`, `ny` and `nz` with type `float`
///   or `double` are found for vertices, a "v:normal" vertex
///   property map is added;
/// - if three PLY properties `red`, `green` and `blue` with type
///   `uchar` are found for vertices, a "v:color" vertex property
///   map is added;
/// - if three PLY properties `red`, `green` and `blue` with type
///   `uchar` are found for faces, a "f:color" face property map is
///   added;
/// - if any other PLY property is found, a "[s]:[name]" property map is
///   added, where `[s]` is `v` for vertex and `f` for face, and
///   `[name]` is the name of the PLY property.
///
/// \tparam Point The type of the \em point property of a vertex. There is no requirement on `P`,
///               besides being default constructible and assignable.
///               In typical use cases it will be a 2D or 3D point type.
///
/// \param is the input stream
/// \param sm the surface mesh to be constructed
/// \param comments a string used to store the potential comments found in the PLY header.
///        Each line starting by "comment " in the header is appended to the `comments` string
///        (without the "comment " word).
/// \param verbose whether extra information is printed when an incident occurs during reading
///
/// \pre The data in the stream must represent a two-manifold. If this is not the case
///      the `failbit` of `is` is set and the mesh cleared.
///
/// \returns `true` if reading was successful, `false` otherwise.
///
template <typename P>
bool read_PLY(std::istream& is,
              Surface_mesh<P>& sm,
              std::string& comments,
              bool verbose = true)
{
  typedef typename Surface_mesh<P>::size_type size_type;

  if(!is.good())
  {
    if(verbose)
      std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  internal::PLY_reader reader(verbose);
  internal::Surface_mesh_filler<P> filler(sm);

  if(!(reader.init(is)))
  {
    is.setstate(std::ios::failbit);
    return false;
  }

  comments = reader.comments();

  for(std::size_t i = 0; i < reader.number_of_elements(); ++ i)
  {
    internal::PLY_element& element = reader.element(i);

    bool is_vertex =(element.name() == "vertex" || element.name() == "vertices");
    bool is_face = false;
    bool is_edge = false;
    bool is_halfedge = false;
    if(is_vertex)
    {
      sm.reserve(sm.number_of_vertices() + size_type(element.number_of_items()),
                 sm.number_of_edges(),
                 sm.number_of_faces());
      filler.instantiate_vertex_properties(element);
    }
    else
      is_face =(element.name() == "face" || element.name() == "faces");

    if(is_face)
    {
      sm.reserve(sm.number_of_vertices(),
                 sm.number_of_edges(),
                 sm.number_of_faces() + size_type(element.number_of_items()));
      filler.instantiate_face_properties(element);
    }
    else
      is_edge =(element.name() == "edge");

    if(is_edge)
      filler.instantiate_edge_properties(element);
    else
      is_halfedge =(element.name() == "halfedge");

    if(is_halfedge)
      filler.instantiate_halfedge_properties(element);

    for(std::size_t j = 0; j < element.number_of_items(); ++ j)
    {
      for(std::size_t k = 0; k < element.number_of_properties(); ++ k)
      {
        internal::PLY_read_number* property = element.property(k);
        property->get(is);
        if(is.fail())
          return false;
      }

      if(is_vertex)
        filler.process_vertex_line(element);
      else if(is_face)
      {
        if(!filler.process_face_line(element))
        {
          is.setstate(std::ios::failbit);
          return false;
        }
      }
      else if(is_edge)
        filler.process_edge_line(element);
      else if(is_halfedge)
        filler.process_halfedge_line(element);
    }
  }

  return true;
}

/// \cond SKIP_IN_MANUAL

template <typename P>
bool read_PLY(std::istream& is, Surface_mesh<P>& sm)
{
  std::string dummy;
  return read_PLY(is, sm, dummy);
}

/// \endcond

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

/*!
  \ingroup PkgSurfaceMeshIOFuncDeprecated
  \deprecated This function is deprecated since \cgal 5.3, `CGAL::IO::read_PLY(std::ostream&, const Surface_mesh<Point>&)` should be used instead.
*/
template <typename P>
CGAL_DEPRECATED bool read_ply(std::istream& is, Surface_mesh<P>& sm, std::string& comments)
{
  return IO::read_PLY(is, sm, comments);
}

#endif // CGAL_NO_DEPRECATED_CODE

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

namespace IO {

/// \ingroup PkgSurfaceMeshIOFuncPLY
///
/// \brief inserts the surface mesh in an output stream in the \ref IOStreamPLY.
///
/// If found, internal property maps with names "v:normal", "v:color" and "f:color" are inserted in the stream.
///
/// All other vertex and face properties with simple types are inserted in the stream.
/// Edges are only inserted in the stream if they have at least one
/// property with simple type: if they do, all edge properties with
/// simple types are inserted in the stream. The halfedges follow
/// the same behavior.
///
/// \attention When writing a binary file, the flag `std::ios::binary` flag must be set during the creation of the `ofstream`.
///
/// \tparam Point The type of the \em point property of a vertex. There is no requirement on `P`,
///               besides being default constructible and assignable.
///               In typical use cases it will be a 2D or 3D point type.
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param os the output stream
/// \param sm the surface mesh to be output
/// \param comments a string included line by line in the header of the PLY stream (each line will be precedeed by "comment ")
/// \param np optional \ref bgl_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{stream_precision}
///     \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
///     \cgalParamType{int}
///     \cgalParamDefault{`the precision of the stream `os``}
///     \cgalParamExtra{This parameter is only meaningful while using ASCII encoding.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \returns `true` if writing was successful, `false` otherwise.
template <typename P,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_PLY(std::ostream& os,
               const Surface_mesh<P>& sm,
               const std::string& comments,
               const CGAL_BGL_NP_CLASS& np)
{
  typedef Surface_mesh<P> SMesh;
  typedef typename SMesh::Vertex_index VIndex;
  typedef typename SMesh::Face_index FIndex;
  typedef typename SMesh::Edge_index EIndex;
  typedef typename SMesh::Halfedge_index HIndex;

  if(!os.good())
    return false;

  set_stream_precision_from_NP(os, np);

  os << "ply" << std::endl
     << ((get_mode(os) == BINARY) ? "format binary_little_endian 1.0" : "format ascii 1.0") << std::endl
     << "comment Generated by the CGAL library" << std::endl;

  if(comments != std::string())
  {
    std::istringstream iss(comments);
    std::string line;
    while(getline(iss, line))
    {
      if(line != "Generated by the CGAL library") // Avoid repeating the line if multiple savings
        os << "comment " << line << std::endl;
    }
  }

  os << "element vertex " << sm.number_of_vertices() << std::endl;

  std::vector<internal::Abstract_property_printer<VIndex>*> vprinters;
  internal::fill_header(os, sm, vprinters);

  os << "element face " << sm.number_of_faces() << std::endl;
  os << "property list uchar int vertex_indices" << std::endl;
  std::vector<internal::Abstract_property_printer<FIndex>*> fprinters;
  internal::fill_header(os, sm, fprinters);


  std::vector<internal::Abstract_property_printer<EIndex>*> eprinters;
  if(sm.template properties<EIndex>().size() > 1)
  {
    std::ostringstream oss;
    internal::fill_header(oss, sm, eprinters);

    if(!eprinters.empty())
    {
      os << "element edge " << sm.number_of_edges() << std::endl;
      os << "property int v0" << std::endl;
      os << "property int v1" << std::endl;
      os << oss.str();
    }
  }

  std::vector<internal::Abstract_property_printer<HIndex>*> hprinters;
  if(sm.template properties<HIndex>().size() > 1)
  {
    std::ostringstream oss;
    internal::fill_header(oss, sm, hprinters);

    if(!hprinters.empty())
    {
      os << "element halfedge " << sm.number_of_halfedges() << std::endl;
      os << "property int source" << std::endl;
      os << "property int target" << std::endl;
      os << oss.str();
    }
  }

  os << "end_header" << std::endl;

  std::vector<int> reindex;
  reindex.resize (sm.num_vertices());
  int n = 0;
  for(VIndex vi : sm.vertices())
  {
    for(std::size_t i = 0; i < vprinters.size(); ++ i)
    {
      vprinters[i]->print(os, vi);
      if(get_mode(os) == ASCII)
        os << " ";
    }
    if(get_mode(os) == ASCII)
      os << std::endl;

    reindex[std::size_t(vi)] = n++;
  }

  std::vector<int> polygon;

  for(FIndex fi : sm.faces())
  {
    // Get list of vertex indices
    polygon.clear();
    for(VIndex vi : CGAL::vertices_around_face(sm.halfedge(fi), sm))
      polygon.push_back(reindex[std::size_t(vi)]);

    if(get_mode(os) == ASCII)
    {
      os << polygon.size() << " ";
      for(std::size_t i = 0; i < polygon.size(); ++ i)
        os << polygon[i] << " ";
    }
    else
    {
      unsigned char size =(unsigned char)(polygon.size());
      os.write(reinterpret_cast<char*>(&size), sizeof(size));
      for(std::size_t i = 0; i < polygon.size(); ++ i)
      {
        int idx = polygon[i];
        os.write(reinterpret_cast<char*>(&idx), sizeof(idx));
      }
    }

    for(std::size_t i = 0; i < fprinters.size(); ++ i)
    {
      fprinters[i]->print(os, fi);
      if(get_mode(os) == ASCII)
        os << " ";
    }

    if(get_mode(os) == ASCII)
      os << std::endl;
  }

  if(!eprinters.empty())
  {
    for(EIndex ei : sm.edges())
    {
      int v0 = reindex[std::size_t(sm.vertex(ei, 0))];
      int v1 = reindex[std::size_t(sm.vertex(ei, 1))];
      if(get_mode(os) == ASCII)
      {
        os << v0 << " " << v1 << " ";
      }
      else
      {
        os.write(reinterpret_cast<char*>(&v0), sizeof(v0));
        os.write(reinterpret_cast<char*>(&v1), sizeof(v1));
      }

      for(std::size_t i = 0; i < eprinters.size(); ++ i)
      {
        eprinters[i]->print(os, ei);
        if(get_mode(os) == ASCII)
          os << " ";
      }

      if(get_mode(os) == ASCII)
        os << std::endl;
    }
  }

  if(!hprinters.empty())
  {
    for(HIndex hi : sm.halfedges())
    {
      int source = reindex[std::size_t(sm.source(hi))];
      int target = reindex[std::size_t(sm.target(hi))];
      if(get_mode(os) == ASCII)
      {
        os << source << " " << target << " ";
      }
      else
      {
        os.write(reinterpret_cast<char*>(&source), sizeof(source));
        os.write(reinterpret_cast<char*>(&target), sizeof(target));
      }

      for(std::size_t i = 0; i < hprinters.size(); ++ i)
      {
        hprinters[i]->print(os, hi);
        if(get_mode(os) == ASCII)
          os << " ";
      }

      if(get_mode(os) == ASCII)
        os << std::endl;
    }
  }

  for(std::size_t i = 0; i < vprinters.size(); ++ i)
    delete vprinters[i];
  for(std::size_t i = 0; i < fprinters.size(); ++ i)
    delete fprinters[i];
  for(std::size_t i = 0; i < eprinters.size(); ++ i)
    delete eprinters[i];
  for(std::size_t i = 0; i < hprinters.size(); ++ i)
    delete hprinters[i];

  return true;
}

/// \cond SKIP_IN_MANUAL

template <typename P>
bool write_PLY(std::ostream& os, const Surface_mesh<P>& sm, const std::string& comments)
{
  return write_PLY(os, sm, comments, parameters::all_default());
}

template <typename P, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_PLY(std::ostream& os, const Surface_mesh<P>& sm, const CGAL_BGL_NP_CLASS& np)
{
  std::string unused_comment;
  return write_PLY(os, sm, unused_comment, np);
}

template <typename P>
bool write_PLY(std::ostream& os, const Surface_mesh<P>& sm)
{
  std::string unused_comment;
  return write_PLY(os, sm, unused_comment, parameters::all_default());
}

/// \endcond

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

/*!
  \ingroup PkgSurfaceMeshIOFuncDeprecated
  \deprecated This function is deprecated since \cgal 5.3, `CGAL::IO::write_PLY(std::ostream&, const Surface_mesh<Point>&)` should be used instead.
*/

template <typename P>
CGAL_DEPRECATED bool write_ply(std::ostream& os, const Surface_mesh<P>& sm, const std::string& comments)
{
  return IO::write_PLY(os, sm, comments);
}

template <typename P>
CGAL_DEPRECATED bool write_ply(std::ostream& os, const Surface_mesh<P>& sm)
{
  return write_PLY(os, sm, "");
}
#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // !defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE) && !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES)

#endif // CGAL_SURFACE_MESH_IO_PLY_H
