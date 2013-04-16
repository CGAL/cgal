// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************


#ifndef CGAL_MESH_VERTEX_BASE_3_H
#define CGAL_MESH_VERTEX_BASE_3_H

#include <CGAL/Mesh_3/config.h>

#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Surface_mesh_vertex_base_3.h>
#include <CGAL/Mesh_3/Has_features.h>
#include <CGAL/internal/Mesh_3/get_index.h>
#include <CGAL/Mesh_3/io_signature.h>

namespace CGAL {
  
// Without erase counter
template <typename Use_erase_counter>
class Mesh_vertex_base_3_base
{
};

#ifdef CGAL_LINKED_WITH_TBB
// Specialized version (with erase counter)
template <>
class Mesh_vertex_base_3_base<Tag_true>
{
public:
  
  // Erase counter (cf. Compact_container)
  unsigned int get_erase_counter() const
  {
    return this->m_erase_counter;
  }
  void set_erase_counter(unsigned int c)
  {
    this->m_erase_counter = c;
  }
  void increment_erase_counter()
  {
    ++this->m_erase_counter;
  }
  
protected:
  typedef tbb::atomic<unsigned int> Erase_counter_type;
  Erase_counter_type                m_erase_counter;

};
#endif // CGAL_LINKED_WITH_TBB

// Class Mesh_vertex_base_3
// Vertex base class used in 3D meshing process.
// Adds information to Vb about the localization of the vertex in regards
// to the 3D input complex.
template<class GT,
         class MT,
         class Vb = Triangulation_vertex_base_3<GT> >
class Mesh_vertex_base_3
: public Surface_mesh_vertex_base_3<GT, Vb>,
  public Mesh_vertex_base_3_base<
    typename Surface_mesh_vertex_base_3<GT, Vb>::Tds::Cell_container_strategy::Uses_erase_counter>
{
public:
  typedef Surface_mesh_vertex_base_3<GT, Vb> Mvb3_base;
  typedef typename Mvb3_base::Vertex_handle Vertex_handle;

  // To get correct vertex type in TDS
  template < class TDS3 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS3>::Other Vb3;
    typedef Mesh_vertex_base_3 <GT, MT, Vb3> Other;
  };

  // Types
  typedef typename MT::Index                      Index;
  typedef typename GT::FT                         FT;

  // Constructor
  Mesh_vertex_base_3() 
    : Surface_mesh_vertex_base_3<GT, Vb>()
    , index_()
    , dimension_(-1)
    , meshing_info_(0)
#ifdef CGAL_INTRUSIVE_LIST
    , next_intrusive_()
    , previous_intrusive_()
#endif //CGAL_INTRUSIVE_LIST
  {}

  // Default copy constructor and assignment operator are ok

  // Returns the dimension of the lowest dimensional face of the input 3D
  // complex that contains the vertex
  int in_dimension() const {
    if(dimension_ < -1) return -2-dimension_;
    else return dimension_; 
  }

  // Sets the dimension of the lowest dimensional face of the input 3D complex
  // that contains the vertex
  void set_dimension(const int dimension) { dimension_ = dimension; }

  // Tells if the vertex is marked as a special protecting ball
  bool is_special() const { return dimension_ < -1; }

  // Marks or unmarks the vertex as a special protecting ball
  void set_special(bool special = true) {
    if(special != (dimension_ < -1) )
      dimension_ = -2-dimension_;
  }

  // Returns the index of the lowest dimensional face of the input 3D complex
  // that contains the vertex
  Index index() const { return index_; }

  // Sets the index of the lowest dimensional face of the input 3D complex
  // that contains the vertex
  void set_index(const Index& index) { index_ = index; }

  // Accessors to meshing_info private data
  const FT& meshing_info() const { return meshing_info_; }
  void set_meshing_info(const FT& value) { meshing_info_ = value; }

  static
  std::string io_signature()
  {
    return
      Get_io_signature<Vb>()() + "+" +
      Get_io_signature<int>()() + "+" +
      Get_io_signature<Index>()();
  }
  
#ifdef CGAL_INTRUSIVE_LIST
public:
  Vertex_handle next_intrusive() const { return next_intrusive_; }
  Vertex_handle& next_intrusive() { return next_intrusive_; }

  Vertex_handle previous_intrusive() const { return previous_intrusive_; }
  Vertex_handle& previous_intrusive() { return previous_intrusive_; }
#endif

private:
  // Index of the lowest dimensional face of the input 3D complex
  // that contains me
  Index index_;
  // Dimension of the lowest dimensional face of the input 3D complex
  // that contains me. Negative values are a marker for special vertices.
  int dimension_;
  // Stores info needed by optimizers
  FT meshing_info_;

#ifdef CGAL_INTRUSIVE_LIST
  Vertex_handle next_intrusive_;
  Vertex_handle previous_intrusive_;
#endif
};  // end class Mesh_vertex_base_3

template<class GT,
         class MT,
         class Vb>
inline
std::istream&
operator>>(std::istream &is, Mesh_vertex_base_3<GT,MT,Vb>& v)
{
  typedef Mesh_vertex_base_3<GT,MT,Vb> Vertex;
  typedef typename Vertex::Mvb3_base Mvb3_base;
  is >> static_cast<Mvb3_base&>(v);
  int dimension;
  if(is_ascii(is)) {
    is >> dimension;

  } else {
    CGAL::read(is, dimension);
  }
  CGAL_assertion(dimension >= 0);
  CGAL_assertion(dimension < 4);
  typename Vertex::Index index = 
    internal::Mesh_3::Read_mesh_domain_index<MT>()(dimension, is);
  v.set_dimension(dimension);
  v.set_index(index);
  return is;
}

template<class GT,
         class MT,
         class Vb>
inline
std::ostream&
operator<<(std::ostream &os, const Mesh_vertex_base_3<GT,MT,Vb>& v)
{
  typedef Mesh_vertex_base_3<GT,MT,Vb> Vertex;
  typedef typename Vertex::Mvb3_base Mvb3_base;
  os << static_cast<const Mvb3_base&>(v);
  if(is_ascii(os)) {
    os << " " << v.in_dimension()
       << " ";
  } else {
    CGAL::write(os, v.in_dimension());
  }
  internal::Mesh_3::Write_mesh_domain_index<MT>()(os, 
                                                  v.in_dimension(),
                                                  v.index());
  return os;
}

}  // end namespace CGAL



#endif // CGAL_MESH_VERTEX_BASE_3_H
