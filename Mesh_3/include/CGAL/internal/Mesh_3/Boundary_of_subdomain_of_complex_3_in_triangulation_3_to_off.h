// Copyright (c) 2008-2014  GeometryFactory (France).
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
// Author(s)     : Laurent Rineau and Sebastien Loriot

#ifndef CGAL_INTERNAL_MESH_3_BOUNDARY_OF_SUDDOMAIN_OF_COMPLEX_3_IN_TRIANGULATION_3_TO_OFF_H
#define CGAL_INTERNAL_MESH_3_BOUNDARY_OF_SUDDOMAIN_OF_COMPLEX_3_IN_TRIANGULATION_3_TO_OFF_H

#include <CGAL/license/Mesh_3.h>


#include <map>
#include <sstream>

#include <CGAL/array.h>

namespace CGAL {

namespace internal{

namespace mesh_3_export{
template <class Vertex_handle>
std::size_t get_vertex_index(Vertex_handle v,std::map<Vertex_handle, std::size_t>& V,std::size_t& inum,std::stringstream& vertex_buffer){
  std::pair<typename std::map<Vertex_handle, std::size_t>::iterator,bool> res=
    V.insert(std::make_pair(v,inum));
  if (res.second){
    ++inum;
    vertex_buffer <<   res.first->first->point().point() <<"\n"; //point is weighted!
  }
  return res.first->second;
}
} // end of namespace mesh_3_export

template <typename C3T3>
std::ostream&
output_boundary_of_c3t3_to_off(const C3T3& c3t3, 
                               typename C3T3::Subdomain_index sd_index,
                               std::ostream& output,
                               bool normals_point_outside_of_the_subdomain=true)
{
  typedef typename C3T3::Triangulation Triangulation;
  typedef typename Triangulation::Vertex_handle Vertex_handle;

  std::map<Vertex_handle, std::size_t> V;
  
  std::size_t inum = 0; 
  std::size_t nfacets = 0;
  cpp0x::array<std::size_t,3> indices={{0,0,0}};
  std::stringstream facet_buffer,vertex_buffer;
  for(typename C3T3::Facets_in_complex_iterator 
        fit = c3t3.facets_in_complex_begin(),
        end = c3t3.facets_in_complex_end();
      fit != end; ++fit) 
  {
    typename C3T3::Subdomain_index cell_sd=c3t3.subdomain_index(fit->first);
    typename C3T3::Subdomain_index opp_sd=c3t3.subdomain_index(fit->first->neighbor(fit->second));
    
    if (cell_sd!=sd_index && opp_sd!=sd_index) continue;

    ++nfacets;
    int j=-1;
    
    
    for (int i = 0; i < 4; ++i)
      if (i != fit->second)
          indices[++j]=mesh_3_export::get_vertex_index((*fit).first->vertex(i), V, inum,vertex_buffer);
    if ( ( (cell_sd==sd_index) == (fit->second%2 == 1) ) == normals_point_outside_of_the_subdomain )
      std::swap(indices[0],indices[1]);
    facet_buffer << "3" << " " << indices[0] <<" " << indices[1] <<" " << indices[2] << "\n";
  }
  
  output << "OFF " << inum << " " << nfacets << " 0\n";
  output << vertex_buffer.str();
  output << facet_buffer.str();
  
  
  return output;
}

template <typename C3T3>
std::ostream&
output_facets_in_complex_to_off(const C3T3& c3t3,
                                std::ostream& output)
{
  typedef typename C3T3::Triangulation Triangulation;
  typedef typename Triangulation::Vertex_handle Vertex_handle;

  std::map<Vertex_handle, std::size_t> V;

  std::size_t inum = 0;
  std::size_t nfacets = 0;
  cpp0x::array<std::size_t,3> indices={{0,0,0}};
  std::stringstream facet_buffer,vertex_buffer;
  for(typename C3T3::Facets_in_complex_iterator
        fit = c3t3.facets_in_complex_begin(),
        end = c3t3.facets_in_complex_end();
      fit != end; ++fit)
  {
    typename C3T3::Subdomain_index cell_sd=c3t3.subdomain_index(fit->first);
    typename C3T3::Subdomain_index opp_sd=c3t3.subdomain_index(fit->first->neighbor(fit->second));

    ++nfacets;
    int j=-1;


    for (int i = 0; i < 4; ++i)
      if (i != fit->second)
          indices[++j]=mesh_3_export::get_vertex_index((*fit).first->vertex(i), V, inum,vertex_buffer);
    if ( (cell_sd > opp_sd) == (fit->second%2 == 1) ) std::swap(indices[0],indices[1]);
    facet_buffer << "3" << " " << indices[0] <<" " << indices[1] <<" " << indices[2] << "\n";
  }

  output << "OFF " << inum << " " << nfacets << " 0\n";
  output << vertex_buffer.str();
  output << facet_buffer.str();


  return output;
}

} } // end of namespace CGAL::internal


#endif // CGAL_INTERNAL_MESH_3_BOUNDARY_OF_SUDDOMAIN_OF_COMPLEX_3_IN_TRIANGULATION_3_TO_OFF_H
