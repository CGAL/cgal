// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Clement Jamin


#ifndef CGAL_TRIANGULATION_IO_H
#define CGAL_TRIANGULATION_IO_H

#include <CGAL/Epick_d.h>
#include <CGAL/Triangulation.h>
#include <sstream>
#include <iostream>

namespace CGAL {

namespace Triangulation_IO
{
// TODO: test if the stream is binary or text?
template<typename Traits, typename P>
int
output_point(std::ostream & os, const Traits &traits, const P & p)
{
  typedef typename Traits::Compute_coordinate_d Ccd;
  const Ccd ccd = traits.compute_coordinate_d_object();
  const int dim = traits.point_dimension_d_object()(p);
  if (dim > 0)
  {
    os << ccd(p, 0);
    for (int i = 1 ; i < dim ; ++i)
      os << " " << CGAL::to_double(ccd(p, i));
  }
  return dim;
}

// TODO: test if the stream is binary or text?
template<typename Traits, typename P>
int
output_weighted_point(std::ostream & os, const Traits &traits, const P & p,
                      bool output_weight = true)
{
  typedef typename Traits::Compute_coordinate_d Ccd;
  typename Traits::Construct_point_d cp =
    traits.construct_point_d_object();
  typename Traits::Compute_weight_d pt_weight = traits.compute_weight_d_object();
  const Ccd ccd = traits.compute_coordinate_d_object();
  const int dim = traits.point_dimension_d_object()(p);
  if (dim > 0)
  {
    output_point(os, traits, p);
    if (output_weight)
      os << " " << pt_weight(p);
  }
  return dim;
}

// TODO: test if the stream is binary or text?
template<typename Traits, typename FCH>
void
output_full_cell(std::ostream & os, const Traits &traits, const FCH & fch,
                      bool output_weights = false)
{
  typename FCH::value_type::Vertex_handle_iterator vit = fch->vertices_begin();
  for( ; vit != fch->vertices_end(); ++vit )
  {
    int dim;
    if (output_weights)
      dim = output_weighted_point(os, traits, (*vit)->point());
    else
      dim = output_point(os, traits, (*vit)->point());
    if (dim > 0)
      os << std::endl;
  }
}

// TODO: test if the stream is binary or text?
/*template<typename Traits, typename P>
void
input_point(std::istream & is, const Traits &traits, P & p)
{
  typedef typename Traits::FT FT;
  std::vector<FT> coords;

  std::string line;
  for(;;)
  {
    if (!std::getline(is, line))
      return is;
    if (line != "")
      break;
  }
  std::stringstream line_sstr(line);
  FT temp;
  while (line_sstr >> temp)
    coords.push_back(temp);

  p = traits.construct_point_d_object()(coords.begin(), coords.end());
}*/

} // namespace Triangulation_IO

///////////////////////////////////////////////////////////////
// TODO: replace these operator>> by an "input_point" function
///////////////////////////////////////////////////////////////

template < class GT, class TDS >
std::ostream &
export_triangulation_to_off(std::ostream & os,
                            const Triangulation<GT,TDS> & tr,
                            bool in_3D_export_surface_only = false)
{
  typedef Triangulation<GT,TDS>                         Tr;
  typedef typename Tr::Vertex_const_handle              Vertex_handle;
  typedef typename Tr::Finite_vertex_const_iterator     Finite_vertex_iterator;
  typedef typename Tr::Finite_full_cell_const_iterator  Finite_full_cell_iterator;
  typedef typename Tr::Full_cell_const_iterator         Full_cell_iterator;
  typedef typename Tr::Full_cell                        Full_cell;
  typedef typename Full_cell::Vertex_handle_const_iterator Full_cell_vertex_iterator;

  if (tr.maximal_dimension() < 2 || tr.maximal_dimension() > 3)
  {
    std::cerr << "Warning: export_tds_to_off => dimension should be 2 or 3.";
    os << "Warning: export_tds_to_off => dimension should be 2 or 3.";
    return os;
  }

  std::size_t n = tr.number_of_vertices();

  std::stringstream output;

  // write the vertices
  std::map<Vertex_handle, int> index_of_vertex;
  int i = 0;
  for(Finite_vertex_iterator it = tr.finite_vertices_begin();
      it != tr.finite_vertices_end(); ++it, ++i)
  {
    Triangulation_IO::output_point(output, tr.geom_traits(), it->point());
    if (tr.maximal_dimension() == 2)
      output << " 0";
    output << std::endl;
    index_of_vertex[it.base()] = i;
  }
  CGAL_assertion( static_cast<std::size_t>(i) == n );

  std::size_t number_of_triangles = 0;
  if (tr.maximal_dimension() == 2)
  {
    for (Finite_full_cell_iterator fch = tr.finite_full_cells_begin() ;
         fch != tr.finite_full_cells_end() ; ++fch)
    {
      output << "3 ";
      for (Full_cell_vertex_iterator vit = fch->vertices_begin() ;
           vit != fch->vertices_end() ; ++vit)
      {
        output << index_of_vertex[*vit] << " ";
      }
      output << std::endl;
      ++number_of_triangles;
    }
  }
  else if (tr.maximal_dimension() == 3)
  {
    if (in_3D_export_surface_only)
    {
      // Parse boundary facets
      for (Full_cell_iterator fch = tr.full_cells_begin() ;
           fch != tr.full_cells_end() ; ++fch)
      {
        if (tr.is_infinite(fch))
        {
          output << "3 ";
          for (Full_cell_vertex_iterator vit = fch->vertices_begin() ;
               vit != fch->vertices_end() ; ++vit)
          {
            if (!tr.is_infinite(*vit))
              output << index_of_vertex[*vit] << " ";
          }
          output << std::endl;
          ++number_of_triangles;
        }
      }
    }
    else
    {
      // Parse finite cells
      for (Finite_full_cell_iterator fch = tr.finite_full_cells_begin() ;
           fch != tr.finite_full_cells_end() ; ++fch)
      {
        output << "3 "
               << index_of_vertex[fch->vertex(0)] << " "
               << index_of_vertex[fch->vertex(1)] << " "
               << index_of_vertex[fch->vertex(2)]
               << std::endl;
        output << "3 "
               << index_of_vertex[fch->vertex(0)] << " "
               << index_of_vertex[fch->vertex(2)] << " "
               << index_of_vertex[fch->vertex(3)]
               << std::endl;
        output << "3 "
               << index_of_vertex[fch->vertex(1)] << " "
               << index_of_vertex[fch->vertex(2)] << " "
               << index_of_vertex[fch->vertex(3)]
               << std::endl;
        output << "3 "
               << index_of_vertex[fch->vertex(0)] << " "
               << index_of_vertex[fch->vertex(1)] << " "
               << index_of_vertex[fch->vertex(3)]
               << std::endl;
        number_of_triangles += 4;
      }
    }
  }

  os << "OFF \n"
     << n << " "
     << number_of_triangles << " 0\n"
     << output.str();

  return os;
}

} //namespace CGAL

#endif // CGAL_TRIANGULATION_IO_H
