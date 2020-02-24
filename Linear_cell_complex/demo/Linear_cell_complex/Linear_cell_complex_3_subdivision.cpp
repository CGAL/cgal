// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#include "typedefs.h"

#define PI 3.1415926535897932

// Smoth a vertex depending on the vertices of its incident facet.
class Smooth_old_vertex
{
public:
  /**  Constructor.
   * @param alcc is the lcc to smooth
   * @param amark is a mark designing old darts (i.e. darts not created during
   *        the triangulation step)
   */
  Smooth_old_vertex (LCC & alcc, LCC::size_type /* TODO amark*/):mlcc (alcc)
  {
  }

  std::pair<Point_3, Dart_handle> operator  () (Vertex & v) const
  {
    Dart_handle d = v.dart ();

    int degree = 0;
    bool open = false;

    LCC::One_dart_per_incident_cell_range<1,0>::iterator it (mlcc, d),
      itend(mlcc.one_dart_per_incident_cell<1,0>(d).end());
    for (; it != itend; ++it)
    {
      ++degree;
      if (mlcc.is_free(it, 2)) open = true;
    }

    if (open)
      return make_pair(v.point(), d);

    LCC::FT alpha = (4.0f - 2.0f *
                     (LCC::FT) cos (2.0f * PI / (LCC::FT) degree)) / 9.0f;
    LCC::Vector vec =
      LCC::Traits::Construct_scaled_vector()
      ( LCC::Traits::Construct_vector() (CGAL::ORIGIN, v.point()), (1.0f - alpha));

    for (it.rewind (); it != itend; ++it)
    {
      CGAL_assertion (!mlcc.is_free(it,2));
      vec = vec + (mlcc.point(mlcc.beta(it,2)) - CGAL::ORIGIN)
        * alpha / degree;
    }

    std::pair<Point_3, Dart_handle> res=std::make_pair
      (LCC::Traits::Construct_translated_point() (CGAL::ORIGIN, vec), d);

    return res;
  }
private:
  LCC & mlcc;
};

// Flip an edge, work only in 2D and 3D
Dart_handle
flip_edge (LCC & m, Dart_handle d)
{
  CGAL_assertion ( !m.is_free(d,2) );
  CGAL_assertion ( !m.is_free(d,1) && !m.is_free(d,0) );
  CGAL_assertion ( !m.is_free(m.beta(d,2), 0) && !m.is_free(m.beta(d, 2), 1) );
  
  if (!m.is_removable<1>(d)) return LCC::null_handle;

  Dart_handle d1 = m.beta(d,1);
  Dart_handle d2 = m.beta(d,2,0);

  CGAL_assertion ( !m.is_free(d1,1) && !m.is_free(d2,0) );

  Dart_handle d3 = m.beta(d1,1);
  Dart_handle d4 = m.beta(d2, 0);

  // We isolated the edge
  m.basic_link_beta_1(m.beta(d,0), m.beta(d,2,1));
  m.basic_link_beta_0(m.beta(d,1), m.beta(d,2,0));
  if ( !m.is_free(d,3) )
  {
    m.basic_link_beta_0(m.beta(d,0,3), m.beta(d,2,1,3));
    m.basic_link_beta_1(m.beta(d,1,3), m.beta(d,2,0,3));
  }

  // Then we push the two extremities.
  m.basic_link_beta_0(d3, d);
  m.basic_link_beta_0(d2, m.beta(d,2));
  m.link_beta_1(d4, d);
  m.link_beta_1(d1, m.beta(d,2));

  if ( !m.is_free(d,3) )
  {
    m.basic_link_beta_0(m.beta(d4,3), m.beta(d,3));
    m.basic_link_beta_0(m.beta(d1,3), m.beta(d,2,3));
    m.link_beta_1(m.beta(d3,3), m.beta(d,3));
    m.link_beta_1(m.beta(d2,3), m.beta(d,2,3));
  }
  
  // CGAL::remove_cell<LCC,1>(m, d);
  // insert_cell_1_in_cell_2(m, d1, d1->beta(1)->beta(1));

  return d;
}

// Subdivide each facet of the lcc by using sqrt(3)-subdivision.
void
subdivide_lcc_3 (LCC & m)
{
  if (m.number_of_darts () == 0)
    return;

  LCC::size_type mark = m.get_new_mark ();
  LCC::size_type treated = m.get_new_mark ();
  m.negate_mark (mark);  // All the old darts are marked in O(1).

  // 1) We smoth the old vertices.
  std::vector <std::pair<Point_3, Dart_handle> > vertices;  // smooth the old vertices
  vertices.reserve (m.number_of_attributes<0> ());  // get intermediate space
  std::transform (m.vertex_attributes().begin (), 
		  m.vertex_attributes().end (),
		  std::back_inserter (vertices), 
		  Smooth_old_vertex (m, mark));

  // 2) We subdivide each facet.
  m.negate_mark (treated);  // All the darts are marked in O(1).
  unsigned int nb = 0;
  for (LCC::Dart_range::iterator it (m.darts().begin ());
       m.number_of_marked_darts (treated) > 0; ++it)
  {
    ++nb;
    if (m.is_marked (it, treated))
    {
      // We unmark the darts of the facet to process only once dart/facet.
      CGAL::unmark_cell < LCC, 2 > (m, it, treated);
      // We triangulate the facet.
      m.insert_barycenter_in_cell<2>(it);
    }
  }

  CGAL_assertion (m.is_whole_map_unmarked (treated));
  CGAL_assertion (m.is_valid ());
  m.free_mark (treated);

  // 3) We update the coordinates of old vertices.
  for (std::vector<std::pair<Point_3, Dart_handle> >::iterator
         vit=vertices.begin(); vit!=vertices.end(); ++vit)
  {
    m.point(vit->second)=vit->first;
  }

  // 4) We flip all the old edges.
  m.negate_mark (mark);  // Now only new darts are marked.
  Dart_handle d2 =LCC::null_handle;
  for (LCC::Dart_range::iterator it (m.darts().begin ());
       it != m.darts().end ();)
  {
    d2 = it++;
    if (!m.is_marked (d2, mark))  // This is an old dart.
    {
      // We process only the last dart of a same edge.
      if (!m.is_free(d2,2) && (m.beta(d2,2,3)==m.beta(d2,3,2)))
      {
        if (m.is_marked(m.beta(d2,2), mark) &&
            (m.is_free(d2,3) ||
             (m.is_marked(m.beta(d2,3), mark) &&
              m.is_marked(m.beta(d2,2,3), mark))))
        {
          flip_edge (m, d2);
          m.mark(d2, mark);
        }
        else
          m.mark (d2, mark);
      }
      else
        m.mark (d2, mark);
    }
  }

  CGAL_assertion (m.is_whole_map_marked (mark));
  m.free_mark (mark);
  
  CGAL_postcondition ( m.is_valid ());
}
