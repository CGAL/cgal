// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
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
  Smooth_old_vertex (LCC & alcc, unsigned int /* TODO amark*/):mlcc (alcc)
  {
  }

  Vertex operator  () (Vertex & v) const
  {
    Dart_handle d = v.dart ();
    CGAL_assertion (d != NULL);

    int degree = 0;
    bool open = false;

    LCC::One_dart_per_incident_cell_range<1,0>::iterator it (mlcc, d),
      itend(mlcc.one_dart_per_incident_cell<1,0>(d).end());
    for (; it != itend; ++it)
    {
      ++degree;
      if (it->is_free (2)) open = true;
    }

    if (open)
      return v;

    LCC::FT alpha = (4.0f - 2.0f *
                     (LCC::FT) cos (2.0f * PI / (LCC::FT) degree)) / 9.0f;
    LCC::Vector vec =
      LCC::Traits::Construct_scaled_vector()
      ( LCC::Traits::Construct_vector() (CGAL::ORIGIN, v.point()), (1.0f - alpha));

    for (it.rewind (); it != itend; ++it)
    {
      CGAL_assertion (!it->is_free (2));
      vec = vec + (mlcc.point(it->beta(2)) - CGAL::ORIGIN)
        * alpha / degree;
    }

    Vertex res= LCC::Traits::Construct_translated_point() (CGAL::ORIGIN, vec);
    res.set_dart (d);

    return res;
  }
private:
  LCC & mlcc;
};

// Flip an edge, work only in 2D and 3D
Dart_handle
flip_edge (LCC & m, Dart_handle d)
{
  CGAL_assertion ( d!=NULL && !d->is_free(2) );
  CGAL_assertion ( !d->is_free(1) && !d->is_free(0) );
  CGAL_assertion ( !d->beta(2)->is_free(0) && !d->beta(2)->is_free(1) );  
  
  if (!CGAL::is_removable<LCC,1>(m,d)) return NULL;

  Dart_handle d1 = d->beta(1);
  Dart_handle d2 = d->beta(2)->beta(0);

  CGAL_assertion ( !d1->is_free(1) && !d2->is_free(0) );

  Dart_handle d3 = d1->beta(1);
  Dart_handle d4 = d2->beta(0);

  // We isolated the edge
  m.basic_link_beta_1(d->beta(0), d->beta(2)->beta(1));
  m.basic_link_beta_0(d->beta(1), d->beta(2)->beta(0));
  if ( !d->is_free(3) )
  {
    m.basic_link_beta_0(d->beta(0)->beta(3), d->beta(2)->beta(1)->beta(3));
    m.basic_link_beta_1(d->beta(1)->beta(3), d->beta(2)->beta(0)->beta(3));
  }

  // Then we push the two extremities.
  m.basic_link_beta_0(d3, d);
  m.basic_link_beta_0(d2, d->beta(2));
  m.link_beta_1(d4, d);
  m.link_beta_1(d1, d->beta(2));

  if ( !d->is_free(3) )
  {
    m.basic_link_beta_0(d4->beta(3), d->beta(3));
    m.basic_link_beta_0(d1->beta(3), d->beta(2)->beta(3));
    m.link_beta_1(d3->beta(3), d->beta(3));
    m.link_beta_1(d2->beta(3), d->beta(2)->beta(3));
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

  unsigned int mark = m.get_new_mark ();
  unsigned int treated = m.get_new_mark ();
  m.negate_mark (mark);  // All the old darts are marked in O(1).

  // 1) We smoth the old vertices.
  std::vector < Vertex > vertices;  // smooth the old vertices
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
  for (std::vector < Vertex >::iterator vit = vertices.begin ();
       vit != vertices.end (); ++vit)
  {
    LCC::point(vit->dart())=vit->point();
  }

  // 4) We flip all the old edges.
  m.negate_mark (mark);  // Now only new darts are marked.
  Dart_handle d2 = NULL;
  for (LCC::Dart_range::iterator it (m.darts().begin ());
       it != m.darts().end ();)
  {
    d2 = it++;
    CGAL_assertion (d2 != NULL);
    if (!m.is_marked (d2, mark))  // This is an old dart.
    {
      // We process only the last dart of a same edge.
      if (!d2->is_free(2) && (d2->beta(2)->beta(3)==d2->beta(3)->beta(2)))
      {
        if (m.is_marked(d2->beta(2), mark) &&
            (d2->is_free(3) ||
             (m.is_marked(d2->beta(3), mark) &&
              m.is_marked(d2->beta(2)->beta(3), mark))))
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
