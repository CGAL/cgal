// Copyright (c) 2010 CNRS, LIRIS, http://liris.cnrs.fr/, All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
   * @param amap is the map to smooth
   * @param amark is a mark designing old darts (i.e. darts not created during
   *        the triangulation step)
   */
  Smooth_old_vertex (Map & amap, unsigned int amark):mmap (amap)
  {
  }

  Vertex operator  () (Vertex & v) const
  {
    Dart_handle d = v.dart ();
    CGAL_assertion (d != NULL);

    int degree = 0;
    bool open = false;

    Map::One_dart_per_incident_cell_range<1,0>::iterator it (mmap, d),
      itend(mmap.one_dart_per_incident_cell<1,0>(d).end());
    for (; it != itend; ++it)
      {
	++degree;
	if (it->is_free (2)) open = true;
      }

    if (open)
      return v;

    Map::FT alpha = (4.0f - 2.0f *
           (Map::FT) cos (2.0f * PI / (Map::FT) degree)) / 9.0f;
    Map::Vector vec = (v - CGAL::ORIGIN) * (1.0f - alpha);

    for (it.rewind (); it != itend; ++it)
      {
	CGAL_assertion (!it->is_free (2));
	vec = vec + (mmap.point(it->beta(2)) - CGAL::ORIGIN)
	  * alpha / degree;
      }

    Vertex res (CGAL::ORIGIN + vec);
    res.set_dart (d);

    //  std::cout<<"operator() "<<v.point()<<" -> "<<res.point()<<std::endl;

    return res;
  }
private:
  Map & mmap;
};

// Flip an edge, work in 2D and in 3D.
Dart_handle
flip_edge (Map & m, Dart_handle d)
{
  CGAL_assertion (d != NULL && !d->is_free (2));

  if (!CGAL::is_removable<Map,1>(m,d))
    return NULL;

  Dart_handle d2 = d->beta(1)->beta(1);
  CGAL::remove_cell<Map,1>(m, d);

  insert_cell_1_in_cell_2(m, d2, d2->beta(1)->beta(1));

  return d2->beta (0);
}

// Subdivide each facet of the map by using sqrt(3)-subdivision.
void
subdivide_map_3 (Map & m)
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
  for (Map::Dart_range::iterator it (m.darts().begin ());
       m.number_of_marked_darts (treated) > 0; ++it)
    {
      ++nb;
      if (m.is_marked (it, treated))
   {
     // We unmark the darts of the facet to process only once dart/facet.
     CGAL::unmark_cell < Map, 2 > (m, it, treated);
     // We triangulate the facet.
     CGAL::insert_center_cell_0_in_cell_2(m, it);
   }
    }

  CGAL_assertion (m.is_whole_map_unmarked (treated));
  CGAL_assertion (m.is_valid ());
  m.free_mark (treated);

  // 3) We update the coordinates of old vertices.
  for (std::vector < Vertex >::iterator vit = vertices.begin ();
       vit != vertices.end (); ++vit)
    {
      m.point(vit->dart())=*vit;
    }

  // 4) We flip all the old edges.
  m.negate_mark (mark);  // Now only new darts are marked.
  Dart_handle d2 = NULL;
  for (Map::Dart_range::iterator it (m.darts().begin ()); it != m.darts().end ();)
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
              m.negate_mark (mark);  // thus new darts will be marked
	      flip_edge (m, d2);
              m.negate_mark (mark);
           }
           else
              m.mark (d2, mark);
        }
        else
           m.mark (d2, mark);
     }
  }

  /*  CGAL::display_darts(m,std::cout)<<std::endl;
  for (Map::Vertex_attribute_iterator it = m.vertex_attributes_begin();
       it!=m.vertex_attributes_end(); ++it)
    {
      std::cout<<it->point()<<", ";
    }
    std::cout<<std::endl;*/

  CGAL_assertion (m.is_whole_map_marked (mark));
  m.free_mark (mark);
  
  CGAL_postcondition ( m.is_valid ());
}
