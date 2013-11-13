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
//		   Jérémy Girerd-Rey <jeremy.girerd-rey@etu.univ-lyon1.fr>
//

#include "typedefs.h"

// Smooth a vertex depending on the vertices of its incidents facets.
class Smooth_edge_pqq
{
public:
  /**  Constructor.
   * @param alcc is the lcc to smooth
   * @param amark is a mark designing old darts (i.e. darts not created during
   *        the subdivision step)
   */
  Smooth_edge_pqq (LCC & alcc, unsigned int o):mlcc (alcc), old(o)
  {
  }

  Vertex operator  () (Vertex & v) const
  {
    Dart_handle d = v.dart ();

    // Old points aren't concerned.
    if (mlcc.is_marked(d, old))
    {       
    return v;
    }

    std::vector<LCC::Point> facetsPoints;
    facetsPoints.resize(0);

    // We search barycenter point of incidents facets.
    for (LCC::One_dart_per_incident_cell_range<1,0>::iterator it =
           mlcc.one_dart_per_incident_cell<1,0>(d).begin();
         it != mlcc.one_dart_per_incident_cell<1,0>(d).end(); ++it)
    {
      // If the vertex is on a border.
      if (mlcc.is_free(it,2))
      { 
        return v;
      }
      // If we found barycenter of a facet.
      if (!mlcc.is_marked(mlcc.opposite(it), old))
      {       
        facetsPoints.push_back(mlcc.point(mlcc.opposite(it)));
      }
    }

    // If we found more than two points we are on a vertice barycenter of a facet.
    // They aren't concerned.
    if (facetsPoints.size() > 2 || facetsPoints.size() < 2)
    {       
      return v;
    }

    // Average.
    LCC::Vector averageFacetsV = LCC::Traits::Construct_vector() 
			          ( CGAL::ORIGIN, LCC::Point(0,0,0) );                                                   

    for (unsigned int i=0; i < facetsPoints.size(); i++)
    {
      averageFacetsV = LCC::Traits::Construct_sum_of_vectors()
      		       (averageFacetsV, 
       		       LCC::Traits::Construct_vector() (CGAL::ORIGIN, facetsPoints[i]));
     }

    averageFacetsV = LCC::Traits::Construct_scaled_vector()
                     ( averageFacetsV, (1.0f/facetsPoints.size()) );

    // Barycenter point.
    LCC::Vector barycenterV = LCC::Traits::Construct_sum_of_vectors()
      		              ( LCC::Traits::Construct_vector() (CGAL::ORIGIN, v.point()), averageFacetsV );

    barycenterV = LCC::Traits::Construct_scaled_vector()
                  ( barycenterV, (1.0f/2.0f) );

    Vertex res = LCC::Traits::Construct_translated_point() (CGAL::ORIGIN, barycenterV);
    res.set_dart (d);

    return res;
  }
private:
  LCC & mlcc;
  unsigned int old;
};

// Smooth an old vertex depending on the vertices of its incident facet and edge.
class Smooth_vertex_pqq
{
public:
  /**  Constructor.
   * @param alcc is the lcc to smooth
   * @param amark is a mark designing old darts (i.e. darts not created during
   *        the triangulation step)
   */
  Smooth_vertex_pqq (LCC & alcc, unsigned int o):mlcc (alcc), old(o)
  {
  }

  Vertex operator  () (Vertex & v) const
  {
    Dart_handle d = v.dart ();

    // Just old points are concerned.
    if (!mlcc.is_marked(d, old))
    {       
      return v;
    }

    unsigned int degree = 0;
    std::vector<LCC::Point> edgesPoints;
    edgesPoints.resize(0);    
    
    std::vector<LCC::Point> facetsPoints;
    facetsPoints.resize(0);    

    // We search barycenter point of incidents facets, and incidents edges points.
    for (LCC::One_dart_per_incident_cell_range<1,0>::iterator it =
           mlcc.one_dart_per_incident_cell<1,0>(d).begin();
        it != mlcc.one_dart_per_incident_cell<1,0>(d).end(); ++it)
    {
      // If the vertex is on a border
      if (mlcc.is_free(it,2))
      { 
        return v;
      }
      // If incident isn't an old point, it's an edge point.
      if (!mlcc.is_marked(mlcc.opposite(it), old))
      {       
        edgesPoints.push_back (mlcc.point(mlcc.opposite(it)));
      }
      // We go find the "facet point" of incidents facet (barycenter of a facet).
      facetsPoints.push_back (mlcc.point(mlcc.beta(mlcc.opposite(it), 0)));
      ++degree;
    }

    CGAL_assertion (facetsPoints.size() != 0 && edgesPoints.size() != 0);

    if (facetsPoints.size() < 3 || edgesPoints.size() < 3 )
    {       
      return v;
    }

    // Average of incidents "edge points".
    LCC::Vector averageEdgesV = LCC::Traits::Construct_vector() 
			        (CGAL::ORIGIN, LCC::Point(0,0,0));                                                   

    for (unsigned int i=0; i < edgesPoints.size(); i++)
    {
      averageEdgesV = LCC::Traits::Construct_sum_of_vectors()
                      ( averageEdgesV, 
       		      LCC::Traits::Construct_vector() (CGAL::ORIGIN, edgesPoints[i]));
    }

    averageEdgesV = LCC::Traits::Construct_scaled_vector()
                    ( averageEdgesV, 1.0f/ (LCC::FT) edgesPoints.size() );

    // Average of incidents "facet points".
    LCC::Vector averageFacetsV = LCC::Traits::Construct_vector() 
			         (CGAL::ORIGIN, LCC::Point(0,0,0));                                                   

    for (unsigned int i=0; i < facetsPoints.size(); i++)
    {
      averageFacetsV = LCC::Traits::Construct_sum_of_vectors()
      		       ( averageFacetsV, 
       		       LCC::Traits::Construct_vector() (CGAL::ORIGIN, facetsPoints[i]));
    }

    averageFacetsV = LCC::Traits::Construct_scaled_vector()
                     ( averageFacetsV, (1.0f/ (LCC::FT) facetsPoints.size()) );

    // COEFFICIENTS of PQQ - Catmull–Clark subdivision : 
    // point = ( averageFacets + 2*averageEdges + point*(degree-3) )/degree

    averageFacetsV = LCC::Traits::Construct_scaled_vector()
      		     ( averageFacetsV, 1.0f/ (LCC::FT) degree);

    averageEdgesV = LCC::Traits::Construct_scaled_vector()
      		    ( averageEdgesV, 2.0f/ (LCC::FT) degree);

    LCC::Vector pointV = LCC::Traits::Construct_scaled_vector()
      		         ( LCC::Traits::Construct_vector() (CGAL::ORIGIN, v.point() ), 
                         (LCC::FT) (degree-3)/ (LCC::FT) degree );

    // New position of the old point.
    LCC::Vector	newPosition = LCC::Traits::Construct_sum_of_vectors()
      		              ( averageFacetsV, averageEdgesV);
    
    newPosition = LCC::Traits::Construct_sum_of_vectors()
      		  ( newPosition, pointV);

    Vertex res = LCC::Traits::Construct_translated_point() (CGAL::ORIGIN, newPosition);
    res.set_dart (d);

    return res;
  }
private:
  LCC & mlcc;
  unsigned int old;
};


// Subdivide each facet of the lcc by using pqq-subdivision.
void
subdivide_lcc_pqq (LCC & m)
{
  if (m.number_of_darts () == 0)
    return;

  unsigned int old = m.get_new_mark ();
  unsigned int treated = m.get_new_mark ();
  m.negate_mark (old);  // All the old darts are marked in O(1).

  // 1) We subdivide each edge.
  m.negate_mark (treated);  // All the darts are marked in O(1).
  
  for (LCC::Dart_range::iterator it (m.darts().begin ());
       m.number_of_marked_darts (treated) > 0; ++it)
  {
    if (m.is_marked (it, treated))
    {
      // We unmark the darts of the facet to process only once dart/facet.
      CGAL::unmark_cell < LCC, 1 > (m, it, treated);

      // We insert barycenter in the middle of the edge.
      m.insert_barycenter_in_cell<1>(it);
    }
  }

  // 2) We create a barycenter point for each facets.
  Dart_handle dc;
  std::vector<Dart_handle> remove;
  remove.resize(0);

  m.negate_mark (treated);  // All the darts are marked in O(1).
  for (LCC::Dart_range::iterator it (m.darts().begin ());
       m.number_of_marked_darts (treated) > 0; ++it)
  {
    if (m.is_marked (it, treated))
    {
      // We unmark the darts of the facet to process only once dart/facet.
      CGAL::unmark_cell < LCC, 2 > (m, it, treated);
      // We insert barycenter of the facet.
      dc =  m.insert_barycenter_in_cell<2>(it);

      // We remove useless edges.
      for (LCC::One_dart_per_incident_cell_range<1,0>::iterator it2 =
             m.one_dart_per_incident_cell<1,0>(dc).begin(); 
          it2 != m.one_dart_per_incident_cell<1,0>(dc).end(); ++it2)
      {
        // If the edge join the center and a corner.
        // We remove the edge.
        if( m.is_marked(m.beta(it2,1), old) )
        {
          remove.push_back(it2);
        }
      }

      // Remove edges. 
      for (std::vector <Dart_handle>::iterator dit = remove.begin ();
          dit != remove.end (); ++dit)
      {
        CGAL_assertion( (CGAL::is_removable<LCC, 1>(m,*dit)) );
        CGAL::remove_cell<LCC, 1>(m, *dit);
      }
      remove.resize(0); 
      // CGAL_assertion( m.is_valid() );
    }
  }

  m.negate_mark (treated);
  CGAL_assertion (m.is_whole_map_marked (treated));
  m.free_mark (treated);

  // 3) Smooth old points. 
  std::vector < Vertex > old_vertices; // smooth the old vertices.
  old_vertices.reserve (m.number_of_attributes<0> ()); // get intermediate space.
  std::transform (m.vertex_attributes().begin (), 
		  m.vertex_attributes().end (),
		  std::back_inserter (old_vertices), 
		  Smooth_vertex_pqq (m, old));

  // Update.
  for (std::vector < Vertex >::iterator vit = old_vertices.begin ();
      vit != old_vertices.end (); ++vit)
  {
    m.point(vit->dart())=vit->point();
  }

  // 4) Smooth new edges points.	  
  std::vector < Vertex > vertices; // smooth the old vertices.
  vertices.reserve (m.number_of_attributes<0> ()); // get intermediate space.
  std::transform (m.vertex_attributes().begin (), 
		  m.vertex_attributes().end (),
		  std::back_inserter (vertices), 
		  Smooth_edge_pqq (m, old));

  // Update.
  for (std::vector < Vertex >::iterator vit = vertices.begin ();
      vit != vertices.end (); ++vit)
  {
    m.point(vit->dart())=vit->point();
  }

  m.unmark_all (old);
  m.free_mark (old);
  CGAL_postcondition ( m.is_valid ());
}

