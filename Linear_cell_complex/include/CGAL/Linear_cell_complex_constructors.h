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
#ifndef CGAL_LINEAR_CELL_COMPLEX_CONSTRUCTORS_H
#define CGAL_LINEAR_CELL_COMPLEX_CONSTRUCTORS_H 1

#include <CGAL/Combinatorial_map_constructors.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/File_header_OFF.h>
#include <CGAL/IO/File_scanner_OFF.h>
#include <CGAL/Linear_cell_complex_incremental_builder.h>
#include <iostream>
#include <map>
#include <vector>
#include <list>

namespace CGAL {

  /** @file Linear_cell_complex_constructors.h
   * Some construction operations for a linear cell complex from other
   * CGAL data structures.
   */

  /** Import an embedded plane graph read into a flux into a
   *  linear cell complex.
   * @param alcc the linear cell complex where the graph will be imported.
   * @param ais the istream where read the graph.
   * @return A dart created during the convertion.
   */
  template< class LCC >
  typename LCC::Dart_handle import_from_plane_graph(LCC& alcc,
                                                    std::istream& ais)
  {
    CGAL_static_assertion( LCC::dimension>=2 && LCC::ambient_dimension==2 );
  
    typedef typename LCC::Dart_handle Dart_handle;
    typedef typename LCC::Traits::Direction_2 Direction;
    typedef typename std::list<Dart_handle>::iterator List_iterator;
    typedef typename std::map<Direction, Dart_handle>::iterator LCC_iterator;
  
    // Arrays of vertices
    std::vector< typename LCC::Vertex_attribute_handle > initVertices;
    std::vector< std::list<Dart_handle> > testVertices;

    std::string txt;
    typename LCC::FT x, y;
    Dart_handle d1 = alcc.null_handle, d2 = alcc.null_handle;
    unsigned int v1, v2;
  
    unsigned int nbSommets = 0;
    unsigned int nbAretes = 0;
  
    ais >> nbSommets >> nbAretes;
    while (nbSommets > 0)
    {
      if (!ais.good())
      {
        std::cout << "Problem: file does not contain enough vertices."
                  << std::endl;
        return alcc.null_handle;
      }

      ais >> iformat(x) >> iformat(y);
      initVertices.push_back(alcc.create_vertex_attribute
                             (typename LCC::Point(x, y)));
      testVertices.push_back(std::list<Dart_handle>());
      --nbSommets;
    }

    while (nbAretes > 0)
    {
      if (!ais.good())
      {
        std::cout << "Problem: file does not contain enough edges."
                  << std::endl;
        return alcc.null_handle;
      }

      // We read an egde (given by the number of its two vertices).
      ais >> v1 >> v2;
      --nbAretes;

      CGAL_assertion(v1 < initVertices.size());
      CGAL_assertion(v2 < initVertices.size());

      d1 = alcc.create_dart(initVertices[v1]);
      d2 = alcc.create_dart(initVertices[v2]);
      alcc.template link_beta<2>(d1, d2);

      testVertices[v1].push_back(d1);
      testVertices[v2].push_back(d2);
    }

    // LCC associating directions and darts.
    std::map<Direction, Dart_handle> tabDart;
    List_iterator it;
    LCC_iterator  it2;

    Dart_handle first = alcc.null_handle;
    Dart_handle prec = alcc.null_handle;
    typename LCC::Point sommet1, sommet2;
  
    for (unsigned int i = 0; i < initVertices.size(); ++i)
    {
      it = testVertices[i].begin();
      if (it != testVertices[i].end()) // Si la liste n'est pas vide.
      {
        // 1. We insert all the darts and sort them depending on the direction
        tabDart.clear();
      
        sommet1 = alcc.point(*it);
        sommet2 = alcc.point(alcc.beta(*it,2));
      
        tabDart.insert(std::pair<Direction, Dart_handle>
                       (typename LCC::Traits::Construct_direction_2()
                        (typename LCC::Traits::Construct_vector()
                         (sommet1,sommet2)), *it));
      
        ++it;
        while (it != testVertices[i].end())
        {
          sommet2 = alcc.point(alcc.beta(*it,2));
          tabDart.insert(std::pair<Direction, Dart_handle>
                         (typename LCC::Traits::Construct_direction_2()
                          (typename LCC::Traits::Construct_vector()
                           (sommet1,sommet2)), *it));
          ++it;
        }
      
        // 2. We run through the array of darts and 1 links darts.
        it2 = tabDart.begin();
        first = it2->second;
        prec = first;
        ++it2;

        while (it2 != tabDart.end())
        {
          alcc.template link_beta<0>(prec, alcc.beta(it2->second,2));
          prec = it2->second;
          ++it2;
        }
        alcc.template link_beta<0>(prec, alcc.beta(first,2));
      }
    }

    // We return a dart from the imported object.
    return first;
  }

  /** Convert a given Triangulation_2 into a 2D linear cell complex.
   * @param alcc the used linear cell complex.
   * @param atr the Triangulation_2.
   * @param aface_to_dart a pointer to a std::map associating to each
   *        triangle of atr a corresponding dart in alcc. Not used if NULL.
   * @return A dart incident to the infinite vertex.
   */
  template < class LCC, class Triangulation >
  typename LCC::Dart_handle import_from_triangulation_2
  (LCC& alcc, const Triangulation &atr,
   std::map<typename Triangulation::Face_handle,
            typename LCC::Dart_handle >* aface_to_dart=NULL)
  {
    CGAL_static_assertion( LCC::dimension>=2 && LCC::ambient_dimension==2 );
    
    // Case of empty triangulations.
    if (atr.number_of_vertices() == 0) return LCC::null_handle;

    // Check the dimension.
    if (atr.dimension() != 2) return LCC::null_handle;
    CGAL_assertion(atr.is_valid());

    typedef typename Triangulation::Vertex_handle         TVertex_handle;
    typedef typename Triangulation::All_vertices_iterator TVertex_iterator;
    typedef typename Triangulation::All_faces_iterator    TFace_iterator;
    typedef typename std::map
      < TFace_iterator, typename LCC::Dart_handle >::iterator itmap_tcell;

    // Create vertices in the map and associate in a map
    // TVertex_handle and vertices in the map.
    std::map< TVertex_handle, typename LCC::Vertex_attribute_handle > TV;
    for (TVertex_iterator itv = atr.all_vertices_begin();
         itv != atr.all_vertices_end(); ++itv)
    {
      TV[itv] = alcc.create_vertex_attribute(itv->point());
    }

    // Create the triangles and create a map to link Cell_iterator
    // and triangles.
    TFace_iterator it;

    std::map<typename Triangulation::Face_handle, typename LCC::Dart_handle> TC;
    std::map<typename Triangulation::Face_handle, typename LCC::Dart_handle>*
      mytc = (aface_to_dart==NULL?&TC:aface_to_dart);
    
    itmap_tcell maptcell_it;

    typename LCC::Dart_handle res=LCC::null_handle, dart=LCC::null_handle;
    typename LCC::Dart_handle cur=LCC::null_handle, neighbor=LCC::null_handle;

    for (it = atr.all_faces_begin(); it != atr.all_faces_end(); ++it)
    {
      /*     if (it->vertex(0) != atr.infinite_vertex() &&
             it->vertex(1) != atr.infinite_vertex() &&
             it->vertex(2) != atr.infinite_vertex() &&
             it->vertex(3) != atr.infinite_vertex())
      */
      {
        res = alcc.make_triangle(TV[it->vertex(0)],
                                 TV[it->vertex(1)],
                                 TV[it->vertex(2)]);

        if ( dart==LCC::null_handle )
        {
          if ( it->vertex(0) == atr.infinite_vertex() )
            dart = res;
          else if ( it->vertex(1) == atr.infinite_vertex() )
            dart = alcc.beta(res,1);
          else if ( it->vertex(2) == atr.infinite_vertex() )
            dart = alcc.beta(res,0);
        }
        
        for (unsigned int i=0; i<3; ++i)
        {
          switch (i)
          {
          case 0: cur = alcc.beta(res,1); break;
          case 1: cur = alcc.beta(res,0); break;
          case 2: cur = res; break;
          }

          maptcell_it = mytc->find(it->neighbor(i));
          if (maptcell_it != mytc->end())
          {
            switch (atr.mirror_index(it,i) )
            {
            case 0: neighbor = alcc.beta(maptcell_it->second,1);
              break;
            case 1: neighbor = alcc.beta(maptcell_it->second,0);
              break;
            case 2: neighbor = maptcell_it->second; break;
            }
            alcc.template topo_sew<2>(cur, neighbor);
          }
        }
        (*mytc)[it] = res;
      }
    }

    CGAL_assertion(dart!=LCC::null_handle);
    return dart;
  }
  
  /** Convert a given Triangulation_3 into a 3D linear cell complex.
   * @param alcc the used linear cell complex.
   * @param atr the Triangulation_3.
   * @param avol_to_dart a pointer to a std::map associating to each
   *        tetrahedron of atr a corresponding dart in alcc. Not used if NULL.
   * @return A dart incident to the infinite vertex.
   */
  template < class LCC, class Triangulation >
  typename LCC::Dart_handle import_from_triangulation_3
  (LCC& alcc, const Triangulation &atr,
   std::map<typename Triangulation::Cell_handle,
            typename LCC::Dart_handle >* avol_to_dart=NULL)
  {
    CGAL_static_assertion( LCC::dimension>=3 && LCC::ambient_dimension==3 );
    
    // Case of empty triangulations.
    if (atr.number_of_vertices() == 0) return LCC::null_handle;

    // Check the dimension.
    if (atr.dimension() != 3) return LCC::null_handle;
    CGAL_assertion(atr.is_valid());

    typedef typename Triangulation::Vertex_handle    TVertex_handle;
    typedef typename Triangulation::Vertex_iterator  TVertex_iterator;
    typedef typename Triangulation::Cell_iterator    TCell_iterator;
    typedef typename std::map
      < TCell_iterator, typename LCC::Dart_handle >::iterator itmap_tcell;

    // Create vertices in the map and associate in a map
    // TVertex_handle and vertices in the map.
    std::map< TVertex_handle, typename LCC::Vertex_attribute_handle > TV;
    for (TVertex_iterator itv = atr.vertices_begin();
         itv != atr.vertices_end(); ++itv)
    {
      TV[itv] = alcc.create_vertex_attribute(itv->point());
    }

    // Create the tetrahedron and create a map to link Cell_iterator
    // and tetrahedron.
    TCell_iterator it;

    std::map<typename Triangulation::Cell_handle, typename LCC::Dart_handle> TC;
    std::map<typename Triangulation::Cell_handle, typename LCC::Dart_handle>*
      mytc = (avol_to_dart==NULL?&TC:avol_to_dart);
    
    itmap_tcell maptcell_it;

    typename LCC::Dart_handle res=LCC::null_handle, dart=LCC::null_handle;
    typename LCC::Dart_handle cur=LCC::null_handle, neighbor=LCC::null_handle;

    for (it = atr.cells_begin(); it != atr.cells_end(); ++it)
    {
      /*     if (it->vertex(0) != atr.infinite_vertex() &&
             it->vertex(1) != atr.infinite_vertex() &&
             it->vertex(2) != atr.infinite_vertex() &&
             it->vertex(3) != atr.infinite_vertex())
      */
      {
        res = alcc.make_tetrahedron(TV[it->vertex(0)],
                                    TV[it->vertex(1)],
                                    TV[it->vertex(2)],
                                    TV[it->vertex(3)]);

        if ( dart==LCC::null_handle )
        {
          if ( it->vertex(0) == atr.infinite_vertex() )
            dart = res;
          else if ( it->vertex(1) == atr.infinite_vertex() )
            dart = alcc.beta(res, 1);
          else if ( it->vertex(2) == atr.infinite_vertex() )
            dart = alcc.beta(res, 1, 1);
          else if ( it->vertex(3) == atr.infinite_vertex() )
            dart = alcc.beta(res, 2, 0);
        }
        
        for (unsigned int i = 0; i < 4; ++i)
        {
          switch (i)
          {
          case 0: cur = alcc.beta(res, 1, 2); break;
          case 1: cur = alcc.beta(res, 0, 2); break;
          case 2: cur = alcc.beta(res, 2); break;
          case 3: cur = res; break;
          }

          maptcell_it = mytc->find(it->neighbor(i));
          if (maptcell_it != mytc->end())
          {
            switch (atr.mirror_index(it,i) )
            {
            case 0: neighbor = alcc.beta(maptcell_it->second, 1, 2);
              break;
            case 1: neighbor = alcc.beta(maptcell_it->second, 0, 2);
              break;
            case 2: neighbor = alcc.beta(maptcell_it->second, 2); break;
            case 3: neighbor = maptcell_it->second; break;
            }
            while (alcc.temp_vertex_attribute(neighbor) !=
                   alcc.temp_vertex_attribute(alcc.other_extremity(cur)) )
              neighbor = alcc.beta(neighbor,1);
            alcc.template topo_sew<3>(cur, neighbor);
          }
        }
        (*mytc)[it] = res;
      }
    }
    CGAL_assertion(dart!=LCC::null_handle);
    return dart;
  }

  /** Import a given Polyhedron_3 into a Linear_cell_complex.
   * @param alcc the linear cell complex where Polyhedron_3 will be converted.
   * @param apoly the Polyhedron.
   * @return A dart created during the convertion.
   */
  template< class LCC, class Polyhedron >
  typename LCC::Dart_handle import_from_polyhedron_3(LCC& alcc, 
                                                     const Polyhedron &apoly)
  {
    CGAL_static_assertion( LCC::dimension>=2 && LCC::ambient_dimension==3 );

    typedef typename Polyhedron::Halfedge_const_handle  Halfedge_handle;
    typedef typename Polyhedron::Facet_const_iterator   Facet_iterator;
    typedef typename Polyhedron::Halfedge_around_facet_const_circulator
      HF_circulator;

    typedef std::map < Halfedge_handle, typename LCC::Dart_handle> 
      Halfedge_handle_map;
    typedef typename Halfedge_handle_map::iterator itmap_hds;
    Halfedge_handle_map TC;

    itmap_hds it;
    typename LCC::Dart_handle d = LCC::null_handle, prev = LCC::null_handle;
    typename LCC::Dart_handle firstFacet = LCC::null_handle, firstAll = LCC::null_handle;

    // First traversal to build the darts and link them.
    for (Facet_iterator i = apoly.facets_begin(); i != apoly.facets_end(); ++i)
    {
      HF_circulator j = i->facet_begin();
      prev = LCC::null_handle;
      do
      {
        d = alcc.create_dart();
        TC[j] = d;
      
        if (prev != LCC::null_handle) alcc.template link_beta<1>(prev, d);
        else firstFacet = d;
        it = TC.find(j->opposite());
        if (it != TC.end())
          alcc.template link_beta<2>(d, it->second);
        prev = d;
      }
      while (++j != i->facet_begin());
      alcc.template link_beta<1>(prev, firstFacet);
      if (firstAll == LCC::null_handle) firstAll = firstFacet;
    }

    // Second traversal to update the geometry.
    // We run one again through the facets of the HDS.
    for (Facet_iterator i = apoly.facets_begin(); i != apoly.facets_end(); ++i)
    {
      HF_circulator j = i->facet_begin();
      do
      {
        d = TC[j]; // Get the dart associated to the Halfedge
        if (alcc.temp_vertex_attribute(d) == LCC::null_handle)
        {
          alcc.set_vertex_attribute
            (d, alcc.create_vertex_attribute(j->opposite()->vertex()->point()));
        }
      }
      while (++j != i->facet_begin());
    }
    return firstAll;
  }

  template < class LCC >
  void load_off(LCC& alcc, std::istream& in)
  {
    File_header_OFF  m_file_header;
    File_scanner_OFF scanner( in, m_file_header.verbose());
    if ( ! in) return;
    m_file_header = scanner;  // Remember file header after return.

    Linear_cell_complex_incremental_builder_3<LCC> B( alcc);
    B.begin_surface( scanner.size_of_vertices(),
                     scanner.size_of_facets(),
                     scanner.size_of_halfedges());

    typedef typename LCC::Point Point;

    // read in all vertices
    std::size_t  i;
    for ( i = 0; i < scanner.size_of_vertices(); i++) {
      Point p;
      file_scan_vertex( scanner, p);
      B.add_vertex( p);
      scanner.skip_to_next_vertex( i);
    }
    /* TODO rollback
       if ( ! in  || B.error()) {
       B.rollback();
       in.clear( std::ios::badbit);
       return;
       }
    */

    // read in all facets
    for ( i = 0; i < scanner.size_of_facets(); i++)
    {
      B.begin_facet();
      std::size_t no;
      scanner.scan_facet( no, i);
      /* TODO manage errors
         if( ! in || B.error() || no < 3) {
         if ( scanner.verbose()) {
         std::cerr << " " << std::endl;
         std::cerr << "Polyhedron_scan_OFF<Traits>::" << std::endl;
         std::cerr << "operator()(): input error: facet " << i
         << " has less than 3 vertices." << std::endl;
         }
         B.rollback();
         in.clear( std::ios::badbit);
         return;
         } */
      for ( std::size_t j = 0; j < no; j++) {
        std::size_t index;
        scanner.scan_facet_vertex_index( index, i);
        B.add_vertex_to_facet( index);
      }
      B.end_facet();
      scanner.skip_to_next_facet( i);
    }
    /* TODO manage errors
       if ( ! in  || B.error()) {
       B.rollback();
       in.clear( std::ios::badbit);
       return;
       }
       if ( B.check_unconnected_vertices()) {
       if ( ! B.remove_unconnected_vertices()) {
       if ( scanner.verbose()) {
       std::cerr << " " << std::endl;
       std::cerr << "Polyhedron_scan_OFF<Traits>::" << std::endl;
       std::cerr << "operator()(): input error: cannot "
       "succesfully remove isolated vertices."
       << std::endl;
       }
       B.rollback();
       in.clear( std::ios::badbit);
       return;
       }
       }*/
    B.end_surface();
  }

  /** Convert a Polyhedron_3 read into a flux into 3D linear cell complex.
   * @param alcc the linear cell complex where Polyhedron_3 will be converted.
   * @param ais the istream where read the Polyhedron_3.
   * @return A dart created during the convertion.
   */
  template < class LCC >
  typename LCC::Dart_handle
  import_from_polyhedron_3_flux(LCC& alcc, std::istream& ais)
  {
    if (!ais.good())
    {
      std::cout << "Error reading flux." << std::endl;
      return LCC::null_handle;
    }
    CGAL::Polyhedron_3<typename LCC::Traits> P;
    ais >> P;
    return import_from_polyhedron_3<LCC, CGAL::Polyhedron_3
                                    <typename LCC::Traits> > (alcc, P);
  }

  /** Export the alcc in off file format. If dimension>2, export all faces but only once.
   */
  template < class LCC >
  void write_off(LCC& alcc, std::ostream& out)
  {
    File_header_OFF header(false);
    header.set_binary(is_binary( out));
    header.set_no_comments(!is_pretty( out));
    File_writer_OFF writer( header);
    writer.header().set_polyhedral_surface(true);
    writer.header().set_halfedges( alcc.number_of_darts());

    // Print header.
    writer.write_header( out,
                         alcc.number_of_vertex_attributes(),
                         alcc.number_of_darts(),
                         alcc.template one_dart_per_cell<2>().size() );

    typedef typename LCC::Vertex_attribute_range::iterator VCI;
    VCI vit, vend = alcc.vertex_attributes().end();
    for ( vit = alcc.vertex_attributes().begin(); vit!=vend; ++vit )
    {
      writer.write_vertex( ::CGAL::to_double( vit->point().x()),
                           ::CGAL::to_double( vit->point().y()),
                           ::CGAL::to_double( vit->point().z()));
    }

    typedef Inverse_index< VCI > Index;
    Index index( alcc.vertex_attributes().begin(),
                 alcc.vertex_attributes().end());
    writer.write_facet_header();

    int m = alcc.get_new_mark();

    for ( typename LCC::Dart_range::iterator itall = alcc.darts().begin(),
            itallend = alcc.darts().end(); itall!=itallend; ++itall )
    {
      if ( !alcc.is_marked(itall, m) )
      {
        std::size_t n = 0;
        // First we count the number of vertices of the face.
        for ( typename LCC::template Dart_of_orbit_range<1>::iterator
                itf=alcc.template darts_of_orbit<1>(itall).begin(),
                itfend=alcc.template darts_of_orbit<1>(itall).end();
              itf!=itfend; ++itf, ++n );

        CGAL_assertion( n>=3 );
        writer.write_facet_begin(n);

        // Second we write the indices of vertices.
        for ( typename LCC::template Dart_of_orbit_range<1>::iterator
                itf=alcc.template darts_of_orbit<1>(itall).begin(),
                itfend=alcc.template darts_of_orbit<1>(itall).end();
              itf!=itfend; ++itf )
        {
          // TODO case with index
          writer.write_facet_vertex_index(index[VCI(alcc.vertex_attribute(itf))]);

          for ( typename LCC::template Dart_of_involution_basic_range<1>::iterator
                  itinv=alcc.template darts_of_involution_basic<1>(itf, m).begin(),
                  itinvend=alcc.template darts_of_involution_basic<1>(itf, m).end();
                itinv!=itinvend; ++itinv )
            alcc.mark(itinv, m);
        }
        writer.write_facet_end();
      }
    }
    writer.write_footer();
    alcc.free_mark(m);
  }

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_CONSTRUCTORS_H //
// EOF //
