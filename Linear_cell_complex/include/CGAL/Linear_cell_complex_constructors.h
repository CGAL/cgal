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
#ifndef CGAL_COMBINATORIAL_MAP_WITH_POINTS_CONSTRUCTORS_H
#define CGAL_COMBINATORIAL_MAP_WITH_POINTS_CONSTRUCTORS_H 1

#include <CGAL/Combinatorial_map_constructors.h>
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

/** @file Combinatorial_map_with_embedding_constructors.h
 * Basic construction operations for an embedded combinatorial map.
 * create edge, triangle, quadrilateral, tetrahedron,hexahedron, plus
 * contruction from other CGAL data structures.
 */

/** Create an edge given 2 Vertex_attribute_handle.
 * @param amap the used combinatorial map.
 * @param ah1 the first vertex handle.
 * @param ah2 the second vertex handle.
 * @return the dart of the new edge incident to ah1.
 */
template <class CMap>
typename CMap::Dart_handle 
make_segment(CMap& amap,
	     typename CMap::Vertex_attribute_handle ah1,
	     typename CMap::Vertex_attribute_handle ah2)
{
   typename CMap::Dart_handle d1 = make_edge(amap);

   amap.set_vertex_attribute_of_dart(d1,ah1);
   amap.set_vertex_attribute_of_dart(d1->beta(2),ah2);

   return d1;
}

/** Create a segment given 2 points.
 * @param amap the used combinatorial map.
 * @param ap1 the first point.
 * @param ap2 the second point.
 * @return the dart of the new segment incident to ap1.
 */
template <class CMap>
typename CMap::Dart_handle make_segment(CMap& amap,
					const typename CMap::Point& ap1,
					const typename CMap::Point& ap2)
{
  return make_segment(amap,
		      amap.create_vertex_attribute(ap1),
		      amap.create_vertex_attribute(ap2));
}

/** Create a triangle given 3 Vertex_attribute_handle.
 * @param amap the used combinatorial map.
 * @param ah1 the first vertex handle.
 * @param ah2 the second vertex handle.
 * @param ah3 the third vertex handle.
 * @return the dart of the new triangle incident to ah1.
 */
template <class CMap>
typename CMap::Dart_handle 
make_triangle(CMap& amap,
	      typename CMap::Vertex_attribute_handle ah1,
	      typename CMap::Vertex_attribute_handle ah2,
	      typename CMap::Vertex_attribute_handle ah3)
{
  typename CMap::Dart_handle d1 = make_combinatorial_polygon(amap,3);
  
  amap.set_vertex_attribute_of_dart(d1,ah1);
  amap.set_vertex_attribute_of_dart(d1->beta(1),ah2);
  amap.set_vertex_attribute_of_dart(d1->beta(0),ah3);
  
  return d1;
}

/** Create a triangle given 3 points.
 * @param amap the used combinatorial map.
 * @param ap1 the first point.
 * @param ap2 the second point.
 * @param ap3 the third point.
 * @return the dart of the new triangle incident to ap1.
 */
template <class CMap>
typename CMap::Dart_handle make_triangle(CMap& amap,
					 const typename CMap::Point& ap1,
					 const typename CMap::Point& ap2,
					 const typename CMap::Point& ap3)
{
   return make_triangle(amap,
                        amap.create_vertex_attribute(ap1),
                        amap.create_vertex_attribute(ap2),
                        amap.create_vertex_attribute(ap3));
}

/** Create a quadrilateral given 4 Vertex_attribute_handle.
 * @param amap the used combinatorial map.
 * @param ah1 the first vertex handle.
 * @param ah2 the second vertex handle.
 * @param ah3 the third vertex handle.
 * @param ah4 the fourth vertex handle.
 * @return the dart of the new quadrilateral incident to ah1.
 */
template <class CMap>
typename CMap::Dart_handle 
make_quadrilateral(CMap& amap,
		   typename CMap::Vertex_attribute_handle ah1,
		   typename CMap::Vertex_attribute_handle ah2,
		   typename CMap::Vertex_attribute_handle ah3,
		   typename CMap::Vertex_attribute_handle ah4)
{
  typename CMap::Dart_handle d1 = make_combinatorial_polygon(amap,4);

   amap.set_vertex_attribute_of_dart(d1,ah1);
   amap.set_vertex_attribute_of_dart(d1->beta(1),ah2);
   amap.set_vertex_attribute_of_dart(d1->beta(1)->beta(1),ah3);
   amap.set_vertex_attribute_of_dart(d1->beta(0),ah4);

   return d1;
}

/** Create a quadrilateral given 4 points.
 * @param amap the used combinatorial map.
 * @param ap1 the first point.
 * @param ap2 the second point.
 * @param ap3 the third point.
 * @param ap4 the fourth point.
 * @return the dart of the new quadrilateral incident to ap1.
 */
template <class CMap>
typename CMap::Dart_handle make_quadrilateral(CMap& amap,
					      const typename CMap::Point& ap1,
					      const typename CMap::Point& ap2,
					      const typename CMap::Point& ap3,
					      const typename CMap::Point& ap4)
{
   return make_quadrilateral(amap,
                             amap.create_vertex_attribute(ap1),
                             amap.create_vertex_attribute(ap2),
                             amap.create_vertex_attribute(ap3),
                             amap.create_vertex_attribute(ap4));
}

/** Create a rectangle given a Iso_rectangle_2
 * @param amap the used combinatorial map.
 * @param r the Iso_rectangle_2
 * @return the dart of the new rectangle incident to the first point of r.
 */
/* template <class CMap>
typename CMap::Dart_handle make_rectangle
(CMap& amap, const typename CMap::Iso_rectangle& r)
{ return make_quadrilateral(amap, r[0], r[1], r[2], r[3]); }*/

/** Create a rectangle given its two extreme points.
 * @param amap the used combinatorial map.
 * @param ap1 the first point.
 * @param ap1 the second point.
 * @return the dart of the new rectangle incident to ap1.
 */
 /*template <class CMap>
typename CMap::Dart_handle make_rectangle(CMap& amap,
					  const typename CMap::Point& ap1,
					  const typename CMap::Point& ap2)
{
  return make_rectangle<CMap>
    (amap, typename CMap::Iso_rectangle()(ap1, ap2));
		}*/

/** Create a square given one point and one length.
 * @param amap the used combinatorial map.
 * @param ap1 the first point.
 * @param al the length.
 * @return the dart of the new square incident to ap1.
 */
/*template < class Map >
typename Map::Dart_handle make_square(Map& amap,
                                      const typename Map::Point& ap1,
                                      typename Map::FT al)
				      { return make_rectangle(amap, ap1, al, al); }*/

/** Create a tetrahedron given 4 Vertex_attribute_handle.
 * @param amap the used combinatorial map.
 * @param ah1 the first vertex handle.
 * @param ah2 the second vertex handle.
 * @param ah3 the third vertex handle.
 * @param ah4 the fourth vertex handle.
 * @return the dart of the new tetrahedron incident to ah1 and to
 *         facet ah1,ah2,ah3.
 */
template < class Map >
typename Map::Dart_handle 
make_tetrahedron(Map& amap,
		 typename Map::Vertex_attribute_handle ah1,
		 typename Map::Vertex_attribute_handle ah2,
		 typename Map::Vertex_attribute_handle ah3,
		 typename Map::Vertex_attribute_handle ah4)
{
   typename Map::Dart_handle d1 = make_triangle(amap, ah1, ah2, ah3);
   typename Map::Dart_handle d2 = make_triangle(amap, ah2, ah1, ah4);
   typename Map::Dart_handle d3 = make_triangle(amap, ah2, ah4, ah3);
   typename Map::Dart_handle d4 = make_triangle(amap, ah4, ah1, ah3);

   return make_combinatorial_tetrahedron(amap, d1, d2, d3, d4);
}

/** Create a tetrahedron given 4 points.
 * @param amap the used combinatorial map.
 * @param ap1 the first point.
 * @param ap2 the second point.
 * @param ap3 the third point.
 * @param ap4 the fourth point.
 * @return the dart of the new tetrahedron incident to ap1 and to
 *         facet ap1,ap2,ap3.
 */
template < class Map >
typename Map::Dart_handle make_tetrahedron(Map& amap,
					   const typename Map::Point& ap1,
					   const typename Map::Point& ap2,
					   const typename Map::Point& ap3,
					   const typename Map::Point& ap4)
{
   return make_tetrahedron(amap,
                           amap.create_vertex_attribute(ap1),
                           amap.create_vertex_attribute(ap2),
                           amap.create_vertex_attribute(ap3),
                           amap.create_vertex_attribute(ap4));
}

/** Create an hexahedron given 8 Vertex_attribute_handle.
 *    (8 vertices, 12 edges and 6 facets)
 * \verbatim
 *       4----7
 *      /|   /|
 *     5----6 |
 *     | 3--|-2
 *     |/   |/
 *     0----1
 * \endverbatim
 * @param amap the used combinatorial map.
 * @param ah0 the first vertex handle.
 * @param ah1 the second vertex handle.
 * @param ah2 the third vertex handle.
 * @param ah3 the fourth vertex handle.
 * @param ah4 the fifth vertex handle.
 * @param ah5 the sixth vertex handle.
 * @param ah6 the seventh vertex handle.
 * @param ah7 the height vertex handle.
 * @return the dart of the new hexahedron incident to ah0 and to
 *         the facet (ah0,ah5,ah6,ah1).
 */
template < class Map >
typename Map::Dart_handle make_hexahedron(Map& amap,
      typename Map::Vertex_attribute_handle ah0,
      typename Map::Vertex_attribute_handle ah1,
      typename Map::Vertex_attribute_handle ah2,
      typename Map::Vertex_attribute_handle ah3,
      typename Map::Vertex_attribute_handle ah4,
      typename Map::Vertex_attribute_handle ah5,
      typename Map::Vertex_attribute_handle ah6,
      typename Map::Vertex_attribute_handle ah7)
{
  typename Map::Dart_handle d1 = make_quadrilateral(amap, ah0, ah5, ah6, ah1);
  typename Map::Dart_handle d2 = make_quadrilateral(amap, ah1, ah6, ah7, ah2);
  typename Map::Dart_handle d3 = make_quadrilateral(amap, ah2, ah7, ah4, ah3);
  typename Map::Dart_handle d4 = make_quadrilateral(amap, ah3, ah4, ah5, ah0);
  typename Map::Dart_handle d5 = make_quadrilateral(amap, ah0, ah1, ah2, ah3);
  typename Map::Dart_handle d6 = make_quadrilateral(amap, ah5, ah4, ah7, ah6);
  
  return make_combinatorial_hexahedron(amap, d1, d2, d3, d4, d5, d6);
}

/** Create an hexahedron given 8 points.
 * \verbatim
 *       4----7
 *      /|   /|
 *     5----6 |
 *     | 3--|-2
 *     |/   |/
 *     0----1
 * \endverbatim
 * @param amap the used combinatorial map.
 * @param ap0 the first point.
 * @param ap1 the second point.
 * @param ap2 the third point.
 * @param ap3 the fourth point.
 * @param ap4 the fifth point.
 * @param ap5 the sixth point.
 * @param ap6 the seventh point.
 * @param ap7 the height point.
 * @return the dart of the new hexahedron incident to ap0
 *         and to the facet (ap0,ap5,ap6,ap1).
 */
template <class Map>
typename Map::Dart_handle make_hexahedron(Map& amap,
      const typename Map::Point& ap0,
      const typename Map::Point& ap1,
      const typename Map::Point& ap2,
      const typename Map::Point& ap3,
      const typename Map::Point& ap4,
      const typename Map::Point& ap5,
      const typename Map::Point& ap6,
      const typename Map::Point& ap7)
{
   return make_hexahedron(amap,
                          amap.create_vertex_attribute(ap0),
			  amap.create_vertex_attribute(ap1),
                          amap.create_vertex_attribute(ap2),
                          amap.create_vertex_attribute(ap3),
                          amap.create_vertex_attribute(ap4),
                          amap.create_vertex_attribute(ap5),
                          amap.create_vertex_attribute(ap6),
                          amap.create_vertex_attribute(ap7));
}

/** Create an iso cuboid given an Iso_cuboid_3.
 * @param amap the used combinatorial map.
 * @param c the iso cuboid.
 * @return the dart of the new cuboid incident to the first vertex of c.
 */
/*template <class CMap>
typename CMap::Dart_handle make_iso_cuboid
(CMap& amap, const typename CMap::Iso_cuboid& r)
{
  return make_hexahedron<CMap>(amap, r[0], r[1], r[2], r[3],
			       r[4], r[5], r[6], r[7]);
						 }*/

/** Create an iso cuboid given its two extreme points.
 * @param amap the used combinatorial map.
 * @param ap1 the first point.
 * @param ap2 the second.
 * @return the dart of the new iso cuboid incident to ap1.
 */
 /*template <class CMap>
typename CMap::Dart_handle make_iso_cuboid(CMap& amap,
					   const typename CMap::Point& ap1,
					   const typename CMap::Point& ap2)
{
  return make_iso_cuboid<CMap>
    (amap, //typename CMap::Construct_iso_cuboid()(ap1, ap2));
     typename CMap::Iso_cuboid(ap1, ap2));
		 }*/

/** Create a cube given one point and one length.
 * @param amap the used combinatorial map.
 * @param ap1 the first point.
 * @param al the length.
 * @return the dart of the new cube incident to ap1.
 */
/*template < class Map >
typename Map::Dart_handle make_cube(Map& amap,
                                    const typename Map::Point& ap1,
                                    typename Map::FT al)
				    { return make_cuboid(amap, ap1, al, al, al); }*/


/** Convert an embedded plane graph read into a flux into combinatorial map.
 * @param amap the combinatorial map where the graph will be converted.
 * @param ais the istream where read the graph.
 * @return A dart created during the convertion.
 */
template< class Map >
typename Map::Dart_handle import_from_plane_graph(Map& amap,
      std::istream& ais)
{
   typedef typename Map::Dart_handle Dart_handle;
   typedef typename Map::Traits::Direction_2 Direction;
   typedef typename std::list<Dart_handle>::iterator List_iterator;
   typedef typename std::map<Direction, Dart_handle>::iterator Map_iterator;

   // Arrays of vertices
   std::vector< typename Map::Vertex_attribute_handle > initVertices;
   std::vector< std::list<Dart_handle> > testVertices;

   std::string txt;
   typename Map::FT x, y;
   Dart_handle d1 = NULL, d2 = NULL;
   unsigned int v1, v2;

   ais >> txt;
   if (txt != "OFF2D")
   {
      std::cout << "Problem: file is not 2D OFF." << std::endl;
      return NULL;
   }

   unsigned int nbSommets = 0;
   unsigned int nbAretes = 0;

   ais >> nbSommets >> nbAretes;
   while (nbSommets > 0)
   {
      if (!ais.good())
      {
         std::cout << "Problem: file does not contain enough vertices."
         << std::endl;
         return NULL;
      }

      ais >> x >> y;
      initVertices.push_back(amap.create_vertex_attribute(typename Map::Point(x, y)));
      testVertices.push_back(std::list<Dart_handle>());
      --nbSommets;
   }

   while (nbAretes > 0)
   {
      if (!ais.good())
      {
         std::cout << "Problem: file does not contain enough edges."
         << std::endl;
         return NULL;
      }

      // We read an egde (given by the number of its two vertices).
      ais >> v1 >> v2; ais.ignore(256, '\n');
      --nbAretes;

      CGAL_assertion(v1 < initVertices.size());
      CGAL_assertion(v2 < initVertices.size());

      d1 = amap.create_dart(initVertices[v1]);
      d2 = amap.create_dart(initVertices[v2]);
      amap.link_beta(d1, d2, 2);

      testVertices[v1].push_back(d1);
      testVertices[v2].push_back(d2);
   }

   // Map associating directions and darts.
   std::map<Direction, Dart_handle> tabDart;
   List_iterator it;
   Map_iterator  it2;

   Dart_handle first = NULL;
   Dart_handle prec = NULL;
   typename Map::Point sommet1, sommet2;

   for (unsigned int i = 0; i < initVertices.size(); ++i)
   {
      it = testVertices[i].begin();
      if (it != testVertices[i].end()) // Si la liste n'est pas vide.
      {
         // 1. We insert all the darts and sort them depending on the direction
         tabDart.clear();

         sommet1 = Map::point(*it);
         sommet2 = Map::point((*it)->beta(2));

         tabDart.insert(std::pair<Direction, Dart_handle>
                        (typename Map::Construct_direction()
			 (typename Map::Construct_vector()
			  (sommet1,sommet2)), *it));

         ++it;
         while (it != testVertices[i].end())
         {
	   sommet2 = Map::point((*it)->beta(2));
            tabDart.insert(std::pair<Direction, Dart_handle>
                        (typename Map::Construct_direction()
			 (typename Map::Construct_vector()
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
	   amap.template link_beta<0>(prec, it2->second->beta(2));
	   prec = it2->second;
	   ++it2;
         }
         amap.template link_beta<0>(prec, first->beta(2));
      }
   }

   // We return a dart from the imported object.
   return first;

}

/** Convert a given Triangulation_3 into the 3D combinatorial map.
 * @param amap the used combinatorial map.
 * @param atr the Triangulation_3.
 * @return A dart incident to the infinite vertex.
 */
template < class Map, class Triangulation >
typename Map::Dart_handle import_from_triangulation_3(Map& amap,
      const Triangulation &atr)
{
   // Case of empty triangulations.
   if (atr.number_of_vertices() == 0) return NULL;

   // Check the dimension.
   if (atr.dimension() != 3) return NULL;
   CGAL_assertion(atr.is_valid());

   typedef typename Triangulation::Vertex_handle    TVertex_handle;
   typedef typename Triangulation::Vertex_iterator  TVertex_iterator;
   typedef typename Triangulation::Cell_iterator    TCell_iterator;
   typedef typename
   std::map< TCell_iterator, typename Map::Dart_handle >::iterator itmap_tcell;

   // Create vertices in the map and associate in a map
   // TVertex_handle and vertices in the map.
   std::map< TVertex_handle, typename Map::Vertex_attribute_handle > TV;
   for (TVertex_iterator it = atr.vertices_begin();
         it != atr.vertices_end(); ++it)
   {
     //  if (it != atr.infinite_vertex())
      {
         TV[it] = amap.create_vertex_attribute(it->point());
      }
   }

   // Create the tetrahedron and create a map to link Cell_iterator
   // and tetrahedron.
   TCell_iterator it;

   std::map< TCell_iterator, typename Map::Dart_handle > TC;
   itmap_tcell maptcell_it;

   typename Map::Dart_handle res=NULL, dart=NULL;
   typename Map::Dart_handle init=NULL, cur=NULL, neighbor=NULL;

   for (it = atr.cells_begin(); it != atr.cells_end(); ++it)
   {
     /*     if (it->vertex(0) != atr.infinite_vertex() &&
            it->vertex(1) != atr.infinite_vertex() &&
            it->vertex(2) != atr.infinite_vertex() &&
            it->vertex(3) != atr.infinite_vertex())
     */
      {
         res = make_tetrahedron(amap,
                                TV[it->vertex(0)],
                                TV[it->vertex(1)],
                                TV[it->vertex(2)],
                                TV[it->vertex(3)]);

	 if ( it->vertex(0) == atr.infinite_vertex() && dart==NULL )
	   dart = res;

         for (unsigned int i = 0; i < 4; ++i)
         {
            switch (i)
            {
               case 0: cur = res->beta(1)->beta(2); break;
               case 1: cur = res->beta(0)->beta(2); break;
               case 2: cur = res->beta(2); break;
               case 3: cur = res; break;
            }

            maptcell_it = TC.find(it->neighbor(i));
            if (maptcell_it != TC.end())
            {
               switch (it->neighbor(i)->index(it))
               {
                  case 0: neighbor =
                        maptcell_it->second->beta(1)->beta(2);
                     break;
                  case 1: neighbor =
                        maptcell_it->second->beta(0)->beta(2);
                     break;
                  case 2: neighbor =
                        maptcell_it->second->beta(2); break;
                  case 3: neighbor = maptcell_it->second; break;
               }
               while (Map::vertex_attribute(neighbor) !=
		      Map::vertex_attribute(cur->other_extremity()) )
		 neighbor = neighbor->beta(1);
               amap.template topo_sew<3>(cur, neighbor);
               if (!neighbor->beta(2)->is_free(3) &&
		   !neighbor->beta(0)->beta(2)->is_free(3) &&
		   !neighbor->beta(1)->beta(2)->is_free(3))
                  TC.erase(maptcell_it);
            }
         }
         if (res->is_free(3) ||
	     res->beta(2)->is_free(3) ||
	     res->beta(0)->beta(2)->is_free(3) ||
	     res->beta(1)->beta(2)->is_free(3))
            TC[it] = res;
      }
   }
   return dart;
}

/// Struct for comparision of two Halfedge_handles (required for stl::map)
namespace internal {
template < class Polyhedron>
struct Hedge_cmp
{
   /// Halfedge handle.
   typedef typename Polyhedron::Halfedge_const_handle  Halfedge_handle;
   /// Operator() to apply the comparison.
   bool operator()(Halfedge_handle ah1,  Halfedge_handle ah2) const
   {
      return &*ah1 < &*ah2;
   }
};
}

/** Convert a given Polyhedron_3 into 3D combinatorial map.
 * @param amap the combinatorial map where Polyhedron_3 will be converted.
 * @param apoly the Polyhedron.
 * @return A dart created during the convertion.
 */
template< class Map, class Polyhedron >
typename Map::Dart_handle import_from_polyhedron(Map& amap, 
						 const Polyhedron &apoly)
{
   typedef typename Polyhedron::Halfedge_const_handle  Halfedge_handle;
   typedef typename Polyhedron::Facet_const_iterator   Facet_iterator;
   typedef typename Polyhedron::Halfedge_around_facet_const_circulator
   HF_circulator;

   typedef std::map < Halfedge_handle, typename Map::Dart_handle> 
   Halfedge_handle_map;
   typedef typename Halfedge_handle_map::iterator itmap_hds;
   Halfedge_handle_map TC;

   itmap_hds it;
   typename Map::Dart_handle d = NULL, prev = NULL,
     firstFacet = NULL, firstAll = NULL;

   // First traversal to build the darts and link them.
   for (Facet_iterator i = apoly.facets_begin(); i != apoly.facets_end(); ++i)
   {
      HF_circulator j = i->facet_begin();
      prev = NULL;
      do
      {
         d = amap.create_dart();
         TC[j] = d;
               
         if (prev != NULL) amap.template link_beta<1>(prev, d);
         else firstFacet = d;
         it = TC.find(j->opposite());
         if (it != TC.end())
	   amap.link_beta(d, it->second, 2);
         prev = d;
      }
      while (++j != i->facet_begin());
      amap.template link_beta<1>(prev, firstFacet);
      if (firstAll == NULL) firstAll = firstFacet;
   }

   // Second traversal to update the geometry.
   // We run one again through the facets of the HDS.
   for (Facet_iterator i = apoly.facets_begin(); i != apoly.facets_end(); ++i)
   {
      HF_circulator j = i->facet_begin();
      do
      {
	d = TC[j]; // Get the dart associated to the Halfedge
	if (Map::vertex_attribute(d) == NULL)
	  {	    
	    amap.set_vertex_attribute(d,
               amap.create_vertex_attribute(j->opposite()->vertex()->point()));
	  }
      }
      while (++j != i->facet_begin());
   }
   return firstAll;
}

template < class Map >
     //	   class Polyhedron=CGAL::Polyhedron_3<typename Map::Kernel> >
void
load_off(Map& amap, std::istream& in)
{
  File_header_OFF  m_file_header;
  File_scanner_OFF scanner( in, m_file_header.verbose());
  if ( ! in) return;
  m_file_header = scanner;  // Remember file header after return.

  Linear_cell_complex_incremental_builder_3<Map> B( amap);
  B.begin_surface( scanner.size_of_vertices(),
                   scanner.size_of_facets(),
                   scanner.size_of_halfedges());

  typedef typename Map::Point Point;

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

/** Convert a Polyhedron_3 read into a flux into 3D combinatorial map.
 * @param amap the combinatorial map where Polyhedron_3 will be converted.
 * @param ais the istream where read the Polyhedron_3.
 * @return A dart created during the convertion.
 */
template < class Map >
	   //	   class Polyhedron=CGAL::Polyhedron_3<typename Map::Kernel> >
typename Map::Dart_handle
import_from_polyhedron_flux(Map& amap, std::istream& ais)
{
   if (!ais.good())
   {
      std::cout << "Error reading flux." << std::endl;
      return NULL;
   }
   CGAL::Polyhedron_3<typename Map::Traits> P;
   ais >> P;
   return import_from_polyhedron<Map, 
    CGAL::Polyhedron_3<typename Map::Traits> > (amap, P);
}

} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_WITH_POINTS_CONSTRUCTORS_H //
// EOF //
