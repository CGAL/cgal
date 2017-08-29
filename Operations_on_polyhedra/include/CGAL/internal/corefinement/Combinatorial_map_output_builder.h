// Copyright (c) 2011 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot

#ifndef CGAL_INTERNAL_COREFINEMENT_COMBINATORIAL_MAP_OUTPUT_BUILDER_H
#define CGAL_INTERNAL_COREFINEMENT_COMBINATORIAL_MAP_OUTPUT_BUILDER_H

#include <CGAL/license/Polygon_mesh_processing.h>


#include <CGAL/internal/corefinement/Combinatorial_map_for_corefinement.h>
#include <CGAL/internal/corefinement/connected_components.h>
#include <CGAL/internal/corefinement/predicates.h>
#include <CGAL/internal/corefinement/utils.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Default.h>
#include <boost/shared_ptr.hpp>

/// \todo right now patches made of coplanar patches are handled by considering
///       each coplanar triangle while it should be sufficient to have one
///       procedure for all coplanar triangles of the patch together with
///       a sew-handling of the boundary of the patch.

namespace CGAL {

namespace internal_IOP {
//turn around the target vertex of dart to find a marked dart
template <class Combinatorial_map_3>
boost::optional<typename Combinatorial_map_3::Dart_handle>
next_marked_dart_around_target_vertex(
  Combinatorial_map_3& final_map,
  typename Combinatorial_map_3::Dart_handle dart,
  std::size_t mark_index)
{
  CGAL_precondition(final_map.is_marked(dart,mark_index));
  typename Combinatorial_map_3::Dart_handle next=final_map.beta(dart,1);
  while ( ! final_map.is_marked(next,mark_index) ) {
    if (final_map.is_free(next,2) )//we reach a boundary
      return  boost::optional<typename Combinatorial_map_3::Dart_handle>();
    next=final_map.beta(next,2,1);
  }
  if (next == dart) //no new dart have been found
    return  boost::optional<typename Combinatorial_map_3::Dart_handle>();
  CGAL_precondition(&final_map.template attribute<0>(final_map.beta(dart,1))->point() ==
                    &final_map.template attribute<0>(next)->point());
  return boost::optional<typename Combinatorial_map_3::Dart_handle> (next);
}

//turn around the target vertex of dart to find a marked dart
//with expected_target as target vertex
template <class Combinatorial_map_3>
typename Combinatorial_map_3::Dart_handle
get_next_marked_dart_around_target_vertex(
  Combinatorial_map_3& final_map,
  typename Combinatorial_map_3::Dart_handle dart,
  std::size_t mark_index)
{
  CGAL_precondition(final_map.is_marked(dart,mark_index));
  typename Combinatorial_map_3::Dart_handle next=final_map.beta(dart,1);
  while ( !final_map.is_marked(next,mark_index) ) {
    CGAL_assertion( !final_map.is_free(next,2) );
    next=final_map.beta(next,2,1);
    CGAL_assertion(next != dart);
  }
  CGAL_precondition(&final_map.template attribute<0>(final_map.beta(dart,1))->point() ==
                    &final_map.template attribute<0>(next)->point());
  return next;
}

//turn around the source vertex of dart to find a marked dart
//with expected_source as source vertex
template <class Combinatorial_map_3>
typename Combinatorial_map_3::Dart_handle
get_next_marked_dart_around_source_vertex(
  Combinatorial_map_3& final_map,
  typename Combinatorial_map_3::Dart_handle dart,
  std::size_t mark_index)
{
  CGAL_precondition(final_map.is_marked(dart,mark_index));
  typename Combinatorial_map_3::Dart_handle next=final_map.beta(dart,0);
  while ( ! final_map.is_marked(next,mark_index) ) {
    CGAL_assertion( !final_map.is_free(next,2) );
    next=final_map.beta(next,2,0);
    CGAL_assertion(next != dart);
  }
  CGAL_precondition(&final_map.template attribute<0>(dart)->point() ==
                    &final_map.template attribute<0>(final_map.beta(next,1))->point());
  return next;
}

//given two marked darts, this function links these two darts with beta<2>
//but in addition it follows the marked darts connected to the same vertex
//(there should be only one) to connect them all together
//( this function is a kind of zipper ;) )
template <class Combinatorial_map_3,class Nodes_vector>
void sew_2_marked_darts( Combinatorial_map_3& final_map,
                         typename Combinatorial_map_3::Dart_handle dart_1,
                         typename Combinatorial_map_3::Dart_handle dart_2,
                         std::size_t mark_index,
                         const Nodes_vector& nodes,
                         const std::pair<int,int>& indices,
                         const std::pair<bool,int>& polyline_info)
{
  CGAL_precondition( final_map.is_free(dart_1, 2) );
  CGAL_precondition( final_map.is_free(dart_2, 2) );
  CGAL_precondition( final_map.is_marked(dart_1,mark_index) );
  CGAL_precondition( final_map.is_marked(dart_2,mark_index) );
  CGAL_precondition( final_map.template attribute<0>(dart_1)->point() ==
                     final_map.template attribute<0>(final_map.beta(dart_2,1))->point() );
  CGAL_precondition( final_map.template attribute<0>(final_map.beta(dart_1,1))->point() ==
                     final_map.template attribute<0>(dart_2)->point() );

  int src_index = ( ( indices.first < indices.second) ==  polyline_info.first )
                  ? indices.second:indices.first;

  if ( final_map.template attribute<0>(dart_1)->point() != nodes[ src_index ] )
    std::swap(dart_1,dart_2);

  int nb_segs=polyline_info.second-1,k=1;

  do {
    CGAL_precondition( final_map.template is_sewable<2>(dart_1,dart_2) );
    final_map.template sew<2>(dart_1,dart_2);

    if (k==nb_segs) break;

    dart_1 = get_next_marked_dart_around_target_vertex(
               final_map, dart_1, mark_index);
    dart_2 = get_next_marked_dart_around_source_vertex(
               final_map, dart_2, mark_index);
  } while(++k);
}

//not_top and not_down are two darts from volumes that get merged with an
//existing other one because of a set of identical coplanar triangles.
//top and down is the dart of the volumes "replacing" that of not_top and
//not down respectively, The function is considering all triangles that
//are bounded by a cycle of marked edges. The volume not_top and not_down
//are part of are those that will disappear at the end of the main algorithm.
template <class Combinatorial_map_3>
void sew_3_marked_darts( Combinatorial_map_3& final_map,
                         typename Combinatorial_map_3::Dart_handle not_top,
                         typename Combinatorial_map_3::Dart_handle not_down,
                         typename Combinatorial_map_3::Dart_handle top,
                         typename Combinatorial_map_3::Dart_handle down,
                         std::size_t mark_index,
                         std::set<typename Combinatorial_map_3::Dart_handle>&
                            darts_to_remove
){
  typedef boost::optional<typename Combinatorial_map_3::Dart_handle>
    O_Dart_handle;

  if ( final_map.template attribute<3>(not_top)->info().is_empty ) {
    CGAL_assertion(final_map.template attribute<3>(not_down)->info().is_empty);
    return;
  }

  CGAL_assertion(!final_map.template attribute<3>(not_down)->info().is_empty);

  //merge attribute of the two volumes:
  internal_IOP::Volume_on_merge merge_attributes;
  merge_attributes(*final_map.template attribute<3>(top),
                   *final_map.template attribute<3>(not_top));
  merge_attributes(*final_map.template attribute<3>(down),
                   *final_map.template attribute<3>(not_down));

  //set volume attributes as empty to avoid double sew_3
  //of the same topological disk of triangles
  final_map.template attribute<3>(not_top)->info().is_empty=true;
  final_map.template attribute<3>(not_down)->info().is_empty=true;

  CGAL_precondition( final_map.is_marked(not_top,mark_index)
                     && final_map.is_marked(top,mark_index) );
  CGAL_precondition( final_map.is_marked(not_down,mark_index)
                     && final_map.is_marked(down,mark_index) );
  CGAL_precondition( final_map.template attribute<0>(not_top)->point() ==
                     final_map.template attribute<0>(final_map.beta(not_down,1))->point() );
  CGAL_precondition( final_map.template attribute<0>(final_map.beta(not_top,1))->point() ==
                     final_map.template attribute<0>(not_down)->point() );
  CGAL_precondition( final_map.template attribute<0>(not_top)->point() ==
                     final_map.template attribute<0>(top)->point() );
  CGAL_precondition( final_map.template attribute<0>(not_down)->point() ==
                     final_map.template attribute<0>(down)->point() );

  CGAL_assertion( final_map.beta(top,3)==down );

  //set to be removed the darts of the two no longer used volumes
  typename Combinatorial_map_3::Dart_handle start=not_top;
  do {
    CGAL_assertion(!final_map.is_free(not_top, 3));
    darts_to_remove.insert(not_top);
    darts_to_remove.insert(final_map.beta(not_top,1));
    darts_to_remove.insert(final_map.beta(not_top,1,1));
    darts_to_remove.insert(final_map.beta(not_top,3));
    darts_to_remove.insert(final_map.beta(not_top,3,1));
    darts_to_remove.insert(final_map.beta(not_top,3,1,1));
    O_Dart_handle current_1 =
      next_marked_dart_around_target_vertex(final_map, not_top, mark_index);
    CGAL_precondition(bool(current_1));
    not_top=*current_1;
  } while(not_top!=start);
}

} //end of internal_IOP namespace

namespace Corefinement
{

//import into the combinatorial map facets in the given range.
//they are supposed to be in the same connected component.
//two volume are created (each facets gives two opposite orientation 2-cell in the map)
template<class Polyhedron, class Map, class Face_iterator,
         class Non_special_edge_predicate,class Halfedge_to_dart_map_,
         class PolyhedronPointPMap >
typename Map::Dart_handle import_from_polyhedron_subset(
    Map& amap,
    Face_iterator faces_begin,
    Face_iterator faces_end,
    const Non_special_edge_predicate& is_non_special_edge,
    Halfedge_to_dart_map_& selected_hedge_to_dart,
    std::size_t mark_index,
    PolyhedronPointPMap ppmap
                                                       )
{
  typedef typename Polyhedron::Halfedge_const_handle  Halfedge_const_handle;
  typedef std::map< Halfedge_const_handle,
                    typename Map::Dart_handle,
                    internal_IOP::Compare_address<Polyhedron> >
    Halfedge_to_dart_map;

  Halfedge_to_dart_map hedge_to_dart;
  typename Map::Dart_handle first_dart = NULL;
  // First traversal to build the darts and link them.
  for (Face_iterator it_face = faces_begin; it_face != faces_end; ++it_face) {
    Halfedge_const_handle start=(*it_face)->halfedge();

    CGAL_precondition(start->next()!=start);

    Halfedge_const_handle current=start;
    typename Map::Dart_handle prev = NULL;
    typename Map::Dart_handle first_dart_of_face = NULL;
    do {
      typename Map::Dart_handle d = amap.create_dart();
      amap.template link_beta<3>(d,amap.create_dart()); //for opposite volume
      hedge_to_dart[current] = d;

      if (prev != NULL) {
        amap.template link_beta<1>(prev, d);
        amap.template link_beta<1>(amap.beta(d,3),amap.beta(prev,3));//for opposite volume
      } else {
        first_dart_of_face = d;
        if (first_dart==NULL) first_dart=d;
      }

      if ( is_non_special_edge (current) ) {
        if ( !current->is_border_edge() ) {
          CGAL_assertion(current != current->opposite());
          typename Halfedge_to_dart_map::iterator it =
            hedge_to_dart.find(current->opposite());
          if (it != hedge_to_dart.end()) {
            //link the opposites halfedges only when both
            //corresponding darts have been created
            amap.template link_beta<2>(d, it->second);
            amap.template link_beta<2>(amap.beta(d,3),
                                       amap.beta(it->second,3));//for opposite volume
          }
        }
      } else {
        typename Halfedge_to_dart_map_::iterator it_hedge_map=
          selected_hedge_to_dart.find(current);
        //all marked hedges are not the selected one for its polyline
        if ( it_hedge_map!=selected_hedge_to_dart.end() ) it_hedge_map->second=d;
        //darts d and d->beta(3) are special edges
        amap.mark(d,mark_index);
        amap.mark(amap.beta(d,3),mark_index);
      }
      prev = d;
      current=current->next();
    } while (current != start);
    amap.template link_beta<1>(prev, first_dart_of_face);
    amap.template link_beta<1>(amap.beta(first_dart_of_face,3),
                               amap.beta(prev,3));//for opposite volume
  }

  // Second traversal to update the geometry.
  // We run one again through the facets of the HDS.
  for (Face_iterator it_face = faces_begin; it_face != faces_end; ++it_face) {
    Halfedge_const_handle start=(*it_face)->halfedge();
    Halfedge_const_handle current=start;
    do {
      typename Map::Dart_handle d =
        hedge_to_dart[current]; // Get the dart associated to the Halfedge
      if (amap.template attribute<0>(d) == NULL) {
        amap.template set_attribute<0>(d,
          amap.template create_attribute<0>(
            get(ppmap, current->opposite()->vertex()))
        );
      }
      current=current->next();
    } while (current != start);
  }

  return first_dart;
}

template <class Polyhedron,
         class Kernel_=Default,
         class PolyhedronPointPMap_=Default >
class Combinatorial_map_output_builder
{
//Default typedefs
  typedef typename  Default::Get<
    PolyhedronPointPMap_,
    Default_polyhedron_ppmap<Polyhedron> >::type     PolyhedronPointPMap;
  typedef typename Default::Get<
    Kernel_,
    typename Kernel_traits<
      typename boost::property_traits<PolyhedronPointPMap>::value_type >::Kernel
  >::type Kernel;
//Polyhedron typedefs
  typedef typename Polyhedron::Halfedge_handle               Halfedge_handle;
  typedef typename Polyhedron::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Polyhedron::Vertex_handle                   Vertex_handle;
//Other typedefs
  typedef internal_IOP::Compare_unik_address<Polyhedron> Cmp_unik_ad;
  typedef std::map< std::pair<int,int>,
                    std::pair< std::map<Polyhedron*,Halfedge_handle>,
                               std::pair<bool,int> > > An_edge_per_polyline_map;
  typedef std::map<int,Vertex_handle> Node_to_polyhedron_vertex_map;
  typedef std::map<Polyhedron*, Node_to_polyhedron_vertex_map> Poly_to_map_node;
public:
//Combinatorial map typedefs
  typedef internal_IOP::Item_with_points_and_volume_info<Kernel,Polyhedron> Items;
  typedef CGAL::Combinatorial_map<3,Items> Combinatorial_map_3;
  typedef typename Combinatorial_map_3::Dart_handle Dart_handle;
  typedef internal_IOP::Volume_info<Polyhedron> Volume_info;
private:
//Data members
  PolyhedronPointPMap ppmap;
  boost::shared_ptr<Combinatorial_map_3>                 final_map_ptr;
private:
  template <class Halfedge_to_dart_map>
  inline Dart_handle get_associated_dart(Halfedge_handle hedge,
                                         Halfedge_to_dart_map& selected_hedge_to_dart) {
    typename Halfedge_to_dart_map::iterator it_saved_dart=
      selected_hedge_to_dart.find(hedge);
    CGAL_assertion(it_saved_dart!=selected_hedge_to_dart.end());
    return it_saved_dart->second;
  }

  //first_hedge defines four volumes, second_hedge only two
  //first_poly is not needed as inside/outside volume is update during the merge
  //of the sew. Only second_poly is needed
  template <class Nodes_vector,class Border_halfedges_map,class Halfedge_to_dart_map>
  void sew_2_three_volumes_case(  Halfedge_handle first_hedge,
                                  Halfedge_handle second_hedge,
                                  const std::pair<int,int>& indices,
                                  const Nodes_vector& nodes,
                                  Border_halfedges_map& border_halfedges,
                                  Halfedge_to_dart_map& selected_hedge_to_dart,
                                  Polyhedron* /*first_poly*/,
                                  Polyhedron* second_poly,
                                  std::size_t mark_index,
                                  std::set<Dart_handle>& darts_to_remove,
                                  const std::pair<bool,int>& polyline_info) {
    bool took_opposite=second_hedge->is_border();
    if (took_opposite) second_hedge=second_hedge->opposite();

    Vertex_handle P1=first_hedge->opposite()->next()->vertex();
    Vertex_handle P2=first_hedge->next()->vertex();
    //    when looking from the side of indices.second, the interior of the first polyhedron is described
    //    by turning counterclockwise from P1 to P2

    Vertex_handle Q = second_hedge->next()->vertex();

    //check if the third point of each triangular face is an original point (stay -1)
    //or a intersection point (in that case we need the index of the corresponding node to
    //have the exact value of the point)
    int index_p1=node_index_of_incident_vertex(first_hedge->opposite()->next(),
                 border_halfedges);
    int index_p2=node_index_of_incident_vertex(first_hedge->next(),
                 border_halfedges);
    int index_q =node_index_of_incident_vertex(second_hedge->next(),
                 border_halfedges);

    //Recover the dart that will be the start point of the different sewing
    //  dof_X_outside = dart of face of, meaning the triangle containing the
    //  point X and part of the volume outside of the corresponding polyhedron
    //-----first polyhedron
    Dart_handle dof_P1_outside = get_associated_dart(first_hedge->opposite(),
                                 selected_hedge_to_dart);
    Dart_handle dof_P2_outside = get_associated_dart(first_hedge,
                                 selected_hedge_to_dart);
    //-----second polyhedron
    Dart_handle dof_Q_outside = get_associated_dart(second_hedge,
                                selected_hedge_to_dart);

    if (index_p1!=-1 && index_p1==index_q) {
      Dart_handle top=final_map().beta(dof_P1_outside,3),
                  not_top=took_opposite?final_map().beta(dof_Q_outside,3):dof_Q_outside;
      Dart_handle down=dof_P1_outside,
                  not_down=took_opposite?dof_Q_outside:final_map().beta(dof_Q_outside,3);

      if ( final_map().template attribute<3>(top)->info().is_empty )
        std::swap(not_top,top);
      if ( final_map().template attribute<3>(down)->info().is_empty )
        std::swap(not_down,down);
      CGAL_assertion( !final_map().template attribute<3>(top)->info().is_empty );
      CGAL_assertion( !final_map().template attribute<3>(down)->info().is_empty );

      //P1P2 or QP2
      sew_2_marked_darts( final_map(), top, final_map().beta(dof_P2_outside,3),
                          mark_index, nodes, indices, polyline_info);
      //P2Q or P2P1
      sew_2_marked_darts( final_map(),dof_P2_outside, down,
                          mark_index, nodes, indices, polyline_info);
      sew_3_marked_darts( final_map(),not_top, not_down, top, down, mark_index,
                          darts_to_remove);

      return;
    }

    if (index_p2!=-1 && index_p2==index_q) {
      Dart_handle top=final_map().beta(dof_P2_outside,3),
                  not_top=took_opposite?dof_Q_outside:final_map().beta(dof_Q_outside,3);
      Dart_handle down=dof_P2_outside,
                  not_down=took_opposite?final_map().beta(dof_Q_outside,3):dof_Q_outside;

      if ( final_map().template attribute<3>(top)->info().is_empty )
        std::swap(not_top,top);
      if ( final_map().template attribute<3>(down)->info().is_empty )
        std::swap(not_down,down);
      CGAL_assertion( !final_map().template attribute<3>(top)->info().is_empty );
      CGAL_assertion( !final_map().template attribute<3>(down)->info().is_empty );

      //P1Q or P1P2
      sew_2_marked_darts( final_map(), final_map().beta(dof_P1_outside,3), top,
                          mark_index, nodes, indices, polyline_info);
      //QP1 or P2P1
      sew_2_marked_darts( final_map(), down, dof_P1_outside,
                          mark_index, nodes, indices, polyline_info);
      sew_3_marked_darts( final_map(), not_top, not_down, top, down, mark_index,
                          darts_to_remove);

      return;
    }

    bool Q_is_between_P1P2 = OOP::sorted_around_edge_filtered(indices.first,
                                                              indices.second,
                                                              index_p1,index_p2,
                                                              index_q,P1,P2,Q,
                                                              nodes,ppmap);

    if (Q_is_between_P1P2) {
      // poly_first  - poly_second            = took_opposite?P1Q:QP2
      // poly_second - poly_first             = {0}
      // poly_first \cap poly_second          = took_opposite?QP2:P1Q
      // opposite( poly_first U poly_second ) = P2P1
      sew_2_marked_darts( final_map(),
                          final_map().beta(dof_P1_outside,3),
                          took_opposite?dof_Q_outside:final_map().beta(dof_Q_outside,3),
                          mark_index, nodes, indices, polyline_info); //P1Q
      sew_2_marked_darts( final_map(),
                          took_opposite?final_map().beta(dof_Q_outside,3):dof_Q_outside,
                          final_map().beta(dof_P2_outside,3),mark_index, nodes, indices,
                          polyline_info); //QP2
      sew_2_marked_darts( final_map(),
                          dof_P2_outside, dof_P1_outside, mark_index,
                          nodes, indices, polyline_info); //P2P1
      final_map().template attribute<3>(dof_P1_outside)->info().outside
        .insert(second_poly); //update P2P1 outside poly
    } else {
      // poly_first  - poly_second            = P1P2
      // poly_second - poly_first             = took_opposite?QP1:P2Q
      // poly_first \cap poly_second          = {0}
      // opposite( poly_first U poly_second ) = took_opposite?P2Q:QP1
      sew_2_marked_darts( final_map(),
                          dof_P2_outside,
                          took_opposite?dof_Q_outside:final_map().beta(dof_Q_outside,3),
                          mark_index, nodes, indices, polyline_info); //P2Q
      sew_2_marked_darts( final_map(),
                          took_opposite?final_map().beta(dof_Q_outside,3):dof_Q_outside,
                          dof_P1_outside,mark_index, nodes,
                          indices, polyline_info); //QP1
      sew_2_marked_darts( final_map(),
                          final_map().beta(dof_P1_outside,3), final_map().beta(dof_P2_outside,3),
                          mark_index, nodes, indices, polyline_info); //P1P2
      final_map().template attribute<3>(final_map().beta(dof_P1_outside,3))->info().outside
        .insert(second_poly); //update P1P2 outside poly
    }
  }

//first_hedge defines two volumes, second_hedge only two
  template <class Halfedge_to_dart_map,
            class Border_halfedges_map,
            class Nodes_vector>
  void sew_2_two_volumes_case(  Halfedge_handle first_hedge,
                                Halfedge_handle second_hedge,
                                Border_halfedges_map& border_halfedges,
                                Halfedge_to_dart_map& selected_hedge_to_dart,
                                std::size_t mark_index,
                                std::set<Dart_handle>& darts_to_remove,
                                const Nodes_vector& nodes,
                                const std::pair<int,int>& indices,
                                const std::pair<bool,int>& polyline_info)
  {
    bool first_took_opposite=first_hedge->is_border();
    if (first_took_opposite) first_hedge=first_hedge->opposite();
    bool second_took_opposite=second_hedge->is_border();
    if (second_took_opposite) second_hedge=second_hedge->opposite();

    //-----first polyhedron
    Dart_handle dof_P_outside = get_associated_dart(first_hedge,
                                selected_hedge_to_dart);
    //-----second polyhedron
    Dart_handle dof_Q_outside = get_associated_dart(second_hedge,
                                selected_hedge_to_dart);




    int index_p =node_index_of_incident_vertex(first_hedge->next(),
                 border_halfedges);
    int index_q =node_index_of_incident_vertex(second_hedge->next(),
                 border_halfedges);

    if (index_p!=-1 && index_q!=-1 && index_p==index_q) {
      Dart_handle top=dof_P_outside, not_top=final_map().beta(dof_Q_outside,3);
      Dart_handle down=final_map().beta(dof_P_outside,3), not_down=dof_Q_outside;

      if (first_took_opposite==second_took_opposite) {
        top=final_map().beta(dof_P_outside,3);
        not_top=final_map().beta(dof_Q_outside,3);
        down=dof_P_outside;
        not_down=dof_Q_outside;
      }

      if ( final_map().template attribute<3>(top)->info().is_empty )
        std::swap(not_top,top);
      if ( final_map().template attribute<3>(down)->info().is_empty )
        std::swap(not_down,down);
      CGAL_assertion( !final_map().template attribute<3>(top)->info().is_empty );
      CGAL_assertion( !final_map().template attribute<3>(down)->info().is_empty );

      sew_3_marked_darts( final_map(),not_top,not_down,top,down,mark_index,
                          darts_to_remove);

      return;
    }



    //since the edge is shared, the inside of each polyhedron must be on opposite orientation halfedges
    if (first_took_opposite==second_took_opposite) {
      //sew out with in
      sew_2_marked_darts( final_map(),final_map().beta(dof_P_outside,3), dof_Q_outside,
                          mark_index, nodes, indices, polyline_info); //PQ
      sew_2_marked_darts( final_map(),final_map().beta(dof_Q_outside,3), dof_P_outside,
                          mark_index, nodes, indices, polyline_info); //QP
    } else {
      //sew in with in
      sew_2_marked_darts( final_map(), dof_P_outside, dof_Q_outside,
                          mark_index, nodes, indices, polyline_info); //PQ
      sew_2_marked_darts( final_map(),final_map().beta(dof_Q_outside,3),
                          final_map().beta(dof_P_outside,3), mark_index, nodes, indices,
                          polyline_info); //QP
    }
  }

//4 volume case with 2 identical volume
//Q2 is supposed to be identical to P2
  template <class Nodes_vector,class Halfedge_to_dart_map>
  void sew_2_four_volumes_case_1(  Halfedge_handle first_hedge,
                                   Halfedge_handle second_hedge,
                                   const std::pair<int,int>& indices,
                                   const Nodes_vector& nodes,
                                   int index_p1, int index_p2, int index_q1,
                                   Halfedge_to_dart_map& selected_hedge_to_dart,
                                   std::size_t mark_index,
                                   std::set<Dart_handle>& darts_to_remove,
                                   const std::pair<bool,int>& polyline_info,
                                   bool swap_in_out_Q=false) {
    Vertex_handle P1=first_hedge->opposite()->next()->vertex();
    Vertex_handle P2=first_hedge->next()->vertex();
    // when looking from the side of indices.second,
    // the interior of the first polyhedron is described
    // by turning counterclockwise from P1 to P2
    Vertex_handle Q1=second_hedge->opposite()->next()->vertex();
    // Vertex_handle Q2=second_hedge->next()->vertex();
    bool Q1_is_between_P1P2 =
      OOP::sorted_around_edge_filtered( indices.first, indices.second, index_p1,
                                        index_p2, index_q1, P1, P2, Q1, nodes, ppmap);


    //Recover the dart that will be the start point of the different sewing
    //  dof_X_outside = dart of face of, meaning the triangle containing the
    //  point X and part of the volume outside of the corresponding polyhedron
    //-----first polyhedron
    Dart_handle dof_P1_outside = get_associated_dart(first_hedge->opposite(),
                                 selected_hedge_to_dart);
    Dart_handle dof_P2_outside = get_associated_dart(first_hedge,
                                 selected_hedge_to_dart);
    //-----second polyhedron
    Dart_handle dof_Q1_outside = get_associated_dart(second_hedge->opposite(),
                                 selected_hedge_to_dart);
    Dart_handle dof_Q2_outside = get_associated_dart(second_hedge,
                                 selected_hedge_to_dart);

    if( swap_in_out_Q ) {
      dof_Q1_outside=final_map().beta(dof_Q1_outside,3);
      dof_Q2_outside=final_map().beta(dof_Q2_outside,3);
    }

    if (Q1_is_between_P1P2) {
      Dart_handle top=final_map().beta(dof_Q2_outside,3), not_top=final_map().beta(dof_P2_outside,3);
      Dart_handle down=dof_Q2_outside, not_down=dof_P2_outside;
      if ( final_map().template attribute<3>(top)->info().is_empty )
        std::swap(not_top,top);
      if ( final_map().template attribute<3>(down)->info().is_empty )
        std::swap(not_down,down);
      CGAL_assertion( !final_map().template attribute<3>(top)->info().is_empty );
      CGAL_assertion( !final_map().template attribute<3>(down)->info().is_empty );

      // poly_first  - poly_second            = P1Q1
      // poly_second - poly_first             = {0}
      // poly_first \cap poly_second          = Q1P2 or Q1Q2
      // opposite( poly_first U poly_second ) = Q2P1 or P2P1
      sew_2_marked_darts( final_map(), final_map().beta(dof_P1_outside,3), dof_Q1_outside,
                          mark_index, nodes, indices, polyline_info); //P1Q1
      sew_2_marked_darts( final_map(),final_map().beta(dof_Q1_outside,3),
                          top, mark_index, nodes, indices,
                          polyline_info); //Q1P2 or Q1Q2
      sew_2_marked_darts( final_map(),down,
                          dof_P1_outside, mark_index, nodes, indices,
                          polyline_info); //Q2P1 or P2P1
      sew_3_marked_darts( final_map(),not_top,not_down,top,down,mark_index,
                          darts_to_remove);

    } else {
      Dart_handle top=final_map().beta(dof_Q2_outside,3), not_top=final_map().beta(dof_P2_outside,3);
      Dart_handle down=dof_Q2_outside, not_down=dof_P2_outside;
      if ( final_map().template attribute<3>(top)->info().is_empty )
        std::swap(not_top,top);
      if ( final_map().template attribute<3>(down)->info().is_empty )
        std::swap(not_down,down);
      CGAL_assertion( !final_map().template attribute<3>(top)->info().is_empty );
      CGAL_assertion( !final_map().template attribute<3>(down)->info().is_empty );

      // poly_first  - poly_second            = {0}
      // poly_second - poly_first             = Q1P1
      // poly_first \cap poly_second          = P1P2 or P1Q2
      // opposite( poly_first U poly_second ) = Q2Q1 or P2Q1
      sew_2_marked_darts( final_map(),final_map().beta(dof_Q1_outside,3), dof_P1_outside,
                          mark_index, nodes, indices, polyline_info); //Q1P1
      sew_2_marked_darts( final_map(),final_map().beta(dof_P1_outside,3),
                          top, mark_index, nodes, indices,
                          polyline_info); //P1P2 or P1Q2
      sew_2_marked_darts( final_map(),down,
                          dof_Q1_outside, mark_index, nodes, indices,
                          polyline_info); //Q2Q1 or P2Q1
      sew_3_marked_darts( final_map(),not_top,not_down,top,down,mark_index,
                          darts_to_remove);
    }
  }

  template <class Nodes_vector,class Halfedge_to_dart_map>
  bool coplanar_triangles_case_handled(
    Halfedge_handle first_hedge, Halfedge_handle second_hedge,
    const std::pair<int,int>& indices,
    const Nodes_vector& nodes,
    int index_p1, int index_p2,
    int index_q1,int index_q2,
    Halfedge_to_dart_map& selected_hedge_to_dart,
    std::size_t mark_index,
    std::set<Dart_handle>& darts_to_remove,
    const std::pair<bool,int>& polyline_info
  ){
    if( index_p1!=-1 ) {
      if (index_p1==index_q1) {
        if(index_p2!=-1) {
          CGAL_assertion(index_p2!=index_q1);
          if(index_p2==index_q2) {
            //-----first polyhedron
            Dart_handle dof_P1_outside = get_associated_dart(first_hedge->opposite(),
                                         selected_hedge_to_dart);
            Dart_handle dof_P2_outside = get_associated_dart(first_hedge,
                                         selected_hedge_to_dart);
            //-----second polyhedron
            Dart_handle dof_Q1_outside = get_associated_dart(second_hedge->opposite(),
                                         selected_hedge_to_dart);
            Dart_handle dof_Q2_outside = get_associated_dart(second_hedge,
                                         selected_hedge_to_dart);

            Dart_handle top_1=final_map().beta(dof_P1_outside,3),
                        not_top_1=final_map().beta(dof_Q1_outside,3);
            Dart_handle top_2=final_map().beta(dof_P2_outside,3),
                        not_top_2=final_map().beta(dof_Q2_outside,3);
            Dart_handle down_1=dof_P1_outside, not_down_1=dof_Q1_outside;
            Dart_handle down_2=dof_P2_outside, not_down_2=dof_Q2_outside;
            if ( final_map().template attribute<3>(top_1)->info().is_empty )
              std::swap(top_1,not_top_1);
            if ( final_map().template attribute<3>(top_2)->info().is_empty )
              std::swap(top_2,not_top_2);
            if ( final_map().template attribute<3>(down_1)->info().is_empty )
              std::swap(down_1,not_down_1);
            if ( final_map().template attribute<3>(down_2)->info().is_empty )
              std::swap(down_2,not_down_2);
            CGAL_assertion( !final_map().template attribute<3>(top_1)->info().is_empty );
            CGAL_assertion( !final_map().template attribute<3>(top_2)->info().is_empty );
            CGAL_assertion( !final_map().template attribute<3>(down_1)->info().is_empty );
            CGAL_assertion( !final_map().template attribute<3>(down_2)->info().is_empty );

            // poly_first  - poly_second            = {0}
            // poly_second - poly_first             = {0}
            // poly_first \cap poly_second          = P1P2 or Q1Q2 or P1Q1 or Q1P2
            // opposite( poly_first U poly_second ) = P2P1 or Q2Q1 or P2Q1 or Q2P1
            sew_2_marked_darts( final_map(),top_1,  top_2,mark_index, nodes,
                                indices, polyline_info); //P1P2 or Q1Q2 or P1Q1 or Q1P2
            sew_2_marked_darts( final_map(),down_2,  down_1, mark_index, nodes,
                                indices, polyline_info); //P2P1 or Q2Q1 or P2Q1 or Q2P1
            sew_3_marked_darts( final_map(),not_top_1, not_down_1, top_1,
                                down_1,mark_index, darts_to_remove);
            sew_3_marked_darts( final_map(),not_top_2, not_down_2, top_2,
                                down_2,mark_index, darts_to_remove);
            return true;
          }
        }
        sew_2_four_volumes_case_1(first_hedge->opposite(),second_hedge->opposite(),
                                  std::make_pair(indices.second,indices.first),
                                  nodes,index_p2,index_p1,index_q2,
                                  selected_hedge_to_dart,mark_index,
                                  darts_to_remove,polyline_info);
        return true;
      }
      if (index_p1==index_q2) {
        if(index_p2!=-1) {
          CGAL_assertion(index_p2!=index_q2);
          if(index_p2==index_q1) {
            //-----first polyhedron
            Dart_handle dof_P1_outside = get_associated_dart(first_hedge->opposite(),
                                         selected_hedge_to_dart);
            Dart_handle dof_P2_outside = get_associated_dart(first_hedge,
                                         selected_hedge_to_dart);
            //-----second polyhedron
            Dart_handle dof_Q1_outside = get_associated_dart(second_hedge->opposite(),
                                         selected_hedge_to_dart);
            Dart_handle dof_Q2_outside = get_associated_dart(second_hedge,
                                         selected_hedge_to_dart);

            Dart_handle top_1=final_map().beta(dof_P1_outside,3), not_top_1=dof_Q2_outside;
            Dart_handle top_2=final_map().beta(dof_P2_outside,3), not_top_2=dof_Q1_outside;
            Dart_handle down_1=dof_P1_outside, not_down_1=final_map().beta(dof_Q2_outside,3);
            Dart_handle down_2=dof_P2_outside, not_down_2=final_map().beta(dof_Q1_outside,3);
            if ( final_map().template attribute<3>(top_1)->info().is_empty )
              std::swap(top_1,not_top_1);
            if ( final_map().template attribute<3>(top_2)->info().is_empty )
              std::swap(top_2,not_top_2);
            if ( final_map().template attribute<3>(down_1)->info().is_empty )
              std::swap(down_1,not_down_1);
            if ( final_map().template attribute<3>(down_2)->info().is_empty )
              std::swap(down_2,not_down_2);
            CGAL_assertion( !final_map().template attribute<3>(top_1)->info().is_empty );
            CGAL_assertion( !final_map().template attribute<3>(top_2)->info().is_empty );
            CGAL_assertion( !final_map().template attribute<3>(down_1)->info().is_empty );
            CGAL_assertion( !final_map().template attribute<3>(down_2)->info().is_empty );

            // poly_first  - poly_second            = P1P2 or P1Q1 or Q2P2 or Q2Q1
            // poly_second - poly_first             = Q1Q2 or Q1P1 or P2P1 or P2Q2
            // poly_first \cap poly_second          = {0}
            // opposite( poly_first U poly_second ) = all space
            sew_2_marked_darts( final_map(),top_1,  top_2, mark_index, nodes,
                                indices, polyline_info); //P1P2 or Q1Q2 or P1Q1 or Q1P2
            sew_2_marked_darts( final_map(),down_2, down_1, mark_index, nodes,
                                indices, polyline_info); //P2P1 or Q2Q1 or P2Q1 or Q2P1
            sew_3_marked_darts( final_map(),not_top_1, not_down_1, top_1,down_1,mark_index,
                                darts_to_remove);
            sew_3_marked_darts( final_map(),not_top_2, not_down_2, top_2,down_2,mark_index,
                                darts_to_remove);
            return true;
          }
        }
        sew_2_four_volumes_case_1(first_hedge->opposite(), second_hedge,
                                  std::make_pair(indices.second,indices.first),
                                  nodes,index_p2, index_p1, index_q1,
                                  selected_hedge_to_dart, mark_index,
                                  darts_to_remove, polyline_info, true);
        return true;
      }
    }

    if(index_p2!=-1) {
      if (index_p2==index_q1) {
        sew_2_four_volumes_case_1(first_hedge,second_hedge->opposite(),
                                  indices,nodes, index_p1,index_p2,index_q2,
                                  selected_hedge_to_dart,mark_index,
                                  darts_to_remove, polyline_info,true);
        return true;
      }
      if(index_p2==index_q2) {
        sew_2_four_volumes_case_1(first_hedge, second_hedge, indices, nodes,
                                  index_p1, index_p2, index_q1,
                                  selected_hedge_to_dart, mark_index,
                                  darts_to_remove, polyline_info);
        return true;
      }
    }


    return false;
  }

public:

  Combinatorial_map_3&                 final_map() {
    return *final_map_ptr;
  }
  Combinatorial_map_3& combinatorial_map() {
    return *final_map_ptr;
  }
  boost::shared_ptr<Combinatorial_map_3> combinatorial_map_shared_ptr() {
    return final_map_ptr;
  }

  Combinatorial_map_output_builder(boost::shared_ptr<Combinatorial_map_3> map_ptr)
    : final_map_ptr(map_ptr) {}

  Combinatorial_map_output_builder(
    PolyhedronPointPMap point_pmap = PolyhedronPointPMap()
  ) : ppmap(point_pmap)
    , final_map_ptr(new  Combinatorial_map_3())
  {}

  void input_have_coplanar_facets(){}

  template <class Nodes_vector, class Bitset>
  void operator()(
    const std::map<Halfedge_const_handle,
                   std::pair<int,int>,Cmp_unik_ad >& border_halfedges,
    const Nodes_vector& nodes,
    const An_edge_per_polyline_map& an_edge_per_polyline,
    const Bitset& /*is_node_of_degree_one*/,
    const Poly_to_map_node& polyhedron_to_map_node_to_polyhedron_vertex) {
    //4) create one output polyhedron per connected component of polyhedron,
    //   connected by an edge which is not an intersection edge
    //5) import into a Combinatorial map
#ifdef CGAL_COREFINEMENT_DEBUG
    std::cout << "Nb marked edges " << border_halfedges.size() << std::endl;
//    for (typename Border_halfedges_map::iterator it=border_halfedges.begin();it!=border_halfedges.end();++it)
//      std::cout << it->first->get(ppmap,opposite()->vertex()) << " " << get(ppmap,it->first->vertex()) << " is constrained " << std::endl;
    std::cout << "Nb polylines " << an_edge_per_polyline.size() << std::endl;
#endif

    internal_IOP::Non_intersection_halfedge<Polyhedron> criterium(border_halfedges);

    std::size_t mark_index=
      final_map().get_new_mark(); //mark used to tag dart that are on an intersection

    //define a map that will contain the correspondance between selected halfedges of the boundary and
    //their corresponding Dart_handle in the cmap.
    typedef std::map< Halfedge_const_handle,
                      Dart_handle,
                      internal_IOP::Compare_address<Polyhedron> >
      Halfedge_to_dart_map;
    Halfedge_to_dart_map selected_hedge_to_dart;
    for (typename An_edge_per_polyline_map::const_iterator it=
           an_edge_per_polyline.begin(); it!=an_edge_per_polyline.end(); ++it) {
      CGAL_assertion(it->second.first.size()==2);
      //orientation of faces around the edge (to be sure we can do it)
      Halfedge_handle first_hedge=it->second.first.begin()->second;
      Halfedge_handle second_hedge=boost::next(it->second.first.begin())->second;

      if (!first_hedge->is_border())
        selected_hedge_to_dart
          .insert(std::make_pair(first_hedge,Dart_handle(NULL)));
      if (!first_hedge->opposite()->is_border())
        selected_hedge_to_dart
          .insert(std::make_pair(first_hedge->opposite(),Dart_handle(NULL)));
      if (!second_hedge->is_border())
        selected_hedge_to_dart
          .insert(std::make_pair(second_hedge,Dart_handle(NULL)));
      if (!second_hedge->opposite()->is_border())
        selected_hedge_to_dart
          .insert(std::make_pair(second_hedge->opposite(),Dart_handle(NULL)));
    }

#ifdef CGAL_COREFINEMENT_DEBUG
    int polynb=0;
#endif
    for (typename Poly_to_map_node::const_iterator
         it=polyhedron_to_map_node_to_polyhedron_vertex.begin();
         it!=polyhedron_to_map_node_to_polyhedron_vertex.end();
         ++it
        ) {
      typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
      typedef ::CGAL::Union_find<Facet_const_handle> UF;
      typedef typename UF::handle UF_handle;
      typedef std::map< Facet_const_handle,
                        std::list<Facet_const_handle>,
                        internal::corefinement::Compare_handle_ptr<Polyhedron> > Result;
      typedef std::map< Facet_const_handle,
                        UF_handle,
                        internal::corefinement::Compare_handle_ptr<Polyhedron> >
          Facet_to_handle_map;

      UF uf;
      Facet_to_handle_map map_f2h;
      Result result;
      Polyhedron* current_poly=it->first;

#ifdef CGAL_COREFINEMENT_DEBUG
      std::cout << "writing poly debug"<< std::endl;
      std::stringstream ss;
      ss << "output_debug-" << ++polynb << ".off";
      std::ofstream output_debug(ss.str().c_str());
      output_debug << *current_poly;
#endif

      internal::corefinement::extract_connected_components(*(static_cast<Polyhedron const *> (current_poly) ),
                                                           criterium,uf,map_f2h,result);

      //add each connected component in the map with 2 volumes per component.
      for (typename Result::iterator it_res=result.begin(); it_res!=result.end();
           ++it_res) {
        //create in the final Cmap a 2D component containing faces of a connected component
        //(twice: one with same orientation and one with the opposite orientation to model the other volume)
        Dart_handle d =
          import_from_polyhedron_subset<Polyhedron>(
            final_map(), it_res->second.begin(), it_res->second.end(),
            criterium, selected_hedge_to_dart, mark_index,ppmap);
        //set an attribute to one volume represented by this component
        //to indicate a part outside of the polyhedron current_poly
        typename Combinatorial_map_3::template Attribute_range<3>::type::iterator
          attrib=final_map().template create_attribute<3>();
        attrib->info().outside.insert(current_poly);
        final_map().template set_attribute<3>(d,attrib);
        //set the attribute for the opposite volume: represent a part inside current_poly
        attrib=final_map().template create_attribute<3>();
        attrib->info().inside.insert(current_poly);
        final_map().template set_attribute<3>(final_map().beta(d, 3),attrib);

#ifdef CGAL_COREFINEMENT_DEBUG
        final_map().display_characteristics(std::cout);
        std::cout << std::endl;
#endif
      }
    }
#ifndef NDEBUG
    for(typename Halfedge_to_dart_map::iterator it=selected_hedge_to_dart.begin();
        it!=selected_hedge_to_dart.end(); ++it)
      CGAL_assertion(it->second!=Dart_handle(NULL));
#endif

    CGAL_assertion(final_map().is_valid());

    std::set<Dart_handle> darts_to_remove;

    //6) Glue pieces together
    //   using one edge per intersection polyline, we merge the different volumes
    for (typename An_edge_per_polyline_map::const_iterator it=
           an_edge_per_polyline.begin(); it!=an_edge_per_polyline.end(); ++it) {
      CGAL_assertion(it->second.first.size()==2);
      //orientation of faces around the edge (to be sure we can do it)
      std::pair<int,int> indices = it->first;
      const std::pair<bool,int>& polyline_info=it->second.second;

      //get the two halfedges incident to the edge [indices.first,indices.second]
      Halfedge_handle first_hedge=it->second.first.begin()->second;
      Halfedge_handle second_hedge=boost::next(it->second.first.begin())->second;

      CGAL_assertion(nodes[indices.second]==get(ppmap,first_hedge->vertex()));
      CGAL_assertion(nodes[indices.first]==get(ppmap,
                     first_hedge->opposite()->vertex()));
      CGAL_assertion(nodes[indices.second]==get(ppmap,second_hedge->vertex()));
      CGAL_assertion(nodes[indices.first]==get(ppmap,
                     second_hedge->opposite()->vertex()));

      Polyhedron* first_poly  = it->second.first.begin()->first;
      Polyhedron* second_poly = boost::next(it->second.first.begin())->first;

      //different handling depending on the number of incident triangles to the edge.
      //After sewing there are two,three or four volumes if there are two,three or four incident triangles respectively
      if ( first_hedge->is_border() || first_hedge->opposite()->is_border() ) {
        if (second_hedge->is_border() || second_hedge->opposite()->is_border())
          sew_2_two_volumes_case(first_hedge, second_hedge, border_halfedges,
                                 selected_hedge_to_dart, mark_index,
                                 darts_to_remove, nodes, indices, polyline_info);
        else
          sew_2_three_volumes_case(second_hedge, first_hedge,indices,nodes,
                                   border_halfedges,selected_hedge_to_dart,
                                   second_poly,first_poly,mark_index,
                                   darts_to_remove,polyline_info);
      } else
        if (second_hedge->is_border()  || second_hedge->opposite()->is_border())
          sew_2_three_volumes_case(first_hedge, second_hedge, indices, nodes,
                                   border_halfedges, selected_hedge_to_dart,
                                   first_poly,second_poly,mark_index,
                                   darts_to_remove,polyline_info);
      else {
        //Sort the four triangle facets around their common edge
        //  we suppose that the exterior of the polyhedron is indicated by
        //  counterclockwise oriented facets.
        Vertex_handle P1=first_hedge->opposite()->next()->vertex();
        Vertex_handle P2=first_hedge->next()->vertex();
        //    when looking from the side of indices.second, the interior of the first polyhedron is described
        //    by turning counterclockwise from P1 to P2
        Vertex_handle Q1=second_hedge->opposite()->next()->vertex();
        Vertex_handle Q2=second_hedge->next()->vertex();
        //    when looking from the side of indices.second, the interior of the second polyhedron is described
        //    by turning from Q1 to Q2

        //check if the third point of each triangular face is an original point (stay -1)
        //or a intersection point (in that case we need the index of the corresponding node to
        //have the exact value of the point)
        int index_p1 = node_index_of_incident_vertex(
          first_hedge->opposite()->next(),border_halfedges);
        int index_p2=node_index_of_incident_vertex(
          first_hedge->next(), border_halfedges);
        int index_q1=node_index_of_incident_vertex(
          second_hedge->opposite()->next(), border_halfedges);
        int index_q2=node_index_of_incident_vertex(
          second_hedge->next(), border_halfedges);

#ifdef CGAL_COREFINEMENT_DEBUG
        std::cout << index_p1 << " " << index_p2 << " " << index_q1 << " " <<index_q2 <<
                  std::endl;
        std::cout << nodes[indices.first] << " | " << nodes[indices.second] <<
                  std::endl;
        std::cout << get(ppmap,P1) << " | " << get(ppmap,P2) << " | " << get(ppmap,
                  Q1) << " | " <<get(ppmap,Q2) << std::endl;
#endif

        if ( coplanar_triangles_case_handled( first_hedge,second_hedge,indices,
                                              nodes, index_p1,index_p2,index_q1,
                                              index_q2,selected_hedge_to_dart,
                                              mark_index, darts_to_remove,
                                              polyline_info) ) continue;

        CGAL_assertion(get(ppmap,P1) !=get(ppmap,Q1) && get(ppmap,P1)!=get(ppmap,Q2)
                       && get(ppmap,P2) !=get(ppmap,Q1) && get(ppmap,P2)!=get(ppmap,Q2));

        bool Q1_is_between_P1P2 =
          OOP::sorted_around_edge_filtered( indices.first, indices.second,
                                            index_p1, index_p2, index_q1,
                                            P1, P2, Q1, nodes, ppmap);
        bool Q2_is_between_P1P2 =
          OOP::sorted_around_edge_filtered( indices.first, indices.second,
                                            index_p1, index_p2, index_q2,
                                            P1, P2, Q2, nodes, ppmap);

        //Recover the dart that will be the start point of the different sewing
        //  dof_X_outside = dart of face of, meaning the triangle containing the
        //  point X and part of the volume outside of the corresponding polyhedron
        //-----first polyhedron
        Dart_handle dof_P1_outside =
          get_associated_dart(first_hedge->opposite(), selected_hedge_to_dart);
        Dart_handle dof_P2_outside =
          get_associated_dart(first_hedge, selected_hedge_to_dart);
        //-----second polyhedron
        Dart_handle dof_Q1_outside =
          get_associated_dart(second_hedge->opposite(),selected_hedge_to_dart);
        Dart_handle dof_Q2_outside =
          get_associated_dart(second_hedge, selected_hedge_to_dart);

        if ( Q1_is_between_P1P2 ) {
          if( Q2_is_between_P1P2 ) {
            bool P1_is_between_Q1Q2 =
              OOP::sorted_around_edge_filtered(indices.first, indices.second,
                                               index_q1, index_q2, index_p1,
                                               Q1, Q2, P1, nodes, ppmap);
            if (!P1_is_between_Q1Q2) {
              // poly_first  - poly_second            = P1Q1 U Q2P2
              // poly_second - poly_first             = {0}
              // poly_first \cap poly_second          = Q1Q2
              // opposite( poly_first U poly_second ) = P2P1
              sew_2_marked_darts( final_map(),final_map().beta(dof_P1_outside, 3),
                                  dof_Q1_outside, mark_index, nodes,
                                  indices, polyline_info); //P1Q1
              sew_2_marked_darts( final_map(),dof_Q2_outside,
                                  final_map().beta(dof_P2_outside,3),mark_index, nodes,
                                  indices, polyline_info); //Q2P2
              sew_2_marked_darts( final_map(),final_map().beta(dof_Q1_outside,3),
                                  final_map().beta(dof_Q2_outside,3), mark_index, nodes,
                                  indices, polyline_info); //Q1Q2
              sew_2_marked_darts( final_map(),dof_P2_outside,
                                  dof_P1_outside, mark_index,
                                  nodes, indices, polyline_info); //P2P1
              //update inside outside info (because darts from the same volume have been merged)
              final_map().template attribute<3>(final_map().beta(dof_Q1_outside,3))->info().inside
                .insert(first_poly); //update Q1Q2 inside poly
              final_map().template attribute<3>(dof_P2_outside)->info().outside
                .insert(second_poly);//update P2P1 outside poly
            } else {
              // poly_first  - poly_second            = Q2Q1
              // poly_second - poly_first             = P2P1
              // poly_first \cap poly_second          = P1Q2 U Q1P2
              // opposite( poly_first U poly_second ) = {O}
              sew_2_marked_darts( final_map(),dof_Q2_outside,
                                  dof_Q1_outside, mark_index, nodes,
                                  indices, polyline_info); //Q2Q1
              sew_2_marked_darts( final_map(),dof_P2_outside,
                                  dof_P1_outside, mark_index, nodes, indices,
                                  polyline_info); //P2P1
              sew_2_marked_darts( final_map(),final_map().beta(dof_Q1_outside,3),
                                  final_map().beta(dof_P2_outside,3), mark_index, nodes,
                                  indices, polyline_info); //Q1P2
              sew_2_marked_darts( final_map(),final_map().beta(dof_P1_outside,3),
                                  final_map().beta(dof_Q2_outside,3),mark_index, nodes,
                                  indices, polyline_info); //P1Q2
              //update inside outside info (because darts from the same volume have been merged)
              final_map().template attribute<3>(dof_Q2_outside)->info().inside.insert(
                first_poly); //update Q2Q1 inside poly
              final_map().template attribute<3>(dof_P2_outside)->info().inside.insert(
                second_poly);//update P2P1 inside poly
            }
          } else {
            // poly_first  - poly_second            = P1Q1
            // poly_second - poly_first             = P2Q2
            // poly_first \cap poly_second          = Q1P2
            // opposite( poly_first U poly_second ) = Q2P1
            sew_2_marked_darts( final_map(),final_map().beta(dof_P1_outside,3),
                                dof_Q1_outside, mark_index, nodes,
                                indices, polyline_info); //P1Q1
            sew_2_marked_darts( final_map(),dof_P2_outside,
                                final_map().beta(dof_Q2_outside,3), mark_index, nodes,
                                indices, polyline_info); //P2Q2
            sew_2_marked_darts( final_map(),final_map().beta(dof_Q1_outside,3),
                                final_map().beta(dof_P2_outside,3), mark_index, nodes,
                                indices, polyline_info); //Q1P2
            sew_2_marked_darts( final_map(),dof_Q2_outside,
                                dof_P1_outside, mark_index, nodes,
                                indices, polyline_info); //Q2P1
          }
        } else {
          if( Q2_is_between_P1P2 ) {
            // poly_first  - poly_second            = Q2P2
            // poly_second - poly_first             = Q1P1
            // poly_first \cap poly_second          = P1Q2
            // opposite( poly_first U poly_second ) = P2Q1
            sew_2_marked_darts( final_map(),dof_Q2_outside,
                                final_map().beta(dof_P2_outside,3), mark_index, nodes,
                                indices, polyline_info); //Q2P2
            sew_2_marked_darts( final_map(),final_map().beta(dof_Q1_outside,3),
                                dof_P1_outside, mark_index, nodes,
                                indices, polyline_info); //Q1P1
            sew_2_marked_darts( final_map(),final_map().beta(dof_P1_outside,3),
                                final_map().beta(dof_Q2_outside,3),mark_index, nodes,
                                indices, polyline_info); //P1Q2
            sew_2_marked_darts( final_map(),dof_P2_outside,
                                dof_Q1_outside, mark_index, nodes,
                                indices, polyline_info); //P2Q1
          } else {
            bool P1_is_between_Q1Q2 =
              OOP::sorted_around_edge_filtered(indices.first, indices.second,
                                               index_q1,index_q2,index_p1,Q1,Q2,P1,
                                               nodes,ppmap);
            if (!P1_is_between_Q1Q2) {
              // poly_first  - poly_second            = P1P2
              // poly_second - poly_first             = Q1Q2
              // poly_first \cap poly_second          = {0}
              // opposite( poly_first U poly_second ) = P2Q1 U Q2P1
              sew_2_marked_darts( final_map(),final_map().beta(dof_P1_outside,3),
                                  final_map().beta(dof_P2_outside,3), mark_index, nodes,
                                  indices, polyline_info); //P1P2
              sew_2_marked_darts( final_map(),final_map().beta(dof_Q1_outside,3),
                                  final_map().beta(dof_Q2_outside,3), mark_index, nodes,
                                  indices, polyline_info); //Q1Q2
              sew_2_marked_darts( final_map(),dof_P2_outside,
                                  dof_Q1_outside, mark_index, nodes,
                                  indices, polyline_info); //P2Q1
              sew_2_marked_darts( final_map(),dof_Q2_outside,
                                  dof_P1_outside, mark_index, nodes,
                                  indices, polyline_info); //Q2P1
              //update inside outside info (because darts from the same volume have been merged)
              final_map().template attribute<3>(final_map().beta(dof_Q1_outside,3))->info().outside
                .insert(first_poly); //update Q1Q2 outside poly
              final_map().template attribute<3>(final_map().beta(dof_P1_outside,3))->info().outside
                .insert(second_poly);//update P2P1 outside poly
            } else {
              // poly_first  - poly_second            = {0}
              // poly_second - poly_first             = Q1P1 U P2Q2
              // poly_first \cap poly_second          = P1P2
              // opposite( poly_first U poly_second ) = Q2Q1
              sew_2_marked_darts( final_map(),final_map().beta(dof_Q1_outside,3),
                                  dof_P1_outside, mark_index, nodes,
                                  indices, polyline_info); //Q1P1
              sew_2_marked_darts( final_map(),dof_P2_outside,
                                  final_map().beta(dof_Q2_outside,3), mark_index, nodes,
                                  indices, polyline_info); //P2Q2
              sew_2_marked_darts( final_map(),final_map().beta(dof_P1_outside,3),
                                  final_map().beta(dof_P2_outside,3), mark_index, nodes,
                                  indices, polyline_info); //P1P2
              sew_2_marked_darts( final_map(),dof_Q2_outside,
                                  dof_Q1_outside, mark_index, nodes, indices,
                                  polyline_info); //Q2Q1
              //update inside outside info (because darts from the same volume have been merged)
              final_map().template attribute<3>(final_map().beta(dof_P1_outside,3))->info().inside
                .insert(second_poly); //update P1P2 inside poly
              final_map().template attribute<3>(dof_Q2_outside)->info().outside
                .insert(first_poly);//update Q2Q1 outside poly
            }
          }
        }
      }
    }

#ifdef CGAL_COREFINEMENT_DEBUG
    std::cout << "number of darts to remove: " << darts_to_remove.size()
              <<std::endl;
#endif
    //remove darts from empty volumes
    for (typename std::set<Dart_handle>::iterator itdart=darts_to_remove.begin(),
         end=darts_to_remove.end(); itdart!=end; ++itdart) {
      final_map().erase_dart(*itdart);
    }

    //remove empty volumes
    typedef typename Combinatorial_map_3::template Attribute_range<3>::type
    Volume_attribute_range;
    Volume_attribute_range& ratrib=final_map().template attributes<3>();
    typename Volume_attribute_range::iterator curr=ratrib.begin(),end=ratrib.end();
    do {
      if (curr->info().is_empty)
        final_map().template erase_attribute<3>(curr++);
      else
        ++curr;
    } while(curr!=end);

    CGAL_assertion(final_map().is_valid());

    //update the info of each volume knowing about only one polyhedron:
    //this happens when one polyhedron has a connected component
    //that do not intersect the other polyhedron

    typedef Side_of_triangle_mesh<Polyhedron, Kernel, PolyhedronPointPMap>
      Inside_poly_test;

    CGAL_precondition(polyhedron_to_map_node_to_polyhedron_vertex.size()==2);
    Polyhedron* Poly_A =
      polyhedron_to_map_node_to_polyhedron_vertex.begin()->first;
    Polyhedron* Poly_B =
      boost::next(polyhedron_to_map_node_to_polyhedron_vertex.begin())->first;
    Inside_poly_test* inside_A_test_ptr=NULL;
    Inside_poly_test* inside_B_test_ptr=NULL;
    bool Poly_A_is_closed = Poly_A->is_closed();
    bool Poly_B_is_closed = Poly_B->is_closed();

#ifdef CGAL_COREFINEMENT_DEBUG
    final_map().display_characteristics(std::cout);
    std::cout << "\n";
#endif

    typename Combinatorial_map_3::template  One_dart_per_cell_range<3> cell_range=
      final_map().template one_dart_per_cell<3>();
    for (typename Combinatorial_map_3::template One_dart_per_cell_range<3>
         ::iterator it = cell_range.begin(), it_end=cell_range.end();
         it_end!= it; ++it )
    {
      internal_IOP::Volume_info<Polyhedron>& info =
        final_map().template attribute<3>(it)->info();
      std::size_t inside_size = info.inside.size();
      std::size_t outside_size = info.outside.size();

      // if a volume is not classified wrt the two polyhedra, it means the component we look at does not
      // is a disjoint (but maybe at a vertex TAG SL001)
      if ( inside_size + outside_size == 1) {
        bool is_inside = (inside_size==1);
        Polyhedron* current_poly= is_inside ? (*info.inside.begin())
                                            : (*info.outside.begin());
        Polyhedron* test_poly;
        Inside_poly_test* inside_test_ptr;
        if ( current_poly==Poly_A) {
          // if the polyhedron is not closed, we set Poly_A to be outside by default
          if (!Poly_B_is_closed){
            info.outside.insert(Poly_B);
            continue;
          }
          test_poly=Poly_B;
          if (inside_B_test_ptr == NULL)
            inside_B_test_ptr = new Inside_poly_test(*Poly_B, ppmap);
          inside_test_ptr=inside_B_test_ptr;
        } else {
          // if the polyhedron is not closed, we set Poly_B to be outside by default
          if (!Poly_A_is_closed){
            info.outside.insert(Poly_A);
            continue;
          }
          test_poly=Poly_A;
          if (inside_A_test_ptr == NULL)
            inside_A_test_ptr = new Inside_poly_test(*Poly_A, ppmap);
          inside_test_ptr=inside_A_test_ptr;
        }

        // We need to find a point of the volume that is not on the boundary of the other volume.
        // Then the position of this point give the position of the volume. If all the points are on
        // the bounday, we take the mid-point of an edge (which must not be on the boundary otherwise
        // it contradicts the fact that volumes are disjoint
        // We first use the dart we have since one_dart_per_incident_cell has a non-negligeable cost.
        typename Kernel::Point_3 query=final_map().template attribute<0>(it)->point();
        CGAL::Bounded_side res = (*inside_test_ptr)(query);
        if (res==ON_BOUNDARY) {
          typedef typename Combinatorial_map_3::
            template One_dart_per_incident_cell_range<0,3> Vertex_range;

          Vertex_range vertex_range =
            final_map().template one_dart_per_incident_cell<0,3>(it);
          typename Vertex_range::iterator vit = vertex_range.begin();

          CGAL_assertion( typename Combinatorial_map_3::Dart_handle(vit) ==
                          typename Combinatorial_map_3::Dart_handle(it) );
          ++vit;
          for ( ; vit!=vertex_range.end(); ++vit) {
            query=final_map().template attribute<0>(vit)->point();
            res = (*inside_test_ptr)(query);
            if ( res != ON_BOUNDARY ) break;
          }

          //take edge midpoint
          if (res == ON_BOUNDARY) {
            /// \todo see do a better test here. At least the construction cannot fail
            ///  but the mid-point can fall outside of the volume...
#ifdef CGAL_COREFINEMENT_DEBUG
#warning this is not exact!!!
#endif
            typename Kernel::Point_3 p1 = final_map().template attribute<0>(it)->point();
            typename Kernel::Point_3 p2 =
              final_map().template attribute<0>(final_map().beta(it,1))->point();
            query=midpoint(p1,p2);
            res = (*inside_test_ptr)(query);
          }

          CGAL_assertion( res!= ON_BOUNDARY );
        }

        if (res  == ON_BOUNDED_SIDE )
          info.inside.insert(test_poly);
        else
          info.outside.insert(test_poly);
      }

#ifdef CGAL_COREFINEMENT_DEBUG
      std::cout << "This volume has inside: ";
      for (typename std::set<Polyhedron*>::iterator itpoly=info.inside.begin();
           itpoly!=info.inside.end(); ++itpoly)
        std::cout << " " << *itpoly;
      std::cout << " and outside: ";
      for (typename std::set<Polyhedron*>::iterator itpoly=info.outside.begin();
           itpoly!=info.outside.end(); ++itpoly)
        std::cout << " " << *itpoly;
      std::cout << std::endl;
#endif
    }
    if (inside_A_test_ptr!=NULL) delete inside_A_test_ptr;
    if (inside_B_test_ptr!=NULL) delete inside_B_test_ptr;
  }

};

}
} //end of namespace CGAL::Corefinement

#endif //CGAL_INTERNAL_COREFINEMENT_COMBINATORIAL_MAP_OUTPUT_BUILDER_H
