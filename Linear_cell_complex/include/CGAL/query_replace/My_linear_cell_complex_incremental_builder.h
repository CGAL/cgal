// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
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
#ifndef MY_LINEAR_CELL_COMPLEX_INCREMENTAL_BUILDER_H
#define MY_LINEAR_CELL_COMPLEX_INCREMENTAL_BUILDER_H 1

#include <vector>
#include <cstddef>
#include <unordered_map>
#include <initializer_list>
#include <CGAL/Linear_cell_complex_base.h>

///////////////////////////////////////////////////////////////////////////////
template<class LCC, class Combinatorial_data_structure=
         typename LCC::Combinatorial_data_structure>
struct Add_vertex_to_face
{
  static typename LCC::Dart_handle run(LCC&,
                                       typename LCC::Vertex_attribute_handle,
                                       typename LCC::Dart_handle)
  {}
};
template<class LCC>
struct Add_vertex_to_face<LCC, CGAL::Combinatorial_map_tag>
{
  static typename LCC::Dart_handle run(LCC& lcc,
                                       typename LCC::Vertex_attribute_handle vh,
                                       typename LCC::Dart_handle prev_dart)
  {
    typename LCC::Dart_handle res=lcc.create_dart(vh);
    if (prev_dart!=lcc.null_handle)
    { 
       lcc.template link_beta<1>(prev_dart, res); 
    }
    return res;
  }
  static void run_for_last(LCC&,
                           typename LCC::Vertex_attribute_handle,
                           typename LCC::Dart_handle)
  { // here nothing to do, all darts were already created.
  }
};
template<class LCC>
struct Add_vertex_to_face<LCC, CGAL::Generalized_map_tag>
{
  static typename LCC::Dart_handle run(LCC& lcc,
                                       typename LCC::Vertex_attribute_handle vh,
                                       typename LCC::Dart_handle prev_dart)
  {
    typename LCC::Dart_handle res=lcc.create_dart(vh);
    if (prev_dart!=lcc.null_handle)
    {
      lcc.template link_alpha<0>(prev_dart, res);
      lcc.template link_alpha<1>(res, lcc.create_dart(vh));
      res=lcc.template alpha<1>(res);
    }
    return res;
  }
  static void run_for_last(LCC& lcc,
                           typename LCC::Vertex_attribute_handle vh,
                           typename LCC::Dart_handle prev_dart)
  {
    // here we need to create a last dart and 0-link it
    assert(prev_dart!=lcc.null_handle);
    lcc.template link_alpha<0>(prev_dart, lcc.create_dart(vh));
  }
};
///////////////////////////////////////////////////////////////////////////////
template<class LCC, class Combinatorial_data_structure=
         typename LCC::Combinatorial_data_structure>
struct Find_opposite_2_no_control // No difference for CMap and GMap
{
  typedef typename LCC::Dart_handle             DH;
  typedef typename LCC::Vertex_attribute_handle VAH;
  static DH run(LCC&,
                std::unordered_map<VAH, std::unordered_map<VAH, DH>>&
                vertex_to_dart_map_in_surface,
                VAH vah1, VAH vah2)
  {
    // We are searching edge vah2->vah1 (the opposite of edge vah1->vah2)
    auto it2=vertex_to_dart_map_in_surface.find(vah2);
    if (it2!=vertex_to_dart_map_in_surface.end())
    {
      auto it1=it2->second.find(vah1);
      if (it1!=it2->second.end())
      { return it1->second; }
    }
    return nullptr;
  }
};
///////////////////////////////////////////////////////////////////////////////
template<class LCC, class Combinatorial_data_structure=
         typename LCC::Combinatorial_data_structure>
struct Find_opposite_2_with_control
{
  typedef typename LCC::Dart_handle             DH;
  typedef typename LCC::Vertex_attribute_handle VAH;
  static DH run(LCC&,
                std::unordered_map<VAH, std::unordered_map<VAH, DH>>&,
                VAH, VAH)
  { return nullptr; }
};
template<class LCC>
struct Find_opposite_2_with_control<LCC, CGAL::Combinatorial_map_tag>
{
  typedef typename LCC::Dart_handle             DH;
  typedef typename LCC::Vertex_attribute_handle VAH;
  static DH run(LCC& lcc,
                std::unordered_map<VAH, std::unordered_map<VAH, DH>>&
                vertex_to_dart_map_in_surface,
                VAH vah1, VAH vah2)
  {
    DH res=Find_opposite_2_no_control<LCC>::run(lcc,
                                                vertex_to_dart_map_in_surface,
                                                vah1, vah2);
    if (res!=nullptr)
    {
      if (!lcc.template is_free<2>(res))
      { // Here a dart vah1->vah2 already exists, and it was already 2-sewn.
        std::cerr<<"ERROR in My_linear_cell_complex_incremental_builder_3: try to use a same oriented edge twice."<<std::endl;
        return nullptr;
      }
    }
    if (Find_opposite_2_no_control<LCC>::run(lcc,
                                             vertex_to_dart_map_in_surface,
                                             vah2, vah1)!=nullptr)
    { // Here a dart vah1->vah2 already exists (but it was not already 2-sewn).
      std::cerr<<"ERROR in My_linear_cell_complex_incremental_builder_3: try to use a same oriented edge twice."<<std::endl;
      return nullptr;
    }

    return res;
  }
};
template<class LCC>
struct Find_opposite_2_with_control<LCC, CGAL::Generalized_map_tag>
{
  typedef typename LCC::Dart_handle             DH;
  typedef typename LCC::Vertex_attribute_handle VAH;
  static DH run(LCC& lcc,
                std::unordered_map<VAH, std::unordered_map<VAH, DH>>&
                vertex_to_dart_map_in_surface,
                VAH vah1, VAH vah2)
  {
    DH res=Find_opposite_2_no_control<LCC>::run(lcc,
                                                vertex_to_dart_map_in_surface,
                                                vah1, vah2);
    if (res!=nullptr)
    {
      if (!lcc.template is_free<2>(res))
      { // Here a dart vah1->vah2 already exists, and it was already 2-sewn.
        std::cerr<<"ERROR in My_linear_cell_complex_incremental_builder_3: try to use a same oriented edge twice."<<std::endl;
        return nullptr;
      }
    }
    return res;
  }
};
///////////////////////////////////////////////////////////////////////////////
template<class LCC, class Combinatorial_data_structure=
         typename LCC::Combinatorial_data_structure>
struct Add_edge_in_associative_array
{
  typedef typename LCC::Dart_handle             DH;
  typedef typename LCC::Vertex_attribute_handle VAH;
  static void run(LCC&, DH,
                  std::unordered_map<VAH, std::unordered_map<VAH, DH>>&)
   {}
};
template<class LCC>
struct Add_edge_in_associative_array<LCC, CGAL::Combinatorial_map_tag>
{
  typedef typename LCC::Dart_handle             DH;
  typedef typename LCC::Vertex_attribute_handle VAH;
  static void run(LCC& lcc, DH dh,
                  std::unordered_map<VAH, std::unordered_map<VAH, DH>>&
                  vertex_to_dart_map_in_surface)
  {
    vertex_to_dart_map_in_surface[lcc.vertex_attribute(dh)].insert
        (std::make_pair(lcc.vertex_attribute(lcc.next(dh)), dh));
  }
};
template<class LCC>
struct Add_edge_in_associative_array<LCC, CGAL::Generalized_map_tag>
{
  typedef typename LCC::Dart_handle             DH;
  typedef typename LCC::Vertex_attribute_handle VAH;
  static void run(LCC& lcc, DH dh,
                  std::unordered_map<VAH, std::unordered_map<VAH, DH>>&
                  vertex_to_dart_map_in_surface)
  {
    vertex_to_dart_map_in_surface[lcc.vertex_attribute(dh)].insert
        (std::make_pair(lcc.vertex_attribute(lcc.template alpha<0>(dh)), dh));

    vertex_to_dart_map_in_surface
        [lcc.vertex_attribute(lcc.template alpha<0>(dh))].insert
        (std::make_pair(lcc.vertex_attribute(dh), lcc.template alpha<0>(dh)));
  }
};
///////////////////////////////////////////////////////////////////////////////
template<class LCC_, unsigned int dim=LCC_::dimension>
struct Sew3_for_LCC_incremental_builder
{
  static void run(LCC_& lcc,
                  typename LCC_::Dart_handle dh1, typename LCC_::Dart_handle dh2)
  {
    if(dh1!=nullptr)
    {
      if(!lcc.template is_free<3>(dh1))
      {
        std::cerr<<"ERROR in My_linear_cell_complex_incremental_builder_3: "
                 <<"it exists more than 2 faces with same indices."<<std::endl;
      }
      else
      { lcc.template sew<3>(lcc.other_orientation(dh1), dh2); }
    }
  }
};
template<class LCC_>
struct Sew3_for_LCC_incremental_builder<LCC_, 2>
{
  static void run(LCC_&, typename LCC_::Dart_handle, typename LCC_::Dart_handle)
  {}
};
///////////////////////////////////////////////////////////////////////////////
// Incremental builder
template < class LCC_ >
class My_linear_cell_complex_incremental_builder_3
{
public:
  typedef LCC_ LCC;
  typedef typename LCC::Dart_handle             DH;
  typedef typename LCC::Vertex_attribute_handle VAH;
  typedef typename LCC::Point                   Point_3;
  typedef typename LCC::size_type               size_type;

  My_linear_cell_complex_incremental_builder_3(LCC & alcc) :
    lcc(alcc)
  {}

  VAH add_vertex(const Point_3& p)
  {
    VAH res=lcc.create_vertex_attribute(p);
    vertex_map.push_back(res);
    return res;
  }

  void begin_facet()
  { // std::cout<<"Begin facet: "<<std::flush;
    first_dart=lcc.null_handle;
    prev_dart =lcc.null_handle;
  }

  void add_vertex_to_facet(size_type i)
  {
    CGAL_assertion(i<vertex_map.size());
    // std::cout<<i<<"  "<<std::flush;
    DH cur_dart=Add_vertex_to_face<LCC>::run(lcc, vertex_map[i], prev_dart);
    if(prev_dart!=lcc.null_handle)
    {
      DH opposite=Find_opposite_2_with_control<LCC>::
                  run(lcc,
                      vertex_to_dart_map_in_surface,
                      lcc.vertex_attribute(prev_dart),
                      lcc.vertex_attribute(cur_dart));
      if(opposite!=lcc.null_handle)
      {
        CGAL_assertion(lcc.template is_free<2>(opposite));
        lcc.template set_opposite<2>(prev_dart, opposite);
      }

      Add_edge_in_associative_array<LCC>::run(lcc, prev_dart,
                                              vertex_to_dart_map_in_surface);

      if (i<min_vertex) { min_vertex=i; min_dart=cur_dart; }
      if (i>max_vertex) { max_vertex=i; }
    }
    else
    { first_dart=cur_dart; min_vertex=max_vertex=i; min_dart=cur_dart; }

    prev_dart=cur_dart;
  }

  // End of the facet. Return the first dart of this facet.
  DH end_facet()
  {
    CGAL_assertion(first_dart!=lcc.null_handle && prev_dart!=lcc.null_handle);

    Add_vertex_to_face<LCC>::run_for_last(lcc,
                                          lcc.vertex_attribute(first_dart),
                                          prev_dart);

    lcc.set_next(prev_dart, first_dart);

    DH opposite=Find_opposite_2_with_control<LCC>::
                run(lcc,
                    vertex_to_dart_map_in_surface,
                    lcc.vertex_attribute(prev_dart),
                    lcc.vertex_attribute(first_dart));
    if(opposite!=lcc.null_handle)
    {
      CGAL_assertion(lcc.template is_free<2>(opposite));
      lcc.template set_opposite<2>(prev_dart, opposite);
    }

    Add_edge_in_associative_array<LCC>::run(lcc, prev_dart,
                                            vertex_to_dart_map_in_surface);

    if(LCC::dimension>2)
    {
      opposite=opposite_face();
      Sew3_for_LCC_incremental_builder<LCC>::run(lcc, opposite, min_dart);
      add_face_in_array();
    }
    return first_dart;
  }

  DH add_facet(std::initializer_list<size_type> l)
  {
    begin_facet();
    for (std::size_t i:l)
    { add_vertex_to_facet(i); }
    return end_facet();
  }

  void begin_surface()
  {
    vertex_to_dart_map_in_surface.clear();
  }

  // End of the surface. Return one dart of the created surface.
  DH end_surface()
  { return first_dart; }

protected:
  /** test if the two given facets have the same vertex handle but with
   *  opposite orientations. For closed facets.
   * @return true iff the two facets have the same vertex handle with opposite
   *         orientation.
   */
  bool are_facets_opposite_and_same_vertex_handles(DH d1, DH d2) const
  {
    DH s1=d1;
    DH s2=d2;
    do
    {
      assert(lcc.is_next_exist(d1) && lcc.is_previous_exist(d2));
      assert(lcc.other_extremity(d2)!=lcc.null_handle);

      if (lcc.vertex_attribute(d1)!=lcc.vertex_attribute(d2))
      { return false; }
      d1=lcc.next(d1);
      d2=lcc.previous(d2);
    }
    while(d1!=s1);

    if (d2!=s2) { return false; }
    return true;
  }

  DH opposite_face()
  {
    auto it1=faces.find(min_vertex);
    if(it1==faces.end()) { return nullptr; }
    auto it2=it1->second.find(max_vertex);
    if(it2==it1->second.end()) { return nullptr; }
    for(auto it3=it2->second.begin(), it3end=it2->second.end(); it3!=it3end; ++it3)
    {
      if (are_facets_opposite_and_same_vertex_handles(*it3, min_dart))
      { return lcc.previous(*it3); }
    }
    return nullptr;
  }

  void add_face_in_array()
  {
    faces[min_vertex][max_vertex].push_back(min_dart);
  }

private:
  LCC& lcc;
  std::vector<VAH> vertex_map; // Map each index to the corresponding vertex handle

  // A map to associate to each edge of a surface its dart. The edge is given
  // by its two vertex handles (source-target).
  std::unordered_map<VAH, std::unordered_map<VAH, DH>> vertex_to_dart_map_in_surface;
  std::unordered_map<std::size_t, std::unordered_map<std::size_t, std::vector<DH>>> faces;

  DH first_dart; /// First dart of the current face
  DH prev_dart; /// Prev dart of the current face
  DH min_dart; /// dart with the min vertex of the current facet.
  std::size_t min_vertex, max_vertex; /// min and max indices of vertices of the current face
};
///////////////////////////////////////////////////////////////////////////////
/* Create an hexahedron, given the indices of its vertices (in the following
*  order), the vertex must already have been added in the incremental builder.
*      3
*     /|\
*    0-|-2
*     \|/
*      1
*/
template<typename IncrementalBuilder>
typename IncrementalBuilder::LCC::Dart_handle
make_tetrahedron_with_builder(IncrementalBuilder& ib,
                                   std::size_t i0,
                                   std::size_t i1,
                                   std::size_t i2,
                                   std::size_t i3)
{
  ib.begin_surface();
  ib.add_facet({i0,i1,i2});
  ib.add_facet({i1,i0,i3});
  ib.add_facet({i2,i1,i3});
  ib.add_facet({i0,i2,i3});
  return ib.end_surface();
}
///////////////////////////////////////////////////////////////////////////////
/*      4
 *     /|\
 *    0-|-3
 *    | | |
 *    1---2
 */
template<typename IncrementalBuilder>
typename IncrementalBuilder::LCC::Dart_handle
 make_pyramid_with_builder(IncrementalBuilder& ib,
                               std::size_t i0,
                               std::size_t i1,
                               std::size_t i2,
                               std::size_t i3,
                               std::size_t i4)
{
  ib.begin_surface();
  ib.add_facet({i0,i1,i2,i3});
  ib.add_facet({i1,i0,i4});
  ib.add_facet({i2,i1,i4});
  ib.add_facet({i3,i2,i4});
  ib.add_facet({i0,i3,i4});
  return ib.end_surface();
}
///////////////////////////////////////////////////////////////////////////////
/*      3
 *     /|\
 *    4---5
 *    | | |
 *    | 0 |
 *    |/ \|
 *    1---2
 */
template<typename IncrementalBuilder>
typename IncrementalBuilder::LCC::Dart_handle
 make_prism_with_builder(IncrementalBuilder& ib,
                             std::size_t i0,
                             std::size_t i1,
                             std::size_t i2,
                             std::size_t i3,
                             std::size_t i4,
                             std::size_t i5)
{
  ib.begin_surface();
  ib.add_facet({i0,i1,i2});
  ib.add_facet({i1,i0,i3,i4});
  ib.add_facet({i2,i1,i4,i5});
  ib.add_facet({i0,i2,i5,i3});
  ib.add_facet({i5,i4,i3});
  return ib.end_surface();
}
///////////////////////////////////////////////////////////////////////////////
/*      7----6
 *     /|   /|
 *    4----5 |
 *    | 3--|-2
 *    |/   |/
 *    0----1
 */
template<typename IncrementalBuilder>
typename IncrementalBuilder::LCC::Dart_handle
 make_hexahedron_with_builder(IncrementalBuilder& ib,
                                  std::size_t i0,
                                  std::size_t i1,
                                  std::size_t i2,
                                  std::size_t i3,
                                  std::size_t i4,
                                  std::size_t i5,
                                  std::size_t i6,
                                  std::size_t i7)
{
  ib.begin_surface();
  ib.add_facet({i0,i1,i2,i3});
  ib.add_facet({i1,i0,i4,i5});
  ib.add_facet({i2,i1,i5,i6});
  ib.add_facet({i3,i2,i6,i7});
  ib.add_facet({i0,i3,i7,i4});
  ib.add_facet({i7,i6,i5,i4});
  return ib.end_surface();
}
///////////////////////////////////////////////////////////////////////////////
template<typename IncrementalBuilder>
typename IncrementalBuilder::LCC::Dart_handle
 make_pentagonal_prism_with_builder(IncrementalBuilder& ib,
                                        std::size_t i0,
                                        std::size_t i1,
                                        std::size_t i2,
                                        std::size_t i3,
                                        std::size_t i4,
                                        std::size_t i5,
                                        std::size_t i6,
                                        std::size_t i7,
                                        std::size_t i8,
                                        std::size_t i9)
{
  ib.begin_surface();
  ib.add_facet({i0,i1,i2,i3,i4});
  ib.add_facet({i1,i0,i5,i6});
  ib.add_facet({i2,i1,i6,i7});
  ib.add_facet({i3,i2,i7,i8});
  ib.add_facet({i4,i3,i8,i9});
  ib.add_facet({i0,i4,i9,i5});
  ib.add_facet({i9,i8,i7,i6,i5});
  return ib.end_surface();
}
///////////////////////////////////////////////////////////////////////////////
template<typename IncrementalBuilder>
typename IncrementalBuilder::LCC::Dart_handle
 make_hexagonal_prism_with_builder(IncrementalBuilder& ib,
                                       std::size_t i0,
                                       std::size_t i1,
                                       std::size_t i2,
                                       std::size_t i3,
                                       std::size_t i4,
                                       std::size_t i5,
                                       std::size_t i6,
                                       std::size_t i7,
                                       std::size_t i8,
                                       std::size_t i9,
                                       std::size_t i10,
                                       std::size_t i11)
{
  ib.begin_surface();
  ib.add_facet({i0,i1,i2,i3,i4,i5});
  ib.add_facet({i1,i0,i6,i7});
  ib.add_facet({i2,i1,i7,i8});
  ib.add_facet({i3,i2,i8,i9});
  ib.add_facet({i4,i3,i9,i10});
  ib.add_facet({i5,i4,i10,i11});
  ib.add_facet({i0,i5,i11,i6});
  ib.add_facet({i11,i10,i9,i8,i7,i6});
  return ib.end_surface();
}
///////////////////////////////////////////////////////////////////////////////
#endif // MY_LINEAR_CELL_COMPLEX_INCREMENTAL_BUILDER_H //
// EOF //
