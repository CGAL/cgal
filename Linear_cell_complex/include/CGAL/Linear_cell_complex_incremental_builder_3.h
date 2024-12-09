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
#ifndef CGAL_LINEAR_CELL_COMPLEX_INCREMENTAL_BUILDER_3_H
#define CGAL_LINEAR_CELL_COMPLEX_INCREMENTAL_BUILDER_3_H 1

#include <vector>
#include <cstddef>
#include <unordered_map>
#include <initializer_list>
#include <CGAL/Linear_cell_complex_base.h>
#include <CGAL/assertions.h>

namespace CGAL {
///////////////////////////////////////////////////////////////////////////////
  template<class LCC, class Combinatorial_data_structure=
           typename LCC::Combinatorial_data_structure>
  struct Add_vertex_to_face
  {
    static typename LCC::Dart_descriptor run(LCC&,
                                         typename LCC::Vertex_attribute_descriptor,
                                         typename LCC::Dart_descriptor)
    {}
  };
  template<class LCC>
  struct Add_vertex_to_face<LCC, Combinatorial_map_tag>
  {
    static typename LCC::Dart_descriptor run(LCC& lcc,
                                         typename LCC::Vertex_attribute_descriptor vh,
                                         typename LCC::Dart_descriptor prev_dart)
    {
      typename LCC::Dart_descriptor res=lcc.create_dart(vh);
      if (prev_dart!=lcc.null_descriptor)
      {
        lcc.template link_beta<1>(prev_dart, res);
      }
      return res;
    }
    static void run_for_last(LCC&,
                             typename LCC::Vertex_attribute_descriptor,
                             typename LCC::Dart_descriptor)
    { // here nothing to do, all darts were already created.
    }
  };
  template<class LCC>
  struct Add_vertex_to_face<LCC, Generalized_map_tag>
  {
    static typename LCC::Dart_descriptor run(LCC& lcc,
                                         typename LCC::Vertex_attribute_descriptor vh,
                                         typename LCC::Dart_descriptor prev_dart)
    {
      typename LCC::Dart_descriptor res=lcc.create_dart(vh);
      if (prev_dart!=lcc.null_descriptor)
      {
        lcc.template link_alpha<0>(prev_dart, res);
        lcc.template link_alpha<1>(res, lcc.create_dart(vh));
        res=lcc.template alpha<1>(res);
      }
      return res;
    }
    static void run_for_last(LCC& lcc,
                             typename LCC::Vertex_attribute_descriptor vh,
                             typename LCC::Dart_descriptor prev_dart)
    {
      // here we need to create a last dart and 0-link it
      CGAL_assertion(prev_dart!=lcc.null_descriptor);
      lcc.template link_alpha<0>(prev_dart, lcc.create_dart(vh));
    }
  };
///////////////////////////////////////////////////////////////////////////////
template<class LCC, class Combinatorial_data_structure=
         typename LCC::Combinatorial_data_structure>
struct Find_opposite_2_no_control // No difference for CMap and GMap
{
  typedef typename LCC::Dart_descriptor             DH;
  typedef typename LCC::Vertex_attribute_descriptor VAH;
  static DH run(LCC& lcc,
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
    return lcc.null_descriptor;
  }
};
///////////////////////////////////////////////////////////////////////////////
template<class LCC, class Combinatorial_data_structure=
         typename LCC::Combinatorial_data_structure>
struct Find_opposite_2_with_control
{
  typedef typename LCC::Dart_descriptor             DH;
  typedef typename LCC::Vertex_attribute_descriptor VAH;
  static DH run(LCC& lcc,
                std::unordered_map<VAH, std::unordered_map<VAH, DH>>&,
                VAH, VAH)
  { return lcc.null_descriptor; }
};
template<class LCC>
struct Find_opposite_2_with_control<LCC, CGAL::Combinatorial_map_tag>
{
  typedef typename LCC::Dart_descriptor             DH;
  typedef typename LCC::Vertex_attribute_descriptor VAH;
  static DH run(LCC& lcc,
                std::unordered_map<VAH, std::unordered_map<VAH, DH>>&
                vertex_to_dart_map_in_surface,
                VAH vah1, VAH vah2)
  {
    DH res=Find_opposite_2_no_control<LCC>::run(lcc,
                                                vertex_to_dart_map_in_surface,
                                                vah1, vah2);
    if (res!=lcc.null_descriptor)
    {
      if (!lcc.template is_free<2>(res))
      { // Here a dart vah1->vah2 already exists, and it was already 2-sewn.
        std::cerr<<"ERROR in Linear_cell_complex_incremental_builder_3: try to use a same oriented edge twice."<<std::endl;
        return lcc.null_descriptor;
      }
    }
    if (Find_opposite_2_no_control<LCC>::run(lcc,
                                             vertex_to_dart_map_in_surface,
                                             vah2, vah1)!=lcc.null_descriptor)
    { // Here a dart vah1->vah2 already exists (but it was not already 2-sewn).
      std::cerr<<"ERROR in Linear_cell_complex_incremental_builder_3: try to use a same oriented edge twice."<<std::endl;
      return lcc.null_descriptor;
    }

    return res;
  }
};
template<class LCC>
struct Find_opposite_2_with_control<LCC, CGAL::Generalized_map_tag>
{
  typedef typename LCC::Dart_descriptor             DH;
  typedef typename LCC::Vertex_attribute_descriptor VAH;
  static DH run(LCC& lcc,
                std::unordered_map<VAH, std::unordered_map<VAH, DH>>&
                vertex_to_dart_map_in_surface,
                VAH vah1, VAH vah2)
  {
    DH res=Find_opposite_2_no_control<LCC>::run(lcc,
                                                vertex_to_dart_map_in_surface,
                                                vah1, vah2);
    if (res!=lcc.null_descriptor)
    {
      if (!lcc.template is_free<2>(res))
      { // Here a dart vah1->vah2 already exists, and it was already 2-sewn.
        std::cerr<<"ERROR in Linear_cell_complex_incremental_builder_3: try to use a same oriented edge twice."<<std::endl;
        return lcc.null_descriptor;
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
  typedef typename LCC::Dart_descriptor             DH;
  typedef typename LCC::Vertex_attribute_descriptor VAH;
  static void run(LCC&, DH,
                  std::unordered_map<VAH, std::unordered_map<VAH, DH>>&)
   {}
};
template<class LCC>
struct Add_edge_in_associative_array<LCC, CGAL::Combinatorial_map_tag>
{
  typedef typename LCC::Dart_descriptor             DH;
  typedef typename LCC::Vertex_attribute_descriptor VAH;
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
  typedef typename LCC::Dart_descriptor             DH;
  typedef typename LCC::Vertex_attribute_descriptor VAH;
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
struct Sew3_for_LCC_incremental_builder_3
{
  static void run(LCC_& lcc,
                  typename LCC_::Dart_descriptor dh1, typename LCC_::Dart_descriptor dh2)
  {
    if(dh1!=lcc.null_descriptor)
    {
      if(!lcc.template is_free<3>(dh1))
      {
        std::cerr<<"ERROR in Linear_cell_complex_incremental_builder_3: "
                 <<"it exist more than 2 faces with same indices."<<std::endl;
      }
      else
      { lcc.template sew<3>(lcc.other_orientation(dh1), dh2); }
    }
  }
};
template<class LCC_>
struct Sew3_for_LCC_incremental_builder_3<LCC_, 2>
{
  static void run(LCC_&, typename LCC_::Dart_descriptor, typename LCC_::Dart_descriptor)
  {}
};
///////////////////////////////////////////////////////////////////////////////
// Incremental builder
template < class LCC_ >
class Linear_cell_complex_incremental_builder_3
{
public:
  typedef LCC_ LCC;
  typedef typename LCC::Dart_descriptor             DH;
  typedef typename LCC::Vertex_attribute_descriptor VAH;
  typedef typename LCC::Point                   Point_3;
  typedef typename LCC::size_type               size_type;

  Linear_cell_complex_incremental_builder_3(LCC & alcc) :
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
    first_dart=lcc.null_descriptor;
    prev_dart =lcc.null_descriptor;
  }

  void add_vertex_to_facet(size_type i)
  {
    CGAL_assertion(i<vertex_map.size());
    // std::cout<<i<<"  "<<std::flush;
    DH cur_dart=Add_vertex_to_face<LCC>::run(lcc, vertex_map[i], prev_dart);
    if ( prev_dart!=lcc.null_descriptor )
    {
      DH opposite=Find_opposite_2_with_control<LCC>::
        run(lcc,
            vertex_to_dart_map_in_surface,
            lcc.vertex_attribute(prev_dart),
            lcc.vertex_attribute(cur_dart));
      if ( opposite!=lcc.null_descriptor )
      {
        CGAL_assertion( lcc.template is_free<2>(opposite) );
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
    CGAL_assertion( first_dart!=lcc.null_descriptor && prev_dart!=lcc.null_descriptor );

    Add_vertex_to_face<LCC>::run_for_last(lcc,
                                          lcc.vertex_attribute(first_dart),
                                          prev_dart);

    lcc.set_next(prev_dart, first_dart);

    DH opposite=Find_opposite_2_with_control<LCC>::
                run(lcc,
                    vertex_to_dart_map_in_surface,
                    lcc.vertex_attribute(prev_dart),
                    lcc.vertex_attribute(first_dart));
    if ( opposite!=lcc.null_descriptor )
    {
      CGAL_assertion( lcc.template is_free<2>(opposite) );
      lcc.template set_opposite<2>(prev_dart, opposite);
    }

    Add_edge_in_associative_array<LCC>::run(lcc, prev_dart,
                                            vertex_to_dart_map_in_surface);

    if(LCC::dimension>2)
    {
      opposite=opposite_face();
      Sew3_for_LCC_incremental_builder_3<LCC>::run(lcc, opposite, min_dart);
      add_face_in_array();
    }
    return first_dart;
  }

  DH add_facet(std::initializer_list<size_type> l)
  {
    begin_facet();
    for (size_type i:l)
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
  bool are_facets_opposite_and_same_vertex_descriptors(DH d1, DH d2) const
  {
    DH s1=d1;
    DH s2=d2;
    do
    {
      CGAL_assertion(lcc.is_next_exist(d1) && lcc.is_previous_exist(d2));
      CGAL_assertion(lcc.other_extremity(d2)!=lcc.null_descriptor);

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
    if(it1==faces.end()) { return lcc.null_descriptor; }
    auto it2=it1->second.find(max_vertex);
    if(it2==it1->second.end()) { return lcc.null_descriptor; }
    for(auto it3=it2->second.begin(), it3end=it2->second.end(); it3!=it3end; ++it3)
    {
      if (are_facets_opposite_and_same_vertex_descriptors(*it3, min_dart))
      { return lcc.previous(*it3); }
    }
    return lcc.null_descriptor;
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
  size_type min_vertex, max_vertex; /// min and max indices of vertices of the current face
};

} //namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_INCREMENTAL_BUILDER_3_H //
// EOF //
