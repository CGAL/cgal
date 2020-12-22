// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_POLYGONAL_SCHEMA_H
#define CGAL_POLYGONAL_SCHEMA_H 1

#include <CGAL/license/Surface_mesh_topology.h>

#include <CGAL/Polygonal_schema_fwd.h>
#include <vector>
#include <unordered_map>
#include <cstddef>
#include <string>
#include <initializer_list>
#include <algorithm>
#include <random>
#include <CGAL/assertions.h>
#include <CGAL/memory.h>
#include <CGAL/Polygonal_schema_min_items.h>
#include <CGAL/Combinatorial_map.h>
#include <CGAL/Generalized_map.h>
#include <CGAL/Random.h>

namespace CGAL {
namespace Surface_mesh_topology {

  namespace internal
  {
    /// @return opposite label of label s
    ///    (i.e. add/remove - depending if s is positive or negative)
    inline std::string opposite_label(const std::string& s)
    {
      CGAL_assertion(!s.empty());
      if (s[0]=='-')
      { return s.substr(1, std::string::npos); }

      return std::string("-")+s;
    }

    template<class Map, class Combinatorial_data_structure=
             typename Map::Combinatorial_data_structure>
    struct Polygonal_schema_tools
    {};
    template<class CMap>
    struct Polygonal_schema_tools<CMap, Combinatorial_map_tag>
    {
      typedef typename CMap::Dart_handle Dart_handle;

      static Dart_handle
      add_edge_to_face(CMap& cmap, const std::string& s,
                       Dart_handle prev_dart,
                       Dart_handle dart_same_label,
                       Dart_handle dart_opposite_label,
                       std::unordered_map<std::string, Dart_handle>&
                       edge_label_to_dart)
      {
        if (dart_same_label!=CMap::null_handle && dart_opposite_label!=CMap::null_handle)
        {
          std::cerr<<"Polygonal_schema ERROR: "<<"both labels "<<s
                   <<" and "<<internal::opposite_label(s)
                   <<" are already added in the surface."
                   <<" This label can not be use anymore."<<std::endl;
          return CMap::null_handle;
        }

        if (dart_same_label!=CMap::null_handle)
        {
          std::cerr<<"Polygonal_schema ERROR: "<<"label "<<s
                   <<" is already added in the surface."
                   <<" Since the surface is orientable, this label can "
                   <<"not be use anymore."<<std::endl;
          return CMap::null_handle;
        }

        Dart_handle res=cmap.create_dart();
        edge_label_to_dart[s]=res;
        cmap.info(res).m_label=s;

        if (prev_dart!=cmap.null_handle)
        { cmap.template link_beta<1>(prev_dart, res); }

        if (dart_opposite_label!=CMap::null_handle)
        { cmap.template link_beta<2>(res, dart_opposite_label); }

        return res;
      }

      const std::string& get_label(CMap& cmap, Dart_handle dh) const
      { return cmap.info(dh).m_label; }
    };
    template<class GMap>
    struct Polygonal_schema_tools<GMap, Generalized_map_tag>
    {
      typedef typename GMap::Dart_handle Dart_handle;

      // In a GMap, if an edge is 2-free, only one of its two dart has one label.
      // Otherwise, d has one label and alpha<0,2>(d) the opposite label.
      static Dart_handle
      add_edge_to_face(GMap& gmap, const std::string& s,
                       Dart_handle prev_dart,
                       Dart_handle dart_same_label,
                       Dart_handle dart_opposite_label,
                       std::unordered_map<std::string, Dart_handle>&
                       edge_label_to_dart)
      {
        if (dart_same_label!=GMap::null_handle && dart_opposite_label!=GMap::null_handle)
        {
          std::cerr<<"Polygonal_schema ERROR: "<<"both labels "<<s
                   <<" and "<<internal::opposite_label(s)
                   <<" are already added in the surface."
                   <<" This label can not be use anymore."<<std::endl;
          return GMap::null_handle;
        }

        Dart_handle res=gmap.create_dart();
        Dart_handle dh2=gmap.create_dart();

        gmap.template link_alpha<0>(res, dh2);
        if (prev_dart!=gmap.null_handle)
        { gmap.template link_alpha<1>(res, gmap.template alpha<0>(prev_dart)); }

        if (dart_same_label!=GMap::null_handle)
        { // Here dart_same_label!=GMap::null_handle
          std::string s2=internal::opposite_label(s);
          edge_label_to_dart[s2]=dh2;
          gmap.info(dh2).m_label=s2;

          gmap.template sew<2>(res, dart_same_label);
        }
        else
        { // Here either dart_opposite_label!=GMap::null_handle, or both are GMap::null_handle
          edge_label_to_dart[s]=res;
          gmap.info(res).m_label=s;

          if (dart_opposite_label!=GMap::null_handle)
          {
            std::string s2=internal::opposite_label(s);
            edge_label_to_dart[s2]=res;
            gmap.info(res).m_label=s2;

            gmap.template sew<2>(dh2, dart_opposite_label);
          }
        }

        return res;
      }

      std::string get_label(GMap& gmap, Dart_handle dh) const
      {
        if (gmap.info(dh).m_label.empty())
        {
          if (!gmap.template is_free<2>(dh))
          { return gmap.info(gmap.template alpha<2>(dh)).m_label; }
          else
          {
            return internal::opposite_label(gmap.info(gmap.template alpha<0>(dh)));
          }
        }
        return gmap.info(dh).m_label;
      }
    };

  }
  // end namespace internal

  struct Combinatorial_map_tag;
  struct Generalized_map_tag;

  template < class BaseModel >
  class Polygonal_schema_base: public BaseModel
  {
  public:
    typedef BaseModel                       Base;
    typedef Polygonal_schema_base           Self;
    typedef BaseModel                       Map; // Either a GMap or a CMap
    typedef typename Map::Dart_handle       Dart_handle;
    typedef typename Map::Dart_const_handle Dart_const_handle;
    typedef typename Map::size_type         size_type;

    Polygonal_schema_base() : Base(),
      mark_perforated(this->get_new_mark()),
      first_dart(this->null_handle),
      prev_dart(this->null_handle),
      facet_started(false)
    {}

    /// Start a new facet.
    void init_facet()
    {
      if (facet_started)
      {
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"you try to start a facet"
                 <<" but the previous facet is not yet ended."<<std::endl;
        return;
      }

      first_dart = this->null_handle;
      prev_dart  = this->null_handle;
      facet_started=true;
    }

    /// Add one edge to the current facet, given by its label
    ///  (any string, using minus sign for orientation)
    void add_edge_to_facet(const std::string& s)
    {
      if (!facet_started)
      {
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"you try to add an edge to a facet"
                 <<" but the facet is not yet started."<<std::endl;
        return;
      }

      Dart_handle dart_same_label=get_dart_labeled(s);
      Dart_handle dart_opposite_label=get_dart_labeled
                                      (internal::opposite_label(s));

      Dart_handle cur=internal::Polygonal_schema_tools<Map>::
        add_edge_to_face(*this, s, prev_dart, dart_same_label,
                         dart_opposite_label, edge_label_to_dart);

      if (prev_dart==this->null_handle)
      { first_dart=cur; }

      prev_dart=cur;
    }

    /// add all the given edges to the current facet.
    /// @param s the sequence of labels of edges to add.
    void add_edges_to_facet(const std::string& s)
    {
      if (!facet_started)
      {
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"you try to add edges to a facet"
                 <<" but the facet is not yet started."<<std::endl;
        return;
      }
      std::istringstream iss(s);
      for (std::string token; std::getline(iss, token, ' '); )
      { add_edge_to_facet(token); }
    }

    /// add one facet, given a sequence of labels.
    /// @param s the sequence of labels of edges to add.
    void add_facet(const std::string& s)
    {
      if (facet_started)
      {
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"you try to add a new facet"
                 <<" but the previous facet is not yet ended."<<std::endl;
        return;
      }
      init_facet();
      add_edges_to_facet(s);
      finish_facet();
    }

    /// add edges to the current facet,
    ///  given a sequence of labels, as an initializer list.
    void add_edges_to_facet(std::initializer_list<const char*> l)
    {
      if (!facet_started)
      {
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"you try to add edges to a facet"
                 <<" but the facet is not yet started."<<std::endl;
        return;
      }
       for (const char* e : l)
       { add_edge_to_facet(e); }
    }

    /// add a new facet, given a sequence of labels, as an initializer list.
    void add_facet(std::initializer_list<const char*> l)
    {
      if (facet_started)
      {
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"you try to add a new facet"
                 <<" but the previous facet is not yet ended."<<std::endl;
        return;
      }
      init_facet();
      add_edges_to_facet(l);
      finish_facet();
    }

    /// End of the facet. Return the first dart of this facet.
    Dart_handle finish_facet()
    {
      if (!facet_started)
      {
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"you try to end a facet"
                 <<" but the facet is not yet started."<<std::endl;
        return Map::null_handle;
      }
      CGAL_assertion( first_dart!=this->null_handle &&
                                  prev_dart!=this->null_handle );
      this->set_next(prev_dart, first_dart);

      facet_started=false;
      return first_dart;
    }

    /// @return dart with the given label, Map::null_handle if this dart does not exist.
    Dart_handle get_dart_labeled(const std::string& s) const
    {
      auto ite=edge_label_to_dart.find(s);
      if (ite==edge_label_to_dart.end())
      { return Map::null_handle; }

      return ite->second;
    }

    std::string get_label(Dart_handle dh) const
    { return internal::Polygonal_schema_tools<Map>::get_label(dh); }

    /// marks the whole facet containing dh as perforated
    /// @return the number of darts of the marked face
    size_type perforate_facet(Dart_handle dh)
    {
      if (this->is_marked(dh, mark_perforated))
      { return 0; }

      return this->template mark_cell<2>(dh, mark_perforated);
    }

    /// same method but using a label
    size_type perforate_facet(const std::string& s)
    {
      auto ite=edge_label_to_dart.find(s);
      if (ite==edge_label_to_dart.end())
      {// maybe there is no need to put an error message
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"you try to label "<<s<<" to be a border"
                 <<" but this label does not exist yet."<<std::endl;
        return 0;
      }

      return perforate_facet(ite->second);
    }

    /// unmark the facet as being perforated, now the facet is filled
    /// @return the number of darts of the unmarked face
    size_type fill_facet(Dart_handle dh)
    {
      if (!this->is_marked(dh, mark_perforated))
      { return 0; }

      return this->template unmark_cell<2>(dh, mark_perforated);
    }

    /// same fonciton but using a label
    size_type fill_facet(const std::string& s)
    {
      auto ite=edge_label_to_dart.find(s);
      if (ite==edge_label_to_dart.end())
      {// maybe there is no need to put an error message
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"you try to label "<<s<<" to be a non border"
                 <<" but this label does not exist yet."<<std::endl;
        return 0;
      }

      return fill_facet(ite->second);
    }

    /// @return true iff dh is on a perforated facet
    bool is_perforated(Dart_const_handle dh) const
    { return this->is_marked(dh, mark_perforated); }

    /// same thing but using a label instead of a dart
    bool is_perforated(const std::string& s) const
    {
      auto ite=edge_label_to_dart.find(s);
      if (ite==edge_label_to_dart.end())
      {// maybe there is no need to put an error message
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"you ask if label "<<s<<" represents a dart border"
                 <<" but this label does not exist yet."<<std::endl;
        return false;
      }

      return is_perforated(ite->second);
    }

    std::size_t get_perforated_mark() const
    { return mark_perforated; }

    void display_perforated_darts() const
    {
      std::cout<<"labels is_free<2> is_perforated"<<std::endl;
      for (auto it=edge_label_to_dart.begin(), itend=edge_label_to_dart.end();
           it!=itend; ++it)
      {
        std::cout<<it->first<<" "<<Self::template is_free<2>(it->second)
                 <<" "<<is_perforated(it->second)<<std::endl;
      }
    }

  protected:
    // For each edge label, its corresponding dart. Stores both association
    // a -a, to allow users to start to add either a or -a.
    std::unordered_map<std::string, Dart_handle> edge_label_to_dart;
    std::size_t mark_perforated; // mark for perforated facets.

    // Data members used when we create a facet.
    Dart_handle first_dart;
    Dart_handle prev_dart;
    bool        facet_started;
  };

  /// Polygonal schema with combinatorial map.
  template <class Items_, class Alloc_, class Storage_>
  class Polygonal_schema_with_combinatorial_map:
    public Polygonal_schema_base<CGAL::Combinatorial_map_base
      <2,
      Polygonal_schema_with_combinatorial_map<Items_, Alloc_, Storage_>,
      Items_, Alloc_, Storage_> >
  {
  public:
    typedef Polygonal_schema_with_combinatorial_map<Items_, Alloc_, Storage_> Self;
    typedef Combinatorial_map_base<2, Self, Items_, Alloc_, Storage_>         CMap_base;
    typedef Polygonal_schema_base<CMap_base>                                  Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Dart_const_handle Dart_const_handle;

    Polygonal_schema_with_combinatorial_map() : Base()
    {}

    Polygonal_schema_with_combinatorial_map(const Self& amap) : Base(amap)
    {}

    template <unsigned int d2, typename Refs2, typename Items2, typename Alloc2,
              typename Storage2>
    Polygonal_schema_with_combinatorial_map(const Combinatorial_map_base
                                            <d2, Refs2, Items2, Alloc2, Storage2>&
                                            amap) : Base(amap)
    {}

    template <unsigned int d2, typename Refs2, typename Items2, typename Alloc2,
              typename Storage2, typename Converters>
    Polygonal_schema_with_combinatorial_map(const Combinatorial_map_base
                                            <d2, Refs2, Items2, Alloc2, Storage2>&
                                            amap, const Converters& converters) :
      Base(amap, converters)
    {}

    template <unsigned int d2, typename Refs2, typename Items2, typename Alloc2,
              typename Storage2, typename Converters,
              typename DartInfoConverter>
    Polygonal_schema_with_combinatorial_map(const Combinatorial_map_base
                                            <d2, Refs2, Items2, Alloc2, Storage2>&
                                            amap, const Converters& converters,
                                            const DartInfoConverter& dartinfoconverter) :
      Base(amap, converters, dartinfoconverter)
    {}

    template <unsigned int d2, typename Refs2, typename Items2, typename Alloc2,
              typename Storage2, typename Converters,
              typename DartInfoConverter, typename PointConverter >
    Polygonal_schema_with_combinatorial_map(const Combinatorial_map_base
                                            <d2, Refs2, Items2, Alloc2, Storage2>&
                                            amap, const Converters& converters,
                                            const DartInfoConverter& dartinfoconverter,
                                            const PointConverter& pointconverter) :
      Base(amap, converters, dartinfoconverter, pointconverter)
    {}
  };

  /// Polygonal schema with generalized map.
  template <class Items_, class Alloc_, class Storage_>
  class Polygonal_schema_with_generalized_map:
    public Polygonal_schema_base<CGAL::Generalized_map_base
      <2,
      Polygonal_schema_with_generalized_map<Items_, Alloc_, Storage_>,
      Items_, Alloc_, Storage_> >
  {
  public:
    typedef Polygonal_schema_with_generalized_map<Items_, Alloc_, Storage_> Self;
    typedef Generalized_map_base<2, Self, Items_, Alloc_, Storage_>         GMap_base;
    typedef Polygonal_schema_base<GMap_base>                                Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Dart_const_handle Dart_const_handle;

    Polygonal_schema_with_generalized_map() : Base()
    {}

    Polygonal_schema_with_generalized_map(const Self& amap) : Base(amap)
    {}

    template <unsigned int d2, typename Refs2, typename Items2, typename Alloc2,
              typename Storage2>
    Polygonal_schema_with_generalized_map(const Generalized_map_base
                                          <d2, Refs2, Items2, Alloc2, Storage2>&
                                            amap) : Base(amap)
    {}

    template <unsigned int d2, typename Refs2, typename Items2, typename Alloc2,
              typename Storage2, typename Converters>
    Polygonal_schema_with_generalized_map(const Generalized_map_base
                                          <d2, Refs2, Items2, Alloc2, Storage2>&
                                            amap, const Converters& converters) :
      Base(amap, converters)
    {}

    template <unsigned int d2, typename Refs2, typename Items2, typename Alloc2,
              typename Storage2, typename Converters,
              typename DartInfoConverter>
    Polygonal_schema_with_generalized_map(const Generalized_map_base
                                          <d2, Refs2, Items2, Alloc2, Storage2>&
                                            amap, const Converters& converters,
                                            const DartInfoConverter& dartinfoconverter) :
      Base(amap, converters, dartinfoconverter)
    {}

    template <unsigned int d2, typename Refs2, typename Items2, typename Alloc2,
              typename Storage2, typename Converters,
              typename DartInfoConverter, typename PointConverter >
    Polygonal_schema_with_generalized_map(const Generalized_map_base
                                          <d2, Refs2, Items2, Alloc2, Storage2>&
                                            amap, const Converters& converters,
                                            const DartInfoConverter& dartinfoconverter,
                                            const PointConverter& pointconverter) :
      Base(amap, converters, dartinfoconverter, pointconverter)
    {}
  };

  /// Generate a random polygonal schema ps.
  /// @param nb_labels the number of labels used to generate ps.
  /// @param seed the seed used for random
  /// @param max_dart_per_face maximal number of darts per face
  /// @param closed if true generates a closed polygonal schema.
  /// @param percentage_of_perforated percentage of perforated faces. If 0
  ///         no perforated faces.
  template<typename PS>
  void generate_random_polygonal_schema(PS& ps, std::size_t nb_labels,
                                        unsigned int seed
                                        =static_cast<unsigned int>(std::time(nullptr)),
                                        std::size_t max_dart_per_face=30,
                                        bool closed=true,
                                        std::size_t percentage_of_perforated=20)
  {
    CGAL::Random random(seed);
    std::vector<std::string> all_labels(nb_labels*2);
    for (std::size_t i=0; i<nb_labels; ++i)
    {
      all_labels[2*i]=std::to_string(static_cast<int>(i)+1);
      all_labels[(2*i)+1]=std::to_string(-(static_cast<int>(i)+1));
    }

    std::shuffle(all_labels.begin(), all_labels.end(),
                 std::default_random_engine(seed));

    std::size_t endlabel=all_labels.size();
    if (!closed)
    { endlabel-=(all_labels.size()/10); } // We remove 10% of labels.

    for (std::size_t i=0; i<endlabel; )
    {
      ps.init_facet();
      for (std::size_t j=0,
           nb=static_cast<std::size_t>
           (random.get_int(1, static_cast<int>(max_dart_per_face)));
           i<endlabel && j<nb; ++i, ++j)
      { ps.add_edges_to_facet(all_labels[i]); }
      typename PS::Dart_handle dh=ps.finish_facet();

      if (static_cast<std::size_t>(rand()%100)<percentage_of_perforated)
      { ps.perforate_facet(dh); }
    }

    ps.keep_biggest_connected_component(); // We keep only the biggest cc.
  }

} //namespace Surface_mesh_topology
} //namespace CGAL

#endif // CGAL_POLYGONAL_SCHEMA_H //
// EOF //
