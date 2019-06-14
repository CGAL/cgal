// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_POLYGONAL_SCHEMA_H
#define CGAL_POLYGONAL_SCHEMA_H 1

#include <vector>
#include <unordered_map>
#include <cstddef>
#include <string>
//#include <CGAL/Path_on_surface.h>
#include <CGAL/Polygonal_schema_min_items.h>

namespace CGAL {
  /// @return opposite label of label s
  ///    (i.e. add/remove - depending if s is positive or negative)
  inline std::string opposite_label(const std::string & s)
  {
    CGAL_assertion(!s.empty());
    if (s[0]=='-')
    { return s.substr(1, std::string::npos); }
    
    return std::string("-")+s;
  }

  struct Combinatorial_map_tag;
  struct Generalized_map_tag;
  
  template<class Map, class Combinatorial_data_structure=
           typename Map::Combinatorial_data_structure>
  struct Map_incremental_builder_tools
  {};
  template<class CMap>
  struct Map_incremental_builder_tools<CMap, Combinatorial_map_tag>
  {
    typedef typename CMap::Dart_handle Dart_handle;
    
    static Dart_handle
    add_edge_to_face(CMap& cmap, const std::string& s,
                     Dart_handle prev_dart,
                     Dart_handle dart_same_label,
                     Dart_handle dart_opposite_label,
                     std::unordered_map<std::string, Dart_handle>& edge_label_to_dart)
    {
      if (dart_same_label!=NULL && dart_opposite_label!=NULL)
      {
        std::cerr<<"Polygonal_schema ERROR: "<<"both labels "<<s
                 <<" and "<<opposite_label(s)<<" are already added in the surface."
                 <<" This label can not be use anymore."<<std::endl;
        return NULL;
      }

      if (dart_same_label!=NULL)
      {
       std::cerr<<"Polygonal_schema ERROR: "<<"label "<<s
                 <<" is already added in the surface."
                 <<" Since the surface is orientable, this label can not be use anymore."<<std::endl;
        return NULL;
      }
      
      Dart_handle res=cmap.create_dart();
      edge_label_to_dart[s]=res;

      cmap.info(res).m_label=new char[s.size()+1];
      strncpy(cmap.info(res).m_label, s.c_str(), s.size()+1); // +1 to copy also the \0 char

      if (prev_dart!=cmap.null_handle)
      { cmap.template link_beta<1>(prev_dart, res); }
      
      if (dart_opposite_label!=NULL)
      { cmap.template link_beta<2>(res, dart_opposite_label); }
      
      return res;
    }
  };
  template<class GMap>
  struct Map_incremental_builder_tools<GMap, Generalized_map_tag>
  {
    typedef typename GMap::Dart_handle Dart_handle;

    static Dart_handle
    add_edge_to_face(GMap& gmap, const std::string& s,
                     Dart_handle prev_dart,
                     Dart_handle dart_same_label,
                     Dart_handle dart_opposite_label,
                     std::unordered_map<std::string, Dart_handle>& edge_label_to_dart)
    {
       if (dart_same_label!=NULL && dart_opposite_label!=NULL)
      {
        std::cerr<<"Polygonal_schema ERROR: "<<"both labels "<<s
                 <<" and "<<opposite_label(s)<<" are already added in the surface."
                 <<" This label can not be use anymore."<<std::endl;
        return NULL;
      }
      
      Dart_handle res=gmap.create_dart();
      Dart_handle dh2=gmap.create_dart();

      gmap.template link_alpha<0>(res, dh2);
      if (prev_dart!=gmap.null_handle)
      { gmap.template link_alpha<1>(res, gmap.template alpha<0>(prev_dart)); }

      if ((dart_same_label==NULL && dart_opposite_label==NULL) ||
          dart_opposite_label!=NULL)
      {
        edge_label_to_dart[s]=res;
        gmap.info(res).m_label=new char[s.size()+1];
        strncpy(gmap.info(res).m_label, s.c_str(), s.size()+1); // +1 to copy also the \0 char

        if (dart_opposite_label!=NULL)
        { gmap.template sew<2>(dh2, dart_opposite_label); }
      }
      else
      { // Here dart_same_label!=NULL
        std::string s2=opposite_label(s);
        edge_label_to_dart[s2]=dh2;
        gmap.info(dh2).m_label=new char[s2.size()+1];
        strncpy(gmap.info(dh2).m_label, s2.c_str(), s2.size()+1); // +1 to copy also the \0 char

        gmap.template sew<2>(res, dart_same_label);
      }

      return res;
    }
  };

  template < class BaseModel >
  class Polygonal_schema_base: public BaseModel
  {
  public:
    typedef BaseModel                       Base;
    typedef BaseModel                       Map; // Either a GMap or a CMap
    typedef typename Map::Dart_handle       Dart_handle;
    typedef typename Map::Dart_const_handle Dart_const_handle;
    typedef typename Map::size_type         size_type;

    Polygonal_schema_base() : Base(),
      first_dart(this->null_handle),
      prev_dart(this->null_handle),
      facet_started(false)
    {}

    ~Polygonal_schema_base()
    {
      for (auto it=this->darts().begin(), itend=this->darts().end(); it!=itend; ++it)
      {
        if (this->info(it).m_label!=NULL)
        {
          delete []this->info(it).m_label;
          this->info(it).m_label=NULL;
        }
      }
    }
    
    /// Start a new facet.
    void begin_facet()
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
      // std::cout<<"Begin facet: "<<std::flush;
    }

    /// Add one edge to the current facet, given by its label (any string, using minus sign for orientation)
    void add_edge_to_facet(const std::string& s)
    {
      if (!facet_started)
      {
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"you try to add an edge to a facet"
                 <<" but the facet is not yet started."<<std::endl;
        return;
      }

      Dart_handle dart_same_label=find_dart_with_label(s);
      Dart_handle dart_opposite_label=find_dart_with_label(opposite_label(s));
      
      Dart_handle cur=Map_incremental_builder_tools<Map>::
        add_edge_to_face(*this, s, prev_dart, dart_same_label, dart_opposite_label,
                         edge_label_to_dart);

      if (prev_dart==this->null_handle)
      { first_dart=cur; }

      prev_dart=cur;
    }

    /// add the given edges to the current facet
    /// s is a sequence of labels, add all the corresponding edges into the current facet.
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
      
    /// add one facet, s is a sequence of labels, add all the corresponding edges into a new facet.
    void add_facet(const std::string& s)
    {
      if (facet_started)
      {
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"you try to add a new facet"
                 <<" but the previous facet is not yet ended."<<std::endl;
        return;
      }
      begin_facet();
      add_edges_to_facet(s);
      end_facet();
    }

    /// End of the facet. Return the first dart of this facet.
    Dart_handle end_facet()
    {
      if (!facet_started)
      {
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"you try to end a facet"
                 <<" but the facet is not yet started."<<std::endl;
        return NULL;
      }
      CGAL_assertion( first_dart!=this->null_handle && prev_dart!=this->null_handle );
      this->set_next(prev_dart, first_dart);
      
      facet_started=false;
      return first_dart;
    }

    /// Start a new surface
    void begin_surface()
    {
      if (facet_started) { end_facet(); }

      first_dart    = this->null_handle;
      prev_dart     = this->null_handle;
      edge_label_to_dart.clear();      
    }

    /// End of the surface. Return one dart of the created surface.
    Dart_handle end_surface()
    { return first_dart; }

    /// Start a path on the surface
    void begin_path()
    {
      /*      if (path_started)
      {
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"you try to start a path"
                 <<" but the previous path is not yet ended."<<std::endl;
        return;
      }
      path_started=true;
      m_first_path_vertex=true;
      m_cur_path.clear();*/
    }

    /// Add edge labeled e at the end of the current path
    void add_edge_to_path(const std::string& e)
    {
      /*      if (!path_started)
      {
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"you try to add an edge to a path"
                 <<" but the path is not yet started."<<std::endl;
        return;
      }

      Dart_const_handle dh=find_dart_with_label(e);
      if (dh==NULL)
      {
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"edge labeled ("<<e<<") does not exists "
                 <<"and thus cannot be added in the path."<<std::endl;
        return;
      }

      if (!m_cur_path.can_be_pushed(dh))
      {
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"edge labeled ("<<e<<") is not adjacent to the previous "
                 <<"edge of the path and thus cannot be added."<<std::endl;
        return;
      }
      m_cur_path.push_back(dh);*/
    }

    /// End the current path
    /*    CGAL::Path_on_surface<Map> end_path()
    {
            if (!path_started)
      {
        std::cerr<<"Polygonal_schema ERROR: "
                 <<"you try to end a path"
                 <<" but the path is not yet started."<<std::endl;
        return m_cur_path;
      }
      path_started=false;
      return m_cur_path;
    }

    /// A shortcut allowing to create a path directly with a sequence
    /// of vertex ids, if the map was built by adding vertices
    /// or edge labels, if the map was built by adding edges
    CGAL::Path_on_surface<Map> create_path(const std::string& s)
    {
      begin_path();

      std::istringstream iss(s);
      for (std::string token; std::getline(iss, token, ' '); )
      { add_edge_to_path(token); }

      return end_path();
      } */

    /// @return dart with the given label, NULL if this dart does not exist.
    Dart_handle find_dart_with_label(const std::string & s) const
    {
      auto ite=edge_label_to_dart.find(s);
      if (ite==edge_label_to_dart.end())
      { return NULL; }

      return ite->second;
    }

    const char* get_label(Dart_handle dh) const
    {
      // TODO
    }
    
  protected:
    // For each edge label, its corresponding dart. Stores both association a -a, to allow
    // users to start to add either a or -a.
    std::unordered_map<std::string, Dart_handle> edge_label_to_dart;
    
    Dart_handle first_dart;
    Dart_handle prev_dart;
    bool        facet_started;
  };

  template <class Items_=Polygonal_schema_min_items,
            class Alloc_=CGAL_ALLOCATOR(int),
            class Storage_= Combinatorial_map_storage_1<2, Items_, Alloc_> >
  class Polygonal_schema_with_combinatorial_map:
    public Polygonal_schema_base<CGAL::Combinatorial_map_base
                                 <2,
                                  Polygonal_schema_with_combinatorial_map<Items_, Alloc_, Storage_>,
                                  Items_, Alloc_, Storage_> >
  {
  public:
    typedef Polygonal_schema_with_combinatorial_map<Items_, Alloc_, Storage_>  Self;
    typedef Combinatorial_map_base<2, Self, Items_, Alloc_, Storage_>          CMap_base;
    typedef Polygonal_schema_base<CMap_base>                                   Base;
    
    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Dart_const_handle Dart_const_handle;

    Polygonal_schema_with_combinatorial_map() : Base()
    {}

    Polygonal_schema_with_combinatorial_map(const Self & amap) : Base(amap)
    {}

    template <unsigned int d2, typename Refs2, typename Items2, typename Alloc2,
              typename Storage2>
    Polygonal_schema_with_combinatorial_map(const Combinatorial_map_base<d2, Refs2, Items2, Alloc2, Storage2>&
                                            amap) : Base(amap)
    {}

    template <unsigned int d2, typename Refs2, typename Items2, typename Alloc2,
              typename Storage2, typename Converters>
    Polygonal_schema_with_combinatorial_map(const Combinatorial_map_base<d2, Refs2, Items2, Alloc2, Storage2>&
                                            amap, const Converters& converters) :
      Base(amap, converters)
    {}
    
    template <unsigned int d2, typename Refs2, typename Items2, typename Alloc2,
              typename Storage2, typename Converters,
              typename DartInfoConverter>
    Polygonal_schema_with_combinatorial_map(const Combinatorial_map_base<d2, Refs2, Items2, Alloc2, Storage2>&
                                            amap, const Converters& converters,
                                            const DartInfoConverter& dartinfoconverter) :
      Base(amap, converters, dartinfoconverter)
    {}

    template <unsigned int d2, typename Refs2, typename Items2, typename Alloc2,
              typename Storage2, typename Converters,
              typename DartInfoConverter, typename PointConverter >
    Polygonal_schema_with_combinatorial_map(const Combinatorial_map_base<d2, Refs2, Items2, Alloc2, Storage2>&
                                            amap, const Converters& converters,
                                            const DartInfoConverter& dartinfoconverter,
                                            const PointConverter& pointconverter) :
      Base(amap, converters, dartinfoconverter, pointconverter)
    {}
  };
  
  template <class Items_=Polygonal_schema_min_items,
            class Alloc_=CGAL_ALLOCATOR(int),
            class Storage_= Generalized_map_storage_1<2, Items_, Alloc_> >
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

    Polygonal_schema_with_generalized_map(const Self & amap) : Base(amap)
    {}

    template <unsigned int d2, typename Refs2, typename Items2, typename Alloc2,
              typename Storage2>
    Polygonal_schema_with_generalized_map(const Generalized_map_base<d2, Refs2, Items2, Alloc2, Storage2>&
                                            amap) : Base(amap)
    {}

    template <unsigned int d2, typename Refs2, typename Items2, typename Alloc2,
              typename Storage2, typename Converters>
    Polygonal_schema_with_generalized_map(const Generalized_map_base<d2, Refs2, Items2, Alloc2, Storage2>&
                                            amap, const Converters& converters) :
      Base(amap, converters)
    {}
    
    template <unsigned int d2, typename Refs2, typename Items2, typename Alloc2,
              typename Storage2, typename Converters,
              typename DartInfoConverter>
    Polygonal_schema_with_generalized_map(const Generalized_map_base<d2, Refs2, Items2, Alloc2, Storage2>&
                                            amap, const Converters& converters,
                                            const DartInfoConverter& dartinfoconverter) :
      Base(amap, converters, dartinfoconverter)
    {}

    template <unsigned int d2, typename Refs2, typename Items2, typename Alloc2,
              typename Storage2, typename Converters,
              typename DartInfoConverter, typename PointConverter >
    Polygonal_schema_with_generalized_map(const Generalized_map_base<d2, Refs2, Items2, Alloc2, Storage2>&
                                            amap, const Converters& converters,
                                            const DartInfoConverter& dartinfoconverter,
                                            const PointConverter& pointconverter) :
      Base(amap, converters, dartinfoconverter, pointconverter)
    {}
  };
  
} //namespace CGAL

#endif // CGAL_POLYGONAL_SCHEMA_H //
// EOF //
