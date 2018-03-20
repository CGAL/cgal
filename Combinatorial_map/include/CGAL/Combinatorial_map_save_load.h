// Copyright (c) 2010-2011 CNRS and LIRIS' Establishments (France).
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
//                 Guillaume Castano <guillaume.castano@gmail.com>
//                 Pascal Khieu <pascal.khieu@gmail.com>
//
#ifndef CGAL_COMBINATORIAL_MAP_SAVE_LOAD_H
#define CGAL_COMBINATORIAL_MAP_SAVE_LOAD_H

#include <boost/foreach.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Combinatorial_map_functors.h>

#include <algorithm>
#include <map>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <typeinfo>

namespace CGAL {

    typedef Exact_predicates_inexact_constructions_kernel::Point_2 RPoint_2;
    typedef Exact_predicates_exact_constructions_kernel::Point_2 EPoint_2;
    typedef Exact_predicates_inexact_constructions_kernel::Point_3 RPoint_3;
    typedef Exact_predicates_exact_constructions_kernel::Point_3 EPoint_3;

  // Tags used in xml tree:
  // For darts:
  //  <darts>
  //    <d> // new dart
  //        <b i="1"> neighbor dart index for beta1 </b>
  //        ...
  //        <v> value of dart (optional) </v>
  //    </d>
  //  ...
  // </darts>
  // For attributes:
  // <attributes>
  //   <dimension index="1"> // new type of non void attribute
  //   <type>type of the info associated</type>
  //   <a> // new attribute
  //    <d> dart index </d>
  //    <v> value of attribute </v>
  //   </a>
  //   ...
  // </attributes>

  // Here T is a Dart_const_handle so we don't need &
  template<typename T>
  void write_cmap_dart_node(boost::property_tree::ptree & /*node*/, T)
  {}

  template<typename T>
  void write_cmap_attribute_node(boost::property_tree::ptree & /*node*/, const T&)
  {}
  
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                       char val)
  {node.add("v",val);}
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                       unsigned char val)
  {node.add("v",val);}
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                       short int val)
  {node.add("v",val);}
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                       unsigned short int val)
  {node.add("v",val);}
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                       int val)
  {node.add("v",val);}
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                       unsigned int val)
  {node.add("v",val);}
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                       long int val)
  {node.add("v",val);}
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                       unsigned long int val)
  {node.add("v",val);}
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                       float val)
  {node.add("v",val);}
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                       double val)
  {node.add("v",val);}
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                       long double val)
  {node.add("v",val);}
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                       bool val)
  {node.add("v",val);}
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                       const std::string& val)
  {node.add("v",val);}
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                                 const RPoint_2& val)
  {
    node.add("p.x",val.x());
    node.add("p.y",val.y());
  }
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                                 const EPoint_2& val)
  {
    node.add("p.x",CGAL::to_double(val.x()));
    node.add("p.y",CGAL::to_double(val.y()));
  }
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                                 const RPoint_3& val)
  {
    node.add("p.x",val.x());
    node.add("p.y",val.y());
    node.add("p.z",val.z());
  }
  inline
  void write_cmap_attribute_node(boost::property_tree::ptree & node,
                                 const EPoint_3& val)
  {
    node.add("p.x",CGAL::to_double(val.x()));
    node.add("p.y",CGAL::to_double(val.y()));
    node.add("p.z",CGAL::to_double(val.z()));
  }

  template<typename CMap, unsigned int i,
           bool WithInfo=CGAL::Is_attribute_has_non_void_info
                          <typename CMap::template Attribute_type<i>::type>::value,
           bool WithPoint=CGAL::Is_attribute_has_point
                          <typename CMap::template Attribute_type<i>::type >::value >
  struct My_functor_cmap_save_one_attrib;

  // An attrib with point and with info
  template<typename CMap, unsigned int i>
  struct My_functor_cmap_save_one_attrib<CMap, i, true, true>
  {
    static void run(CMap& amap, boost::property_tree::ptree& ptree,
                    std::map<typename CMap::Dart_const_handle,
                               typename CMap::size_type>& myDarts)
    {
      // to check all i-cells of the map
      typename CMap::template Attribute_range<i>::type::const_iterator
        it_attrib, itend_attrib;
      it_attrib=amap.template attributes<i>().begin();
      itend_attrib=amap.template attributes<i>().end();

      // add dimension & type
      boost::property_tree::ptree& ndim = ptree.add("dimension", "");
      ndim.put("<xmlattr>.index", i);
      ndim.add("type", typeid(typename CMap::template Attribute_type<i>::type::Info).name());
      ndim.add("type_point", typeid(typename CMap::Point).name());

      // for every attribute of the dimension
      for (; it_attrib!=itend_attrib; ++it_attrib)
      {
        // make composant, dart and property node
        boost::property_tree::ptree & nattr = ndim.add("a", "");
        /* boost::property_tree::ptree & ndarts = */
          nattr.add("d", myDarts[it_attrib->dart()]);

        // update property node to add a value node (from basic or custom type
        write_cmap_attribute_node(nattr, it_attrib->info());
        write_cmap_attribute_node(nattr, it_attrib->point());
      }
    }
  };

  // An attribute with point and without info
  template<typename CMap, unsigned int i>
  struct My_functor_cmap_save_one_attrib<CMap, i, false, true>
  {
    static void run(CMap& amap, boost::property_tree::ptree& ptree,
                    std::map<typename CMap::Dart_const_handle,
                               typename CMap::size_type>& myDarts)
    {
      // to check all i-cells of the map
      typename CMap::template Attribute_range<i>::type::const_iterator
        it_attrib, itend_attrib;
      it_attrib=amap.template attributes<i>().begin();
      itend_attrib=amap.template attributes<i>().end();

      // add dimension & type
      boost::property_tree::ptree& ndim = ptree.add("dimension", "");
      ndim.put("<xmlattr>.index", i);
      ndim.add("type", "void");
      ndim.add("type_point", typeid(typename CMap::Point).name());

      // for every attribute of the dimension
      for (; it_attrib!=itend_attrib; ++it_attrib)
      {
        // make composant, dart and property node
        boost::property_tree::ptree & nattr = ndim.add("a", "");
        /* boost::property_tree::ptree & ndarts = */
        nattr.add("d", myDarts[it_attrib->dart()]);

        // update property node to add a value node (from basic or custom type
        write_cmap_attribute_node(nattr, it_attrib->point());
      }
    }
  };

  // An attribute without point and with info
  template<typename CMap, unsigned int i>
  struct My_functor_cmap_save_one_attrib<CMap, i, true, false>
  {
    static void run(CMap& amap, boost::property_tree::ptree& ptree,
                    std::map<typename CMap::Dart_const_handle,
                               typename CMap::size_type>& myDarts)
    {
      // to check all i-cells of the map
      typename CMap::template Attribute_range<i>::type::const_iterator
        it_attrib, itend_attrib;
      it_attrib=amap.template attributes<i>().begin();
      itend_attrib=amap.template attributes<i>().end();

      // add dimension & type
      boost::property_tree::ptree& ndim = ptree.add("dimension", "");
      ndim.put("<xmlattr>.index", i);
      ndim.add("type", typeid(typename CMap::template
                              Attribute_type<i>::type::Info).name());
      ndim.add("type_point", "void");

      // for every attribute of the dimension
      for (; it_attrib!=itend_attrib; ++it_attrib)
      {
        // make composant, dart and property node
        boost::property_tree::ptree & nattr = ndim.add("a", "");
        /* boost::property_tree::ptree & ndarts = */
        nattr.add("d", myDarts[it_attrib->dart()]);

        // update property node to add a value node (from basic or custom type
        write_cmap_attribute_node(nattr, it_attrib->info());
      }
    }
  };

  // An attrib without point and without info
  template<typename CMap, unsigned int i>
  struct My_functor_cmap_save_one_attrib<CMap, i, false, false>
  {
    static void run(CMap& amap, boost::property_tree::ptree& ptree,
                    std::map<typename CMap::Dart_const_handle,
                               typename CMap::size_type>& myDarts)
    {
      // to check all i-cells of the map
      typename CMap::template Attribute_range<i>::type::const_iterator
        it_attrib, itend_attrib;
      it_attrib=amap.template attributes<i>().begin();
      itend_attrib=amap.template attributes<i>().end();

      // add dimension & type
      boost::property_tree::ptree& ndim = ptree.add("dimension", "");
      ndim.put("<xmlattr>.index", i);
      ndim.add("type", "void");
      ndim.add("type_point", "void");

      // for every attribute of the dimension
      for (; it_attrib!=itend_attrib; ++it_attrib)
      {
        // make composant, dart and property node
        boost::property_tree::ptree & nattr = ndim.add("a", "");
        /* boost::property_tree::ptree & ndarts = */
        nattr.add("d", myDarts[it_attrib->dart()]);
      }
    }
  };

  template<typename CMap>
  struct My_functor_cmap_save_attrib
  {
    template <unsigned int i>
    static void run(CMap& amap, boost::property_tree::ptree& ptree,
                    std::map<typename CMap::Dart_const_handle,
                               typename CMap::size_type>& myDarts)
    {
      My_functor_cmap_save_one_attrib<CMap, i>::run(amap, ptree, myDarts);
    }
  };

  template < class CMap >
  boost::property_tree::ptree cmap_save_darts
  (CMap& amap, std::map<typename CMap::Dart_const_handle,
                                typename CMap::size_type>& myDarts)
  {
    CGAL_assertion( myDarts.empty() );
    
    // First we numbered each dart by using the std::map.
    typename CMap::Dart_range::const_iterator it(amap.darts().begin());
    for(typename CMap::size_type num=1; num<=amap.number_of_darts();
        ++num, ++it)
    {
      myDarts[it] = num;
    }

    // make a tree
    using boost::property_tree::ptree;
    ptree pt;

    // Now we save each dart, and its neighbors.
    it=amap.darts().begin();
    for(typename CMap::size_type num=0; num<amap.number_of_darts(); ++num, ++it)
    {
      // make a dart node
      ptree& ndart = pt.add("d", "");

      // the beta, only for non free sews
      for(unsigned int dim=1; dim<=amap.dimension; dim++)
      {
        if(!amap.is_free(it, dim))
        {
          ptree& currentNext = ndart.add("b", myDarts[amap.beta(it, dim)]);
          currentNext.put("<xmlattr>.i", dim);
        }
      }

      // update property node to add a value node (if user defined its own
      // function)
      write_cmap_dart_node(ndart, it);
    }
    
    return pt;
  }

  template < class CMap >
  boost::property_tree::ptree cmap_save_attributes
  (const CMap& amap, std::map<typename CMap::Dart_const_handle,
                                typename CMap::size_type>& myDarts)
  {
    using boost::property_tree::ptree;
    ptree pt;

    // update pt adding nodes containing attributes informations
    CMap::Helper::template Foreach_enabled_attributes
      <My_functor_cmap_save_attrib<CMap> >::run(const_cast<CMap&>(amap), pt, myDarts);

    return pt;
  }

  struct EmptyFunctor
  {
    void operator() (boost::property_tree::ptree & /*node*/) const
    {
      // node.add("myinfo.myvalie",15);
    }
  };

  template < class CMap, class Functor >
  bool save_combinatorial_map(const CMap& amap, std::ostream & output,
                              const Functor& f)
  {
    using boost::property_tree::ptree;
    ptree tree;
    tree.put("data", "");

    /** First we save general information of the map (by default nothing,
        the fuction can be specialized by users). */
    f(tree);

    // map dart => number
    std::map<typename CMap::Dart_const_handle, typename CMap::size_type> myDarts;

    // Save darts
    ptree pt_darts=cmap_save_darts(amap, myDarts);
    tree.add_child("data.darts",pt_darts);

    // Save attributes
    ptree pt_attr=cmap_save_attributes(amap, myDarts);
    tree.add_child("data.attributes", pt_attr);

    // save data in output
    write_xml(output, tree);

    return true;
  }

  template < class CMap, class Functor >
  bool save_combinatorial_map(const CMap& amap, const char* filename,
                              const Functor& f)
  {
    std::ofstream output(filename);
    if (!output) return false;
    return save_combinatorial_map(amap, output, f);
  }
  
  template < class CMap >
  bool save_combinatorial_map(const CMap& amap, std::ostream & output)
  {
    EmptyFunctor f;
    return save_combinatorial_map(amap, output, f);
  }

  template < class CMap >
  bool save_combinatorial_map(const CMap& amap, const char* filename)
  {
    EmptyFunctor f;
    return save_combinatorial_map(amap, filename, f);
  }

  // Here T is a Dart_handle so no need of &
  template<typename T>
  void read_cmap_dart_node
  (const boost::property_tree::ptree::value_type &/*v*/, T /*val*/)
  {}
  template<typename T>
  void read_cmap_attribute_node
  (const boost::property_tree::ptree::value_type &/*v*/, T &/*val*/)
  {}
  template<> inline
  void read_cmap_attribute_node
  (const boost::property_tree::ptree::value_type &v,char &val)
  {val=boost::lexical_cast< char >(v.second.data());}
  template<> inline
  void read_cmap_attribute_node
  (const boost::property_tree::ptree::value_type &v,unsigned char &val)
  {val=boost::lexical_cast< unsigned char >(v.second.data());}
  template<> inline
  void read_cmap_attribute_node
  (const boost::property_tree::ptree::value_type &v,short int &val)
  {val=boost::lexical_cast< short int >(v.second.data());}
  template<> inline
  void read_cmap_attribute_node
  (const boost::property_tree::ptree::value_type &v,unsigned short int &val)
  {val=boost::lexical_cast< unsigned short int >(v.second.data());}
  template<> inline
  void read_cmap_attribute_node
  (const boost::property_tree::ptree::value_type &v,int &val)
  {val=boost::lexical_cast< int >(v.second.data());}
  template<> inline
  void read_cmap_attribute_node
  (const boost::property_tree::ptree::value_type &v,unsigned int &val)
  {val=boost::lexical_cast< unsigned int >(v.second.data());}
  template<> inline
  void read_cmap_attribute_node
  (const boost::property_tree::ptree::value_type &v,long int &val)
  {val=boost::lexical_cast< long int >(v.second.data());}
  template<> inline
  void read_cmap_attribute_node
  (const boost::property_tree::ptree::value_type &v,unsigned long int &val)
  {val=boost::lexical_cast< unsigned long int >(v.second.data());}
  template<> inline
  void read_cmap_attribute_node
  (const boost::property_tree::ptree::value_type &v,float &val)
  {val=boost::lexical_cast< float >(v.second.data());}
  template<> inline
  void read_cmap_attribute_node
  (const boost::property_tree::ptree::value_type &v,double &val)
  {val=boost::lexical_cast< double >(v.second.data());}
  template<> inline
  void read_cmap_attribute_node
  (const boost::property_tree::ptree::value_type &v,long double &val)
  {val=boost::lexical_cast< long double >(v.second.data());}
  template<> inline
  void read_cmap_attribute_node
  (const boost::property_tree::ptree::value_type &v,bool &val)
  {val=boost::lexical_cast< bool >(v.second.data());}
  template<> inline
  void read_cmap_attribute_node
  (const boost::property_tree::ptree::value_type &v,std::string &val)
  {val=boost::lexical_cast< std::string >(v.second.data());}
  template<> inline
  void read_cmap_attribute_node
  (const boost::property_tree::ptree::value_type &v,RPoint_2 &val)
  {
    double x=v.second.get<double>("x");
    double y=v.second.get<double>("y");
    val = RPoint_2(x,y);
  }
  template<> inline
  void read_cmap_attribute_node
  (const boost::property_tree::ptree::value_type &v,RPoint_3 &val)
  {
    double x=v.second.get<double>("x");
    double y=v.second.get<double>("y");
    double z=v.second.get<double>("z");
    val = RPoint_3(x,y,z);
  }

  template<typename CMap, unsigned int i,
           bool WithInfo=CGAL::Is_attribute_has_non_void_info
                          <typename CMap::template Attribute_type<i>::type>::value,
           bool WithPoint=CGAL::Is_attribute_has_point
                          <typename CMap::template Attribute_type<i>::type >::value >
  struct My_functor_cmap_load_one_attrib;

  // An attrib with point and with info
  template<typename CMap, unsigned int i>
  struct My_functor_cmap_load_one_attrib<CMap, i, true, true>
  {
    static void run(const boost::property_tree::ptree& pt, CMap& amap,
                    const std::vector<typename CMap::Dart_handle>& myDarts)
    {
      BOOST_FOREACH( const boost::property_tree::ptree::value_type &v0,
                     pt.get_child("data.attributes") )
      {
        // <dimension>
        if (v0.first == "dimension")
        {
          int dimension=v0.second.get("<xmlattr>.index", -1);
          
          // if map.dimension == dimension saved in the xml file
          if (dimension==i)
          {
            unsigned int id_dart_cellule=0;
            std::string type =  v0.second.get<std::string>("type");
            std::string type_map=std::string
              (typeid(typename CMap::template Attribute_type<i>::type::Info).name());
            
            std::string ptype =  v0.second.get<std::string>("type_point");
            std::string ptype_map= std::string
              (typeid(typename CMap::template Attribute_type<i>::type::Point).name());
            
                //  std::cout<<"ptype="<<ptype<<"  and type_map="<<type_map<<std::endl;
                /* if(type!=type_map && ptype!=ptype_map)
                {
                  //  std::cout<<"Not loaded."<<std::endl;
                  return;
                  }*/

            BOOST_FOREACH(const boost::property_tree::ptree::value_type &v1,
                          v0.second )
            {
              if( v1.first == "a" )
              {
                id_dart_cellule=v1.second.get<unsigned int>("d")-1;
                
                BOOST_FOREACH(const boost::property_tree::ptree::value_type &v2,
                              v1.second )
                {
                  if( type==type_map && v2.first == "v" )
                  {
                    if (amap.template attribute<i>(myDarts[id_dart_cellule])
                        ==NULL )
                      amap.template set_attribute<i>
                        (myDarts[id_dart_cellule],
                         amap.template create_attribute<i>());
                    read_cmap_attribute_node
                      (v2,
                       amap.template info<i>(myDarts[id_dart_cellule]));
                  }
                  if( ptype==ptype_map && v2.first == "p" )
                  {
                    if (amap.template attribute<i>(myDarts[id_dart_cellule])
                        ==NULL )
                      amap.template set_attribute<i>
                        (myDarts[id_dart_cellule],
                         amap.template create_attribute<i>());
                    read_cmap_attribute_node
                      (v2,
                       amap.template attribute<i>(myDarts[id_dart_cellule])
                       ->point());
                  }
                }
              }
            }
          }
        }
      }
    }
  };

  // An attribute with point and without info
  template<typename CMap, unsigned int i>
  struct My_functor_cmap_load_one_attrib<CMap, i, false, true>
  {
    static void run(const boost::property_tree::ptree& pt, CMap& amap,
                    const std::vector<typename CMap::Dart_handle>& myDarts)
    {
      BOOST_FOREACH( const boost::property_tree::ptree::value_type &v0,
                     pt.get_child("data.attributes") )
      {
        // <dimension>
        if (v0.first == "dimension")
        {
          int dimension=v0.second.get("<xmlattr>.index", -1);
          
          // if map.dimension == dimension saved in the xml file
          if (dimension==i)
          {
            unsigned int id_dart_cellule=0;
            std::string ptype =  v0.second.get<std::string>("type_point");
            std::string type_map= typeid
              (typename CMap::template Attribute_type<i>::type::Point).name();
                //  std::cout<<"ptype="<<ptype<<"  and type_map="<<type_map<<std::endl;
                /*                if(ptype!=type_map)
                {
                  //  std::cout<<"Not loaded."<<std::endl;
                  return;
                  }*/

            BOOST_FOREACH(const boost::property_tree::ptree::value_type &v1,
                          v0.second )
            {
              if( v1.first == "a" )
              {
                id_dart_cellule=v1.second.get<unsigned int>("d")-1;
                
                BOOST_FOREACH(const boost::property_tree::ptree::value_type &v2,
                              v1.second )
                {
                  if( v2.first == "p" )
                  {
                    if (amap.template attribute<i>
                        (myDarts[id_dart_cellule])==NULL )
                      amap.template set_attribute<i>
                        (myDarts[id_dart_cellule],
                         amap.template create_attribute<i>());
                    
                    read_cmap_attribute_node
                      (v2,
                       (amap.template attribute<i>
                        (myDarts[id_dart_cellule])->point()));
                  }
                }
              }
            }
          }
        }
      }
    }
  };

  // An attribute without point and with info
  template<typename CMap, unsigned int i>
  struct My_functor_cmap_load_one_attrib<CMap, i, true, false>
  {
    static void run(const boost::property_tree::ptree& pt, CMap& amap,
                    const std::vector<typename CMap::Dart_handle>& myDarts)
    {
      BOOST_FOREACH( const boost::property_tree::ptree::value_type &v0,
                     pt.get_child("data.attributes") )
      {
        // <dimension>
        if (v0.first == "dimension")
        {
          int dimension=v0.second.get("<xmlattr>.index", -1);
          
          // if map.dimension == dimension saved in the xml file
          if (dimension==i)
          {
            unsigned int id_dart_cellule=0;
            std::string ptype =  v0.second.get<std::string>("type");
            std::string type_map= typeid
              (typename CMap::template Attribute_type<i>::type::Info).name();
                //  std::cout<<"ptype="<<ptype<<"  and type_map="<<type_map<<std::endl;
                /*      if(ptype!=type_map)
                {
                  //  std::cout<<"Not loaded."<<std::endl;
                  return;
                  } */

            BOOST_FOREACH(const boost::property_tree::ptree::value_type &v1,
                          v0.second )
            {
              if( v1.first == "a" )
              {
                id_dart_cellule=v1.second.get<unsigned int>("d")-1;
                
                BOOST_FOREACH(const boost::property_tree::ptree::value_type &v2,
                              v1.second )
                {
                  if( v2.first == "v" )
                  {
                    if (amap.template attribute<i>
                        (myDarts[id_dart_cellule])==NULL)
                      amap.template set_attribute<i>
                        (myDarts[id_dart_cellule],
                         amap.template create_attribute<i>());
                    read_cmap_attribute_node
                      (v2,
                       amap.template info<i>(myDarts[id_dart_cellule]));
                  }
                }
              }
            }
          }
        }
      }
    }
  };

  // An attribute without point and without info
  template<typename CMap, unsigned int i>
  struct My_functor_cmap_load_one_attrib<CMap, i, false, false>
  {
    static void run(const boost::property_tree::ptree& pt, CMap& amap,
                    const std::vector<typename CMap::Dart_handle>& myDarts)
    {
      BOOST_FOREACH( const boost::property_tree::ptree::value_type &v0,
                     pt.get_child("data.attributes") )
      {
        // <dimension>
        if (v0.first == "dimension")
        {
          int dimension=v0.second.get("<xmlattr>.index", -1);
          
          // if map.dimension == dimension saved in the xml file
          if (dimension==i)
          {
            unsigned int id_dart_cellule=0;
            
            BOOST_FOREACH(const boost::property_tree::ptree::value_type &v1,
                          v0.second )
            {
              if( v1.first == "a" )
              {
                id_dart_cellule=v1.second.get<unsigned int>("d")-1;
                
                if (amap.template attribute<i>(myDarts[id_dart_cellule])==NULL)
                  amap.template set_attribute<i>
                    (myDarts[id_dart_cellule],
                   amap.template create_attribute<i>());
              }
            }
          }
        }
      }
    }
  };

  /** Functor called to load i-attributes.
   *  @param pt a boost::property_tree::ptree load from an xml file
   *  @param amap a pointer to the map to load into
   *  @param myDarts an array of Dart_handle st myDarts[i] is the ith dart.
   */
  template<class CMap>
  struct My_functor_cmap_load_attrib
  {
    template <unsigned int i>
    static void run(const boost::property_tree::ptree& pt, CMap& amap,
                    const std::vector<typename CMap::Dart_handle>& myDarts)
    {
       My_functor_cmap_load_one_attrib<CMap, i>::run(pt, amap, myDarts);
    }
  };

  template < class CMap >
  bool cmap_load_darts(boost::property_tree::ptree &pt, CMap& amap,
                       std::vector<typename CMap::Dart_handle>& myDarts)
  {
    // use a boost::property_tree
    using boost::property_tree::ptree;

    // make darts
    BOOST_FOREACH( const ptree::value_type &v, pt.get_child("data.darts") )
    {
      if( v.first == "d" )
        myDarts.push_back(amap.create_dart());
    }

    // update beta links
    unsigned int index;
    unsigned int currentDartInt = 0;
    unsigned int nextDartInt;

    BOOST_FOREACH( const ptree::value_type &v, pt.get_child("data.darts") )
    {
      if( v.first == "d" )
      {
        BOOST_FOREACH( const ptree::value_type &v2, v.second )
        {
          if (v2.first == "b")
          {
            index = v2.second.get("<xmlattr>.i", 0);
            nextDartInt = boost::lexical_cast< int >(v2.second.data())-1;
            
            if ( index<=amap.dimension )
            {
              amap.basic_link_beta(myDarts[currentDartInt],
                                   myDarts[nextDartInt],
                                   index);
            }
          }
          else if (v2.first=="v")
            read_cmap_dart_node(v2,myDarts[currentDartInt]);
        }
      }
      ++currentDartInt;
    }

    return true;
  }

  template < class CMap >
  void cmap_load_attributes(const boost::property_tree::ptree& pt, CMap& amap,
                            const std::vector<typename CMap::Dart_handle>& myDarts)
  {
    CMap::Helper::template Foreach_enabled_attributes
      <My_functor_cmap_load_attrib<CMap> >::run(pt, amap, myDarts);
  }

  template < class CMap, class Functor >
  bool load_combinatorial_map(std::ifstream & input, CMap& amap,
                              Functor& f)
  {
    using boost::property_tree::ptree;
    ptree pt;
    read_xml(input, pt);

    /** First we load general information of the map (by default nothing,
        the fuction can be specialized by users). */
    f(pt);

    // Then we load darts and attributes.
    std::vector<typename CMap::Dart_handle> myDarts;
    cmap_load_darts(pt,amap,myDarts);
    cmap_load_attributes(pt,amap,myDarts);
    return true;
  }
  
  template < class CMap, class Functor >
  bool load_combinatorial_map(const char* filename, CMap& amap,
                              Functor& f)
  {
    std::ifstream input(filename);
    if (!input) return false;
    return load_combinatorial_map(input, amap, f);
  }

  template < class CMap >
  bool load_combinatorial_map(std::ifstream & input, CMap& amap)
  {
    EmptyFunctor f;
    return load_combinatorial_map(input, amap, f);
  }

  template < class CMap >
  bool load_combinatorial_map(const char* filename, CMap& amap)
  {
    EmptyFunctor f;
    return load_combinatorial_map(filename, amap, f);
  }
} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_SAVE_LOAD_H //
// EOF //
