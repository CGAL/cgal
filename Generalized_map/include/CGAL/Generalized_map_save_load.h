// Copyright (c) 2016 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_GENERALIZED_MAP_SAVE_LOAD_H
#define CGAL_GENERALIZED_MAP_SAVE_LOAD_H

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Combinatorial_map_save_load.h>

#include <algorithm>
#include <map>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <typeinfo>

/* We reuse the following functions from Combinatorial_map_save_load.h:
 *   write_cmap_dart_node
 *   write_cmap_attribute_node
 *   My_functor_cmap_save_one_attrib
 *   My_functor_cmap_save_attrib
 *   cmap_save_attributes
 *
 *   read_cmap_dart_node
 *   read_cmap_attribute_node
 *   My_functor_cmap_load_one_attrib
 *   My_functor_cmap_load_attrib
 *   cmap_load_attributes
 */

namespace CGAL {

    typedef Exact_predicates_inexact_constructions_kernel::Point_3 RPoint_3;
    typedef Exact_predicates_exact_constructions_kernel::Point_3 EPoint_3;

  // Tags used in xml tree:
  // For darts:
  //  <darts>
  //    <d> // new dart
  //        <a i="0"> neighbor dart index for alpha0 </b>
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

  template < class GMap >
  boost::property_tree::ptree gmap_save_darts
  (const GMap& amap,
   std::map<typename GMap::Dart_const_handle,
              typename GMap::size_type>& myDarts)
  {
    CGAL_assertion( myDarts.empty() );

    // First we numbered each dart by using the std::map.
    typename GMap::Dart_range::const_iterator it(amap.darts().begin());
    for(typename GMap::size_type num=1; num<=amap.number_of_darts();
        ++num, ++it)
    {
      myDarts[it] = num;
    }

    // make a tree
    using boost::property_tree::ptree;
    ptree pt;

    // Now we save each dart, and its neighbors.
    it=amap.darts().begin();
    for(typename GMap::size_type num=0; num<amap.number_of_darts(); ++num, ++it)
    {
      // make a dart node
      ptree& ndart = pt.add("d", "");

      // the beta, only for non free sews
      for(unsigned int dim=0; dim<=amap.dimension; dim++)
      {
        if(!amap.is_free(it, dim))
        {
          ptree& currentNext = ndart.add("a", myDarts[amap.alpha(it, dim)]);
          currentNext.put("<xmlattr>.i", dim);
        }
      }

      // update property node to add a value node (if user defined its own
      // function)
      write_cmap_dart_node(ndart, it);
    }

    return pt;
  }

  template < class GMap >
  bool save_generalized_map(const GMap& amap, std::ostream & output)
  {
    using boost::property_tree::ptree;
    ptree tree;

    // map dart => number
    std::map<typename GMap::Dart_const_handle, typename GMap::size_type> myDarts;

    // Save darts
    ptree pt_darts=gmap_save_darts(amap, myDarts);
    tree.add_child("data.darts",pt_darts);

    // Save attributes
    ptree pt_attr=cmap_save_attributes(amap, myDarts);
    tree.add_child("data.attributes", pt_attr);

    // save data in output
    write_xml(output, tree);

    return true;
  }

  template < class GMap >
  bool save_generalized_map(const GMap& amap, const char* filename)
  {
    std::ofstream output(filename);
    if (!output) return false;
    return save_generalized_map(amap, output);
  }

  template < class GMap >
  bool gmap_load_darts(boost::property_tree::ptree &pt, GMap& amap,
                       std::vector<typename GMap::Dart_handle>& myDarts)
  {
    // use a boost::property_tree
    using boost::property_tree::ptree;

    // make darts
    for(const ptree::value_type& v : pt.get_child("data.darts") )
    {
      if( v.first == "d" )
        myDarts.push_back(amap.create_dart());
    }

    // update beta links
    unsigned int index;
    unsigned int currentDartInt = 0;
    unsigned int nextDartInt;

    for(const ptree::value_type& v : pt.get_child("data.darts") )
    {
      if( v.first == "d" ) // d for dart
      {
        for(const ptree::value_type& v2 : v.second )
        {
          if (v2.first == "a") // a for alpha
          {
            index = v2.second.get("<xmlattr>.i", 0);
            nextDartInt = boost::lexical_cast< int >(v2.second.data())-1;

            if ( index<=amap.dimension )
            {
              amap.basic_link_alpha(myDarts[currentDartInt],
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

  template < class GMap >
  bool load_generalized_map(std::ifstream & input, GMap& amap)
  {
    using boost::property_tree::ptree;
    ptree pt;
    read_xml(input, pt);
    std::vector<typename GMap::Dart_handle> myDarts;
    gmap_load_darts(pt,amap,myDarts);
    cmap_load_attributes(pt,amap,myDarts);
    return true;
  }

  template < class GMap >
  bool load_generalized_map(const char* filename, GMap& amap)
  {
    std::ifstream input(filename);
    if (!input) return false;
    return load_generalized_map(input, amap);
  }

} // namespace CGAL

#endif // CGAL_GENERALIZED_MAP_SAVE_LOAD_H //
// EOF //
