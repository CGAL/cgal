// Copyright (c) 2010-2011 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_COMBINATORIAL_MAP_BASIC_OPERATIONS_H
#define CGAL_COMBINATORIAL_MAP_BASIC_OPERATIONS_H 1

#include <CGAL/Combinatorial_map_iterators_base.h>
#include <CGAL/tags.h>

#include <boost/type_traits/is_same.hpp>

namespace CGAL
{
  /** @file Combinatorial_map_basic_operations.h
   * Basic operations on a combinatorial map.
   */

  /** Test if the two given darts belong to the same given orbit.
   * @param amap a combinatorial map.
   * @param adart1 the first dart.
   * @param adart2 the second dart.
   * @return true iff the two darts belong to the same orbit.
   */
  template < class Map, class Iterator >
  bool belong_to_same_orbit(const Map & amap,
                            typename Map::Dart_const_handle adart1,
                            typename Map::Dart_const_handle adart2)
  {
    CGAL_static_assertion( (boost::is_same<typename Iterator::Basic_iterator,
                                           Tag_false>::value) );
    bool found=false;

    for (Iterator it(amap, adart1); !found && it.cont(); ++it)
    {
      if (it == adart2) found = true;
    }

    return found;
  }

  /** Test if all the darts of a given orbit are marked.
   * @param amap a combinatorial map.
   * @param adart a dart of the orbit.
   * @param amark the mark.
   * @return true iff all the darts are marked.
   */
  template < class Map, class Iterator >
  bool is_whole_orbit_marked(const Map & amap,
                             typename Map::Dart_const_handle adart,
                             typename Map::size_type amark)
  {
    CGAL_static_assertion( (boost::is_same<typename Iterator::Basic_iterator,
                                           Tag_false>::value) );
    bool res=true;

    for ( Iterator it(amap, adart); res && it.cont(); ++it )
    {
      if (!amap.is_marked(it, amark)) res = false;
    }

    return res;
  }

  /** Test if all the darts of a given orbit are unmarked.
   * @param amap a combinatorial map.
   * @param adart a dart of the orbit.
   * @param amark the mark.
   * @return true iff all the darts are unmarked.
   */
  template < class Map, class Iterator >
  bool is_whole_orbit_unmarked(const Map & amap,
                               typename Map::Dart_const_handle adart,
                               typename Map::size_type amark)
  {
    amap.negate_mark(amark);
    bool res=CGAL::template is_whole_orbit_marked<Map,Iterator>(amap, adart, amark);
    amap.negate_mark(amark);
    return res;
  }

  /** Mark a given orbit with a given mark.
   * @param amap a combinatorial map.
   * @param adart a dart of the orbit.
   * @param amark the mark.
   * @return the number of darts of the marked orbit.
   * @pre The whole orbit must be unmarked.
   */
  template < class Map, class Iterator >
  typename Map::size_type mark_orbit(const Map & amap,
                                     typename Map::Dart_const_handle adart,
                                     typename Map::size_type amark)
  {
    CGAL_static_assertion( (boost::is_same<typename Iterator::Basic_iterator,
                            Tag_true>::value) );
    CGAL_assertion( (is_whole_orbit_unmarked<Map,
                     CMap_non_basic_iterator<Map,Iterator> >
                     (amap, adart, amark)) );
    typename Map::size_type res=0;

    for (Iterator it(amap, adart, amark); it.cont(); ++it)
    {
      amap.mark(it, amark);
      ++res;
    }

    return res;
  }

  /** Unmark a given orbit with a given mark.
   * @param amap a combinatorial map.
   * @param adart a dart of the orbit.
   * @param amark the mark.
   * @return the number of darts of the unmarked orbit.
   * @pre The whole orbit must be marked.
   */
  template < class Map, class Iterator >
  typename Map::size_type unmark_orbit(const Map & amap,
                                       typename Map::Dart_const_handle adart,
                                       typename Map::size_type amark)
  {
    amap.negate_mark(amark);
    typename Map::size_type
        res=CGAL::template mark_orbit<Map, Iterator>(amap, adart, amark);
    amap.negate_mark(amark);
    return res;
  }

  /** Test if the two given darts belong to the same cell.
   * @param amap a combinatorial map.
   * @param adart1 the first dart.
   * @param adart2 the second dart.
   * @return true iff the two darts belong to the same cell.
   */
  template < class Map, unsigned int i, unsigned int d>
  bool belong_to_same_cell(const Map & amap,
                           typename Map::Dart_const_handle adart1,
                           typename Map::Dart_const_handle adart2)
  {
    return CGAL::belong_to_same_orbit<Map,
        typename Map::template Dart_of_cell_range<i,d>::const_iterator>
        (amap, adart1, adart2);
  }

  template < class Map, unsigned int i>
  bool belong_to_same_cell(const Map & amap,
                           typename Map::Dart_const_handle adart1,
                           typename Map::Dart_const_handle adart2)
  {
    return CGAL::template belong_to_same_cell<Map,i,Map::dimension>(amap,adart1,adart2);
  }


  /** Test if all the darts of a given cell are marked.
   * @param amap a combinatorial map.
   * @param adart a dart of the cell.
   * @param amark the mark.
   * @return true iff all the darts are marked.
   */
  template < class Map, unsigned int i, unsigned int d>
  bool is_whole_cell_marked(const Map & amap,
                            typename Map::Dart_const_handle adart,
                            typename Map::size_type amark)
  {
    return CGAL::template is_whole_orbit_marked<Map,
        typename Map::template Dart_of_cell_range<i,d>::const_iterator>
        (amap, adart, amark);
  }

  template < class Map, unsigned int i>
  bool is_whole_cell_marked(const Map & amap,
                            typename Map::Dart_const_handle adart,
                            typename Map::size_type amark)
  {
    return CGAL::template is_whole_cell_marked<Map,i,Map::dimension>
      (amap,adart,amark);
  }

  /** Test if all the darts of a given cell are unmarked.
   * @param amap a combinatorial map.
   * @param adart a dart of the cell.
   * @param amark the mark.
   * @return true iff all the darts are marked.
   */
  template < class Map, unsigned int i, unsigned int d >
  bool is_whole_cell_unmarked(const Map & amap,
                              typename Map::Dart_const_handle adart,
                              typename Map::size_type amark)
  {
    return CGAL::template is_whole_orbit_unmarked<Map,
        typename Map::template Dart_of_cell_range<i,d>::const_iterator>
        (amap, adart, amark);
  }

  template < class Map, unsigned int i>
  bool is_whole_cell_unmarked(const Map & amap,
                              typename Map::Dart_const_handle adart,
                              typename Map::size_type amark)
  {
    return CGAL::template is_whole_cell_unmarked<Map,i,Map::dimension>
        (amap,adart,amark);
  }

  /** Mark a given cell with a given mark.
   * @param amap a combinatorial map.
   * @param adart a dart of the cell.
   * @param amark the mark.
   * @return the number of darts of the marked cell.
   * @pre The whole cell must be unmarked.
   */
  template < class Map, unsigned int i, unsigned int d >
  typename Map::size_type mark_cell(const Map & amap,
                                    typename Map::Dart_const_handle adart,
                                    typename Map::size_type amark)
  { return CGAL::template mark_orbit<Map,
        typename Map::template Dart_of_cell_basic_range<i,d>::const_iterator>
        (amap, adart, amark); }

  template < class Map, unsigned int i>
  typename Map::size_type mark_cell(const Map & amap,
                                    typename Map::Dart_const_handle adart,
                                    typename Map::size_type amark)
  { return CGAL::template mark_cell<Map,i,Map::dimension>(amap, adart, amark);}

  /** Unmark a given cell with a given mark.
   * @param amap a combinatorial map.
   * @param adart a dart of the cell.
   * @param amark the mark.
   * @return the number of darts of the unmarked cell.
   * @pre The whole cell must be marked.
   */
  template < class Map, unsigned int i, unsigned int d >
  typename Map::size_type unmark_cell(const Map & amap,
                                      typename Map::Dart_const_handle adart,
                                      typename Map::size_type amark)
  { return CGAL::template unmark_orbit<Map,
        typename Map::template Dart_of_cell_basic_range<i,d>::const_iterator>
        (amap, adart, amark);}

  template < class Map, unsigned int i >
  typename Map::size_type unmark_cell(const Map & amap,
                                      typename Map::Dart_const_handle adart,
                                      typename Map::size_type amark)
  { return CGAL::template unmark_cell<Map,i,Map::dimension>(amap, adart, amark); }

  /// Functor to test if the given dart must be marked, when we want
  /// to mark an oriented orbit. By default (for CMap), return true.
  template<class Map, class T=typename Map::Combinatorial_data_structure>
  struct Need_to_mark_for_oriented_cell
  {
    static bool run(const Map & /*amap*/,
                    typename Map::Dart_const_handle /*adart*/,
                    typename Map::size_type /*amark*/,
                    OperationState /*op*/)
    { return true; }
  };

  struct Generalized_map_tag;

  /// Specialization for GMap. We mark one out of two darts.
  template<class Map>
  struct Need_to_mark_for_oriented_cell<Map, CGAL::Generalized_map_tag>
  {
    static bool run(const Map & amap,
                    typename Map::Dart_const_handle adart,
                    typename Map::size_type amark,
                    OperationState op)
    {
      CGAL_assertion(op!=OP_END);

      if (op==OP_NONE) return true;

      bool all_darts_unmarked=true;
      for (int i=0; all_darts_unmarked && i<=(int)Map::dimension; ++i)
      {
        if (!amap.is_free(adart, i) &&
            amap.is_marked(amap.alpha(adart, i), amark))
        { all_darts_unmarked=false; }
      }
      return all_darts_unmarked;
    }
  };

  /** Mark a given oriented orbit with a given mark. For combinatorial map,
   *  equivalent to mark_orbit. For generalized map, mark one out of two dart
   *  if the orbit is orientable. If amark2!=INVALID_MARK, mark also all the darts
   *  of the orbit with amark2.
   * @param amap a combinatorial map / generalized.
   * @param adart a dart of the orbit.
   * @param amark the mark.
   * @param amark2 a second mark.
   * @return the number of darts of the marked orbit.
   * @pre The whole orbit must be unmarked.
   */
  template < class Map, class Iterator >
  typename Map::size_type mark_oriented_orbit(const Map & amap,
                                              typename Map::Dart_const_handle adart,
                                              typename Map::size_type amark,
                                              typename Map::size_type amark2=Map::INVALID_MARK)
  {
    CGAL_static_assertion( (boost::is_same<typename Iterator::Basic_iterator,
                            Tag_true>::value) );
    CGAL_assertion( (is_whole_orbit_unmarked<Map,
                     CMap_non_basic_iterator<Map,Iterator> >
                     (amap, adart, amark)) );
    typename Map::size_type res=0;
    typename Map::size_type amark3=(amark2==Map::INVALID_MARK?amap.get_new_mark():amark2);
    for (Iterator it(amap, adart, amark3); it.cont(); ++it)
    {
      amap.mark(it, amark3);
      if (CGAL::template Need_to_mark_for_oriented_cell<Map>::
          run(amap, it, amark, it.prev_operation()))
      {
        amap.mark(it, amark);
        ++res;
      }
    }

    if (amark2==Map::INVALID_MARK)
    {
      CGAL::template unmark_orbit<Map, Iterator>(amap, adart, amark3);
      amap.free_mark(amark3);
    }

    return res;
  }

  template < class Map, class Iterator >
  typename Map::size_type unmark_oriented_orbit(const Map & amap,
                                                typename Map::Dart_const_handle adart,
                                                typename Map::size_type amark,
                                                typename Map::size_type amark2=Map::INVALID_MARK)
  {
    amap.negate_mark(amark);
    if (amark2!=Map::INVALID_MARK) { amap.negate_mark(amark2); }
    typename Map::size_type
      res=CGAL::template mark_oriented_orbit<Map, Iterator>(amap, adart, amark, amark2);
    amap.negate_mark(amark);
    if (amark2!=Map::INVALID_MARK) { amap.negate_mark(amark2); }
    return res;
  }

  /** Mark a given oriented cell with a given mark. For combinatorial map,
   *  equivalent to mark_cell. For generalized map, mark one out of two dart
   *  if the cell is orientable. If amark2!=INVALID_MARK, mark also all the darts
   *  of the cell with amark2.
   * @param amap a combinatorial map / generalized.
   * @param adart a dart of the cell.
   * @param amark the mark.
   * @param amark2 a second mark.
   * @return the number of darts of the marked cell.
   * @pre The whole cell must be unmarked.
   */
  template < class Map, unsigned int i, unsigned int d >
  typename Map::size_type mark_oriented_cell(const Map & amap,
                                             typename Map::Dart_const_handle adart,
                                             typename Map::size_type amark,
                                             typename Map::size_type amark2=Map::INVALID_MARK)
  { return CGAL::template mark_oriented_orbit
      <Map, typename Map::template Dart_of_cell_basic_range<i,d>::const_iterator>
      (amap, adart, amark, amark2); }

  template < class Map, unsigned int i>
  typename Map::size_type mark_oriented_cell(const Map & amap,
                                             typename Map::Dart_const_handle adart,
                                             typename Map::size_type amark,
                                             typename Map::size_type amark2=Map::INVALID_MARK)
  { return CGAL::template mark_oriented_cell<Map,i,Map::dimension>
      (amap, adart, amark, amark2);}

  /** Unmark a given oriented cell with a given mark.
   * @param amap a combinatorial map.
   * @param adart a dart of the cell.
   * @param amark the mark.
   * @return the number of darts of the unmarked cell.
   * @pre The whole cell must be marked.
   */
  template < class Map, unsigned int i, unsigned int d >
  typename Map::size_type unmark_oriented_cell(const Map & amap,
                                               typename Map::Dart_const_handle adart,
                                               typename Map::size_type amark,
                                               typename Map::size_type amark2=Map::INVALID_MARK)
  { return CGAL::template unmark_oriented_orbit
      <Map, typename Map::template Dart_of_cell_basic_range<i,d>::const_iterator>
      (amap, adart, amark, amark2);}

  template < class Map, unsigned int i >
  typename Map::size_type unmark_oriented_cell(const Map & amap,
                                               typename Map::Dart_const_handle adart,
                                               typename Map::size_type amark,
                                               typename Map::size_type amark2=Map::INVALID_MARK)
  { return CGAL::template unmark_oriented_cell<Map,i,Map::dimension>
      (amap, adart, amark, amark2); }

  /** Compute the degree of a given i-cell c.
   * The degree is the number of distinct i+1 cells incident to c.
   * @param amap a combinatorial map.
   * @param adart a dart of the cell.
   * @return the degree of the cell.
   */
  template < class Map, unsigned int i >
  typename Map::size_type degree(const Map & amap,
                                 typename Map::Dart_const_handle adart)
  {
    CGAL_assertion(adart != nullptr);

    typename Map::size_type nbIncident = 0;
    typename Map::size_type mark;
    typename Map::size_type treated;
    mark = amap.get_new_mark();
    treated = amap.get_new_mark();

    typename Map::template
      Dart_of_cell_basic_range<i>::const_iterator it(amap, adart, mark);
    for ( ;it.cont(); ++it )
    {
      if (!amap.is_marked(it, treated))
      {
        ++nbIncident;
        CGAL::mark_cell<Map,i+1>(amap, it, treated);
      }
      amap.mark(it,mark);
    }

    amap.negate_mark(mark);
    for (it.rewind(); it.cont(); ++it)
    {
      if (amap.is_marked(it, treated))
      { CGAL::unmark_cell<Map,i+1>(amap, it, treated); }
      amap.mark(it,mark);
    }

    amap.negate_mark(mark);
    CGAL_assertion( amap.is_whole_map_unmarked(mark) );
    CGAL_assertion( amap.is_whole_map_unmarked(treated) );

    amap.free_mark(mark);
    amap.free_mark(treated);

    CGAL_assertion(nbIncident != 0);
    return nbIncident;
  }

  /** Compute the co-degree of a given i-cell c.
   * The co-degree is the number of distinct i-1 cells incident to c.
   * @param amap a combinatorial map.
   * @param adart a dart of the cell.
   * @return the co-degree of the cell.
   */
  template < class Map, unsigned int i >
  typename Map::size_type codegree(const Map & amap,
                                   typename Map::Dart_const_handle adart)
  {
    CGAL_assertion(adart != nullptr);

    typename Map::size_type nbIncident = 0;
    typename Map::size_type mark;
    typename Map::size_type treated;
    mark = amap.get_new_mark();
    treated = amap.get_new_mark();

    typename Map::template
      Dart_of_cell_basic_range<i>::const_iterator it(amap, adart, mark);
    for ( ; it.cont(); ++it)
    {
      if (!amap.is_marked(it, treated))
      {
        ++nbIncident;
        CGAL::mark_cell<Map,i-1>(amap, it, treated);
      }
      amap.mark(it,mark);
    }

    amap.negate_mark(mark);
    for (it.rewind(); it.cont(); ++it)
    {
      if (amap.is_marked(it, treated))
      { CGAL::unmark_cell<Map,i-1>(amap, it, treated); }
      amap.mark(it,mark);
    }

    amap.negate_mark(mark);
    CGAL_assertion( amap.is_whole_map_unmarked(mark) );
    CGAL_assertion( amap.is_whole_map_unmarked(treated) );

    amap.free_mark(mark);
    amap.free_mark(treated);

    CGAL_assertion(nbIncident != 0);
    return nbIncident;
  }

} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_BASIC_OPERATIONS_H //
// EOF //
