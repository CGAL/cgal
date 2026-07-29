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
////////////////////////////////////////////////////////////////////////////////
#ifndef CMAP_COPY_H
#define CMAP_COPY_H

#include <unordered_map>
#include <CGAL/Combinatorial_map/internal/Combinatorial_map_copy_functors.h>

////////////////////////////////////////////////////////////////////////////////
namespace CGAL::CMap_copy
{
////////////////////////////////////////////////////////////////////////////////
/*!
 * \brief Copy marked darts from lcc1 to lcc2.
 * BECAREFUL: depending on which darts are marked, can produce an invalid map!
 * \return the number of new darts.
 */
template<typename CMap1, typename CMap2, typename Converters,
         typename DartInfoConverter, typename PointConverter>
std::size_t partial_copy(
    const CMap1& amap1,
    const typename CMap1::size_type amark,
    CMap2& amap2,
    std::unordered_map<typename CMap1::Dart_const_descriptor,
                       typename CMap2::Dart_descriptor>       *origin_to_copy,
    std::unordered_map<typename CMap2::Dart_descriptor,
                       typename CMap1::Dart_const_descriptor> *copy_to_origin,
    const Converters&         converters,
    const DartInfoConverter&  dartinfoconverter,
    const PointConverter&     pointconverter,
    bool                      copy_perforated_darts=false,
    typename CMap1::size_type mark_perforated=CMap1::INVALID_MARK)
{
  std::size_t res = 0;

  std::unordered_map<typename CMap1::Dart_const_descriptor,
                     typename CMap2::Dart_descriptor> local_dartmap;
  // Use local_dartmap if user does not provides its own unordered_map
  if (origin_to_copy == nullptr)
  { origin_to_copy = &local_dartmap; }

  typename CMap2::Dart_descriptor new_dart;
  for (auto it = amap1.darts().begin(); it != amap1.darts().end(); ++it)
  {
    if (amap1.is_marked(it, amark) &&
       (copy_perforated_darts || !amap1.is_perforated(it)))
    {
      new_dart = amap2.create_dart(); // , amap.get_marks(it));

      if (mark_perforated != CMap1::INVALID_MARK && amap1.is_perforated(it))
      { amap2.mark(new_dart, mark_perforated); }

      (*origin_to_copy)[it] = new_dart;
      if (copy_to_origin != nullptr)
      { (*copy_to_origin)[new_dart] = it; }

      CGAL::internal::Copy_dart_info_functor<typename CMap1::Refs,
                                             typename CMap2::Refs,
                                             DartInfoConverter>::run(
            static_cast<const typename CMap1::Refs&>(amap1),
            static_cast<typename CMap2::Refs&>(amap2),
            it, new_dart, dartinfoconverter);
    }
  }

  unsigned int min_dim = std::min(amap1.dimension, amap2.dimension);

  typename std::unordered_map<
      typename CMap1::Dart_const_descriptor,
      typename CMap2::Dart_descriptor>::iterator dartmap_iter;
  for (dartmap_iter = origin_to_copy->begin();
       dartmap_iter != origin_to_copy->end(); ++dartmap_iter)
  {
    for (unsigned int i = 0; i <= min_dim; ++i)
    {
      if (!amap1.is_free(dartmap_iter->first, i) &&
          amap2.is_free(dartmap_iter->second, i) &&
          amap1.is_marked(amap1.beta(dartmap_iter->first, i), amark))
      {
        amap2.basic_link_beta(
              dartmap_iter->second,
              (*origin_to_copy)[amap1.beta(dartmap_iter->first, i)], i);
      }
    }
  }

  /** Copy attributes */
  for (dartmap_iter = origin_to_copy->begin();
       dartmap_iter != origin_to_copy->end(); ++dartmap_iter)
  {
    CMap2::Helper::template Foreach_enabled_attributes
      <CGAL::internal::Copy_attributes_functor<typename CMap1::Refs,
        typename CMap2::Refs, Converters, PointConverter>>::run(
          static_cast<const typename CMap1::Refs&>(amap1),
          static_cast<typename CMap2::Refs&>(amap2),
          dartmap_iter->first, dartmap_iter->second,
          converters, pointconverter);
  }

  CGAL_assertion (amap2.is_valid());

  return res;
}
////////////////////////////////////////////////////////////////////////////////
template<typename CMap1, typename CMap2>
std::size_t partial_copy(
    const CMap1& amap1,
    const typename CMap1::size_type amark,
    CMap2& amap2,
    std::unordered_map
        <typename CMap1::Dart_const_descriptor,
         typename CMap2::Dart_descriptor>       *origin_to_copy=nullptr,
    std::unordered_map
        <typename CMap2::Dart_descriptor,
         typename CMap1::Dart_const_descriptor> *copy_to_origin=nullptr,
    bool copy_perforated_darts=false,
    typename CMap1::size_type mark_perforated=CMap1::INVALID_MARK)
{
  std::tuple<> converters;
  CGAL::Default_converter_dart_info<
      typename CMap1::Refs, typename CMap2::Refs> dartinfoconverter;
  CGAL::Default_converter_cmap_0attributes_with_point<
      typename CMap1::Refs, typename CMap2::Refs> pointconverter;
  return partial_copy(amap1, amark, amap2, origin_to_copy, copy_to_origin,
                      converters, dartinfoconverter, pointconverter,
                      copy_perforated_darts, mark_perforated);
}
////////////////////////////////////////////////////////////////////////////////
/** Copy each cell<i,dim> incident to the cell<j,dim> containing *dh from
 *  lcc1 to lcc2.
 *  @return the number of new darts.
 */
template<unsigned int i, unsigned int j, unsigned int dim,
         typename CMap1, typename CMap2, typename Converters,
         typename DartInfoConverter, typename PointConverter>
std::size_t copy_incident_cells(
    const CMap1& amap1,
    const typename CMap1::Dart_const_descriptor d,
    CMap2& amap2,
    std::unordered_map
        <typename CMap1::Dart_const_descriptor,
         typename CMap2::Dart_descriptor>       *origin_to_copy,
    std::unordered_map
        <typename CMap2::Dart_descriptor,
         typename CMap1::Dart_const_descriptor> *copy_to_origin,
    const Converters&         converters,
    const DartInfoConverter&  dartinfoconverter,
    const PointConverter&     pointconverter,
    bool                      copy_perforated_darts=false,
    typename CMap1::size_type mark_perforated=CMap1::INVALID_MARK)
{
  typename CMap1::size_type amark = amap1.get_new_mark();
  for (auto it = amap1.template one_dart_per_incident_cell<i,j, dim>(d).begin();
       it != amap2.template one_dart_per_incident_cell<i,j, dim>(d).end(); ++it)
  { amap1.template mark_cell<i, dim>(it, amark); }
  std::size_t res = partial_copy(amap1, amark, amap2,
                                 origin_to_copy, copy_to_origin,
                                 converters, dartinfoconverter, pointconverter,
                                 copy_perforated_darts, mark_perforated);
  amap1.free_mark(amark);
  return res;
}
////////////////////////////////////////////////////////////////////////////////
/** Copy each cell<i> incident to the cell<j> containing *dh from lcc1 to lcc2.
 *  @return the number of new darts.
 */
template<unsigned int i, unsigned int j,
         typename CMap1, typename CMap2, typename Converters,
         typename DartInfoConverter, typename PointConverter>
std::size_t copy_incident_cells(
    const CMap1& amap1,
    const typename CMap1::Dart_const_descriptor d,
    CMap2& amap2,
    std::unordered_map
    <typename CMap1::Dart_const_descriptor,
     typename CMap2::Dart_descriptor>       *origin_to_copy,
    std::unordered_map
    <typename CMap2::Dart_descriptor,
     typename CMap1::Dart_const_descriptor> *copy_to_origin,
    const Converters&         converters,
    const DartInfoConverter&  dartinfoconverter,
    const PointConverter&     pointconverter,
    bool                      copy_perforated_darts=false,
    typename CMap1::size_type mark_perforated=CMap1::INVALID_MARK)
{ return copy_incident_cells<i, j, CMap1::dim, CMap1, CMap2, Converters,
                             DartInfoConverter, PointConverter>
      (amap1, d, amap2, origin_to_copy, copy_to_origin, converters,
       dartinfoconverter, pointconverter, copy_perforated_darts,
       mark_perforated);
}
////////////////////////////////////////////////////////////////////////////////
/** Copy each cell<i,dim> incident to the cell<j,dim> containing *dh from
 *  lcc1 to lcc2.
 *  Version with default dart and point converters.
 *  @return the number of new darts.
 */
template<unsigned int i, unsigned int j, unsigned int dim,
         typename CMap1, typename CMap2>
std::size_t copy_incident_cells(
    const CMap1& amap1,
    const typename CMap1::Dart_const_descriptor d,
    CMap2& amap2,
    std::unordered_map
    <typename CMap1::Dart_const_descriptor,
     typename CMap2::Dart_descriptor>       *origin_to_copy=nullptr,
    std::unordered_map
    <typename CMap2::Dart_descriptor,
     typename CMap1::Dart_const_descriptor> *copy_to_origin=nullptr,
    bool copy_perforated_darts=false,
    typename CMap1::size_type mark_perforated=CMap1::INVALID_MARK)
{
  std::tuple<> converters;
  CGAL::Default_converter_dart_info<
      typename CMap1::Refs, typename CMap2::Refs> dartinfoconverter;
  CGAL::Default_converter_cmap_0attributes_with_point<
      typename CMap1::Refs, typename CMap2::Refs> pointconverter;
  return copy_incident_cells<i,j, dim>(
      amap1, d, amap2, origin_to_copy, copy_to_origin,
      converters, dartinfoconverter, pointconverter,
      copy_perforated_darts, mark_perforated);
}
////////////////////////////////////////////////////////////////////////////////
/** Copy each cell<i> incident to the cell<j> containing *dh from lcc1 to lcc2.
 *  Version with default dart and point converters.
 *  @return the number of new darts.
 */
template<unsigned int i, unsigned int j, typename CMap1, typename CMap2>
std::size_t copy_incident_cells(
    const CMap1& amap1,
    const typename CMap1::Dart_const_descriptor d,
    CMap2& amap2,
    std::unordered_map
    <typename CMap1::Dart_const_descriptor,
     typename CMap2::Dart_descriptor>       *origin_to_copy=nullptr,
    std::unordered_map
    <typename CMap2::Dart_descriptor,
     typename CMap1::Dart_const_descriptor> *copy_to_origin=nullptr,
    bool copy_perforated_darts=false,
    typename CMap1::size_type mark_perforated=CMap1::INVALID_MARK)
{
  return copy_incident_cells<i,j, CMap1::dimension>(
      amap1, d, amap2, origin_to_copy, copy_to_origin,
      copy_perforated_darts, mark_perforated);
}
////////////////////////////////////////////////////////////////////////////////
/** Copy all cells<i,dim>(dh) from lcc1 to lcc2, with converters for darts
 *  and attributes.
 *  @return the number of new darts.
 */
template<unsigned int i, unsigned int dim, typename CMap1, typename CMap2,
         typename Converters, typename DartInfoConverter, typename PointConverter>
std::size_t copy_cells(CMap1& amap1,
                       const std::vector<typename CMap1::Dart_descriptor>& cells,
                       CMap2& amap2,
                       std::unordered_map
                       <typename CMap1::Dart_descriptor, typename CMap2::Dart_descriptor>*
                           origin_to_copy,
                       std::unordered_map
                       <typename CMap1::Dart_descriptor, typename CMap2::Dart_descriptor>*
                           copy_to_origin,
                       const Converters& converters,
                       const DartInfoConverter& dartinfoconverter,
                       const PointConverter& pointconverter,
                       bool copy_perforated_darts=false,
                       typename CMap1::size_type mark_perforated=CMap1::INVALID_MARK)
{
  std::size_t res=0;

  std::unordered_map<typename CMap1::Dart_descriptor,
                     typename CMap2::Dart_descriptor> local_dartmap;
  if (origin_to_copy==nullptr) // Use local_dartmap if user does not provides its own unordered_map
  { origin_to_copy=&local_dartmap; }

  typename CMap2::Dart_descriptor new_dart;
  auto mark_copy=amap1.get_new_mark(); // mark already copied darts

  /// 1) Copy all the darts
  for(typename CMap1::Dart_descriptor dh: cells)
  {
    if(!amap1.is_marked(dh, mark_copy))
    {
      for(auto it=amap1.template darts_of_cell_basic<i, dim>(dh, mark_copy).begin(),
           itend=amap1.template darts_of_cell_basic<i, dim>(dh, mark_copy).end();
           it!=itend; ++it)
      {
        amap1.mark(it, mark_copy); // in case the basic iterator does not mark
        if (copy_perforated_darts || !amap1.is_perforated(it))
        {
          new_dart=amap2.create_dart(); // , amap.get_marks(it));

          if (mark_perforated!=CMap1::INVALID_MARK && amap1.is_perforated(it))
          { amap2.mark(new_dart, mark_perforated); }

          (*origin_to_copy)[it]=new_dart;
          if(copy_to_origin!=nullptr) { (*copy_to_origin)[new_dart]=it; }

          CGAL::internal::Copy_dart_info_functor
              <typename CMap1::Refs, typename CMap2::Refs, DartInfoConverter>::run
              (static_cast<const typename CMap1::Refs&>(amap1),
               static_cast<typename CMap2::Refs&>(amap2),
           it, new_dart, dartinfoconverter);
        }
      }
    }
  }

  /// 2) Link the different new darts between them
  unsigned int min_dim=std::min({amap1.dimension, amap2.dimension, dim});
  typename std::unordered_map<typename CMap1::Dart_descriptor,
                              typename CMap2::Dart_descriptor>::iterator
      dartmap_iter, dartmap_iter_end=origin_to_copy->end();
  for (dartmap_iter=origin_to_copy->begin(); dartmap_iter!=dartmap_iter_end;
       ++dartmap_iter)
  {
    for (unsigned int d=0; d<=min_dim; ++d)
    {
      if (!amap1.is_free(dartmap_iter->first,d) &&
          amap2.is_free(dartmap_iter->second,d) &&
          amap1.is_marked(amap1.beta(dartmap_iter->first,d), mark_copy))
      {
        amap2.basic_link_beta(dartmap_iter->second,
                              (*origin_to_copy)[amap1.beta(dartmap_iter->first,d)],
                              d);
      }
    }
  }

  /// 3) Copy attributes
  for (dartmap_iter=origin_to_copy->begin(); dartmap_iter!=dartmap_iter_end;
       ++dartmap_iter)
  {
    amap1.unmark(dartmap_iter->first, mark_copy);
    CMap2::Helper::template Foreach_enabled_attributes
        <CGAL::internal::Copy_attributes_functor<typename CMap1::Refs,
                                                 typename CMap2::Refs, Converters, PointConverter>>::
        run(static_cast<const typename CMap1::Refs&>(amap1),
            static_cast<typename CMap2::Refs&>(amap2),
            dartmap_iter->first, dartmap_iter->second, converters, pointconverter);
  }
  CGAL_assertion(amap1.is_whole_map_unmarked(mark_copy));
  amap1.free_mark(mark_copy);
  // CGAL_assertion (amap2.is_valid());

  return res;
}
////////////////////////////////////////////////////////////////////////////////
/** Copy all cells<i>(dh) from lcc1 to lcc2, with converters for darts
 *  and attributes.
 *  @return the number of new darts.
 */
template<unsigned int i, typename CMap1, typename CMap2,
         typename Converters, typename DartInfoConverter, typename PointConverter>
std::size_t copy_cells(CMap1& amap1,
    const std::vector<typename CMap1::Dart_descriptor>& cells,
    CMap2& amap2,
    std::unordered_map
    <typename CMap1::Dart_descriptor, typename CMap2::Dart_descriptor>*
        origin_to_copy,
    std::unordered_map
    <typename CMap1::Dart_descriptor, typename CMap2::Dart_descriptor>*
        copy_to_origin,
    const Converters& converters,
    const DartInfoConverter& dartinfoconverter,
    const PointConverter& pointconverter,
    bool copy_perforated_darts=false,
    typename CMap1::size_type mark_perforated=CMap1::INVALID_MARK)
{ return copy_cells<i, CMap1::dimension, CMap1, CMap2, Converters,
                    DartInfoConverter, PointConverter>
      (amap1, cells, amap2, origin_to_copy, copy_to_origin, converters,
       dartinfoconverter, pointconverter, copy_perforated_darts,
       mark_perforated);
}
////////////////////////////////////////////////////////////////////////////////
/** Copy all cells<i,dim>(dh) from lcc1 to lcc2, from lcc1 to lcc2 without
 *  converters.
 *  @return the number of new darts.
 */
template<unsigned int i, unsigned int dim, typename CMap1, typename CMap2>
std::size_t copy_cells(CMap1& amap1,
                       const std::vector<typename CMap1::Dart_descriptor>& cells,
                       CMap2& amap2,
                       std::unordered_map
                       <typename CMap1::Dart_descriptor, typename CMap2::Dart_descriptor>*
                           origin_to_copy=nullptr,
                       std::unordered_map
                       <typename CMap1::Dart_descriptor, typename CMap2::Dart_descriptor>*
                           copy_to_origin=nullptr,
                       bool copy_perforated_darts=false,
                       typename CMap1::size_type mark_perforated=CMap1::INVALID_MARK)
{
  std::tuple<> converters;
  CGAL::Default_converter_dart_info<typename CMap1::Refs,
                                    typename CMap2::Refs> dartinfoconverter;
  CGAL::Default_converter_cmap_0attributes_with_point<typename CMap1::Refs,
                                                      typename CMap2::Refs> pointconverter;
  return copy_cells<i, dim>(amap1, cells, amap2, origin_to_copy, copy_to_origin,
                         converters, dartinfoconverter, pointconverter,
                         copy_perforated_darts, mark_perforated);
}
////////////////////////////////////////////////////////////////////////////////
/** Copy all cells<i>(dh) from lcc1 to lcc2, from lcc1 to lcc2 without
 *  converters.
 *  @return the number of new darts.
 */
template<unsigned int i, typename CMap1, typename CMap2>
std::size_t copy_cells(CMap1& amap1,
    const std::vector<typename CMap1::Dart_descriptor>& cells,
    CMap2& amap2,
    std::unordered_map
    <typename CMap1::Dart_descriptor, typename CMap2::Dart_descriptor>*
        origin_to_copy=nullptr,
    std::unordered_map
    <typename CMap1::Dart_descriptor, typename CMap2::Dart_descriptor>*
        copy_to_origin=nullptr,
    bool copy_perforated_darts=false,
    typename CMap1::size_type mark_perforated=CMap1::INVALID_MARK)
{ return copy_cells<i, CMap1::dimension, CMap1, CMap2>
      (amap1, cells, amap2, origin_to_copy, copy_to_origin,
       copy_perforated_darts, mark_perforated); }
////////////////////////////////////////////////////////////////////////////////
/** Copy cell<i,dim>(dh) from lcc1 to lcc2, with converters for darts and
 *  attributes.
 *  @return the number of new darts.
 */
template<unsigned int i, unsigned int dim, typename CMap1, typename CMap2,
         typename Converters, typename DartInfoConverter, typename PointConverter>
std::size_t copy_cell(CMap1& amap1, typename CMap1::Dart_descriptor dh,
                      CMap2& amap2,
                      std::unordered_map
                      <typename CMap1::Dart_descriptor, typename CMap2::Dart_descriptor>*
                          origin_to_copy,
                      std::unordered_map
                      <typename CMap1::Dart_descriptor, typename CMap2::Dart_descriptor>*
                          copy_to_origin,
                      const Converters& converters,
                      const DartInfoConverter& dartinfoconverter,
                      const PointConverter& pointconverter,
                      bool copy_perforated_darts=false,
                      typename CMap1::size_type mark_perforated=CMap1::INVALID_MARK)
{
  std::vector<typename CMap1::Dart_descriptor> cells;
  cells.push_back(dh);
  return copy_cells<i, dim>(amap1, cells, amap2, origin_to_copy, copy_to_origin,
                         converters, dartinfoconverter, pointconverter,
                         copy_perforated_darts, mark_perforated);
}
////////////////////////////////////////////////////////////////////////////////
/** Copy cell<i>(dh) from lcc1 to lcc2, with converters for darts and
 *  attributes.
 *  @return the number of new darts.
 */
template<unsigned int i, typename CMap1, typename CMap2,
         typename Converters, typename DartInfoConverter, typename PointConverter>
std::size_t copy_cell(CMap1& amap1, typename CMap1::Dart_descriptor dh,
    CMap2& amap2,
    std::unordered_map
    <typename CMap1::Dart_descriptor, typename CMap2::Dart_descriptor>*
        origin_to_copy,
    std::unordered_map
    <typename CMap1::Dart_descriptor, typename CMap2::Dart_descriptor>*
        copy_to_origin,
    const Converters& converters,
    const DartInfoConverter& dartinfoconverter,
    const PointConverter& pointconverter,
    bool copy_perforated_darts=false,
    typename CMap1::size_type mark_perforated=CMap1::INVALID_MARK)
{ return copy_cell<i, CMap1::dimension, CMap1, CMap2, Converters,
                   DartInfoConverter, PointConverter>
      (amap1, dh, amap2, origin_to_copy, copy_to_origin, converters,
       dartinfoconverter, pointconverter, copy_perforated_darts,
       mark_perforated); }
////////////////////////////////////////////////////////////////////////////////
/** Copy cell<i, dim>(dh) from lcc1 to lcc2 without converters.
 *  @return the number of new darts.
 */
template<unsigned int i, unsigned int dim, typename CMap1, typename CMap2>
std::size_t copy_cell(CMap1& amap1, typename CMap1::Dart_descriptor dh,
                      CMap2& amap2,
                      std::unordered_map
                      <typename CMap1::Dart_descriptor, typename CMap2::Dart_descriptor>*
                      origin_to_copy=nullptr,
                      std::unordered_map
                      <typename CMap1::Dart_descriptor, typename CMap2::Dart_descriptor>*
                      copy_to_origin=nullptr,
                      bool copy_perforated_darts=false,
                      typename CMap1::size_type mark_perforated=CMap1::INVALID_MARK)
{
  std::tuple<> converters;
  CGAL::Default_converter_dart_info<typename CMap1::Refs,
      typename CMap2::Refs> dartinfoconverter;
  CGAL::Default_converter_cmap_0attributes_with_point<typename CMap1::Refs,
       typename CMap2::Refs> pointconverter;
  return copy_cell<i>(amap1, dh, amap2, origin_to_copy, copy_to_origin,
                        converters, dartinfoconverter, pointconverter,
                        copy_perforated_darts, mark_perforated);
}
////////////////////////////////////////////////////////////////////////////////
/** Copy cell<i>(dh) from lcc1 to lcc2 without converters.
 *  @return the number of new darts.
 */
template<unsigned int i, typename CMap1, typename CMap2>
std::size_t copy_cell(CMap1& amap1, typename CMap1::Dart_descriptor dh,
    CMap2& amap2,
    std::unordered_map
    <typename CMap1::Dart_descriptor, typename CMap2::Dart_descriptor>*
        origin_to_copy=nullptr,
    std::unordered_map
    <typename CMap1::Dart_descriptor, typename CMap2::Dart_descriptor>*
        copy_to_origin=nullptr,
    bool copy_perforated_darts=false,
    typename CMap1::size_type mark_perforated=CMap1::INVALID_MARK)
{ return copy_cell<i, CMap1::dimension, CMap1, CMap2>
      (amap1, dh, amap2, origin_to_copy, copy_to_origin,
       copy_perforated_darts, mark_perforated); }

////////////////////////////////////////////////////////////////////////////////
} // namespace CGAL::CMap_copy
////////////////////////////////////////////////////////////////////////////////
#endif // CMAP_COPY_H
