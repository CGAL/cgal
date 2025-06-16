// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of 3d-query-replace.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef CMAP_COPY_H
#define CMAP_COPY_H

#include <unordered_map>
#include <CGAL/Combinatorial_map/internal/Combinatorial_map_copy_functors.h>

/** Copy volume(dh) from lcc1 to lcc2.
 *  @return the number of new darts.
 */
template<unsigned int dim, typename CMap1, typename CMap2,
         typename Converters, typename DartInfoConverter, typename PointConverter>
std::size_t copy_cell(CMap1& amap1, typename CMap1::Dart_handle dh,
                      CMap2& amap2,
                      std::unordered_map
                      <typename CMap1::Dart_handle, typename CMap2::Dart_handle>*
                      origin_to_copy,
                      std::unordered_map
                      <typename CMap1::Dart_handle, typename CMap2::Dart_handle>*
                      copy_to_origin,
                      const Converters& converters,
                      const DartInfoConverter& dartinfoconverter,
                      const PointConverter& pointconverter,
                      bool copy_perforated_darts=false,
                      typename CMap1::size_type mark_perforated=CMap1::INVALID_MARK)
{
  std::size_t res=0;

  std::unordered_map<typename CMap1::Dart_handle, typename CMap2::Dart_handle> local_dartmap;
  if (origin_to_copy==nullptr) // Use local_dartmap if user does not provides its own unordered_map
  { origin_to_copy=&local_dartmap; }

  typename CMap2::Dart_handle new_dart;
  for(auto it=amap1.template darts_of_cell<dim>(dh).begin(),
      itend=amap1.template darts_of_cell<dim>(dh).end(); it!=itend; ++it)
  {
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

  // unsigned int min_dim=std::min(amap1.dimension, amap2.dimension);
  unsigned int min_dim=3;

  typename std::unordered_map<typename CMap1::Dart_handle,
      typename CMap2::Dart_handle>::iterator
    dartmap_iter, dartmap_iter_end=origin_to_copy->end();
  for (dartmap_iter=origin_to_copy->begin(); dartmap_iter!=dartmap_iter_end;
       ++dartmap_iter)
  {
    for (unsigned int i=0; i<=min_dim; i++)
    {
      if (i!=dim &&
          !amap1.is_free(dartmap_iter->first,i) &&
          amap2.is_free(dartmap_iter->second,i))
      {
        amap2.basic_link_beta(dartmap_iter->second,
                              (*origin_to_copy)[amap1.beta(dartmap_iter->first,i)], i);
      }
    }
  }

  /** Copy attributes */
  for (dartmap_iter=origin_to_copy->begin(); dartmap_iter!=dartmap_iter_end;
       ++dartmap_iter)
  {
    CMap2::Helper::template Foreach_enabled_attributes
      <CGAL::internal::Copy_attributes_functor<typename CMap1::Refs,
        typename CMap2::Refs, Converters, PointConverter>>::
      run(static_cast<const typename CMap1::Refs&>(amap1),
          static_cast<typename CMap2::Refs&>(amap2),
          dartmap_iter->first, dartmap_iter->second, converters, pointconverter);
  }

  CGAL_assertion (amap2.is_valid());

  return res;
}

template<unsigned int dim, typename CMap1, typename CMap2>
std::size_t copy_cell(CMap1& amap1, typename CMap1::Dart_handle dh,
                      CMap2& amap2,
                      std::unordered_map
                      <typename CMap1::Dart_handle, typename CMap2::Dart_handle>*
                      origin_to_copy=nullptr,
                      std::unordered_map
                      <typename CMap1::Dart_handle, typename CMap2::Dart_handle>*
                      copy_to_origin=nullptr,
                      bool copy_perforated_darts=false,
                      typename CMap1::size_type mark_perforated=CMap1::INVALID_MARK)
{
  std::tuple<> converters;
  CGAL::Default_converter_dart_info<typename CMap1::Refs,
      typename CMap2::Refs> dartinfoconverter;
  CGAL::Default_converter_cmap_0attributes_with_point<typename CMap1::Refs,
       typename CMap2::Refs> pointconverter;
  return copy_cell<dim>(amap1, dh, amap2, origin_to_copy, copy_to_origin,
                        converters, dartinfoconverter, pointconverter,
                        copy_perforated_darts, mark_perforated);
}

#endif // CMAP_COPY_H
