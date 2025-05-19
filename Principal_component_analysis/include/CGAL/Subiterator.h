// Copyright (c) 2020  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Simon Giraudot
//

#ifndef CGAL_PCA_SUBITERATOR_H
#define CGAL_PCA_SUBITERATOR_H

#include <CGAL/license/Principal_component_analysis.h>


#include <boost/iterator/iterator_facade.hpp>

namespace CGAL
{

template <typename InputIterator, typename ValueType, int Size>
class Subiterator
  : public boost::iterator_facade<Subiterator<InputIterator, ValueType, Size>,
                                  ValueType,
                                  std::input_iterator_tag>
{
public:
  using Self = Subiterator<InputIterator, ValueType, Size>;
  using Facade = boost::iterator_facade<Self, ValueType, std::input_iterator_tag>;

  using Input_type = typename std::iterator_traits<InputIterator>::value_type;
  using Output_type = ValueType;

  using Converter = std::function<Output_type(const Input_type&, int)>;

private:

  Converter m_converter;
  int m_next_index;
  InputIterator m_base;
  mutable Output_type m_current;

public:

  Subiterator() { }

  Subiterator(InputIterator begin, const Converter& converter)
    : m_converter(converter), m_next_index(0), m_base(begin)
  { }

  Subiterator(InputIterator end)
    : m_next_index(0), m_base(end)
  { }

private:

  friend class boost::iterator_core_access;

  void increment()
  {
    ++ m_next_index;
    if (m_next_index == Size)
    {
      ++ m_base;
      m_next_index = 0;
    }
  }

  bool equal(const Self& other) const
  {
    return this->m_base == other.m_base && this->m_next_index == other.m_next_index;
  }

  Output_type& dereference() const
  {
    m_current = m_converter (*m_base, m_next_index);
    return const_cast<Output_type&>(m_current);
  }
};

template <typename ValueType, int Size, typename InputIterator>
Subiterator<InputIterator, ValueType, Size>
make_subiterator (InputIterator begin,
                  const typename Subiterator<InputIterator, ValueType, Size>::Converter& converter)
{
  return Subiterator<InputIterator, ValueType, Size>(begin, converter);
}

template <typename ValueType, int Size, typename InputIterator>
Subiterator<InputIterator, ValueType, Size>
make_subiterator (InputIterator end)
{
  return Subiterator<InputIterator, ValueType, Size>(end);
}


} // namespace CGAL

#endif // CGAL_PCA_SUBITERATOR_H
