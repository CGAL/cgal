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
#ifndef ONE_INFO_H
#define ONE_INFO_H

#include <cstddef>
#include <ostream>

/// To store min, max and sum of a given info
template<typename T>
class One_info
{
public:
  One_info(): m_min(static_cast<T>(0)),
              m_max(static_cast<T>(0)),
              m_sum(static_cast<T>(0)),
              m_nb(0)
  {}

  void update(T val)
  {
    if (m_nb==0)
    {
      m_min=val;
      m_max=val;
      m_sum=val;
    }
    else
    {
      if (val<m_min) m_min=val;
      if (val>m_max) m_max=val;
      m_sum+=val;
    }
    ++m_nb;
  }

  T min() const
  { return m_min; }

  T max() const
  { return m_max; }

  double mean() const
  {
    if (m_nb==0) {return 0.; }
    return static_cast<double>(m_sum)/static_cast<double>(m_nb);
  }

  T sum() const
  { return m_sum; }

  std::size_t number_of_elements() const
  { return  m_nb; }

  void reinit()
  {
    m_min=static_cast<T>(0);
    m_max=static_cast<T>(0);
    m_sum=static_cast<T>(0);
    m_nb=0;
  }

  friend std::ostream& operator<<(std::ostream& os, const One_info& info)
  {
    os<<info.mean()<<" ("<<info.min()<<", "<<info.max()<<")";
    return os;
  }

protected:
  T m_min, m_max, m_sum;
  std::size_t m_nb;
};

#endif // ONE_INFO_H
