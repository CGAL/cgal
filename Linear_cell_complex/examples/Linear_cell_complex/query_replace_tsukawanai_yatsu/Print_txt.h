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
#ifndef PRINT_TXT_H
#define PRINT_TXT_H

#include <sstream>
#include <utility>
#include <chrono>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
template<typename T>
void put_one_value(std::ostringstream& txt, const T& v)
{ txt<<v; }
//----------------------------------------------------------------------------
template<> inline
void put_one_value(std::ostringstream& txt,
                   const std::chrono::duration<double>& v)
{
  //std::chrono::seconds s=std::chrono::duration_cast<std::chrono::seconds>(v);
  txt<<v.count();
}
///////////////////////////////////////////////////////////////////////////////
template <typename Arg, typename... Args>
void print_txt(Arg&& arg, Args&&... args)
{
  std::ostringstream txt;
  put_one_value(txt, std::forward<Arg>(arg));
  using expander = int[];
  (void)expander{0, (void(put_one_value(txt, std::forward<Args>(args))),0)...};
  std::cout<<txt.str()<<std::flush;
}
///////////////////////////////////////////////////////////////////////////////
template <typename Arg, typename... Args>
void print_txt_with_endl(Arg&& arg, Args&&... args)
{
  std::ostringstream txt;
  put_one_value(txt, std::forward<Arg>(arg));
  using expander = int[];
  (void)expander{0, (void(put_one_value(txt, std::forward<Args>(args))),0)...};
  txt<<std::endl;
  std::cout<<txt.str()<<std::flush;
}
///////////////////////////////////////////////////////////////////////////////

#endif // PRINT_TXT_H
