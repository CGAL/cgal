//Check if the system is x32 or x64.
//Copyright (C) 2013  INRIA - Sophia Antipolis
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Author(s):      Thijs van Lankveld

    
template<int sz> struct _Env{ typedef void s_ptr_type; typedef void ptr_type; static const int ptr_size = sizeof(ptr_type)*8; static const bool x64 = false; };
template<> struct _Env<4> { typedef int s_ptr_type; typedef unsigned int ptr_type; static const int ptr_size = sizeof(ptr_type)*8; };
template<> struct _Env<8> { typedef long long s_ptr_type; typedef unsigned long long ptr_type; static const int ptr_size = sizeof(ptr_type)*8; static const bool x64 = true;};

typedef _Env<sizeof(void*)> _ENV;