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


#ifndef __CHECK_SYSTEM_PARAMS__
#define __CHECK_SYSTEM_PARAMS__

namespace CGAL {

namespace internal {

/// A general system parameters definition.
template< int sz > struct _EnvDef{
    typedef void s_ptr_type;                        ///< The pointer size type.
    typedef void ptr_type;                          ///< The pointer type.
}; // struct _EnvDef

/// The x32 system parameters definition.
template<> struct _EnvDef<4> {
    typedef int s_ptr_type;                         ///< The pointer size type.
    typedef unsigned int ptr_type;                  ///< The pointer type.
}; // struct _EnvDef<4>

/// The x64 system parameters definition.
template<> struct _EnvDef<8> {
    typedef long long s_ptr_type;                   ///< The pointer size type.
    typedef unsigned long long ptr_type;            ///< The pointer type.
}; // struct _EnvDef<8>

/// The current system configuration.
struct _Env {
    typedef _EnvDef< sizeof( void* ) > _DEF;        ///< The system definitions.

    typedef _DEF::s_ptr_type    s_ptr_type;         ///< The pointer size type.
    typedef _DEF::ptr_type      ptr_type;           ///< The pointer type.

    static const int ptr_size = sizeof(ptr_type)*8; ///< The size of a pointer.
    static const bool is_x32 = ( ptr_size == 32 );  ///< Whether the system is 32-bit.
    static const bool is_x64 = ( ptr_size == 64 );  ///< Whether the system is 64-bit.
}; // struct _Env

typedef _Env _ENVIRONMENT;                          ///< The system environment.

} // namespace internal

} // namespace CGAL

#endif // __CHECK_SYSTEM_PARAMS__