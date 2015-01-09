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

//  provides general system environment types.
//  \ingroup PkgScaleSpaceReconstruction3Auxiliary
template< int sz > struct _EnvTypes{
#ifdef DOXYGEN_RUNNING
/// \name Types
/// \{
    typedef unspecified_type s_ptr_type;                        ///< defines the pointer size type.
    typedef unspecified_type ptr_type;                          ///< defines the pointer type.

/// \}
#else // DOXYGEN_RUNNING
    typedef void s_ptr_type;
    typedef void ptr_type;
#endif // DOXYGEN_RUNNING
}; // struct _EnvTypes

// The x32 system environment types.
template<> struct _EnvTypes<4> {
    typedef int s_ptr_type;
    typedef unsigned int ptr_type;
}; // struct _EnvTypes<4>

// The x64 system environment types.
template<> struct _EnvTypes<8> {
    typedef long long s_ptr_type;
    typedef unsigned long long ptr_type;
}; // struct _EnvTypes<8>

//  specifies the current system environment.
//  \ingroup PkgScaleSpaceReconstruction3Auxiliary
struct _Env
: public _EnvTypes< sizeof( void* ) > {
/// \name Types
/// \{
    typedef _EnvTypes< sizeof( void* ) > Types;     ///< defines the system environment types.

    typedef Types::s_ptr_type    s_ptr_type;        ///< defines the pointer size type.
    typedef Types::ptr_type      ptr_type;          ///< defines the pointer type.

/// \}

/// \name Environment Parameters
/// \{
    static const int ptr_size = sizeof(ptr_type)*8; ///< stores the size of a pointer.
    static const bool is_x32 = ( ptr_size == 32 );  ///< stores whether the system is 32-bit.
    static const bool is_x64 = ( ptr_size == 64 );  ///< stores whether the system is 64-bit.

/// \}
}; // struct _Env

typedef _Env _ENVIRONMENT;                          //< specifies the system environment.

} // namespace internal

} // namespace CGAL

#endif // __CHECK_SYSTEM_PARAMS__
