// Copyright (c) 1997-2002  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Radu Ursu


//+--------------------------------------------------------------------------
// The compiler has to provide a Standard Template Library
//+--------------------------------------------------------------------------
// STL test ok

//+--------------------------------------------------------------------------
//| The flag CGAL_CFG_CCTYPE_MACRO_BUG is set, if a compiler defines the
//| standard C library functions in cctype (isdigit etc.) as macros.
//| According to the standard they have to be functions.
//+--------------------------------------------------------------------------
//#define CGAL_CFG_CCTYPE_MACRO_BUG 1

//+--------------------------------------------------------------------------
//| The flag CGAL_CFG_EARLY_INSTANTIATION_BUG is set, if a compiler does not 
//| not how to compile the following code. See the solution bellow.
//| Created to workaround a cl1300 bug
//+--------------------------------------------------------------------------
//#define CGAL_CFG_EARLY_INSTANTIATION_BUG 1

//+--------------------------------------------------------------------------
//| This flag is set, if the compiler does not promote enumeration types
//| (which depend on a template parameter) correctly when they are used 
//| as int template arguments. (e.g. Borland 5.5)
//+--------------------------------------------------------------------------
//#define CGAL_CFG_ENUM_BUG 1

//+--------------------------------------------------------------------------
//| If a compiler (or assembler or linker) has problems with long names
//| CGAL_CFG_LONGNAME_BUG is set.
//+--------------------------------------------------------------------------
#define CGAL_CFG_LONGNAME_BUG 1

//+--------------------------------------------------------------------------
//| This flag is set, if the compiler does not match the most
//| specialized instance of a function template correctly,
//| but complains about multiple matches.
//| (e.g. VC++ 7)
//+--------------------------------------------------------------------------
//#define CGAL_CFG_MATCHING_BUG_2 1

//+--------------------------------------------------------------------------
//| This flag is set, if the compiler does not match function arguments 
//| of pointer type correctly, when the return type depends on 
//| the parameter's type. (e.g. sun C++ 5.3)
//+--------------------------------------------------------------------------
//#define CGAL_CFG_MATCHING_BUG_3 1

//+--------------------------------------------------------------------------
//| This flag is set, if the compiler does not match member 
//| definition to an existing declaration (eg. cl1310 Beta)
//+--------------------------------------------------------------------------
#define CGAL_CFG_NET2003_MATCHING_BUG 1

//+--------------------------------------------------------------------------
//| When template implementation files are not included in the source files,
//| a compiler may attempt to find the unincluded template bodies
//| automatically. For example, suppose that the following conditions are
//| all true.
//|
//| - template entity ABC::f is declared in file xyz.h
//| - an instantiation of ABC::f is required in a compilation
//| - no definition of ABC::f appears in the source code processed by the
//|   compilation
//| 
//| In this case, the compiler may look to see if the source file xyz.n exists,
//| where n is .c, .C, .cpp, .CPP, .cxx, .CXX, or .cc. If this feature is
//| missing, the flag CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION is set.
//+--------------------------------------------------------------------------
#define CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION 1

//+--------------------------------------------------------------------------
//| The byte order of a machine architecture distinguishes into
//| big-endian and little-endian machines.
//| The following definition is set if it is a little-endian machine.
//+--------------------------------------------------------------------------
#define CGAL_CFG_NO_BIG_ENDIAN 1

//+--------------------------------------------------------------------------
//| This flag is set if the compiler doesn't support the operator Koenig
//| lookup. That is, it does not search in the namespace of the arguments for
//| the function.
//+--------------------------------------------------------------------------
//#define CGAL_CFG_NO_KOENIG_LOOKUP 1

//+--------------------------------------------------------------------------
//| If a compiler doesn't know <limits> (g++-2.95)
//| or has a bug in the implementation (Sun CC 5.4, MipsPro CC)
//| CGAL_CFG_NO_LIMITS is set. 
//+--------------------------------------------------------------------------
//#define CGAL_CFG_NO_LIMITS 1

//+--------------------------------------------------------------------------
//| If a compiler doesn't know the locale classic
//| CGAL_CFG_NO_LOCALE is set. 
//+--------------------------------------------------------------------------
//#define CGAL_CFG_NO_LOCALE 1

//+--------------------------------------------------------------------------
//| The long long built-in integral type is not part of the ISO C++ standard,
//| but many compilers support it nevertheless since it's part of the ISO
//| C standard.
//| The following definition is set if it is supported.
//+--------------------------------------------------------------------------
//#define CGAL_CFG_NO_LONG_LONG 1

//+--------------------------------------------------------------------------
//| If a compiler doesn't support partial specialisation of class templates,
//| the flag CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION is set.
//+--------------------------------------------------------------------------
//#define CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION 1

//+--------------------------------------------------------------------------
//| The flag CGAL_CFG_NO_STDC_NAMESPACE is set, if a compiler does not
//| put the parts of the standard library inherited from the standard
//| C library in namespace std. (only tests for the symbols used in CGAL)
//+--------------------------------------------------------------------------
//#define CGAL_CFG_NO_STDC_NAMESPACE 1

//+--------------------------------------------------------------------------
//| G++ 2.95.2 has problems with member functions implemented outside of
//| the class body if this member function has a parameter type that is
//| dependant on a template in the template parameter list of the class. A
//| workaround would be to implement the member function inline in the class.
//| The following definition is set if this error error occurs.
//+--------------------------------------------------------------------------
//#define CGAL_CFG_NO_TMPL_IN_TMPL_DEPENDING_FUNCTION_PARAM 1

//+--------------------------------------------------------------------------
//| Nested templates in template parameter, such as 'template <
//| template <class T> class A>' are not supported by any compiler. 
//| The following definition is set if they are not supported.
//+--------------------------------------------------------------------------
//#define CGAL_CFG_NO_TMPL_IN_TMPL_PARAM 1

//+--------------------------------------------------------------------------
//| The flag CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG is set, if 
//| a compiler does not support the definition of the members templates 
//| out of line. The solution is to put the definition inside the class.
//| This is a feature of cl1200 and cl1300
//+--------------------------------------------------------------------------
//#define CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG 1

//+--------------------------------------------------------------------------
//| If a compiler does not accept the overloading of a template function, when
//| the template returns a reference, while the overloading doesn't.
//| In that case, CGAL_CFG_RETURN_TYPE_BUG is set.
//| This bug shows up on VC++ 6, and VC++ 7 beta 2.
//+--------------------------------------------------------------------------
//#define CGAL_CFG_RETURN_TYPE_BUG 1

//+--------------------------------------------------------------------------
//| The flag CGAL_CFG_USING_NAMESPACE_BUG is set, if a compiler does not 
//| not how to compile the following code.
//| Created to workaround a cl1300 bug
//+--------------------------------------------------------------------------
//#define CGAL_CFG_USING_NAMESPACE_BUG 1

//+--------------------------------------------------------------------------
//| This is a test-case for a bug in VC++ 7.0 beta2 that occurs in the kernel.
//| When the bug is present, CGAL_CFG_VC7_PRIVATE_TYPE_BUG is set.
//+--------------------------------------------------------------------------
//#define CGAL_CFG_VC7_PRIVATE_TYPE_BUG 1

