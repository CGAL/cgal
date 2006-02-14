// Copyright (c) 1999  Martin-Luther-University Halle-Wittenberg (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Matthias Baesken, Algorithmic Solutions


#ifndef CGAL_WINDOW_SYS_CYGWIN32_H
#define CGAL_WINDOW_SYS_CYGWIN32_H

#define __win32__

#define __HAS_EXPLICIT_KEYWORD__
#define __HAS_TYPENAME_KEYWORD__
#define __HAS_MEMBER_TEMPLATES__


#define __explicit          explicit
#define __typename          typename
#define __temp_friend_decl  <>
#define __temp_func_inline


#define LITTLE_ENDIAN_MACHINE


//------------------------------------------------------------------------------
//  DLL definitions
//------------------------------------------------------------------------------

#define __exportC
#define __exportF
#define __exportD

#define _exportC
#define _exportF
#define _exportD



#endif
