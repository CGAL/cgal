// Copyright (c) 2009  GeometryFactory Sarl (France).
// All rights reserved.
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
// $URL$
// $Id$
//
// Author(s)     : Laurent Rineau

//| If a compiler does not support the alternative tokens for logicial
//| operators (section 2.5 Alternative tokens [lex.digraph] of the C++
//| norm, 2003), then CGAL_CFG_NO_LOGICAL_OPERATORS_ALTERNATIVES is set.

int main()
{
  if( true and (not false) )
    if(1 not_eq 2)
      if(false or true)
        return 0;

  return 1;
}
