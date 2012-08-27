// Copyright (c) 1997-2001  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
//
// Author(s)     : Andreas Fabri

//| If a compiler has problems with restoring the rounding mode in a try/catch
//| CGAL_CFG_FPU_ROUNDING_MODE_UNWINDING_VC_BUG


#ifdef  _MSC_VER

#include <float.h>

int
main()
{
  try
    {
      _controlfp( _RC_UP, _MCW_RC );
      throw 1;
    }
  catch ( ... )
    {
      // Sets the rounding mode to 0 and show that it's realy the case.
      _controlfp( 0, _MCW_RC );
    }

  // After the catch, this value will not be 0 in x64
  unsigned int cw = _controlfp( 0, 0 ) & _MCW_RC;
  return (cw != 0);
}

#else

int main()
{
  return 0;
}
#endif

