// Copyright (c) 1997-2001
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
      // Sets the rounding mode to 0 and show that it's really the case.
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

