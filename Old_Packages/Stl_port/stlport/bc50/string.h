/*
 * Copyright (c) 1999 
 * Boris Fomitchev
 *
 * This material is provided "as is", with absolutely no warranty expressed
 * or implied. Any use is at your own risk.
 *
 * Permission to use or copy this software for any purpose is hereby granted 
 * without fee, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is granted,
 * provided the above notices are retained, and a notice that the code was
 * modified is included with the above copyright notice.
 *
 */

/*
 *
 *  This wrapper is needed for Borland C++ 5.0 to get STLport 
 *  header properly included
 */


#  include <..\string.h>

# ifndef __IN_STLPORT_CSTRING
// okay, include STLPort header
#  include  <..\string.>
# endif /* __IN_STLPORT_CSTRING */


// Local Variables:
// mode:C++
// End:
