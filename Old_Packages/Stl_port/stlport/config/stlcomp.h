/*
 * Copyright (c) 1997
 * Moscow Center for SPARC Technology
 *
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
 * Purpose of this file :
 *
 * To hold COMPILER-SPECIFIC portion of STLport settings.
 * In general, user should not edit this file unless 
 * using the compiler not recognized below.
 *
 * If your compiler is not being recognized yet, 
 * please look for definitions of macros in stl_mycomp.h,
 * copy stl_mycomp.h to stl_YOUR_COMPILER_NAME, 
 * adjust flags for your compiler, and add  <include config/stl_YOUR_COMPILER_NAME>
 * to the secton controlled by unique macro defined internaly by your compiler.
 *
 * NOTE : you may adjust these settings automatically for your 
 * compiler with the aid of "configure" script that comes with STLport.
 * 
 * To change user-definable settings, please edit <../stl_user_config.h> 
 *
 */

#ifndef __STLCOMP_H
# define __STLCOMP_H

//==========================================================
// per-version compiler features recognition
//==========================================================

// reporting of incompatibility
#  define __GIVE_UP_WITH_STL(message) void give_up() \
   { upgrade_the_compiler_to_use_STL;}

// distinguish real MSC from Metrowerks and Intel
# if defined(_MSC_VER) && !defined(__MWERKS__) && !defined (__ICL) && !defined (__COMO__)
#  define __STL_MSVC _MSC_VER
# endif


//==========================================================
// below, configuration for each compiler is being kept
// in separate configuration file.
//==========================================================

# if defined(__sgi) && !defined(__GNUC__)

// SGI CC compilers
#  include <config/stl_sgi.h>

# elif (defined(__OS400__))

// AS/400 C++
#  include <config/stl_as400.h>

# elif ( defined (__xlC__) && __xlC__ < 0x400 ) || \
    (defined(__MVS__) && defined ( __IBMCPP__ ) && (__IBMCPP__ < 23000 )) || \
    ( defined (  __IBMCPP__ ) && (  __IBMCPP__ < 400 ) && !defined(__OS400__) )

// AIX xlC 3.1 , 3.0.1 ==0x301
// Visual Age C++ 3.x
// OS-390 C++
#  include <config/stl_ibm.h>

# elif defined(__STL_MSVC)

// Microsoft Visual C++ 4.0, 4.1, 4.2, 5.0
#  include <config/stl_msvc.h>

# elif defined ( __BORLANDC__ )

// Borland C++ ( 4.x - 5.x )
#  include <config/stl_bc.h>

# elif defined(__SUNPRO_CC)

// SUN CC 4.0.1-5.0 
#  include <config/stl_sunpro.h>

# elif defined (__GNUC__ )

// g++ 2.7.x and above 
#  include <config/stl_gcc.h>

# elif defined (__WATCOM_CPLUSPLUS__)

// Watcom C++
#  include <config/stl_watcom.h>

# elif defined(__COMO__) || defined (__COMO_VERSION_)

#  include <config/stl_como.h>

# elif defined (__SC__) && (__SC__ < 0x800)		

// Symantec 7.5
#  include <config/stl_symantec.h>

# elif defined (__MRC__) || (defined (__SC__) && (__SC__ >= 0x882))

// Apple MPW SCpp 8.8.2  
// Apple MPW MrCpp 4.1.0
#  include <config/stl_apple.h>

# elif defined (__MWERKS__)

// Metrowerks CodeWarrior
#  include <config/stl_mwerks.h>

# elif defined(__hpux) && !defined(__GNUC__)

// HP compilers
#  include <config/stl_hpacc.h>

# elif defined(__ICL)

// Intel compiler
#  include <config/stl_intel.h>

// SCO UDK 7 compiler (UnixWare 7x, OSR 5, UnixWare 2x)
#elif defined(__USLC__)

#  include <config/stl_sco.h>

// Apogee 4.x, use "configure" for older versions
#elif defined (__APOGEE__) 

#  include <config/stl_apcc.h>

#elif defined (__KCC)

#  include <config/stl_kai.h>

#elif defined (__DECCXX)

#  include <config/stl_dec.h>
 
#else

// Unable to identify the compiler, issue error diagnostic.
// Edit <config/stl_mycomp.h> to set STLport up for your compiler.

#  include <config/stl_mycomp.h>

# endif /* end of compiler choice */

# undef __GIVE_UP_WITH_STL

#endif
