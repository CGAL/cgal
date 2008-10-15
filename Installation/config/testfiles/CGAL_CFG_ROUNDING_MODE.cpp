// Copyright (c) 1998-2007  Utrecht University (The Netherlands),
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
// $URL$
// $Id$
//
//
// Author(s)     : Sylvain Pion



// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by install_cgal.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler has problems with restoring the rounding mode in a try/catch
//| CGAL_CFG_ROUNDING_MODE 


#if defined __alpha__  && defined __linux__
extern "C" {
#  include <fenv.h>
}
#elif defined __powerpc__
#  include <fpu_control.h>
#elif defined __SUNPRO_CC && defined __sun
#  include <ieeefp.h>
#elif defined __osf || defined __osf__
#  ifdef __GNUG__
     // GCC seems to remove (fixincludes) read_rnd/write_rnd...
#    include "/usr/include/float.h"
#  else
#    include <cfloat>
#  endif
#elif defined __sparc__ || (defined __i386__ && !defined __PGI && !defined __SUNPRO_CC)
   // Nothing to include.
#elif defined _MSC_VER
#include <float.h>

#else
   // By default we use the ISO C99 version.
#  include <fenv.h>
#endif




// Pure and safe SSE2 mode (g++ -mfpmath=sse && (-msse2 || -march=pentium4))
// can be detected by :
// TODO : see what Intel and VC++ have to say about this.
#if defined __FLT_EVAL_METHOD__ && defined __SSE2_MATH__ && \
      (__FLT_EVAL_METHOD__ == 0 || __FLT_EVAL_METHOD__ == 1)
#  define CGAL_SAFE_SSE2
#  include <xmmintrin.h>
#endif

// We do not handle -mfpmath=387,sse yet.
#if defined __SSE2_MATH__ && \
    ! (__FLT_EVAL_METHOD__ == 0 || __FLT_EVAL_METHOD__ == 1)
#  warning Unsafe SSE2 mode : not supported yet.
#endif




#if defined __i386__ && !defined __PGI && !defined __SUNPRO_CC

#  if defined CGAL_SAFE_SSE2

#define CGAL_IA_SETFPCW(CW) _MM_SET_ROUNDING_MODE(CW)
#define CGAL_IA_GETFPCW(CW) CW = _MM_GET_ROUNDING_MODE()
typedef unsigned int FPU_CW_t;
#define CGAL_FE_UPWARD       _MM_ROUND_UP

#  else
// The GNU libc version (cf powerpc) is nicer, but doesn't work on libc 5 :(
// This one also works with CygWin.
// Note that the ISO C99 version may not be enough because of the extended
// mantissa issue on x86 (may be required by some kinds of computation, but
// as far as CGAL::Interval_nt<> is concerned, the double-rounding issues
// are taking care of there).
#define CGAL_IA_SETFPCW(CW) asm volatile ("fldcw %0" : :"m" (CW))
#define CGAL_IA_GETFPCW(CW) asm volatile ("fnstcw %0" : "=m" (CW))
typedef unsigned short FPU_CW_t;
#define CGAL_FE_UPWARD       (0x800 | 0x127f)

#  endif

#elif defined __powerpc__
#define CGAL_IA_SETFPCW(CW) _FPU_SETCW(CW)
#define CGAL_IA_GETFPCW(CW) _FPU_GETCW(CW)
typedef fpu_control_t FPU_CW_t;
#define CGAL_FE_UPWARD       (_FPU_RC_UP      | _FPU_DEFAULT)

#elif defined __SUNPRO_CC && defined __sun
#define CGAL_IA_SETFPCW(CW) fpsetround(fp_rnd(CW))
#define CGAL_IA_GETFPCW(CW) CW = fpgetround()
typedef unsigned int FPU_CW_t;
#define CGAL_FE_UPWARD       FP_RP

#elif defined __sparc__
#define CGAL_IA_SETFPCW(CW) asm volatile ("ld %0,%%fsr" : :"m" (CW))
#define CGAL_IA_GETFPCW(CW) asm volatile ("st %%fsr,%0" : "=m" (CW))
typedef unsigned int FPU_CW_t;
#define CGAL_FE_UPWARD       (0x80000000 | 0x20000000 | 0x1f)

#elif defined __mips__
#define CGAL_IA_SETFPCW(CW) asm volatile ("ctc1 %0,$31" : :"r" (CW))
#define CGAL_IA_GETFPCW(CW) asm volatile ("cfc1 %0,$31" : "=r" (CW))
typedef unsigned int FPU_CW_t;
#define CGAL_FE_UPWARD       (0x2)

#elif defined __osf || defined __osf__  // Not yet supported.
#define CGAL_IA_SETFPCW(CW) write_rnd(CW)
#define CGAL_IA_GETFPCW(CW) CW = read_rnd()
typedef unsigned int FPU_CW_t;
#define CGAL_FE_UPWARD       FP_RND_RP

#elif defined ( _MSC_VER )
#if ( _MSC_VER < 1400)
#define CGAL_IA_SETFPCW(CW) _controlfp (CW, _MCW_RC )
#define CGAL_IA_GETFPCW(CW) CW = _controlfp (0, 0 ) &  _MCW_RC
typedef unsigned short FPU_CW_t;
#else
#define CGAL_IA_SETFPCW(CW) unsigned int dummy; _controlfp_s (&dummy, CW, _MCW_RC )
#define CGAL_IA_GETFPCW(CW)_controlfp_s (&CW, 0, 0 ); CW  &=  _MCW_RC
typedef unsigned int FPU_CW_t;
#endif

#define CGAL_FE_UPWARD       _RC_UP

#else
// This is a version following the ISO C99 standard, which aims at portability.
// The drawbacks are speed on one hand, and also, on x86, it doesn't fix the
// extended mantissa issue (this is not a problem for IA, but it is one for
// some future modular computations).
#define CGAL_IA_SETFPCW(CW)  fesetround(CW)
#define CGAL_IA_GETFPCW(CW)  CW = fegetround()
typedef int FPU_CW_t;
#define CGAL_FE_UPWARD       FE_UPWARD
#endif


inline
FPU_CW_t
FPU_get_cw (void)
{
    FPU_CW_t cw;
    CGAL_IA_GETFPCW(cw);
    return cw;
}

inline
void
FPU_set_cw (FPU_CW_t cw)
{
    CGAL_IA_SETFPCW(cw);
}

inline
FPU_CW_t
FPU_get_and_set_cw (FPU_CW_t cw)
{
    FPU_CW_t old = FPU_get_cw();
    FPU_set_cw(cw);
    return old;
}


// A class whose constructor sets the FPU mode to +inf, saves a backup of it,
// and whose destructor resets it back to the saved state.


struct Protect_FPU_rounding{
  Protect_FPU_rounding(FPU_CW_t r = CGAL_FE_UPWARD)
    : backup( FPU_get_and_set_cw(r) ) {}

  ~Protect_FPU_rounding()
  {
     FPU_set_cw(backup);
  }


  FPU_CW_t backup;
};




int main()
{
  try {
    Protect_FPU_rounding p;
    throw 1;
  } catch(...){}

  if(FPU_get_cw() == CGAL_FE_UPWARD){
    return 1;
  }
  return 0;
}
