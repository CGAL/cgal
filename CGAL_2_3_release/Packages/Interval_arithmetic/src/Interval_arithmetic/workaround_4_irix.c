/*
 * This is to workaround various problems with IRIX 5.x and MipsPro.
 * This file must be compiled with GCC on IRIX, as a C file (not C++),
 * then linked into libCGAL.{so|a}.
 *
 * Sylvain.Pion@sophia.inria.fr, April 1999.
 */

void CGAL_workaround_IRIX_set_FPU_cw (int cw)
{
  asm volatile ("ctc1 %0,$31" : :"r" (cw));
}

int CGAL_workaround_IRIX_get_FPU_cw (void)
{
  int cw;
  asm volatile ("cfc1 %0,$31" : "=r" (cw));
  return cw;
}
