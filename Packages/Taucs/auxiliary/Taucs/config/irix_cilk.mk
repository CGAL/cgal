#########################################################
# Irix Cilk
# Sivan's comments, 4 Sep 2003:
# This file is for use with the current branch of Cilk
# on or.iucc.ac.il
#########################################################
FC        = g77
FFLAGS    = -O3 -fno-second-underscore -Wall
FOUTFLG   = -o ./

CC        = gcc
COUTFLG   = -o ./
CFLAGS    = -O3 -Wall -Werror -std=c89

CILKC      = ../cilk-devel/current/support/cilkclocal
CILKOUTFLG = -o ./
CILKFLAGS  = -O3 -x cilk

LD        = $(CILKC)
LDFLAGS   = 
LOUTFLG   = -o ./

LIBBLAS   = -lcomplib.sgimath
LIBMETIS  = -Lexternal/lib/irix -lmetis32

LIBF77    = -L /usr/local/gcc3/lib/gcc-lib/mips-sgi-irix6.5/3.0.4/mabi=64 -lg2c -lgcc
LIBF77    = 
LIBC      = -lc -lm
