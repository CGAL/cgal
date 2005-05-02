#########################################################
# Solaris Cilk
# Sivan's comments, 4 Sep 2003:
# This file is for use with the current branch of Cilk
# on pisces.math.tau.ac.il
#########################################################
FC        = g77
FFLAGS    = -O3 -fno-second-underscore -Wall
FOUTFLG   = -o ./

CC        = gcc
COUTFLG   = -o ./
CFLAGS    = -O3 -Wall -Werror -std=c89

CILKC      = /export/home/stoledo/cilk-devel/current/support/cilkclocal
CILKOUTFLG = -o ./
CILKFLAGS  = -O3 -x cilk

LD        = $(CILKC)
LDFLAGS   = 
LOUTFLG   = -o ./

LIBMETIS  = -Lexternal/lib/solaris -lmetis

# the -lpthread is to compensate in some bug in cilkclocal
# that Bradley reported to me.

LIBF77    = 
LIBC      = -lpthread -lc -lm

#########################################################








