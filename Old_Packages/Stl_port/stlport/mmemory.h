CAUUSE A SYNTAX ERROR
/*
 *	mmemory.h
 *	
 * Copyright (c) 1997
 * Mark of the Unicorn, Inc.
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Mark of the Unicorn makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 * Adaptation note: this file is for use with the Metrowerks Standard Library only!
 */
#ifndef __MMEMORY_H_
#define __MMEMORY_H_
#ifdef __MWERKS__
#include <mcompile.h>
#include <stddef.h>
#include <pair.h>
#include <stdio.h>
#include <stdlib.h>
#include <iterator.h>
#include <utility.h>
#include <algorithm.h>
#include <new.h>
#include <tempbuf.h>

#ifndef DefAllocator
#define DefAllocator allocator
#endif

#include <defalloc.h>
#include <assert.h>
#include <auto_ptr.h>
# else
#  error "<mmemory.h> is intended for use with MetroWerks CodeWarrior only. 
Just remove <mmemory.h> from the adaptation directory if using other 
compiler"
# endif /* __MWERKS */

#endif //__MMEMORY_H_
