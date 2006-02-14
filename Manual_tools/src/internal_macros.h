/**************************************************************************
 
  internal_macros.h
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Internal macro definitions.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

#ifndef INTERNAL_MACROS_H
#define INTERNAL_MACROS_H 1

#include <stdlib.h>
#include <mstring.h>

void init_internal_macros();


void handleIndex();
void handleIndex2(string main_item, string sub_item, string sub_sub_item, int modifier);
void OpenFileforIndex();
void CloseFileforIndex();
void handleIndexTraitsClass();
void handleIndexRefName();


#endif // INTERNAL_MACROS_H //
// EOF //

