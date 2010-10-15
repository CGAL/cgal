/**************************************************************************
 
  html_lex.h
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : External symbols of the HTML scanner.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

#ifndef HTML_LEX_H
#define HTML_LEX_H 1

#include <stdio.h>
#include <mstring.h>
#include <output.h>

const char SEPARATOR = '\001';

extern int set_old_state;

extern string cc_string;
extern string cc_filename;

extern int number_of_args;

extern char *yytext;

extern "C" {
    int yylex( void);
    void init_scanner( FILE* in);
}

#endif // HTML_LEX_H 1 //
// EOF //

