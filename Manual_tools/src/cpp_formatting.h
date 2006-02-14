/**************************************************************************
 
  cpp_formatting.h
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : C++ Formatting functions.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

#ifndef CPP_FORMATTING_H
#define CPP_FORMATTING_H 1

#include <mstring.h>

extern string class_name;
extern string template_class_name;

extern int HREF_counter;

void handle_two_column_layout( char key, const char* decl);
void handle_three_column_layout( char key, const char* decl, bool empty);
void  make_index(string main_item, string sub_item, string sub_sub_item,
                 const char* signature);



/* Flexibility for HTML class files. */
/* ================================= */
extern bool html_no_class_toc;
extern bool html_no_class_file;
extern bool html_no_class_links;
extern bool html_no_class_index;
extern bool html_no_links;
extern bool html_no_index;

extern bool html_inline_classes;


#endif // CPP_FORMATTING_H 1 //
// EOF //

