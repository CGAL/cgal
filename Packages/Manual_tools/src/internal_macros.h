/**************************************************************************
 
  internal_macros.h
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Internal macro definitions.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#ifndef INTERNAL_MACROS_H
#define INTERNAL_MACROS_H 1

#include <stdlib.h>
#include <mstring.h>
#include <buffer.h>

void init_internal_macros();

void handleChapter(  const Buffer_list& T);
void handlePart(  const Buffer_list& T);
void handleBiblio(  const Buffer_list& T);

void handleIndex();
void handleIndex2(string main_item, string sub_item, string sub_sub_item, int modifier);
void OpenFileforIndex();
void CloseFileforIndex();
void handleIndexTraitsClass();
void handleIndexRefName();


void handleClassEnvironment();
void handleClassEnd();

void handleClassNameEnd( void);
void handleClassFileEnd( void);

void handleHtmlClassFile( const string& filename, const Buffer_list& T);

string handleHtmlCrossLink( string key, bool tmpl_class = false);

// Index sorting.
// ==============
// sort_keys are used to sort the index according to different sections.
// A sort_key is followed by a 0 to indicate the section title.
// A sort_key is followed by a 1 to indicate a normal entry.

extern const string sort_key_class;
extern const string sort_key_nested_type;
extern const string sort_key_struct;
extern const string sort_key_enum;
extern const string sort_key_enum_tags;
extern const string sort_key_typedef;
extern const string sort_key_variable;
extern const string sort_key_function;
extern const string sort_key_member_function;


const string& find_sort_key( string txt);

#endif // INTERNAL_MACROS_H //
// EOF //

