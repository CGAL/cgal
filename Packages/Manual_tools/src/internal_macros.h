/**************************************************************************
 
  internal_macros.h
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Internal macro definitions.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#ifndef INTERNAL_MACROS_H
#define INTERNAL_MACROS_H 1

#include <stdlib.h>
#include <mstring.h>
#include <database.h>

void init_internal_macros();

void handleChapter(  const Text& T);
void handleBiblio(  const Text& T);

void handleClassEnvironment();
void handleClassEnd();

void handleClassNameEnd( void);
void handleClassFileEnd( void);

void handleHtmlClassFile( const string& filename, const Text& T);

string handleHtmlIndexC( const string& category, const string& item);
string handleHtmlIndex( const string& category, const string& item);
string handleHtmlCrossLink( string key, bool tmpl_class = false);


/* Functions to manage footnotes. */
/* ============================== */
extern int  pre_footnote_counter;
extern int  main_footnote_counter;
extern int  class_footnote_counter;
extern int* footnote_counter;

extern PQueue<char*>  pre_footnotes;
extern PQueue<char*>  main_footnotes;
extern PQueue<char*>  class_footnotes;
extern PQueue<char*>* footnotes;

void insertFootnote( char* s);
// increments counter and returns current value.
int nextFootnoteCounter();
// format footnote reference and hyperlink based on actual counter 
// into a Buffer.
TextToken* formattedFootnoteNumber();
// prints footnote reference and hyperlink based on actual counter.
void printFootnoteCounter( ostream& out);
// prints footnotes and resets counter.
void printFootnotes( ostream& out);

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

void write_headers_to_index( ostream& out);

#endif // INTERNAL_MACROS_H //
// EOF //

