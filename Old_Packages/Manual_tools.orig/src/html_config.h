/**************************************************************************
 
  html_config.h
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : Configuration constants and variables
  System    : C++ (g++)
  Author    : (c) 1997 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#if ! defined( MODULE_CONFIG)
#define MODULE_CONFIG 1

#include <buffer.h>
#include <mstring.h>


/* Configurable command line options */
/* ================================= */
typedef char Switch;
 
#define NO_SWITCH    0
#define MINUS_SWITCH 1
#define PLUS_SWITCH  2
 
extern Switch  trace_switch;
extern Switch  line_switch;

extern Switch  config_switch;
extern Switch  quiet_switch;
extern Switch  macro_def_switch;
extern Switch  macro_exp_switch;
extern Switch  stack_trace_switch;

extern Switch  noheader_switch;
extern Switch  onlyheader_switch;


// Global declarations that are implemented in the main module.
// There they can be taylored to the specific application, i.e.
// extraction or checker.

extern const string prog_name;
extern const string prog_release;
extern const string reference_icon;


// outdated
int   text_block_length(     const Buffer_list& T);
char* text_block_to_string(  const Buffer_list& T);
bool  is_text_block_empty(   const Buffer_list& T);
void  print_html_text_block( ostream &out, const Buffer_list& T);


void handleText(   const Buffer_list& T);
void handleBuffer( const Buffer&      TT);
void handleString( const char*        s);
void handleString( const string&      s);
void handleChar(   char               c);

#endif // MODULE_CONFIG //
