/**************************************************************************
 
  html_syntax.y
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : grammatical parser for TeX and C++ code mixed files.
              Taylored for HTML manual generation.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1996 Lutz Kettner
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

%{
#include <html_syntax.h>
#include <html_lex.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <string_conversion.h>
#include <html_error.h>
#include <lex_include.h>
#include <macro_dictionary.h>
#include <cpp_formatting.h>
#include <internal_macros.h>

#include <database.h>
#include <html_config.h>


/* Own prototypes */
/* ============== */
int yyerror( char *s);
Text* blockintroProcessing( const char* text, int len, Text* t);
Buffer* blockintroProcessing( const char* text, int len, Buffer* t);

%}

%union {
    struct {
        const char* text;   /* a chunk of zero terminated text */
        int         len;    /* its length */
    }          string;
    char       character;   /* a character */
    int        number;      /* an integer  */
    Buffer*    pBuffer;     /* Buffer for collected strings and characters */
    TextToken* pTextToken;
    Text*      pText;
}

/* Elementary data types */
/* --------------------- */
%token <string>    STRING
%token <character> CHAR

/* Newcommand and parameter parsing */
/* -------------------------------- */
%token <string>    DEFWITHARGS
%token <string>    GDEFWITHARGS
%token <string>    DEFWITHUNKNOWNARGS
%token <string>    PARAMETER
%token <string>    PARAMETER_OPTION

/* Special Scanning Modes */
/* ---------------------- */
%token             ASCIITOHTML
%token             RAWOUTPUT
%token <string>    RAWOUTPUTN

/* Special TeX Parsing    */
/* ---------------------- */
%token             CHAPTER
%token             SECTION
%token <string>    LABEL

%token             TTBLOCKINTRO
%token             EMBLOCKINTRO
%token             ITBLOCKINTRO
%token             SCBLOCKINTRO
%token             BFBLOCKINTRO
%token             RMBLOCKINTRO
%token             SFBLOCKINTRO

/* C++ Parsing            */
/* ---------------------- */
%token             BEGINCLASS
%token             ENDCLASS

/* HTML Specialties       */
/* ---------------------- */
%token <string>    HTMLINDEX

/* File Management        */
/* ---------------------- */
%token             HTMLBEGINCLASSFILE


/* Type declarations for production rules       */
/* -------------------------------------------- */
%type  <string>     blockintro

%type  <pTextToken> string_token 
%type  <pText>      comment_group  comment_sequence  comment_token



%%
/* Grammar: Top Level */
/* ================== */

input:            /* empty */
                | input stmt
;

stmt:             string_token        {   handleText( * $1); delete $1; }
                | CHAPTER             { pushMacroScope(); }
		  comment_sequence
                  '}'                 {   handleChapter( * $3); 
                                          delete $3;
                                      }
                | SECTION             { pushMacroScope(); }
		  comment_sequence
                  '}'                 {   handleSection( * $3); 
                                       ; //delete $3;
                                      }
		| BEGINCLASS          {   handleClassEnvironment(); }
		| ENDCLASS            {   handleClassEnd(); }
		| HTMLBEGINCLASSFILE comment_group 
		                      {  handleHtmlClassFile( cc_filename,*$2);
		                          delete $2;
                                      }
		| group
;

group:            '{'
                  input
		  '}'
                  | blockintro { pushMacroScope(); handleString( $1.text);}
                  input
                  '}'          { handleString( $1.text + strlen( $1.text)+1);}
;

/* Auxiliary Rules                 */
/* =============================== */
blockintro:       TTBLOCKINTRO { $$.text = "<TT>\0</TT>"; $$.len = 4; }
                | EMBLOCKINTRO { $$.text = "<EM>\0</EM>"; $$.len = 4; }
                | ITBLOCKINTRO { $$.text = "<I>\0</I>";   $$.len = 3; }
                | SCBLOCKINTRO { $$.text = "<TT>\0</TT>"; $$.len = -1; }
                | BFBLOCKINTRO { $$.text = "<B>\0</B>";   $$.len = 3; }
                  /* Sorry: \rm not supported. TT might be fine. */
                | RMBLOCKINTRO { $$.text = "<TT>\0</TT>"; $$.len = 4; }
                | SFBLOCKINTRO { $$.text = "<TT>\0</TT>"; $$.len = 4; }
;


string_token:     STRING        { $$ = new TextToken( $1.text); }
                | CHAR          { $$ = new TextToken(); $$->add( $1); }
                | LABEL         {
                                  handleLabel( $1.text, $1.len);
                                  $$ = new TextToken( "<A NAME=\"");
				  $$->add( $1.text, $1.len);
				  $$->add( "\"></A>");
                                }
                | HTMLINDEX  comment_group     {
		                  char* s = text_block_to_string(* $2);
				  delete $2;
                                  $$ = new TextToken();
				  $$->add( handleHtmlIndex( $1.text, s)
					   .c_str());
				  delete[] s;
                                }
                | ASCIITOHTML PARAMETER  {
				  $$ = new TextToken(
				      convert_ascii_to_html($2.text));
				  delete[] $2.text;
				  set_old_state = 1;
				}
                | RAWOUTPUT PARAMETER  {
				  $$ = new TextToken( $2.text);
				  delete[] $2.text;
				  set_old_state = 1;
				}
                | RAWOUTPUTN   {
				  $$ = new TextToken( $1.text);
				  delete[] $1.text;
				  set_old_state = 1;
				}
                | DEFWITHARGS  PARAMETER {
				  $$ = new TextToken();
				  int m = atoi($1.text + $1.len-1);
				  if ( m < 1 || m > 9)
				      printErrorMessage( NParamRangeError);
				  else
				      insertMacro( $1.text,in_string->name(),
						   in_string->line(), 
						   $2.text, m);
				  delete[] $1.text;
				  delete[] $2.text;
				  set_old_state = 1;
				}
                | GDEFWITHARGS  PARAMETER {
				  $$ = new TextToken();
				  int m = atoi($1.text + $1.len-1);
				  if ( m < 1 || m > 9)
				      printErrorMessage( NParamRangeError);
				  else
				      insertGlobalMacro( $1.text,
							 in_string->name(),
							 in_string->line(), 
							 $2.text, m);
				  delete[] $1.text;
				  delete[] $2.text;
				  set_old_state = 1;
				}
                | DEFWITHUNKNOWNARGS  PARAMETER {
				  $$ = new TextToken();
				  if ( ! quiet_switch) {
				      cerr << endl 
					   << "Unknown macro definition: " 
					   << $1.text << " in `" 
					   << in_string->name() << " in line "
					   << in_string->line() << "'.";
				      printErrorMessage( MacroDefUnknownError);
				  }
				  delete[] $1.text;
				  delete[] $2.text;
				  set_old_state = 1;
				}
;


/* A sequence of words forming a comment */
/* ===================================== */
comment_group:      '{' comment_sequence '}'  { $$ = $2; }
                   | blockintro { pushMacroScope(); }
                     comment_sequence '}'  { 
                              $$ = blockintroProcessing( $1.text, $1.len, $3); 
                           }
;

comment_sequence:        /* empty */  { $$ = new Text( managed); }
                       | comment_sequence comment_token {
		                  $$ = $1;
				  $$->append( * $2);
				  delete $2;
		                }
;

comment_token:      string_token  { $$ = new Text( * $1, managed); }
                  | comment_group { $$ = $1; }
;


/* End of Grammar */
/* ============== */
%%

// support functions
// -----------------
Text* blockintroProcessing( const char* text, int len, Text* t) { 
    if ( len < 0) {  /* Hack! Here we know that t has to get capitalized.*/
        len = strlen(text);
        InListFIter< TextToken> ix( *t);
        ForAll( ix) {
	    if ( ! (*ix).isSpace) {
	        char *s = (*ix).string;
		while ( *s) {
                    *s = toupper( *s);
		    s++;
		}
	    }
        }
    }
    if ( ! t->isEmpty()) {
        t->head().prepend( TextToken( text, len));
	/* Hack! ptr arithmetic points to the closing tag text */
	t->append( * new TextToken( text + len + 1));
    }
    return t;
}

Buffer* blockintroProcessing( const char* text, int len, Buffer* t) { 
    if ( len < 0) {  /* Hack! Here we know that t has to get capitalized.*/
        len = strlen(text);
        t->capitalize();
    }
    t->prepend( text, len);
    /* Hack! ptr arithmetic points to the closing tag text */
    t->add( text + len + 1);
    return t;
}

