/**************************************************************************
 
  syntax.y
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : grammatical parser for TeX and C++ code mixed files.
              Taylored for HTML manual generation.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1996 Lutz Kettner
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

%{
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

/* Declarations from lex.yy */
/* ======================== */
extern int set_CCMode;
extern int set_NestingMode;
extern int set_INITIAL;
extern int set_HTMLMODE;
extern int set_MMODE;
extern int line_number;

extern const char* in_filename;
extern       char* creationvariable;

extern char *yytext;

extern "C" {
int yylex( void);
void init_scanner( FILE* in);
}

/* Declarations for the parser */
/* =========================== */
/* This variable flags for bison that we are in the CCMode */
int CCMode = 0;

/* Datastructures for the parser */
/* ============================= */
#include <database.h>
#include <confightml.h>

/* Declarations from the cgal_extract_html.cc file */
/* =============================================== */
extern char* class_name;
extern char* formatted_class_name;
char* text_block_to_string( const Text& T);

/* for the bibliography */
/* ==================== */
extern bool first_bibitem;

/* Own prototypes */
/* ============== */
int yyerror( char *s);
Text* blockintroProcessing( const char* text, int len, Text* t);

%}

%union {
    struct {
        const char* text;   /* a (meaningless) chunk of zero terminated text */
        int         len;    /* its length */
    }          string;
    char       character;   /* a (meaningless) character */
    Buffer*    pBuffer;     /* Buffer for collected strings and characters */
    TextToken* pTextToken;
    Text*      pText;
}

/* Elementary data types */
/* --------------------- */
%token <string>    STRING
%token <string>    SPACE
%token <character> CHAR
%token             NEWLINE

/* Keywords to trigger on */
/* ---------------------- */
%token             CHAPTER
%token             SECTION
%token             SUBSECTION
%token             SUBSUBSECTION
%token             BEGINBIBLIO
%token             ENDBIBLIO
%token <string>    BIBITEM
%token <string>    CITE
%token <string>    LABEL
%token             BEGINCLASS
%token             ENDCLASS
%token             BEGINCLASSTEMPLATE
%token             ENDCLASSTEMPLATE
%token <string>    CREATIONVARIABLE
%token             CONSTRUCTOR
%token             METHOD
%token             FUNCTION
%token             FUNCTIONTEMPLATE
%token             VARIABLE
%token             TYPEDEF
%token             ENUM
%token             GLOBALFUNCTION
%token             GLOBALFUNCTIONTEMPLATE
%token             GLOBALVARIABLE
%token             GLOBALTYPEDEF
%token             GLOBALENUM
%token             DECLARATION

/* Special action keywords */
/* ----------------------- */
%token             CPROGBEGIN
%token             CPROGEND
%token             HIDDEN

%token             TEXONLYBEGIN
%token             TEXONLYEND
%token             LATEXHTML
%token             ANCHOR
%token             HTMLBEGIN
%token             HTMLEND

%token             CCSTYLE
%token             CCSECTION
%token             GOBBLEONEPARAM
%token             GOBBLETWOPARAMS
%token             IGNOREBLOCK
%token             IGNORETWOBLOCKS
%token             TTBLOCKINTRO
%token             EMBLOCKINTRO
%token             ITBLOCKINTRO
%token             SCBLOCKINTRO
%token             BFBLOCKINTRO

%token             BEGINMATH
%token             ENDMATH
%token <character> SINGLESUBSCRIPT
%token <character> SINGLESUPERSCRIPT
%token             BEGINSUBSCRIPT
%token             BEGINSUPERSCRIPT
%token             FRACTION

/* handle LALR(1) restriction */
/* -------------------------- */
%token             LALRRESTRICTION

%type  <string>     blockintro

%type  <pBuffer>    string  string_token string_with_nl string_with_nl_token
%type  <pBuffer>    string_with_nl_or_mt string_with_spcnl_token
%type  <pBuffer>    verbatim_style
%type  <pBuffer>    math_sequence   math_token
%type  <pBuffer>    declaration classname template_params
%type  <pBuffer>    cc_stmts cc_stmt cc_stmts_skip_space
%type  <pText>      comment_group  comment_sequence  
%type  <pText>      nested_token_sequence nested_token
%type  <pText>      compound_comment  full_comment_sequence  
%type  <pText>      non_empty_comment_sequence

%type  <pText>      whitespaces  optional_whitespaces

%type  <pTextToken> comment_token  non_empty_token  whitespace



%%
/* Grammar: Top Level */
/* ================== */

input:            /* empty */
                | input stmt
;

stmt:             string              {   handleBuffer( * $1); delete $1;}
                | whitespaces         {   handleText(   * $1, true); 
		                          delete $1;
		                      }
		| verbatim_style         {   handleBuffer( * $1); delete $1;}
                | CHAPTER
		  comment_sequence
                  '}'                 {   handleChapter( * $2); delete $2;}
		| BEGINCLASS
		  classname           {   handleClass( $2->string());
		                          delete $2;}
		    decl_sequence 
		  ENDCLASS
                                      {
					  handleClassEnd();
                                          free( creationvariable);
                                          creationvariable = NULL;
				      }
		| BEGINCLASSTEMPLATE
		  classname           {   handleClassTemplate( $2->string());
		                          delete $2;}
		  decl_sequence
		  ENDCLASSTEMPLATE
                                      {
					  handleClassTemplateEnd();
                                          free( creationvariable);
                                          creationvariable = NULL;
				      }
		| CREATIONVARIABLE    {}
		| CCSTYLE  '{' nested_token_sequence '}'  {
					  set_INITIAL = 1;
					  handleString( "<VAR>");
					  handleText( * $3);
					  handleString( "</VAR>");
					  delete $3;
		                      }
                | global_tagged_declarator
		| group
;

group:            '{'
                  input
		  '}'
                | blockintro
                  input
                  '}'         {}
;

/* Auxiliary Rules                 */
/* =============================== */
blockintro:       TTBLOCKINTRO { $$.text = "<TT>\0</TT>"; $$.len = 4; }
                | EMBLOCKINTRO { $$.text = "<EM>\0</EM>"; $$.len = 4; }
                | ITBLOCKINTRO { $$.text = "<I>\0</I>";   $$.len = 3; }
                | SCBLOCKINTRO { $$.text = "<TT>\0</TT>"; $$.len = -1; }
                | BFBLOCKINTRO { $$.text = "<B>\0</B>";   $$.len = 3; }
;


string_with_nl_or_mt:   { $$ = new Buffer; /* Empty */ }
                | string_with_nl
;

string_with_nl:   string_with_nl_token
                | string_with_nl  string_with_nl_token {
				  $$ = $1;
		                  $$->add( $2);
				  delete $2;
		                }
;

string:           string_token
                | string string_token {
				  $$ = $1;
		                  $$->add( $2);
				  delete $2;
		                }
;

string_with_spcnl_token: 
		  string_token
                | whitespace    {
                                  $$ = new Buffer;
				  $$->add( $1->string, $1->len);
				  delete $1;
                                }
;
string_with_nl_token: 
		  string_token
                | NEWLINE       {
                                  $$ = new Buffer;
				  $$->add( '\n');
                                }
;
string_token:     STRING       {
                                  $$ = new Buffer;
				  $$->add( $1.text, $1.len);
                                }
                | LABEL         {
                                  handleLabel( $1.text);
                                  $$ = new Buffer;
				  $$->add( "<A NAME=\"");
				  $$->add( $1.text, $1.len);
				  $$->add( "\"></A>");
                                }
                | CITE          { $$ = handleCite( $1.text); }
                | BIBITEM  '{' nested_token_sequence '}'     {
		                  set_INITIAL = 1;
		                  char* s = text_block_to_string(* $3);
                                  $$ = handleBibItem( s);
				  delete[] s;
				  delete $3;
				}
                | BIBITEM  '[' nested_token_sequence ']'
		           '{' nested_token_sequence '}'     {
		                  set_INITIAL = 1;
		                  char* s = text_block_to_string(* $3);
		                  char* p = text_block_to_string(* $6);
                                  $$ = handleBibItem( p, s);
				  delete[] s;
				  delete[] p;
				  delete $3;
				  delete $6;
				}
                | BEGINMATH math_sequence ENDMATH {
		                  $$ = $2;
				  $$->prepend( "<MATH>", 6);
				  $$->add( "</MATH>", 7);
		                }
                | CHAR          {
                                  $$ = new Buffer;
				  $$->add( $1);
                                }
;

non_empty_token:    string      { $$ = new TextToken( 
					       $1->string(), $1->length());
                                  delete $1;
                                }
;

optional_whitespaces:  /* empty */    { $$ = new Text( managed); }
                  | whitespaces       { $$ = $1; }
;

whitespaces:        whitespace  { $$ = new Text( * $1, managed); }
                  | whitespaces whitespace {
                                  $$ = $1;
                                  $$->append( * $2);
                                }
;

whitespace:         SPACE       { $$ = new TextToken( $1.text, $1.len, true); }
                  | NEWLINE     { $$ = new TextToken( "\n", 1, true); }
		  | GOBBLETWOPARAMS comment_group comment_group
                                { $$ = new TextToken( " ", 1, true);
				  delete $2; 
				  delete $3; 
				}
		  | GOBBLEONEPARAM  comment_group
                                { $$ = new TextToken( " ", 1, true); 
				  delete $2;
				}
                  | IGNOREBLOCK comment_sequence '}'  { 
                                  $$ = new TextToken( " ", 1, true);
				  delete $2;
				  set_INITIAL = 1;
                                }
                  | IGNORETWOBLOCKS comment_sequence '}'
		                '{' comment_sequence '}'  { 
                                  $$ = new TextToken( " ", 1, true);
				  delete $2;
				  delete $5;
				  set_INITIAL = 1;
                                }
                  | texonly_style
                                { $$ = new TextToken( " ", 1, true); 
				}
;



/* Class Declaration with Comments */
/* =============================== */
decl_sequence:    comment_sequence  {
				  handleMainComment( * $1);
				  delete $1;
		                }
		| decl_sequence
		  tagged_declarator 
		  comment_sequence {
				  handleMainComment( * $3);
				  delete $3;
		                }
;

tagged_declarator:
		  global_tagged_declarator
		| CONSTRUCTOR   declaration   comment_group {
		                  handleConstructorDeclaration( $2->string(),
								* $3);
				  delete $2;
				  delete $3;
		                }
		| METHOD        declaration   comment_group {
		                  handleMethodDeclaration( $2->string(), * $3);
				  delete $2;
				  delete $3;
		                }
;

global_tagged_declarator:
                  SECTION
		  comment_sequence
                  '}'                 {   handleSection( * $2); delete $2;}
                | SUBSECTION
		  comment_sequence
                  '}'                 {   
		                          handleString( "<H3>");
		                          handleText( * $2);
		                          handleString( "</H3>\n");
					  delete $2;
		                      }
                | SUBSUBSECTION
		  comment_sequence
                  '}'                 {   
		                          handleString( "<H4>");
		                          handleText( * $2);
		                          handleString( "</H4>\n");
					  delete $2;
		                      }
                | BEGINBIBLIO
		  comment_sequence
                  ENDBIBLIO           {   handleBiblio( * $2); delete $2;}
		| FUNCTION      declaration   comment_group {
		                  handleFunctionDeclaration( $2->string(),
							     * $3);
				  delete $2;
				  delete $3;
		                }
		| FUNCTIONTEMPLATE 
                      template_params
                      optional_whitespaces
		      declaration 
		      comment_group {
		                  handleFunctionTemplateDeclaration(
					  $2->string(),
					  $4->string(),
					  * $5);
				  delete $2;
				  delete $3;
				  delete $4;
				  delete $5;
		                }
 		| VARIABLE      declaration   comment_group {
		                  handleVariableDeclaration( $2->string(),
							     * $3);
				  delete $2;
				  delete $3;
		                }
 		| TYPEDEF       declaration   comment_group {
		                  handleVariableDeclaration( $2->string(), * $3);
				  delete $2;
				  delete $3;
		                }
 		| ENUM          declaration   comment_group {
		                  handleEnumDeclaration( $2->string(), * $3);
				  delete $2;
				  delete $3;
		                }
		| GLOBALFUNCTION      declaration   {
		                  handleFunctionDeclaration( $2->string());
				  delete $2;
		                }
		| GLOBALFUNCTIONTEMPLATE 
		      template_params 
		      optional_whitespaces 
		      declaration  {
		                  handleFunctionTemplateDeclaration( 
					  $2->string(),
					  $4->string());
				  delete $2;
				  delete $3;
				  delete $4;
		                }
 		| GLOBALVARIABLE      declaration   {
		                  handleVariableDeclaration( $2->string());
				  delete $2;
		                }
 		| GLOBALTYPEDEF       declaration   {
		                  handleVariableDeclaration( $2->string());
				  delete $2;
		                }
 		| GLOBALENUM          declaration   {
		                  handleEnumDeclaration( $2->string());
				  delete $2;
		                }
 		| DECLARATION   declaration {
		                  handleDeclaration( $2->string());
				  delete $2;
		                }
                | HIDDEN 
		  optional_whitespaces 
		  hidden_keys      
		  declaration   
		  comment_group {
		                  delete $2;
		                  delete $4;
		                  delete $5;
		                }
;

hidden_keys:      CONSTRUCTOR
                | METHOD
                | FUNCTION
                | VARIABLE
                | TYPEDEF
                | ENUM
;

/* A sequence of words forming a comment */
/* ===================================== */
comment_group:      optional_whitespaces '{' comment_sequence '}'  { 
                                  $$ = $3; 
				  delete $1;
                                }
                  | optional_whitespaces blockintro comment_sequence '}'  { 
                                  $$ = blockintroProcessing( $2.text, 
							     $2.len,
							     $3); 
				  delete $1;
                                }
;

comment_sequence:   optional_whitespaces { $$ = new Text( managed); }
                  | optional_whitespaces
                    non_empty_comment_sequence
                    optional_whitespaces { $$ = $2; }
;

full_comment_sequence:   /* empty */  { $$ = new Text( managed); }
                  | whitespaces       { $$ = $1; }
                  | optional_whitespaces
                    non_empty_comment_sequence
                    optional_whitespaces  { 
		                  $$ = $1;
				  $$->append( * $2);
				  $$->append( * $3);
		                }
;

non_empty_comment_sequence:
		    comment_token     { $$ = new Text( * $1, managed); }
                  | compound_comment  { $$ = $1; }
                  | non_empty_comment_sequence optional_whitespaces comment_token {
		                  $$ = $1;
				  $$->append( * $2);
				  $$->append( * $3);
		                }
                  | non_empty_comment_sequence 
		    optional_whitespaces
		    compound_comment {
		                  $$ = $1;
				  $$->append( * $2);
				  $$->append( * $3);
		                }
;

comment_token:      non_empty_token   { $$ = $1; }
                  | comment_token non_empty_token {
		                  $$ = $1;
				  $$->add( * $2);
				  delete $2;
		                }
;

compound_comment:   '{' full_comment_sequence '}' {
				  $$ = $2;
		                  /* $$->cons(   *new TextToken( "{", 1)); */
		                  /* $$->append( *new TextToken( "}", 1)); */
		                }
                  | '(' full_comment_sequence ')' {
				  $$ = $2;
		                  $$->cons(   *new TextToken( "(", 1));
		                  $$->append( *new TextToken( ")", 1));
		                }
                  | '[' full_comment_sequence ']' {
				  $$ = $2;
		                  $$->cons(   *new TextToken( "[", 1));
		                  $$->append( *new TextToken( "]", 1));
		                }
                  | blockintro full_comment_sequence '}' {
                                  $$ = blockintroProcessing( $1.text,
							     $1.len, 
							     $2); 
		                }
                  | CCSTYLE '{' nested_token_sequence '}'  { 
				  $$ = $3;
		                  set_INITIAL = 1;
                                  if ( $$->isEmpty() || 
				       $$->head().isSpace)  // should not
		                      $$->cons(   *new TextToken( "<VAR>", 1));
                                  else
		                      $$->head().prepend( "<VAR>");
                                  InListFIter< TextToken> ix( * $$);
				  ForAll( ix) {
				      if ( ix.isLast())
					  if ( ix->isSpace)
					      $$->append( *new TextToken(
                                                              "</VAR>", 1));
					  else
					      ix->add( "</VAR>");
				  }
                                }
                  | verbatim_style {
		                  $$ = new Text( managed);
				  $$->cons(   *new TextToken( $1->string()));
		               }
                  | CCSECTION '{' comment_sequence '}'  {
				  $$ = $3;
		                  $$->cons(   *new TextToken( " ", 1, true));
		                  $$->cons(   *new TextToken( "<H1>"));
		                  $$->cons(   *new TextToken( "\n", 1, true));
		                  $$->append( *new TextToken( " ("));
		                  $$->append( *new TextToken( 
						       formatted_class_name));
		                  $$->append( *new TextToken( ")</H1>"));
		                  $$->append( *new TextToken( "\n", 1, true));
		                }
                  | CREATIONVARIABLE   { $$ = new Text( managed);}
;

/* Parsing of a C++ expression/statement with nested expressions */
/* ============================================================= */
nested_token_sequence:
		    /* empty */ {
		                  $$ = new Text(managed);
		                }
		  | nested_token_sequence nested_token
		                {
				  $1->append( * $2);
				  $$ = $1;
				}
;

nested_token:       string      {
                                  $$ = new Text(*new TextToken( 
						    $1->string(),
						    $1->length()),
						managed);
				  delete $1;
                                }
                  | SPACE       {
                                  $$ = new Text(*new TextToken( 
						    $1.text,
						    $1.len,
						    true),
						managed);
		  }
                  | NEWLINE     {
                                  $$ = new Text(*new TextToken( "\n", 1, true),
						managed);
                                }
		  | '{' nested_token_sequence '}' {
		                  $2->cons(   *new TextToken( "{", 1));
		                  $2->append( *new TextToken( "}", 1));
				  $$ = $2;
		                }
		  | blockintro nested_token_sequence '}' {
                                  $$ = blockintroProcessing( $1.text,
							     $1.len, 
							     $2); 
		                }
		  | '[' nested_token_sequence ']' {
		                  $2->cons(   *new TextToken( "[", 1));
		                  $2->append( *new TextToken( "]", 1));
				  $$ = $2;
		                }
		  | '(' nested_token_sequence ')' {
		                  $2->cons(   *new TextToken( "(", 1));
		                  $2->append( *new TextToken( ")", 1));
				  $$ = $2;
		                }
;

/* Parsing of a C++ Declaration (function, method ..., not class) */
/* ============================================================== */
declaration:      '{'           { 
                                  CCMode = 1; 
                                }
                  cc_stmts_skip_space
		  '}'           { 
		                  set_INITIAL = 1;
				  CCMode = 0;
				  $$ = $3;
		                }
;

classname:        '{'           {
                                  CCMode = 1; 
                                }
                  cc_stmts_skip_space
		  '}'           { 
		                  set_INITIAL = 1;
				  CCMode = 0;
				  $$ = $3;
		                }
;

template_params: '{'           {
                                  CCMode = 1; 
                                }
                  cc_stmts_skip_space
		  '}'           { 
		                  /* set_INITIAL = 1; */
				  CCMode = 0;
				  $$ = $3;
		                }
;

cc_stmts:         /* empty */
                                { $$ = new Buffer;}
		  | cc_stmts cc_stmt {
		                  $$ = $1;
				  $$->add( $2);
				  delete $2;
		                }
;

cc_stmt:          string        { $$ = $1;
                                }
                | SPACE         { $$ = new Buffer;
                                  $$->add( ' ');
		                }
                | NEWLINE       { $$ = new Buffer;
                                  $$->add( ' ');
		                }
		| '{'
                  cc_stmts
		  '}'           {
		                  $$ = $2;
				  $$->prepend( '{');
				  $$->add( '}');
		                }
;

cc_stmts_skip_space: 
		  /* empty */
                                { $$ = new Buffer;}
		| string
		  cc_stmts      { $$ = $1;
				  $$->add( $2);
				  delete $2;
		                }
                | SPACE         cc_stmts { $$ = $2;}
                | NEWLINE       cc_stmts { $$ = $2;}
		| '{'
                  cc_stmts
		  '}'
                  cc_stmts      { 
		                  $$ = $2;
				  $$->prepend( '{');
				  $$->add( '}');
				  $$->add( $4);
				  delete $4;
		                }
;

/* Parsing of the CPROG environment and other verbatim environments */
/* ================================================================ */
verbatim_style:   CPROGBEGIN string_with_nl_or_mt CPROGEND {
                                  $$ = $2;
				  $$->prepend( "<PRE>" , 5);
				  $$->add(     "</PRE>", 6);
                                }
                | HTMLBEGIN string_with_nl_or_mt HTMLEND {
                                  $$ = $2;
                                }
                | LATEXHTML 
                  comment_sequence 
                  '}' '{'       { 
                                  delete $2;
				  set_HTMLMODE = 1;
                                }
                  string_with_nl_or_mt '}'
                                {
				  $$ = $6;
				}
                | ANCHOR
                  string_with_nl_or_mt 
		  '}'
		  comment_group
                    {
		        $$ = $2;
			$$->prepend( "<A HREF=\"");
			$$->add(     "\">");
			char* s = text_block_to_string( * $4);
			$$->add( s);
			$$->add(     "</A>");
			delete s;
			delete $4;
		    }
;
texonly_style:    TEXONLYBEGIN string_with_nl TEXONLYEND {
                                  delete $2;
                                }
;

/* Parsing of mathematical formulas from TeX */
/* ========================================= */
math_sequence:
      /* empty */
        {
            $$ = new Buffer;
	}
    | math_sequence
      math_token
        {
            $$ = $1;
	    $$->add( $2);
	    delete $2;
	}
;

math_token:
      string_with_spcnl_token
    | '{'
      math_sequence
      '}'
        {
	    $$ = $2;
	}
    | SINGLESUBSCRIPT
        {
            $$ = new Buffer;
	    $$->add( "<SUB>", 5);
	    $$->add( $1);
	    $$->add( "</SUB>", 6);
	}
    | SINGLESUPERSCRIPT
        {
            $$ = new Buffer;
	    $$->add( "<SUP>", 5);
	    $$->add( $1);
	    $$->add( "</SUP>", 6);
	}
    | BEGINSUBSCRIPT  math_sequence  '}'
        {
            $$ = $2;
	    $$->prepend( "<SUB>", 5);
	    $$->add( "</SUB>", 6);
	}
    |  BEGINSUPERSCRIPT  math_sequence  '}'
        {
            $$ = $2;
	    $$->prepend( "<SUP>", 5);
	    $$->add( "</SUP>", 6);
	}
    | FRACTION  '{'  math_sequence  '}'  '{'  math_sequence  '}'
        {
	    $$ = $3;
	    $$->prepend( "<BOX>", 5);
	    $$->add( "<OVER>", 6);
	    $$->add( $6);
	    $$->add( "</BOX>", 6);
	    delete $6;
	}
;

/* End if Grammar */
/* ============== */
%%

int yyerror( char *s) {
    fprintf( stderr,
	     "error 1 in line %d in %s: in %s-code: %s.\n", 
	     line_number,
	     in_filename,
	     (CCMode ? "CC" : "TeX"),
	     s);
    return 0;
}

// Functions belonging to the Error messages
// -----------------------------------------
// See their implementations in parser.y

const char* errorMessage( ErrorNumber n) {
    switch ( n) {
    case NoError:
	return "NO ERROR";
    case ParseError:
	return "parse error";
    case VariableUsedError:
	return "The creationvariable was used but not defined";
    case ClassnameUsedError:
	return "The classname was used out of scope of any class";
    case TemplateParamExpectedError:
        return "A template parameter is missing";
    case MalformedTemplateParamError:
        return "The template parameter is malformed (<> nesting ..)";
    case MalformedFunctionDeclaration:
        return "The function declaration is malformed";
    case SemicolonMissingError:
        return "The declaration does not end in a semicolon";
    case IncludeNestingTooDeepError:
        return "Includes nested too deeply";
    case IncludeOpenError:
        return "Cannot open include file";
    case ChapterStructureError:
        return "Malformed chapter structure: one chapter per file";
    }
    return "UNKNOWN ERROR MESSAGE NUMBER";
}

void  printErrorMessage( ErrorNumber n){
    cerr << "error " << n << " in line " << line_number << " in `" 
         << in_filename << "': " << errorMessage( n) << "." << endl;
}


// support functions
// -----------------
Text* blockintroProcessing( const char* text, int len, Text* t) { 
    if ( len < 0) {  /* Hack! Here we know that t has to get capitalized.*/
        len = 4;
        InListFIter< TextToken> ix( *t);
        ForAll( ix) {
	    if ( ! (*ix).isSpace) {
	        char *s = (*ix).string;
		while ( *s) {
                    *s++ = toupper( *s);
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
