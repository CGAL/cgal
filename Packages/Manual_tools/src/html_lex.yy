/**************************************************************************
 
  html_lex.yy
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : lexical scanner for TeX and C++ code mixed files.
              Taylored for HTML manual generation.
  System    : flex, bison, C++ (g++)
  Author    : (c) 1996 Lutz Kettner
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

%{
extern "C" int yylex( void );

#include <stdlib.h>
#include <stdio.h>
#include <stream.h>
#include <string.h>
#include <database.h>
#include <html_config.h>
#include <html_syntax.tab.h>

// This flag is true if an \mbox encounters in a math environment.
bool mbox_within_math = false;

// This flag is true for the first macro in a newcommand. It should not
// get replaced by previous definitions.
bool actual_defining = false;

// This flag indicates whether we are in a tabbing or tabular environment.
bool tab_tag = false;

// This flag is true if we are inside a definition. here, \input and 
// \include should not open a file.
bool ignore_input_tag = false;

/* Set this flag to 1 to switch immediately to CCMode. */
int set_CCMode      = 0;
/* Set this flag to 1 to switch immediately to NestingMode. */
int set_NestingMode = 0;
/* Set this flag to 1 to switch back to INITIAL. */
int set_INITIAL     = 0;
int set_HTMLMODE    = 0;
int set_MMODE       = 0;

/* Store an old state like MMODE before treating macros like ccAnchor */
/* Some constants for this purpose. */
const int state_INITIAL     = 0;
const int state_MMODE       = 1;
const int state_NestingMode = 2;

int set_old_state = 0;
int old_state     = state_INITIAL;

/* Tag to mark whenever the unchecked keyword occurs. */
int   unchecked_tag = 0;

/* Count the linenumber for better errormessages. */
int   line_number  = 1;

/* store the creationvariable */
char* creationvariable = NULL;
char* formatted_creationvariable = NULL;
extern char* class_name;
extern char* template_class_name;
extern char* formatted_class_name;
extern char* formatted_template_class_name;

/* match math mode delimiters */
bool math_mode_toggle = false;

/* remember the necessary stop character for \verb"..." */
char stop_character;

/* prototypes */
bool is_html_multi_character( char c);
const char* html_multi_character( char c);
char* addSuffix( const char* name, const char* suffix);

/* Enable line number printing to cerr. */
extern char  line_switch;

/* Enable warnings for undefined macros. */
extern char  warn_switch;

/* Hack, to get rid of the yywrap. */
#define YY_SKIP_YYWRAP 1
#define yywrap() 1

void skipspaces( void);
void skipoptionalparam( void);

/* Include file handling */
FILE* current_in_stream;
extern const char* in_filename;
#define MAX_INCLUDE_DEPTH 32
YY_BUFFER_STATE buffer_stack[ MAX_INCLUDE_DEPTH];
const char*     file_name_stack[ MAX_INCLUDE_DEPTH];
FILE*           in_stream_stack[ MAX_INCLUDE_DEPTH];
int             line_number_stack[ MAX_INCLUDE_DEPTH];
int stack_ptr = 0;


#define SET( s) ((yylval.string.text = (s)), ( yylval.string.len = -1))
%}

/* The normal scanning mode parses TeX conventions.      */
/* In CCMode, it parses C++ conventions.                 */
/* The NestingMode parses only (){}[] nested expressions */
/* The VerbMode parses LaTeX \verb"..." statements as    */
/*     a sequence of characters                          */
/* ccStyleMode parses only ccStyle expressions.          */
%x CCMode
%x ccStyleMode
%x NestingMode
%x VerbMode
%x CPROGMode
%x ITEMMODE
%x INCLUDEMODE
%x TEXONLYMODE
%x HTMLMODE
%x HTMLGROUPMode
%x MMODE

letter          [a-zA-Z]
noletter        [^a-zA-Z]
digit           [0-9]
CCletter        [a-zA-Z_]
idfier          {letter}+
texmacro        [\\]{idfier}
CCidfier        ({CCletter}({CCletter}|{digit})*)
filename        [^ \t\n\\\{\}\[\]()]+
space           [\t ]
w               {space}*
ws              {space}+
escchar         [\\]
sign            [+-]
exp             [eE]
number          {digit}+
signNumber      ({sign}?{number})
floatNumber     ({signNumber}\.|{signNumber}\.{number})
expNumber       ({floatNumber}|{signNumber}){exp}{signNumber}
No              ({signNumber}|{floatNumber}|{expNumber})
operator        [^a-zA-Z_0-9\n\r\t \\]
measure         (({signNumber})|({floatNumber})){letter}{letter}
rmblockintro    ([\{][\\](rm))|([\\]((text)|(math))rm[\{])
ttblockintro    ([\{][\\](tt))|([\\]((text)|(math))tt[\{])
emblockintro    ([\{][\\](em))|([\\]emph[\{])
itblockintro    ([\{][\\]((it)|(sl)))|([\\]((text)|(math))((it)|(sl))[\{])
scblockintro    ([\{][\\](sc))|([\\]textsc[\{])
sfblockintro    ([\{][\\](sf))|([\\]((text)|(math))sf[\{])
bfblockintro    ([\{][\\]((bf)|(mathbold)))|([\\]((text)|(math))bf[\{])
calblockintro   ([\{][\\](cal))|([\\]mathcal[\{])

%%
 /* Mode switching can be triggered from the parser */
 /* ----------------------------------------------- */
	if (set_CCMode) {
	    BEGIN( CCMode);
	    set_CCMode = 0;
	} 
	if (set_NestingMode) {
	    BEGIN( NestingMode);
	    set_NestingMode = 0;
	} 
	if (set_INITIAL) {
	    BEGIN( INITIAL);
	    set_INITIAL = 0;
	}
	if (set_HTMLMODE) {
	    BEGIN( HTMLMODE);
	    set_HTMLMODE = 0;
	}
	if (set_MMODE) {
	    BEGIN( MMODE);
	    set_MMODE = 0;
	}
	if (set_old_state) {
	    switch ( old_state) {
	    case state_INITIAL:
		BEGIN( INITIAL);
		break;
	    case state_MMODE:
		BEGIN( MMODE);
		break;
	    case state_NestingMode:
		BEGIN( NestingMode);
		break;
	    }
	    set_old_state = 0;
	}

 /* Count line numbers in all modes for better error messages */
 /* --------------------------------------------------------- */
<INITIAL,CCMode,NestingMode,ccStyleMode,CPROGMode,MMODE,ITEMMODE,TEXONLYMODE,HTMLMODE,HTMLGROUPMode>[\n]	{
		    line_number++;
		    if ( line_switch)
		        cerr << "src-line " << line_number << endl;
		    return NEWLINE;
		}
<INITIAL,MMODE,NestingMode,ccStyleMode>[\\]"\n"      {
		    line_number++;
		    if ( line_switch)
		        cerr << "src-line " << line_number << endl;
	            yylval.string.text = " ";
		    yylval.string.len  = 1;
	  	    return SPACE;
}

 /* Handle include files      */
 /* ------------------------- */

[\\]((include)|(input))[\{]{w}   {  
	if (ignore_input_tag) {
	    int c = yyinput();
	    while( c && c != '}')
	        c = yyinput();
            yylval.string.text = " ";
            yylval.string.len  = 0;
            return SPACE;
        }
        BEGIN ( INCLUDEMODE); 
}
<INCLUDEMODE>{filename}          {
        /* remove trailing characters from the input/include statement */
        int c = yyinput();
	while( c && c != '}')
	    c = yyinput();
        if ( stack_ptr >= MAX_INCLUDE_DEPTH)
	    printErrorMessage( IncludeNestingTooDeepError);
	else {
	    buffer_stack[ stack_ptr]      = YY_CURRENT_BUFFER;
	    file_name_stack[ stack_ptr]   = in_filename;
	    in_stream_stack[ stack_ptr]   = current_in_stream;
	    line_number_stack[ stack_ptr] = line_number;
	    ++stack_ptr;
	    /* check whether any suffix is already there or not */
	    int i = 0;
	    while ( i < yyleng && yytext[i] != '.')
	        ++i;
	    char* tmp_name;
	    if ( i < yyleng)
	        /* this is done to get a 'delete'able copy */
	        tmp_name = addSuffix( yytext, "");
	    else
	        tmp_name = addSuffix( yytext, ".tex");
	    current_in_stream = fopen( tmp_name, "r");
	    if ( ! current_in_stream) {
	        int j;
	        printErrorMessage( IncludeOpenError);
		delete[] tmp_name;
		--stack_ptr;
		for ( j = stack_ptr; j >= 0; j--) {
		    cerr << "file include from `" << file_name_stack[ j]
			 << "'" << endl;
		}
	    } else {
	        yyin = current_in_stream;
	        in_filename = tmp_name;
	        yy_switch_to_buffer( yy_create_buffer( yyin, YY_BUF_SIZE));
		line_number = 1;
	    }
	}
	BEGIN( INITIAL);
}
<<EOF>> {
        if ( stack_ptr <= 0) {
	    yyterminate();
	} else {
	    -- stack_ptr;
	    delete[] (char*)in_filename;
	    line_number       = line_number_stack[ stack_ptr];
	    in_filename       = file_name_stack[ stack_ptr];
	    current_in_stream = in_stream_stack[ stack_ptr];
	    yy_switch_to_buffer( buffer_stack[ stack_ptr]);
	}
}


 /* Rules for TeX conventions */
 /* ------------------------- */
<INITIAL,MMODE>[\\]"%"  {   /* Avoid the quoted comment symbol */
		    yylval.character   = '%';
		    return CHAR;
		}
<INITIAL,MMODE>"%".*[\n]{w}  { /* Match one line TeX comments */
		    /* remove spaces in next line  */
		    unput( '\n');
		}
<INITIAL,MMODE>"%".*  { /* Match one line TeX comments */
                        /* at the last line in file */
		}
[\\]verb{noletter}   {   /* match LaTeX \verb"..." constructs */
		    BEGIN( VerbMode);
		    stop_character = yytext[ yyleng-1];
		    yylval.string.text = "<TT>";
		    yylval.string.len  = 4;
	  	    return STRING;
                }
<VerbMode>{ws}	{
	            yylval.string.text = yytext;
		    yylval.string.len  = yyleng;
	  	    return SPACE;
		}
<VerbMode>.	{
		    if ( yytext[0] == stop_character) {
		        BEGIN( INITIAL);
	                yylval.string.text = "</TT>";
		        yylval.string.len  = 5;
	  	        return STRING;
                    }
		    if ( yytext[0] == '\n') {
		        line_number++; 
			if ( line_switch)
			    cerr << "src-line " << line_number << endl;
                    }
	            yylval.character = yytext[0];
	  	    return CHAR;
		}

 /* Chapter and labels triggering new file and linking */
 /* -------------------------------------------------- */
[\\]chapter{w}[\{]  {
		    skipspaces();
		    return CHAPTER;
}
[\\]section[*]?{w}[\{]  {
		    return SECTION;
}
[\\]subsection[*]?{w}[\{]  {
		    return SUBSECTION;
}
[\\]subsubsection[*]?{w}[\{]  {
		    return SUBSUBSECTION;
}
[\\]label{w}[\{][^\}]+/[\}]  {
                    yyinput();
		    char* s = yytext;
		    while( *s != '{')
		        ++s;
		    ++s;
		    while( *s && *s <= ' ')
		        ++s;
	            yylval.string.text = s;
		    yylval.string.len  = -1;
		    return LABEL;
}


 /* Different keywords from the manual style triggering C++ formatting */
 /* ------------------------------------------------------------------ */
[\\]begin{w}[\{]ccClass[\}]{w}   {
		    BEGIN( CCMode);
		    current_font = it_font;
		    return BEGINCLASS;
		}
[\\]end{w}[\{]ccClass[\}]   {
		    return ENDCLASS;
		}
[\\]begin{w}[\{]ccClassTemplate[\}]{w}   {
		    BEGIN( CCMode);
		    current_font = it_font;
		    return BEGINCLASSTEMPLATE;
		}
[\\]end{w}[\{]ccClassTemplate[\}]   {
		    return ENDCLASSTEMPLATE;
		}
[\\]ccCreationVariable{w}[\{]{w}[^\}]*{w}[\}]   {
		    char *s = yytext + yyleng - 2;
		    while (( *s == ' ') || ( *s == '\t'))
		        s--;
		    char *r = yytext;
		    while ( *r != '{')
		        r++;
		    r++;
		    while (( *r == ' ') || ( *r == '\t'))
		        r++;
		    s[1]  = 0;
		    if ( creationvariable)
	                delete[] creationvariable;
		    if ( formatted_creationvariable)
		        delete[] formatted_creationvariable;
		    creationvariable = newstr( r);
		    formatted_creationvariable = new char[ strlen( 
		            creationvariable) + 8];
		    strcpy( formatted_creationvariable, "<I>");
		    strcat( formatted_creationvariable, creationvariable);
		    strcat( formatted_creationvariable, "</I>");
	            yylval.string.text = r;
		    yylval.string.len  = s - r + 1;
		    return CREATIONVARIABLE;
		}		    
[\\]ccConstructor/{noletter} { /* constructor declaration: change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return CONSTRUCTOR;
		}
[\\]ccMemberFunction/{noletter}  { /* method declaration: change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return METHOD;
		}
[\\]ccMethod/{noletter}   {   /* method declaration: change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return METHOD;
		}
[\\]ccFunction/{noletter} {   /* function declaration: change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return FUNCTION;
		}
[\\]ccFunctionTemplate/{noletter} {   /* function template declaration: 
			                 change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return FUNCTIONTEMPLATE;
		}
[\\]ccVariable/{noletter} {   /* variable declaration: change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return VARIABLE;
		}
[\\]ccTypedef/{noletter} {   /* typedef declaration: change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return TYPEDEF;
		}
[\\]ccNestedType/{noletter} {   /* nested type declaration: change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return NESTEDTYPE;
		}
[\\]ccEnum/{noletter} {   /* enum declaration: change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return ENUM;
		}
[\\]ccStruct/{noletter} {   /* struct declaration: change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return STRUCT;
		}
[\\]ccGlobalFunction/{noletter} {  /* function declaration: change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return GLOBALFUNCTION;
		}
[\\]ccGlobalFunctionTemplate/{noletter} {   /* function template declaration: 
			       change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return GLOBALFUNCTIONTEMPLATE;
		}
[\\]ccGlobalVariable/{noletter} {  /* variable declaration: change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return GLOBALVARIABLE;
		}
[\\]ccGlobalTypedef/{noletter} {   /* typedef declaration: change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return GLOBALTYPEDEF;
		}
[\\]ccGlobalEnum/{noletter} {   /* enum declaration: change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return GLOBALENUM;
		}
[\\]ccGlobalStruct/{noletter} {   /* struct declaration: change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return GLOBALSTRUCT;
		}
[\\]ccDeclaration/{noletter} {   /* general declaration: change to CCMode */
		    skipspaces();
		    current_font = it_font;
		    BEGIN( CCMode);
		    return DECLARATION;
		}
[\\]ccHidden/{noletter}   {
		    skipspaces();
	  	    return HIDDEN;
		}
[\\]ccUnchecked/{noletter}  { 
	            /* trigger a global boolean and treat it like a space */
		    skipspaces();
		    unchecked_tag = 1;
		    yylval.string.text = "          ";
		    yylval.string.len  = 10;
	  	    return SPACE;
		}

[\\]"begin{verbatim}" |
[\\]"begin{cprog}"    {
		    BEGIN( CPROGMode);
		    return CPROGBEGIN;
                }
<CPROGMode>[\\]"end{verbatim}"  |
<CPROGMode>[\\]"end{cprog}"     {
		    BEGIN( INITIAL);
		    return CPROGEND;
                }
[\\]"cprogfile{"[^\}]*"}"    {
	            yylval.string.text = yytext + 11;
		    yylval.string.len  = yyleng - 12;
	  	    return CPROGFILE;
                 }

[\\]"begin{ccHtmlOnly}" {
                    old_state = state_INITIAL;
		    BEGIN( HTMLGROUPMode);
		    return HTMLBEGIN;
                }
<MMODE>[\\]"begin{ccHtmlOnly}" {
                    old_state = state_MMODE;
		    BEGIN( HTMLGROUPMode);
		    return HTMLBEGIN;
                }
<NestingMode>[\\]"begin{ccHtmlOnly}" {
                    old_state = state_NestingMode;
		    BEGIN( HTMLGROUPMode);
		    return HTMLBEGIN;
                }
<HTMLGROUPMode>[\\]"end{ccHtmlOnly}"  {
		    set_old_state = 1;
		    return HTMLEND;
                }

[\\]"begin{ccTexOnly}" {
                    old_state = state_INITIAL;
		    BEGIN( TEXONLYMODE);
		    return TEXONLYBEGIN;
                }
<MMODE>[\\]"begin{ccTexOnly}" {
                    old_state = state_MMODE;
		    BEGIN( TEXONLYMODE);
		    return TEXONLYBEGIN;
                }
<NestingMode>[\\]"begin{ccTexOnly}" {
                    old_state = state_NestingMode;
		    BEGIN( TEXONLYMODE);
		    return TEXONLYBEGIN;
                }
<TEXONLYMODE>[\\]"end{ccTexOnly}"  {
		    set_old_state = 1;
		    return TEXONLYEND;
                }
<TEXONLYMODE>.    {
		    yylval.character   = yytext[0];
		    return CHAR;
		}
[\\]ccTexHtml{w}[\{]  {
                    old_state = state_INITIAL;
		    return LATEXHTML;
                }
<MMODE>[\\]ccTexHtml{w}[\{]  {
                    old_state = state_MMODE;
		    return LATEXHTML;
                }
<NestingMode>[\\]ccTexHtml{w}[\{]  {
                    old_state = state_NestingMode;
		    return LATEXHTML;
                }
[\\]ccAnchor{w}[\{]  {
		    /* The first parameter is the URL, the second is the */
                    /* message that will be highlighted */
                    old_state = state_INITIAL;
		    BEGIN( HTMLMODE);
		    return ANCHOR;
                }
<MMODE>[\\]ccAnchor{w}[\{]  {
                    old_state = state_MMODE;
		    BEGIN( HTMLMODE);
		    return ANCHOR;
                }
<NestingMode>[\\]ccAnchor{w}[\{]  {
                    old_state = state_NestingMode;
		    BEGIN( HTMLMODE);
		    return ANCHOR;
                }
<HTMLMODE>[\}]  {
		    set_old_state = 1;
		    return '}';
                }
<HTMLMODE,HTMLGROUPMode>.     {
		    yylval.character   = yytext[0];
		    return CHAR;
		}
<INITIAL,MMODE,NestingMode>[\\]path[|][^|]*[|] {
	            yylval.string.text = yytext + 6;
		    yylval.string.len  = -1;
                    return HTMLPATH;
}

 /* Flexibility for HTML class files. */
 /* -------------------------------------------------------------- */
<INITIAL,MMODE,NestingMode>[\\]ccHtmlNoClassLinks/{noletter}   {
		    skipspaces();
		    html_no_class_links = true;
}
<INITIAL,MMODE,NestingMode>[\\]ccHtmlNoClassFile/{noletter}    {
		    skipspaces();
		    html_no_class_file = true;
}
<INITIAL,MMODE,NestingMode>[\\]ccHtmlNoClassIndex/{noletter}   {
		    skipspaces();
		    html_no_class_index = true;
}
[\\]begin{w}[\{]ccHtmlClassFile[\}]{w}   {
		    BEGIN( CCMode);
		    return HTMLBEGINCLASSFILE;
}
[\\]end{w}[\{]ccHtmlClassFile[\}]   {
		    skipspaces();
		    return HTMLENDCLASSFILE;
}
[\\]ccHtmlIndex/{noletter}                                     {
		    skipspaces();
		    yylval.string.text = sort_key_class;
		    return HTMLINDEX;
}
[\\]ccHtmlIndex[\[][^\]][\]]/{noletter}                        {
		    skipspaces();
		    yylval.string.text = find_sort_key( yytext + 13);
		    return HTMLINDEX;
}
[\\]ccHtmlIndexC/{noletter}                                    {
		    skipspaces();
		    yylval.string.text = sort_key_class;
		    BEGIN( CCMode);
		    return HTMLINDEXC;
}
[\\]ccHtmlIndexC[\[][^\]][\]]/{noletter}                       {
		    skipspaces();
		    yylval.string.text = find_sort_key( yytext + 14);
		    BEGIN( CCMode);
		    return HTMLINDEXC;
}
[\\]ccHtmlCrossLink/{noletter}                                    {
		    skipspaces();
		    BEGIN( CCMode);
		    return HTMLCROSSLINK;
}

 /* Specialized keywords from the manual style */
 /* -------------------------------------------------------------- */
[\\]cc((Style)|(c))/{noletter}  {
                    /* CCstyle formatting: change to ccStyleMode */
		    skipspaces();
		    BEGIN( ccStyleMode);
		    current_font = it_font;
		    return CCSTYLE;
		}
<INITIAL,MMODE,NestingMode>[\\]ccVar/{noletter}      {
		    skipspaces();
		    if ( creationvariable) {
	                yylval.string.text = formatted_creationvariable;
		        yylval.string.len  = strlen( yylval.string.text);
		    } else {
		        printErrorMessage( VariableUsedError);
	                yylval.string.text = "Unknown creationvariable";
		        yylval.string.len  = strlen( yylval.string.text);
		    }
		    return STRING;
                 }
<INITIAL,MMODE,NestingMode>[\\]ccPureVar/{noletter}      {
		    skipspaces();
		    if ( creationvariable) {
	                yylval.string.text = creationvariable;
		        yylval.string.len  = strlen( yylval.string.text);
		    } else {
		        printErrorMessage( VariableUsedError);
	                yylval.string.text = "Unknown creationvariable";
		        yylval.string.len  = strlen( yylval.string.text);
		    }
		    return STRING;
                 }
<INITIAL,MMODE,NestingMode>[\\]ccClassName/{noletter} {
		    skipspaces();
		    if ( formatted_class_name) {
	                yylval.string.text = formatted_class_name;
		        yylval.string.len  = strlen( yylval.string.text);
		    } else {
		        printErrorMessage( ClassnameUsedError);
	                yylval.string.text = "Unknown classname";
		        yylval.string.len  = strlen( yylval.string.text);
		    }
		    return STRING;
                 }
<INITIAL,MMODE,NestingMode>[\\]ccPureClassName/{noletter} {
		    skipspaces();
		    if ( class_name) {
	                yylval.string.text = class_name;
		        yylval.string.len  = strlen( yylval.string.text);
		    } else {
		        printErrorMessage( ClassnameUsedError);
	                yylval.string.text = "Unknown classname";
		        yylval.string.len  = strlen( yylval.string.text);
		    }
		    return STRING;
                 }
<INITIAL,MMODE,NestingMode>[\\]ccClassTemplateName/{noletter} {
		    skipspaces();
		    if ( formatted_template_class_name) {
	                yylval.string.text = formatted_template_class_name;
		        yylval.string.len  = strlen( yylval.string.text);
		    } else {
		        printErrorMessage( ClassnameUsedError);
	                yylval.string.text = "Unknown template_classname";
		        yylval.string.len  = strlen( yylval.string.text);
		    }
		    return STRING;
                 }
<INITIAL,MMODE,NestingMode>[\\]ccPureClassTemplateName/{noletter} {
		    skipspaces();
		    if ( template_class_name) {
	                yylval.string.text = template_class_name;
		        yylval.string.len  = strlen( yylval.string.text);
		    } else {
		        printErrorMessage( ClassnameUsedError);
	                yylval.string.text = "Unknown template_classname";
		        yylval.string.len  = strlen( yylval.string.text);
		    }
		    return STRING;
                 }
[\\]ccSection/{noletter} {  
		    skipspaces();
		    return CCSECTION; 
                 }
[\\]ccSubsection/{noletter} {  
		    skipspaces();
		    return CCSUBSECTION; 
                 }
[\\]RCSdef\{{texmacro}    {
	            yylval.string.text = newstr(yytext + 8);
		    yylval.string.len  = -1;
		    return RCSDEF;
                 }
[\\]RCSdefDate\{{texmacro}    {
	            yylval.string.text = newstr(yytext + 12);
		    yylval.string.len  = -1;
		    return RCSDEFDATE;
                 }
[\\]CC/{noletter}        {
		    skipspaces();
	            yylval.string.text = "C++";
		    yylval.string.len  = 3;
		    return STRING;
                 }
[\\]gcc/{noletter}        {
		    skipspaces();
	            yylval.string.text = "g++";
		    yylval.string.len  = 3;
		    return STRING;
                 }
<INITIAL,MMODE>[\\]nat/{noletter}       {
		    skipspaces();
	            yylval.string.text = "<B>N</B>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
<INITIAL,MMODE>[\\]real/{noletter}      {
		    skipspaces();
	            yylval.string.text = "<B>R</B>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
<INITIAL,MMODE>[\\]N/{noletter}         {
		    skipspaces();
	            yylval.string.text = "<B>N</B>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
<INITIAL,MMODE>[\\]R/{noletter}         {
		    skipspaces();
	            yylval.string.text = "<B>R</B>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
<INITIAL,MMODE>[\\]Q/{noletter}         {
		    skipspaces();
	            yylval.string.text = "<B>Q</B>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
<INITIAL,MMODE>[\\]Z/{noletter}         {
		    skipspaces();
	            yylval.string.text = "<B>Z</B>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
<INITIAL,MMODE>[\\]E/{noletter}         {
		    skipspaces();
	            yylval.string.text = "<B>E</B>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]leda/{noletter}        {
		    skipspaces();
	            /* yylval.string.text = "<TT>LEDA</TT>"; */
	            yylval.string.text = "LEDA";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]cgal/{noletter}      {
		    skipspaces();
		    /* yylval.string.text = "<TT>CGAL</TT>"; */
	            yylval.string.text = "CGAL";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]protocgal/{noletter} {
		    skipspaces();
	            /* yylval.string.text = "<TT>C++GAL</TT>"; */
	            yylval.string.text = "C++GAL";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]plageo/{noletter}    {
		    skipspaces();
	            /* yylval.string.text = "<TT>PLAGEO</TT>"; */
	            yylval.string.text = "PLAGEO";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccDefinition/{noletter} {
		    skipspaces();
	            yylval.string.text = "<H3>Definition</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccInheritsFrom/{noletter} {
		    skipspaces();
	            yylval.string.text = "<H3>Inherits From</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccParameters/{noletter} {
		    skipspaces();
	            yylval.string.text = "<H3>Parameters</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccTypes/{noletter} {
		    skipspaces();
	            yylval.string.text = "<H3>Types</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccConstants/{noletter} {
		    skipspaces();
	            yylval.string.text = "<H3>Constants</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccCreation/{noletter}  {
		    skipspaces();
	            yylval.string.text = "<H3>Creation</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccOperations/{noletter} {
		    skipspaces();
	            yylval.string.text = "<H3>Operations</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccAccessFunctions/{noletter} {
		    skipspaces();
	            yylval.string.text = "<H3>Access Functions</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccPredicates/{noletter} {
		    skipspaces();
	            yylval.string.text = "<H3>Predicates</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccModfiers/{noletter} {
		    skipspaces();
	            yylval.string.text = "<H3>Modifiers</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccImplementation/{noletter} {
		    skipspaces();
	            yylval.string.text = "<H3>Implementation</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccExample/{noletter}   {
		    skipspaces();
	            yylval.string.text = "<H3>Example</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]begin[\{]ccAdvanced[\}]   {
		    skipspaces();
	            yylval.string.text = "<BR><IMG BORDER=0 SRC=\""
                        "cc_advanced_begin.gif\" ALT=\"begin of advanced "
                        "section\"><BR>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]end[\{]ccAdvanced[\}]   {
		    skipspaces();
	            yylval.string.text = "<BR><IMG BORDER=0 SRC=\""
                        "cc_advanced_end.gif\" ALT=\"end of advanced "
			"section\"><BR>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccInclude/{noletter}   {
 		    current_font = it_font;
                    return INCLUDE;
                 }
[\\]ccHeading/{noletter}   {
                    return HEADING;
                 }
[\\]ccPrecond/{noletter}      {
		    skipspaces();
	            yylval.string.text = "<BR><STRONG>Precondition: </STRONG>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccPostcond/{noletter}      {
		    skipspaces();
	            yylval.string.text ="<BR><STRONG>Postcondition: </STRONG>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccCommentHeading/{noletter}   {
                    return COMMENTHEADING;
                 }
[\\]ccSetTwoOfThreeColumns/{noletter} {
		    skipspaces();
		    return GOBBLETWOPARAMS;
                 }
[\\]cc((SetThreeColumns)|(Three))/{noletter} {
		    skipspaces();
		    return GOBBLETHREEPARAMS;
                 }
[\\]ccSetOneOfTwoColumns/{noletter} {
		    skipspaces();
		    return GOBBLEONEPARAM;
                 }
[\\]cc((SetTwoColumns)|(Two))/{noletter} {
		    skipspaces();
		    return GOBBLETWOPARAMS;
                 }
[\\]cc((PropagateThreeToTwoColumns)|(ThreeToTwo))/{noletter}  {
                    yylval.string.text = " ";
                    yylval.string.len  = 0;
                    return SPACE;
}
[\\]g?def{w}[\\]{idfier}[^\{]*    {
                    ignore_input_tag = true;
                    BEGIN( NestingMode);
	            yylval.string.text = yytext;
		    yylval.string.len  = yyleng;
                    return TEXDEF;
}
[\\](re)?newcommand{w}/{noletter}    {
                    ignore_input_tag = true;
                    actual_defining = true;
                    BEGIN( NestingMode);
                    return NEWCOMMAND;
}
[\\]ccTagDefaults{w}/{noletter}  {
                    tag_defaults();
                    yylval.string.text = " ";
                    yylval.string.len  = 0;
                    return SPACE;
}
[\\]ccTagFullDeclarations{w}/{noletter}  {
                    tag_full_declarations();
                    yylval.string.text = " ";
                    yylval.string.len  = 0;
                    return SPACE;
}
[\\]ccChapterAuthor{w}/{noletter}   {
                    return CHAPTERAUTHOR;
}
[\\]ccChapterSubTitle{w}/{noletter}   {
                    return CHAPTERSUBTITLE;
}
[\\]cc((Glue((Begin)|(End)|(Declarations))?)|(ParDims))/{noletter}  {
                    yylval.string.text = " ";
                    yylval.string.len  = 0;
                    return SPACE;
}



 /* Try tabular and tabbing environments                           */
 /* -------------------------------------------------------------- */

<INITIAL,MMODE,NestingMode>[\\]begin[\{]((tabbing)|(eqnarray"*"?))[\}]     |
<INITIAL,MMODE,NestingMode>[\\]begin[\{]((tabular)|(array))[\}][\{][^\}]*[\}] {
	            tab_tag = true;
	            yylval.string.text = 
                        "<TABLE><TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
<INITIAL,MMODE,NestingMode>[\\]end[\{]((tabbing)|(tabular)|(eqnarray"*"?)|(array))[\}]   {
	            tab_tag = false;
	            yylval.string.text = "</TD></TR></TABLE>";
		    yylval.string.len  = -1;
		    return STRING;		    
                 }
<INITIAL,MMODE,NestingMode>[\\][=>]         {
	            yylval.string.text = 
                        "</TD><TD ALIGN=LEFT VALIGN=TOP NOWRAP>";
		    yylval.string.len  = -1;
		    return STRING;		    
                 }
	

 /* keywords from TeX/LaTeX that have an easy HTML counterpart     */
 /* -------------------------------------------------------------- */

[\\]begin[\{]itemize[\}]     {
	            yylval.string.text = "<UL>";
		    yylval.string.len  = 4;
		    return STRING;		    
                 }
[\\]end[\{]itemize[\}]   {
	            yylval.string.text = "</UL>";
		    yylval.string.len  = 5;
		    return STRING;		    
                 }
[\\]item[\[]     {
		    BEGIN( ITEMMODE);
	            yylval.string.text = "<DT><B>";
		    yylval.string.len  = 7;
		    return STRING;		    
                 }
<ITEMMODE>[\]]   {
		    BEGIN( INITIAL);
	            yylval.string.text = "</B><DD>";
		    yylval.string.len  = 8;
		    return STRING;		    
                 }
<ITEMMODE>[^\]\n]+ {
	            yylval.string.text = yytext;
		    yylval.string.len  = yyleng;
	  	    return STRING;
                 }
[\\]item/{noletter}     {
		    skipspaces();
	            yylval.string.text = "<LI>";
		    yylval.string.len  = 4;
		    return STRING;		    
                 }
[\\]begin[\{]enumerate[\}]     {
	            yylval.string.text = "<OL>";
		    yylval.string.len  = 4;
		    return STRING;		    
                 }
[\\]end[\{]enumerate[\}]   {
	            yylval.string.text = "</OL>";
		    yylval.string.len  = 5;
		    return STRING;		    
                 }
[\\]begin[\{]description[\}]     {
	            yylval.string.text = "<DL>";
		    yylval.string.len  = 4;
		    return STRING;		    
                 }
[\\]end[\{]description[\}]   {
	            yylval.string.text = "</DL>";
		    yylval.string.len  = 5;
		    return STRING;		    
                 }
[\\]begin[\{]quote[\}]     {
	            yylval.string.text = "<BLOCKQUOTE>";
		    yylval.string.len  = -1;
		    return STRING;		    
                 }
[\\]end[\{]quote[\}]    {
	            yylval.string.text = "</BLOCKQUOTE>";
		    yylval.string.len  = -1;
		    return STRING;		    
                 }
[\\]((pagebreak)|(newpage)|(clearpage)|(cleardoublepage))/{noletter}     {
		    skipoptionalparam();
	            yylval.string.text = "<P>";
		    yylval.string.len  = 3;
		    return STRING;		    
                 }
<INITIAL,MMODE,NestingMode,ccStyleMode>[\\](l?)dots/{noletter}     {
		    skipspaces();
	            yylval.string.text = "...";
		    yylval.string.len  = 3;
		    return STRING;
                 }
<MMODE>[\\]le[q]?/{noletter}          {
		    skipspaces();
	            yylval.string.text = "&lt;=";  // &le;  not yet supported
		    yylval.string.len  = 5;
		    return STRING;
                 }
<MMODE>[\\]ge[q]?/{noletter}          {
		    skipspaces();
	            yylval.string.text = "&gt;=";  // &ge;  not yet supported
		    yylval.string.len  = 5;
		    return STRING;
                 }
<MMODE>[\\]neq/{noletter}          {
		    skipspaces();
	            yylval.string.text = "&ne;";
		    yylval.string.len  = 4;
		    return STRING;
                 }
[\\]((big)|(med))skip/{noletter}          {
		    return NEWLINE;
                 }
<INITIAL,MMODE,NestingMode,ccStyleMode>[\\]"&"          {
		     yylval.string.text = "&amp;";
		     yylval.string.len  = 5;
		     return STRING;
		 }
<INITIAL,MMODE,NestingMode,ccStyleMode>[\\][_^#$~%]  {
		    yylval.character = yytext[1];
		    return CHAR;
                 }
<INITIAL,MMODE,NestingMode>[~]              |
<INITIAL,MMODE,NestingMode,ccStyleMode>[\\]{space}      {
	            yylval.string.text = " ";
		    yylval.string.len  = 1;
	  	    return SPACE;
		 }
<INITIAL,MMODE,NestingMode>[\\][\\]         {
		    skipoptionalparam();
		    if ( tab_tag) {
	                yylval.string.text = 
                            "</TD></TR><TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>";
		        yylval.string.len  = -1;
		        return STRING;
                    }
	            yylval.string.text = "<BR>";
		    yylval.string.len  = 4;
	  	    return STRING;
		 }

 /* Mathmode                                                       */
 /* -------------------------------------------------------------- */
[\\][(\[]        |
[$]              {
		    BEGIN (MMODE);
		    return BEGINMATH;
                 }
<MMODE>[\\][)\]] |
<MMODE>[$]       {
		    BEGIN (INITIAL);
		    return ENDMATH;
                 }
<MMODE>"_"[^\{]  {
		    yylval.character = yytext[1];
		    return SINGLESUBSCRIPT;
                 }
<MMODE>"^"[^\{]  {
		    yylval.character = yytext[1];
		    return SINGLESUPERSCRIPT;
                 }
<MMODE>"_"[\{]   {
		    return BEGINSUBSCRIPT;
                 }
<MMODE>"^"[\{]   {
		    return BEGINSUPERSCRIPT;
                 }
<MMODE>\frac/{noletter}   {
		    return FRACTION;
                 }
<MMODE>[\\][\{\}]   {
		    yylval.character = yytext[1];
		    return CHAR;
                 }
<MMODE>[\{\}]    {
		    return yytext[0];
                 }

  /* yet not supported characters ...
  <MMODE>[\\]delta/{noletter}       { SET( "&delta;");    return STRING;}
  <MMODE>[\\]epsilon/{noletter}     { SET( "&epsi;");     return STRING;}
  <MMODE>[\\]varepsilon/{noletter}  { SET( "&epsi;");     return STRING;}
  <MMODE>[\\]lambda/{noletter}      { SET( "&lambda;");   return STRING;}
  <MMODE>[\\]pi/{noletter}          { SET( "&pi;");       return STRING;}
  <MMODE>[\\]varpi/{noletter}       { SET( "&piv;");      return STRING;}
  ... */

<INITIAL,NestingMode>[\\]["]a     { SET( "&auml;");     return STRING;}
<INITIAL,NestingMode>[\\]["]o     { SET( "&ouml;");     return STRING;}
<INITIAL,NestingMode>[\\]["]u     { SET( "&uuml;");     return STRING;}
<INITIAL,NestingMode>[\\]["]A     { SET( "&Auml;");     return STRING;}
<INITIAL,NestingMode>[\\]["]O     { SET( "&Ouml;");     return STRING;}
<INITIAL,NestingMode>[\\]["]U     { SET( "&Uuml;");     return STRING;}
<INITIAL,NestingMode>[\\][']a     { SET( "&aacute;");   return STRING;}
<INITIAL,NestingMode>[\\][']e     { SET( "&eacute;");   return STRING;}
<INITIAL,NestingMode>[\\][`]a     { SET( "&agrave;");   return STRING;}
<INITIAL,NestingMode>[\\][`]e     { SET( "&egrave;");   return STRING;}
<INITIAL,NestingMode>[\\]"^"a     { SET( "&acirc;");    return STRING;}
<INITIAL,NestingMode>[\\]"^"e     { SET( "&ecirc;");    return STRING;}
<INITIAL,NestingMode>[\\]ss[\{][\}]  { SET( "&szlig;");    return STRING;}
<MMODE>[\\]times/{noletter}       { SET( "&times;");    return STRING;}
<MMODE>[\\]in/{noletter}          { SET( " is in ");    return STRING;}
<MMODE>[\\]alpha/{noletter}       { SET( "&alpha;");    return STRING;}
<MMODE>[\\]beta/{noletter}        { SET( "&beta;");     return STRING;}
<MMODE>[\\]gamma/{noletter}       { SET( "&gamma;");    return STRING;}
<MMODE>[\\]delta/{noletter}       { SET( "delta");    return STRING;}
<MMODE>[\\]epsilon/{noletter}     { SET( "eps");     return STRING;}
<MMODE>[\\]varepsilon/{noletter}  { SET( "eps");     return STRING;}
<MMODE>[\\]zeta/{noletter}        { SET( "&zeta;");     return STRING;}
<MMODE>[\\]eta/{noletter}         { SET( "&eta;");      return STRING;}
<MMODE>[\\]theta/{noletter}       { SET( "&theta;");    return STRING;}
<MMODE>[\\]vartheta/{noletter}    { SET( "&thetav;");   return STRING;}
<MMODE>[\\]iota/{noletter}        { SET( "&iota;");     return STRING;}
<MMODE>[\\]kappa/{noletter}       { SET( "&kappa;");    return STRING;}
<MMODE>[\\]lambda/{noletter}      { SET( "lambda");   return STRING;}
<MMODE>[\\]mu/{noletter}          { SET( "&mu;");       return STRING;}
<MMODE>[\\]nu/{noletter}          { SET( "&nu;");       return STRING;}
<MMODE>[\\]xi/{noletter}          { SET( "&xi;");       return STRING;}
<MMODE>[\\]pi/{noletter}          { SET( "pi");       return STRING;}
<MMODE>[\\]varpi/{noletter}       { SET( "pi");      return STRING;}
<MMODE>[\\]rho/{noletter}         { SET( "&rho;");      return STRING;}
<MMODE>[\\]varrho/{noletter}      { SET( "&rho;");      return STRING;}
<MMODE>[\\]sigma/{noletter}       { SET( "&sigma;");    return STRING;}
<MMODE>[\\]varsigma/{noletter}    { SET( "&sigmav;");   return STRING;}
<MMODE>[\\]tau/{noletter}         { SET( "&tau;");      return STRING;}
<MMODE>[\\]upsilon/{noletter}     { SET( "&upsi;");     return STRING;}
<MMODE>[\\]phi/{noletter}         { SET( "&phi;");      return STRING;}
<MMODE>[\\]varphi/{noletter}      { SET( "&phiv;");     return STRING;}
<MMODE>[\\]chi/{noletter}         { SET( "&chi;");      return STRING;}
<MMODE>[\\]psi/{noletter}         { SET( "&psi;");      return STRING;}
<MMODE>[\\]omega/{noletter}       { SET( "&omega;");    return STRING;}
<MMODE>[\\]Gamma/{noletter}       { SET( "&Gamma;");    return STRING;}
<MMODE>[\\]Delta/{noletter}       { SET( "&Delta;");    return STRING;}
<MMODE>[\\]Theta/{noletter}       { SET( "&Theta;");    return STRING;}
<MMODE>[\\]Lambda/{noletter}      { SET( "&Lambda;");   return STRING;}
<MMODE>[\\]Xi/{noletter}          { SET( "&Xi;");       return STRING;}
<MMODE>[\\]Pi/{noletter}          { SET( "&Pi;");       return STRING;}
<MMODE>[\\]Sigma/{noletter}       { SET( "&Sigma;");    return STRING;}
<MMODE>[\\]Upsilon/{noletter}     { SET( "&Upsi;");     return STRING;}
<MMODE>[\\]Phi/{noletter}         { SET( "&Phi;");      return STRING;}
<MMODE>[\\]Psi/{noletter}         { SET( "&Psi;");      return STRING;}
<MMODE>[\\]Omega/{noletter}       { SET( "&Omega;");    return STRING;}

  /* math symbols */
  /* ------------ */

<MMODE>[\\]((arc)?)|((tan)|(sin)|(cos))/{noletter}  { 
	            yylval.string.text = yytext+1;
		    yylval.string.len  = yyleng-1;
	  	    return STRING;
}
<MMODE>[\\]((arg)|(cosh)|(sinh)|(cot)|(coth)|(csc)|(deg)|(det))/{noletter} { 
	            yylval.string.text = yytext+1;
		    yylval.string.len  = yyleng-1;
	  	    return STRING;
}
<MMODE>[\\]((dim)|(exp)|(gcd)|(hom)|(inf)|(ker)|(lg)|(lim))/{noletter}  { 
	            yylval.string.text = yytext+1;
		    yylval.string.len  = yyleng-1;
	  	    return STRING;
}
<MMODE>[\\]((liminf)|(limsup)|(ln)|(log)|(max)|(min)|(Pr)|(sec))/{noletter}  { 
	            yylval.string.text = yytext+1;
		    yylval.string.len  = yyleng-1;
	  	    return STRING;
}
<MMODE>[\\]((sinh)|(sup)|(tanh)|(bmod)|(pmod))/{noletter}  { 
	            yylval.string.text = yytext+1;
		    yylval.string.len  = yyleng-1;
	  	    return STRING;
}

 /* keywords from TeX/LaTeX that should vanish in HTML             */
 /* -------------------------------------------------------------- */
<INITIAL,MMODE,NestingMode>[\\]((smallskip)|(protect)|(sloppy))/{noletter}   {}
[\\]((maketitle)|(tableofcontents))/{noletter}            {}
[\\]((begin)|(end))[\{]document[\}]                       {}
<INITIAL,MMODE,NestingMode>[\\]((tiny)|(scriptsize)|(footnotesize)|(small)|(normalsize)|(large)|(Large))/{noletter}                      {}


[\\]newsavebox{w}[\{]          |
[\\]usebox{w}[\{]              |
[\\]hspace[*]?{w}[\{]          |
[\\]vspace[*]?{w}[\{]          {
	            /* CCstyle formatting: change to NestingMode */
		    BEGIN( NestingMode);
		    return IGNOREBLOCK;
		}
[\\]((documentclass)|(documentstyle)|(usepackage)|(pagestyle)|(pagenumbering)|(bibliographystyle)|(bibliography)|(title)|(author)|(date)){w}/{noletter}  {
		    skipoptionalparam();
		    yyinput();  /* gobble opening brace */
		    BEGIN( NestingMode);
		    return IGNOREBLOCK;
		}
[\\]((textwidth)|(textheight)|(topmargin)|(evensidemargin)|(oddsidemargin)|(headsep)|(parindent)|(parskip)|(beforecprogskip)|(aftercprogskip)){w}(({texmacro})|({measure}))  {}

[\\]((savebox)|(setlength)|(setcounter)){w}[\{]     {
	            /* CCstyle formatting: change to NestingMode */
		    BEGIN( NestingMode);
		    return IGNORETWOBLOCKS;
		}
[\\]begin[\{]minipage[\}]      |
[\\]end[\{]minipage[\}]	       {
		    skipoptionalparam();
	            yylval.string.text = " ";
		    yylval.string.len  = 1;
		    return SPACE;		    
                 }

[\\](([vf]box)|(parbox)|([hv]fill)|(nopagebreak)|(nolinebreak)|(linebreak)|(samepage))/{noletter}        {
		    skipoptionalparam();
	            yylval.string.text = " ";
		    yylval.string.len  = 1;
		    return SPACE;		    
                 }
<INITIAL,NestingMode>[\\][mh]box/{noletter}        {
		    return MBOX;		    
                 }
<MMODE>[\\][mh]box/{noletter}        {
                    mbox_within_math = true;
                    BEGIN( INITIAL);
		    return MBOX;		    
                 }
<INITIAL,MMODE,NestingMode,ccStyleMode>[\\][,;:!]  {}



 /* Footnotemanagement           */
 /* ---------------------------- */
<INITIAL,MMODE,NestingMode>[\\]footnotemark/{noletter}  {
                    return FOOTNOTEMARK;
                 }
<INITIAL,MMODE,NestingMode>[\\]footnotetext/{noletter}  {
                    return FOOTNOTETEXT;
                 }
<INITIAL,MMODE,NestingMode>[\\]footnote/{noletter}  {
                    return FOOTNOTE;
                 }


 /* Support for the Bibliography */
 /* ---------------------------- */
[\\]begin[\{]thebibliography[\}][\{][^\}]*[\}]    {
		    return BEGINBIBLIO;		    
                 }
[\\]end[\{]thebibliography[\}]   {
		    return ENDBIBLIO;		    
                 }
[\\]newblock/{noletter}    {}
[\\]cite{w}/{noletter}     {
		    BEGIN( NestingMode);
		    return CITE;		    
                 }
[\\]bibcite{w}/{noletter}  {
		    BEGIN( NestingMode);
		    return BIBCITE;
                }
[\\]bibitem{w}/{noletter}  {
		    BEGIN( NestingMode);
		    return BIBITEM;
                }


 /* Grouping symbols */
 /* ---------------- */
<INITIAL,MMODE,NestingMode,ccStyleMode>[\\][\{]        {
	            yylval.character = '{';
	  	    return CHAR;
		}
<INITIAL,MMODE,NestingMode,ccStyleMode>[\\][\}]        {
	            yylval.character = '}';
	  	    return CHAR;
		}
<INITIAL,MMODE,NestingMode>[\\]left.       {
	            yylval.character = yytext[5];
	  	    return CHAR;
		}
<INITIAL,MMODE,NestingMode>[\\]right.      {
	            yylval.character = yytext[6];
	  	    return CHAR;
		}
<INITIAL,CCMode,NestingMode,ccStyleMode>[\{]    {
		    return '{';
		}
<INITIAL,CCMode,NestingMode,ccStyleMode>[\}]    {
	  	    return '}';
		}

<INITIAL,MMODE,NestingMode>{ttblockintro}  {  /* TeX styles like {\tt ... */
		    skipspaces();
		    return TTBLOCKINTRO;
		}
<INITIAL,MMODE,NestingMode>{emblockintro}  { skipspaces();return EMBLOCKINTRO;}
<INITIAL,MMODE,NestingMode>{itblockintro}  { skipspaces();return ITBLOCKINTRO;}
<INITIAL,MMODE,NestingMode>{scblockintro}  { skipspaces();return SCBLOCKINTRO;}
<INITIAL,MMODE,NestingMode>{bfblockintro}  { skipspaces();return BFBLOCKINTRO;}
<INITIAL,MMODE,NestingMode>{rmblockintro}  { skipspaces();return RMBLOCKINTRO;}
<INITIAL,MMODE,NestingMode>{sfblockintro}  { skipspaces();return SFBLOCKINTRO;}
<INITIAL,MMODE,NestingMode>{calblockintro} {skipspaces();return CALBLOCKINTRO;}

<CCMode,ccStyleMode>[\\]tt/{noletter}      {
		        skipspaces();
		        yylval.string.text = "\\T\\";
		        yylval.string.len  = -1;
			return STRING;
                }
<CCMode,ccStyleMode>[\\]bf/{noletter}      {
		        skipspaces();
		        yylval.string.text = "\\B\\";
		        yylval.string.len  = -1;
			return STRING;
                }
<CCMode,ccStyleMode>[\\]ccFont/{noletter}  {
		        skipspaces();
		        yylval.string.text = "\\I\\";
		        yylval.string.len  = -1;
			return STRING;
                }
<CCMode,ccStyleMode>[\\](l?)dots/{noletter}  {
		        skipspaces();
		        yylval.string.text = "...";
		        yylval.string.len  = 3;
			return STRING;
                }
<CCMode,ccStyleMode>[\\]ccEndFont/{noletter}   {}

<NestingMode>[\[]    {
		    return '[';
		}
<NestingMode>[\]]    {
	  	    return ']';
		}

<NestingMode>[\(]    {
		    return '(';
		}
<NestingMode>[\)]    {
	  	    return ')';
		}

<INITIAL,MMODE,NestingMode>"---" {
		        skipspaces();
		        yylval.string.text = " - ";
		        yylval.string.len  = 3;
			return STRING;
                }

<INITIAL,MMODE,NestingMode>"--" {
		        yylval.string.text = "-";
		        yylval.string.len  = 1;
			return STRING;
                }

 /* TeX macros                             */
 /* -------------------------------------- */
<INITIAL,MMODE,NestingMode>[\\]((ref)|(ccTrue)|(ccFalse)|(kill)|(parskip)|(parindent))/{noletter}    {  // copy without warning
                    yylval.string.text = yytext;
                    yylval.string.len  = -1;
		    return STRING;    
                }
<INITIAL,MMODE,NestingMode,ccStyleMode>{texmacro}    {
                    if (actual_defining) {
                        yylval.string.text = yytext;
                    } else {
			yylval.string.text = fetchMacro( yytext);
			if ( yylval.string.text == NULL) {
			    yylval.string.text = yytext;
			    if ( warn_switch) {
				cerr << prog_name 
				     << ": warning: unknown macro " << yytext 
				     << " in line " << line_number
				     << " of file `" << in_filename << "'."
				     << endl;
			    }
			}
                    }
                    yylval.string.len  = -1;
                    actual_defining = false;
		    return STRING;
                }


 /* The rest: spaces and single characters */
 /* -------------------------------------- */
[\\]?{ws}	{
		    if ( *yytext == '\\') {
		        yylval.string.text = yytext + 1;
		        yylval.string.len  = yyleng - 1;
		    } else {
		        yylval.string.text = yytext;
		        yylval.string.len  = yyleng;
		    }
	  	    return SPACE;
		}
<CCMode,NestingMode,ccStyleMode>{ws}	{
	            yylval.string.text = yytext;
		    yylval.string.len  = yyleng;
	  	    return SPACE;
		}
<CCMode>.	{
	            yylval.character = yytext[0];
	  	    return CHAR;
		}
<INITIAL,NestingMode,ccStyleMode>[\\][/]  {}
<INITIAL,MMODE,NestingMode>[&]  {
                    if ( tab_tag) {
	                yylval.string.text = 
                            "</TD><TD ALIGN=LEFT VALIGN=TOP NOWRAP>";
		        yylval.string.len  = -1;
		        return STRING;		    
                    }
		    yylval.string.text = html_multi_character(yytext[0]);
		    yylval.string.len  = strlen( yylval.string.text);
		    return STRING;
                }

<INITIAL,NestingMode,MMODE,CPROGMode,ccStyleMode>.	{
	            yylval.character = yytext[0];
		    if ( is_html_multi_character( yylval.character)) {
		        yylval.string.text = html_multi_character(
                                                 yylval.character);
		        yylval.string.len  = strlen( yylval.string.text);
		        return STRING;
		    }
	  	    return CHAR;
		}
%%

void init_scanner( FILE* in){
    line_number  = 1;
    set_CCMode   = 0;
    set_NestingMode   = 0;
    set_INITIAL  = 0;
    unchecked_tag = 0;
    current_in_stream = in;
    stack_ptr = 0;
    yyrestart( in);
}

void skipspaces( void) {
    int c = yyinput();
    while( c && c <= ' ') {
	if ( c == '\n')
	    line_number++;
        c = yyinput();
    }
    unput( c);
}

void skipoptionalparam( void) {
    int c = yyinput();
    while( c && c <= ' ') {
	if ( c == '\n')
	    line_number++;
        c = yyinput();
    }
    if ( c == '[')
        while( c && c != ']') {
	    if ( c == '\n')
	        line_number++;
            c = yyinput();
        }
    else
        unput( c);
}
