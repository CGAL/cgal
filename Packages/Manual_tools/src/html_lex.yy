/**************************************************************************
 
  lexhtml.yy
  =============================================================
  Project   : CGAL merger tool for the specification task
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
#include <string.h>
#include <database.h>
#include <confightml.h>
#include <syntaxhtml.tab.h>

/* Set this flag to 1 to switch immediately to CCMode. */
int set_CCMode      = 0;
/* Set this flag to 1 to switch immediately to NestingMode. */
int set_NestingMode = 0;
/* Set this flag to 1 to switch back to INITIAL. */
int set_INITIAL     = 0;
int set_HTMLMODE    = 0;
int set_MMODE       = 0;

/* Tag to mark the unchecked keyword */
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

/* enable line number printing to cerr */
extern char  line_switch;

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
%x CCMode
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
CCidfier        ({CCletter}({CCletter}|{digit})*)
filename        [^ \t\n/\\\{\}\[\]()]+
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
ttblockintro    [\{][\\](tt)
emblockintro    [\{][\\](em)
itblockintro    [\{][\\]((it)|(sl))
scblockintro    [\{][\\](sc)
bfblockintro    [\{][\\](bf)

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

 /* Count line numbers in all modes for better error messages */
 /* --------------------------------------------------------- */
<INITIAL,CCMode,NestingMode,CPROGMode,MMODE,ITEMMODE,TEXONLYMODE,HTMLMODE,HTMLGROUPMode>[\n]	{
		    line_number++;
		    if ( line_switch)
		        cerr << "src-line " << line_number << endl;
		    return NEWLINE;
		}
<INITIAL,MMODE,NestingMode>[\\]"\n"      {
		    line_number++;
		    if ( line_switch)
		        cerr << "src-line " << line_number << endl;
	            yylval.string.text = " ";
		    yylval.string.len  = 1;
	  	    return SPACE;
}

 /* Handle include files      */
 /* ------------------------- */

[\\]((include)|(input))[\{]{w}   {  BEGIN ( INCLUDEMODE); }
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
<INITIAL,MMODE>[\\]"%"         {   /* Avoid the quoted comment symbol */
		    yylval.character   = '%';
		    return CHAR;
		}
<INITIAL,MMODE>"%".*[\n]{w}    {   /* Match one line TeX comments */
		    /* remove spaces in next line  */
		    unput( '\n');
		}
<INITIAL,MMODE>"%".*  { /* Match one line TeX comments */
                        /* at the last line in file */
		}
[\\]verb{noletter}   {   /* match LaTeX \verb"..." constructs */
		    BEGIN( VerbMode);
		    stop_character = yytext[ yyleng-1];
		    yylval.string.text = "<PRE>";
		    yylval.string.len  = 5;
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
	                yylval.string.text = "</PRE>";
		        yylval.string.len  = 6;
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
[\\]begin{w}[\{]class[\}]{w}   {
		    BEGIN( CCMode);
		    return BEGINCLASS;
		}
[\\]end{w}[\{]class[\}]   {
		    return ENDCLASS;
		}
[\\]begin{w}[\{]classtemplate[\}]{w}   {
		    BEGIN( CCMode);
		    return BEGINCLASSTEMPLATE;
		}
[\\]end{w}[\{]classtemplate[\}]   {
		    return ENDCLASSTEMPLATE;
		}
[\\]creationvariable{w}[\{]{w}[^\}]*{w}[\}]   {
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
		    delete[] creationvariable;
		    if ( formatted_creationvariable)
		        delete[] formatted_creationvariable;
		    creationvariable = newstr( r);
		    formatted_creationvariable = new char[ strlen( 
		            creationvariable) + 12];
		    strcpy( formatted_creationvariable, "<VAR>");
		    strcat( formatted_creationvariable, creationvariable);
		    strcat( formatted_creationvariable, "</VAR>");
	            yylval.string.text = r;
		    yylval.string.len  = s - r + 1;
		    return CREATIONVARIABLE;
		}		    
[\\]constructor/{noletter} {   /* constructor declaration: change to CCMode */
		    skipspaces();
		    BEGIN( CCMode);
		    return CONSTRUCTOR;
		}
[\\]method/{noletter}   {   /* method declaration: change to CCMode */
		    skipspaces();
		    BEGIN( CCMode);
		    return METHOD;
		}
[\\]function/{noletter} {   /* function declaration: change to CCMode */
		    skipspaces();
		    BEGIN( CCMode);
		    return FUNCTION;
		}
[\\]functiontemplate/{noletter} {   /* function template declaration: 
			       change to CCMode */
		    skipspaces();
		    BEGIN( CCMode);
		    return FUNCTIONTEMPLATE;
		}
[\\]variable/{noletter} {   /* variable declaration: change to CCMode */
		    skipspaces();
		    BEGIN( CCMode);
		    return VARIABLE;
		}
[\\]typedef/{noletter} {   /* typedef declaration: change to CCMode */
		    skipspaces();
		    BEGIN( CCMode);
		    return TYPEDEF;
		}
[\\]enum/{noletter} {   /* enum declaration: change to CCMode */
		    skipspaces();
		    BEGIN( CCMode);
		    return ENUM;
		}
[\\]globalfunction/{noletter} {   /* function declaration: change to CCMode */
		    skipspaces();
		    BEGIN( CCMode);
		    return GLOBALFUNCTION;
		}
[\\]globalfunctiontemplate/{noletter} {   /* function template declaration: 
			       change to CCMode */
		    skipspaces();
		    BEGIN( CCMode);
		    return GLOBALFUNCTIONTEMPLATE;
		}
[\\]globalvariable/{noletter} {   /* variable declaration: change to CCMode */
		    skipspaces();
		    BEGIN( CCMode);
		    return GLOBALVARIABLE;
		}
[\\]globaltypedef/{noletter} {   /* typedef declaration: change to CCMode */
		    skipspaces();
		    BEGIN( CCMode);
		    return GLOBALTYPEDEF;
		}
[\\]globalenum/{noletter} {   /* enum declaration: change to CCMode */
		    skipspaces();
		    BEGIN( CCMode);
		    return GLOBALENUM;
		}
[\\]declaration/{noletter} {   /* general declaration: change to CCMode */
		    skipspaces();
		    BEGIN( CCMode);
		    return DECLARATION;
		}
[\\]hidden/{noletter}   {
		    skipspaces();
	  	    return HIDDEN;
		}
[\\]unchecked/{noletter}  { 
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

[\\]"begin{HtmlOnly}" {
		    BEGIN( HTMLGROUPMode);
		    return HTMLBEGIN;
                }
<HTMLGROUPMode>[\\]"end{HtmlOnly}"  {
		    BEGIN( INITIAL);
		    return HTMLEND;
                }

[\\]"begin{TexOnly}" {
		    BEGIN( TEXONLYMODE);
		    return TEXONLYBEGIN;
                }
<TEXONLYMODE>[\\]"end{TexOnly}"  {
		    BEGIN( INITIAL);
		    return TEXONLYEND;
                }
<TEXONLYMODE>.    {
		    yylval.character   = yytext[0];
		    return CHAR;
		}
[\\]LatexHtml{w}[\{]  {
		    return LATEXHTML;
                }
[\\]Anchor{w}[\{]  {
		    /* The first parameter is the URL, the second is the */
                    /* message that will be highlighted */
		    BEGIN( HTMLMODE);
		    return ANCHOR;
                }
<HTMLMODE>[\}]  {
		    BEGIN( INITIAL);
		    return '}';
                }
<HTMLMODE,HTMLGROUPMode>.     {
		    yylval.character   = yytext[0];
		    return CHAR;
		}

 /* Specialized keywords from the manual style */
 /* -------------------------------------------------------------- */
[\\]CCstyle/{noletter}  {   /* CCstyle formatting: change to NestingMode */
		    skipspaces();
		    BEGIN( NestingMode);
		    return CCSTYLE;
		}
<INITIAL,MMODE,NestingMode>[\\]var/{noletter}      {
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
<INITIAL,MMODE,NestingMode>[\\]purevar/{noletter}      {
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
<INITIAL,MMODE,NestingMode>[\\]classname/{noletter} {
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
<INITIAL,MMODE,NestingMode>[\\]pureclassname/{noletter} {
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
<INITIAL,MMODE,NestingMode>[\\]classtemplatename/{noletter} {
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
<INITIAL,MMODE,NestingMode>[\\]puretemplatename/{noletter} {
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
[\\]CCsection/{noletter} {  
		    skipspaces();
		    return CCSECTION; 
                 }
[\\]CC/{noletter}        {
		    skipspaces();
	            yylval.string.text = "C++";
		    yylval.string.len  = 3;
		    return STRING;
                 }
[\\]gg/{noletter}        {
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
[\\]definition/{noletter} {
		    skipspaces();
	            yylval.string.text = "<H3>Definition</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]creation/{noletter}  {
		    skipspaces();
	            yylval.string.text = "<H3>Creation</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]operations/{noletter} {
		    skipspaces();
	            yylval.string.text = "<H3>Operations</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]implementation/{noletter} {
		    skipspaces();
	            yylval.string.text = "<H3>Implementation</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]example/{noletter}   {
		    skipspaces();
	            yylval.string.text = "<H3>Example</H3>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]precond/{noletter}      {
		    skipspaces();
	            yylval.string.text = "<BR><STRONG>Precondition: </STRONG>";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]threecolumns/{noletter} {
		    skipspaces();
		    return GOBBLETWOPARAMS;
                 }
[\\]constructorcolumn/{noletter} {
		    skipspaces();
		    return GOBBLEONEPARAM;
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
<INITIAL,MMODE,NestingMode>[\\]ldots/{noletter}     {
		    skipspaces();
	            yylval.string.text = "...";
		    yylval.string.len  = 3;
		    return STRING;
                 }
<MMODE>[\\]leq/{noletter}          {
		    skipspaces();
	            yylval.string.text = "&lt;=";  // &le;  not yet supported
		    yylval.string.len  = 5;
		    return STRING;
                 }
<MMODE>[\\]geq/{noletter}          {
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
<INITIAL,MMODE,NestingMode>[\\]"&"          {
		     yylval.string.text = "&amp;";
		     yylval.string.len  = 5;
		     return STRING;
		 }
<INITIAL,MMODE,NestingMode>[\\][_^#$]  {
		    yylval.character = yytext[1];
		    return CHAR;
                 }
<INITIAL,MMODE,NestingMode>[~]              |
<INITIAL,MMODE,NestingMode>[\\]{space}      {
	            yylval.string.text = " ";
		    yylval.string.len  = 1;
	  	    return SPACE;
		 }
<INITIAL,MMODE,NestingMode>[\\][\\]         {
		    skipoptionalparam();
	            yylval.string.text = " ";
		    yylval.string.len  = 1;
	  	    return SPACE;
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
  <MMODE>[\\]times/{noletter}       { SET( "&times;");    return STRING;}
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
<MMODE>[\\]times/{noletter}       { SET( "x");    return STRING;}
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


 /* keywords from TeX/LaTeX that should vanish in HTML             */
 /* -------------------------------------------------------------- */
[\\]((smallskip)|(protect)|(sloppy))/{noletter}           {}
[\\]((maketitle)|(tableofcontents))/{noletter}            {}
[\\]((begin)|(end))[\{]document[\}]                       {}

[\\]newsavebox{w}[\{]          |
[\\]usebox{w}[\{]              |
[\\][*]?hspace{w}[\{]          |
[\\][*]?vspace{w}[\{]          |
[\\]g?def{w}[\\]{letter}+[^\{]*[\{]  {
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
[\\]((textwidth)|(textheight)|(topmargin)|(evensidemargin)|(oddsidemargin)|(headsep)|(parindent)|(parskip)){w}[-+0-9.]+{w}..  {}

[\\]newcommand{w}[\{][^\}]*[\}]([\[][^\]]*[\]])[\{]   {
	            /* CCstyle formatting: change to NestingMode */
		    BEGIN( NestingMode);
		    return IGNOREBLOCK;
		}
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

[\\](([mvhf]box)|(parbox)|([hv]fill)|(nopagebreak)|(nolinebreak)|(linebreak)|(samepage))/{noletter}        {
		    skipoptionalparam();
	            yylval.string.text = " ";
		    yylval.string.len  = 1;
		    return SPACE;		    
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
[\\]cite{w}([\[][^\]]*[\]])?[\{][^\}]*[\}]     {
	            yylval.string.text = yytext;
		    yylval.string.len  = yyleng;
		    return CITE;		    
                 }
[\\]bibitem{w}/{noletter}  {
		    BEGIN( NestingMode);
		    return BIBITEM;
                }


 /* Grouping symbols */
 /* ---------------- */
[\\][\{]        {
	            yylval.character = '{';
	  	    return CHAR;
		}
[\\][\}]        {
	            yylval.character = '}';
	  	    return CHAR;
		}
[\\]left.       {
	            yylval.character = yytext[5];
	  	    return CHAR;
		}
[\\]right.      {
	            yylval.character = yytext[6];
	  	    return CHAR;
		}
<INITIAL,CCMode,NestingMode>[\{]    {
		    return '{';
		}
<INITIAL,CCMode,NestingMode>[\}]    {
	  	    return '}';
		}

{ttblockintro}  {  /* A couple of TeX styles like {\tt ... */
		    return TTBLOCKINTRO;
		}
{emblockintro}  {   return EMBLOCKINTRO; }
{itblockintro}  {   return ITBLOCKINTRO; }
{scblockintro}  {   return SCBLOCKINTRO; }
{bfblockintro}  {   return BFBLOCKINTRO; }

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
<CCMode,NestingMode>{ws}	{
	            yylval.string.text = yytext;
		    yylval.string.len  = yyleng;
	  	    return SPACE;
		}
<CCMode>.	{
	            yylval.character = yytext[0];
	  	    return CHAR;
		}
<INITIAL,NestingMode,MMODE,CPROGMode>.	{
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
    while( c && c <= ' ')
        c = yyinput();
    unput( c);
}

void skipoptionalparam( void) {
    int c = yyinput();
    while( c && c <= ' ')
        c = yyinput();
    if ( c == '[')
        while( c && c != ']')
            c = yyinput();
    else
        unput( c);
}
