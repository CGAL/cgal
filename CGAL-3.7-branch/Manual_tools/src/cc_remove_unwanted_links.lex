%{
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>



using namespace std;

/* Hack, to get rid of the yywrap. */
#define YY_SKIP_YYWRAP 1
#define yywrap() 1


bool removeLinks = false;

%}

%x ANCHORMODE
%x HEADINGMODE


ws              [ \t\n\r]*
nonws           [^ \t\n\r]
par             {ws}("<"[pP]">"{ws})*
glue            {ws}("<"[pP]">"{ws})*"<!GLUE>"{ws}("<"[pP]">"{ws})*

%%


"<!-- REMOVE_LINKS_BEGIN -->"            {  removeLinks = true; }
"<!-- REMOVE_LINKS_END -->"              {  removeLinks = false; }
[<][aA]{ws}[hH][rR][eE][fF][^>]*[>]      { if( removeLinks ) 
                                             BEGIN(ANCHORMODE); 
                                           else
                                             ECHO;
                                         }
<ANCHORMODE>"</A>"                       { BEGIN(INITIAL); }

"<"[hH][1-6]                   { BEGIN( HEADINGMODE); ECHO;}
<HEADINGMODE>"</"[hH][1-6]">"  { BEGIN( INITIAL); ECHO;}
<HEADINGMODE>.                 { ECHO;}


"</TABLE><!3>"{glue}"<!3><TABLE"[^>]*">" {
                               ; /* smooth adjacent tables */
}
"</TABLE><!2>"{glue}"<!2><TABLE"[^>]*">" { 
                               ; /* smooth adjacent tables */
}
"</TABLE><!".">"{glue}"<!"."><TABLE"     { 
                               /* smooth other tables   */
                               cout << "</TABLE>" << endl << "        <TABLE";
}
"</TABLE><!3>"{par}"<!3><TABLE"[^>]*">"  {
                               /* smooth adjacent tables */
                               cout  << "<TR><TD><BR></TD></TR>";
}
"</TABLE><!2>"{par}"<!2><TABLE"[^>]*">"  {
                               /* smooth adjacent tables */
                               cout << "<TR><TD><BR></TD></TR>";
}
"</TABLE><!".">"{par}"<!"."><TABLE"      { 
                               /* smooth other tables   */
                               cout << "</TABLE>" << endl << "<P>" << endl << endl << "        <TABLE";
}
{ws}"<"[pP]">"{ws}("<"[pP]">"{ws})*      {
                               /* Reduce <P>'s to one <P> */
                               cout << endl << "<P>" << endl << endl;
}
{ws}"<"[bB][rR]">"{ws}     { /* Pretty print single <BR> */
                               cout << "<BR>" << endl << endl;
}
"<!GLUE>"  { ; /* remove superfluous glues */ }
"<!3>"     { ; /* remove superfluous table hints */ }
"<!2>"     { ; /* remove superfluous table hints */ }

%%

int main( int argc, char** argv ) {
  yylex();
  return 0;
}
