/**************************************************************************
 
  internal_macros.C
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Internal macro definitions.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/


#include <internal_macros.h>

#include <mstring.h>
#include <ctype.h>
#include <strstream.h>
#include <string_conversion.h>
#include <macro_dictionary.h>
#include <cpp_formatting.h>
#include <lex_include.h>
#include <html_lex.h>
#include <html_error.h>
#include <html_config.h>
#include <output.h>

// ======================================================================
// External visible functions
// ======================================================================


// ======================================================================
// Index keys and sorting
// ======================================================================

// Sort keys for the index
// ----------------------------------
// The index is organized in a set of fixed topics. 

// sort_keys are used to sort the index according to different sections.
// A sort_key is followed by a 0 to indicate the section title.
// A sort_key is followed by a 1 to indicate a normal entry.

const string sort_key_class            = "<!sort B";
const string  sort_key_nested_type      = "<!sort C";
const string  sort_key_struct           = "<!sort B";  // structs like classes
const string  sort_key_enum             = "<!sort E";
const string  sort_key_enum_tags        = "<!sort F";
const string  sort_key_typedef          = "<!sort G";
const string  sort_key_variable         = "<!sort H";
const string  sort_key_function         = "<!sort I";
const string  sort_key_member_function  = "<!sort J";

const string& find_sort_key( const string& txt) {
    if ( txt == "class")
	return sort_key_class;
    if ( txt == "nested_type")
	return sort_key_nested_type;
    if ( txt == "struct")
	return sort_key_struct;
    if ( txt == "enum")
	return sort_key_enum;
    if ( txt == "enum_tags")
	return sort_key_enum_tags;
    if ( txt == "typedef")
	return sort_key_typedef;
    if ( txt == "variable")
	return sort_key_variable;
    if ( txt == "function")
	return sort_key_function;
    if ( txt == "member_function")
	return sort_key_member_function;
    printErrorMessage( UnknownIndexCategoryError);
    return sort_key_class;
}

void write_headers_to_index( ostream& out){
    out << sort_key_class       << "0 !><P><LI><B>Classes</B><P>" << endl;
    out << sort_key_nested_type << "0 !><P><LI><B>Nested Types</B><P>" << endl;
    //out << sort_key_struct      << "0 !><P><LI><B>Structs</B><P>" << endl;
    out << sort_key_enum        << "0 !><P><LI><B>Enums</B><P>" << endl;
    out << sort_key_enum_tags   << "0 !><P><LI><B>Enum Tags</B><P>" << endl;
    out << sort_key_typedef     << "0 !><P><LI><B>Typedefs</B><P>" << endl;
    out << sort_key_variable    << "0 !><P><LI><B>Global Variables and "
                                   "Consts</B><P>" << endl;
    out << sort_key_function    << "0 !><P><LI><B>Functions</B><P>" << endl;
    //    out << sort_key_member_function  << "0 !><P><LI><B>Member Functions</B><P>"
    //	                             << endl;
}

static int index_anchor_counter = 0;

string handleHtmlIndex( const string& category, 
			const string& sort_key, 
			const string& formatted_reference) {
    *index_stream << category << '1';
    filter_for_index_comment( *index_stream, sort_key);
    *index_stream << "!><UL><LI><A HREF=\"" << current_filename
		  << "#Index_anchor_" << index_anchor_counter << "\">"
		  << formatted_reference << "</A></UL>" << endl;
    return string("\n<A NAME=\"Index_anchor_") 
	 + int_to_string( index_anchor_counter++)
	 + "\"></A>\n";
}

string handleHtmlIndexC( const string& category, const string& item) {
    return handleHtmlIndex( category, item, convert_C_to_html( item));
}

string handleHtmlIndex( const string& category, const string& item) {
    return handleHtmlIndex( category, item, 
			    convert_fontified_ascii_to_html( item));
}


static int cross_link_anchor_counter = 0;

string handleHtmlCrossLink( string key, bool tmpl_class) {
    crop_string( key);
    if ( key.empty()) {
	printErrorMessage( EmptyCrossLinkError);
	exit( 1);
    }

    char *tmp_name = convert_fontified_ascii_to_html( key);
    *anchor_stream << "[a-zA-Z0-9_]\"" << tmp_name
		   << "\"    { ECHO; }" << endl;
    *anchor_stream << '"' << tmp_name
		   << "\"/{noCCchar}    { wrap_anchor( \"" << current_filename
		   << "#Cross_link_anchor_" << cross_link_anchor_counter 
		   << "\", yytext); }" << endl;
    if ( tmpl_class) {
        *anchor_stream << '"' << tmp_name
		       << "\"{ws}\"&lt;\"{ws}{CCidfier}{ws}\"&gt;\"    {\n"
		       << "        wrap_anchor( \"" << current_filename
		       << "#Cross_link_anchor_" << cross_link_anchor_counter 
		       << "\", yytext); }" << endl;
    }
//     *anchor_stream << '"' << tmp_name
// 		   << "\"/{noCCchar}    { fputs( \"<A HREF=\\\""
// 		   << current_filename << "#Cross_link_anchor_" 
// 		   << cross_link_anchor_counter << "\\\">"
// 		   << tmp_name << "</A>\", stdout); }" 
// 		   << endl;
//     if ( tmpl_class) {
//         *anchor_stream << '"' << tmp_name
// 		       << "\"{ws}\"&lt;\"{ws}{CCidfier}{ws}\"&gt;\"    {\n"
// 		       << "        fputs( \"<A HREF=\\\""
// 		       << current_filename << "#Cross_link_anchor_" 
// 		       << cross_link_anchor_counter << "\\\">\", stdout);\n" 
// 		       << "        ECHO;\n"
// 		       << "        fputs( \"</A>\", stdout); }\n"
// 		       << endl;
//     }
    delete[] tmp_name;

    return string("\n<A NAME=\"Cross_link_anchor_") 
	 + int_to_string( cross_link_anchor_counter++) + "\"></A>\n";
}


// Chapter and section
// =================================================

string chapter_title;

void handleChapter(  const Text& T) {
    string new_main_filename = macroX( "\\lciInputPath")
                             + macroX( "\\lciChapterPrefix")
                             + macroX( "\\lciInputFilenameBase")
	                     + macroX( "\\lciHtmlSuffix");
    if ( new_main_filename == main_filename) {
        printErrorMessage( ChapterStructureError);
	return;
    }
    if ( class_stream != 0) {
        if ( ! chapter_title.empty())
	    // navigation footer
	    *class_stream << "<HR><B> Return to chapter:</B> <A HREF=\"" 
			  << main_filename 
			  << "\">" << chapter_title << "</A>" << endl;
        close_html( *class_stream);
	assert_file_write( *class_stream, class_filename);
	delete   class_stream;
	class_filename = string();
	class_stream = 0;
    }
    chapter_title = string( text_block_to_string( T));
    if ( main_stream != &cout && main_stream != pre_stream) {
        // navigation footer
        *main_stream << "<HR> Next chapter: <A HREF=\"" 
		     << new_main_filename 
		     << "\">" << chapter_title << "</A>" << endl;
        close_html( *main_stream);
	assert_file_write( *main_stream, main_filename);
	delete   main_stream;
	main_stream = 0;
    }
    main_filename = new_main_filename;
    main_stream = open_file_for_write( tmp_path + main_filename);
    current_ostream  = main_stream;
    current_filename = main_filename;
    insertInternalGlobalMacro( "\\lciOutputFilename", current_filename);
    insertInternalGlobalMacro( "\\lciMainFilename",   main_filename);
    open_html( *main_stream);

    *main_stream << "<H1>" << chapter_title << "</H1>" << endl;

    // table of contents
    *contents_stream << "    <LI> <A HREF=\"" << main_filename 
		     << "\">" << chapter_title << "</A>" << endl;
}

void handleBiblio(  const Text& T) {
    ostream* out = open_file_for_write( tmp_path + macroX("\\lciBibFilename"));
    istream* in  = open_config_file( macroX( "\\lciBiblioHeader"));
    filter_config_file( *in, *out);

    print_html_text_block( *out, T);
    *out << "</TD></TR>" << endl;

    delete in;
    in = open_config_file( macroX( "\\lciBiblioFooter"));
    filter_config_file( *in, *out);

    assert_file_write( *out, macroX( "\\lciBibFilename"));
    delete in;
    delete out;
}

void handleSection(  const Text& T) {
    static int section_counter = 1;
    char* section = text_block_to_string( T);
    *current_ostream << "<A NAME=\"Section_" << section_counter << "\"></A>"
		    << endl;
    *current_ostream << "<H2>" << section << "</H2>" << endl;

    // table of contents
    *contents_stream << "        <UL><LI><A HREF=\"" << current_filename
		     << "#Section_" << section_counter << "\">"
		     << section << "</A></UL>" << endl;

    ++ section_counter;
    delete[] section;
}

void handleLabel( const char* l, size_t len) { // trust only len!
    /* The lexical processing has removed the parantheses around */
    /* \ref{...} macros from TeX, so here is the correct pattern match */
    /* to find them in the pre-HTML text */
    char* s = new char[len+1];
    strncpy( s, l, len);
    s[len] = '\0';
    *anchor_stream << "[\\[]ref[:]\"" << s
		   << "\"[\\]]    { fputs( \"<A HREF=\\\""
		   << current_filename 
		   << "#" << s << "\\\">" << reference_icon 
		   << "</A>\", stdout); }" << endl;    
    // There are two special ref commands defined within the manual
    *anchor_stream << "[\\\\]Chapter[ \\t\\n]*\"" << s
		   << "\"    { fputs( \"Chapter <A HREF=\\\""
		   << current_filename 
		   << "#" << s << "\\\">" << reference_icon 
		   << "</A>\", stdout); }" << endl;    
    *anchor_stream << "[\\\\]Section[ \\t\\n]*\"" << s
		   << "\"    { fputs( \" Section <A HREF=\\\""
		   << current_filename 
		   << "#" << s << "\\\">"  << reference_icon 
		   << "</A>\", stdout); }" << endl;    
    delete[] s;
}

// Opens a new classfile. Only a filename and a HTML formatted reference
// text are given.
void handleClassFile( const string& filename, 
		      const string& formatted_reference) {
    if ( class_stream != 0) {
        // navigation footer
        *class_stream << "<HR> <B>Next:</B> " << formatted_reference << endl;
        close_html( *class_stream);
	assert_file_write( *class_stream, class_filename);
	delete   class_stream;
	class_stream = 0;
    }

    class_filename = filename;
    class_stream = open_file_for_write( tmp_path + class_filename);
    open_html( *class_stream);
    // Make a hyperlink in the chapter to the class file.
    if ( main_stream != &cout) {
        *main_stream  << "<UL><LI>\n" << formatted_reference
		      << ".</UL>\n" << endl;
    }

    // table of contents
    if ( macroIsTrue( "\\lciIfHtmlClassToc"))
	*contents_stream << "        <UL><LI> " << formatted_reference
			 << "</UL>" << endl;

    current_ostream  = class_stream;
    current_filename = class_filename;
    insertInternalGlobalMacro( "\\lciOutputFilename", current_filename);
}

void handleClassFileEnd( void) {
    /* ...  Hack to implement the link from one class to the next class
	close_html( *class_stream);
	delete   class_stream;
	delete[] class_filename;
	class_stream = 0;
	class_filename = 0;
    ... */
    current_ostream  = main_stream;
    current_filename = main_filename;
    insertInternalGlobalMacro( "\\lciOutputFilename", current_filename);
}

void handleClassEnvironment() {
    template_class_name = macroX("\\ccPureClassTemplateName");
    string formatted_template_class_name = 
	convert_C_to_html( template_class_name);
    class_name = macroX( "\\ccPureClassName");

    if ( macroIsTrue( "\\lciIfHtmlClassNotInline") &&
	 macroIsTrue( "\\lciIfHtmlClassFile")) {
	// Start a new class file.
	string filename = remove_font_commands( class_name)
	                + macroX( "\\lciHtmlSuffix");
	string contents( " Class declaration of ");
	contents += string("<A HREF=\"") + filename + "\">"
	            + formatted_template_class_name + "</A>";
	handleClassFile( filename, contents);
    }

    if ( macroIsTrue( "\\lciIfHtmlClassIndex") && 
	 macroIsTrue( "\\lciIfHtmlIndex"))    // Index.
	*current_ostream << handleHtmlIndex( sort_key_class, 
					     template_class_name,
					     formatted_template_class_name);

    if ( macroIsTrue( "\\lciIfHtmlClassLinks") && 
	 macroIsTrue( "\\lciIfHtmlLinks")) {   // Cross links.
	// Generate a substitution rule for hyperlinking,
        // but only if classname is longer than 1 character.
	if ( class_name.size() > 1) 
	    *current_ostream << handleHtmlCrossLink( class_name, true);
    }
}

void handleClassNameEnd( void) {
    class_name = string();
    template_class_name = string();
}

void handleClassEnd( void) {
    handleClassNameEnd();
    if ( macroIsTrue( "\\lciIfHtmlClassNotInline") &&
	 macroIsTrue( "\\lciIfHtmlClassFile"))
	handleClassFileEnd();
}

void handleHtmlClassFile( const string& filename, const Text& T) {
    char*  s = text_block_to_string( T);
    string p = string("<A HREF=\"") + filename + "\">" + s + "</A>";
    handleClassFile( filename, p);
    delete[] s;    
}




// ======================================================================
// Functions used internally
// ======================================================================


#define NParamCheck( m, o) \
    if ( n != (m) && opt != (o)) {        \
	printErrorMessage( MacroParamNumberError);           \
	return string("<internal macro expansion failed>");  \
    }

// Major macro definition support
// ======================================================================
string                // macro
newcommand( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    string macroname( param[0]);
    crop_string( macroname);
    remove_separator( macroname);
    insertMacro( macroname, in_string->name(), in_string->line(), param[1], 0);
    return string();
}

string                // macro
newcommand_opt( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 1);
    string macroname( param[0]);
    crop_string( macroname);
    remove_separator( macroname);
    int m = atoi( param[1].c_str());
    if ( m < 1 || m > 9)
	printErrorMessage( NParamRangeError);
    else
	insertMacro( macroname, in_string->name(), in_string->line(), 
		     param[2], m);
    return string();
}

string                // macro
global_newcommand( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    string macroname( param[0]);
    crop_string( macroname);
    remove_separator( macroname);
    insertGlobalMacro( macroname, in_string->name(), in_string->line(), 
		       param[1], 0);
    return string();
}

string
undef( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);
    string macroname( param[0]);
    crop_string( macroname);
    remove_separator( macroname);
    eraseMacro( macroname);
    return string();
}

string
optional_param_at_end( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    string macroname( param[0]);
    crop_string( macroname);
    remove_separator( macroname);
    string s( param[1]);
    crop_string( s);
    setMacroOptEnd( macroname, atoi(s.c_str()));
    return string();
}

string
begin_scope( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);
    pushMacroScope();
    return string();
}

string
end_scope( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);
    popMacroScope();
    return string();
}

// If control structures
// ======================================================================
string                // macro
if_defined( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);
    string macroname( param[0]);
    crop_string( macroname);
    remove_separator( macroname);
    if ( definedMacro( macroname))
	return string("\\lcTrue");
    return string("\\lcFalse");
}

string
if_equal( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    string s0( param[0]);
    string s1( param[1]);
    remove_separator( s0);
    remove_separator( s1);
    if ( s0 == s1)
	return string("\\lcTrue");
    return string("\\lcFalse");
}

string
if_less( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    string s0( param[0]);
    string s1( param[1]);
    remove_separator( s0);
    remove_separator( s1);
    if ( atoi( s0.c_str()) < atoi(s1.c_str()))
	return string("\\lcTrue");
    return string("\\lcFalse");
}

string
if_less_or_equal( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    string s0( param[0]);
    string s1( param[1]);
    remove_separator( s0);
    remove_separator( s1);
    if ( atoi( s0.c_str()) <= atoi(s1.c_str()))
	return string("\\lcTrue");
    return string("\\lcFalse");
}

// Arithmetic
// ======================================================================
string
add_to( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    string macroname( param[0]);
    crop_string( macroname);
    int sum = atoi( expandFirstMacro( macroname).c_str()) 
            + atoi(param[1].c_str());
    insertMacro( macroname, in_string->name(), in_string->line(), 
		 int_to_string( sum));
    return string();
}

string
mult_to( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    string macroname( param[0]);
    crop_string( macroname);
    int prod = atoi( expandFirstMacro( macroname).c_str()) 
             * atoi(param[1].c_str());
    insertMacro( macroname, in_string->name(), in_string->line(),
		 int_to_string( prod));
    return string();
}

string
global_add_to( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    string macroname( param[0]);
    crop_string( macroname);
    int sum = atoi( expandFirstMacro( macroname).c_str()) 
            + atoi(param[1].c_str());
    insertGlobalMacro( macroname, in_string->name(), in_string->line(), 
		       int_to_string( sum));
    return string();
}

string
global_mult_to( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    string macroname( param[0]);
    crop_string( macroname);
    int prod = atoi( expandFirstMacro( macroname).c_str()) 
             * atoi(param[1].c_str());
    insertGlobalMacro( macroname, in_string->name(), in_string->line(),
		       int_to_string( prod));
    return string();
}

// String conversion
// ======================================================================
string 
string_to_upper( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);
    string s( param[0]);
    for ( int i = 0; i < s.size(); i++)
	s[i] = toupper( s[i]);
    return s;
}

// Error and message output
// ======================================================================
string 
html_error( const string&, string param[], size_t n, size_t opt) {
    cerr << endl;
    if ( n != 1 && opt != 0)
	printErrorMessage( MacroParamNumberError);
    else {
	string s( param[0]);
	remove_separator( s);
	cerr << "Error: " << s << '.';
    }
    printErrorMessage( UserDefinedError);
    return string();
}

string
html_message( const string&, string param[], size_t n, size_t opt) {
    if ( n != 1 && opt != 0)
	printErrorMessage( MacroParamNumberError);
    else if ( ! quiet_switch) {
	string s( param[0]);
	remove_separator( s);
	cerr << s;
    }
    return string();
}

// LaTeX commands
// ======================================================================
string
bib_cite( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    string s0( param[0]);
    string s1( param[1]);
    remove_separator( s0);
    remove_separator( s1);
    // A rule to substitute key by item in the bibliography.
    // param[0] == key, s1 == replacement text
    *anchor_stream << "\"<A NAME=\"[\"]\"Biblio_" << s0
		   << "\"[\"]\"></A><B>["  << s0
		   << "]</B>\"        { fputs( \"<A NAME=\\\"Biblio_"
		   << s0 << "\\\"></A><B>[" << s1 
		   << "]</B>\", stdout); }" 
		   << endl;
    // A rule to substitute key by item for cite's in the main text.
    *anchor_stream << "\"<A HREF=\"[\"]\"" 
		   << fetchMacroBody( "\\lciBibFilename") << "#Biblio_" 
		   << s0 << "\"[\"]\">"  << s0
		   << "\"        { fputs( \"<A HREF=\\\"" 
		   << fetchMacroBody( "\\lciBibFilename")
		   << "#Biblio_" << s0 << "\\\">" << s1 
		   << "\", stdout); }" 
		   << endl;
    return string();
}

// C++ Formatting
// ======================================================================
string
cc_style( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);  // uses cc_string as parameter, param[0] is a dummy
                         // to avoid scanning of spaces behind \ccc
    string s( "\\lciRawOutputN{");
    char* p = convert_C_to_html( cc_string);
    s += int_to_string( strlen(p)) + '}' + p;
    delete[] p;
    return s;
}

string
cc_parameter( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);  // uses cc_string as parameter
    return cc_string;
}

string
two_column_layout( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);  // uses cc_string as parameter, param[0] gives
                         // the declaration category.
    crop_string( cc_string);
    compress_spaces_in_string( cc_string);
    handle_two_column_layout( param[0][0], cc_string.c_str());
    return string();
}

string
three_column_layout( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);  // uses cc_string as parameter, param[0] gives
                         // the declaration category, param[1] says whether
                         // the comment is empty.
    crop_string( cc_string);
    compress_spaces_in_string( cc_string);
    handle_three_column_layout( param[0][0], cc_string.c_str(), 
				param[1][0] == '1');
    return string();
}


// HTML Converter Features
// ======================================================================
string
html_index( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);  // param[0] is index category, param[1] is text
    string key = find_sort_key( param[0]);
    string s = string( "\n\\lcRawHtml{<A NAME=\"Index_anchor_") 
	 + int_to_string( index_anchor_counter) + "\"></A>}\n"
	   "\\lciPushOutput{index}\\lcRawHtml{" + key + "1"
         + filter_for_index_comment( param[1]) + "!><UL><LI><A HREF=\""
         + current_filename + "#Index_anchor_" 
	 + int_to_string( index_anchor_counter) + "\">}" + param[1] 
	 + "\\lcRawHtml{</A></UL>\n}\\lciPopOutput";
    ++index_anchor_counter;
    return s;
}

string
html_index_C( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);  // uses cc_string as index item, param[0] is
                         // index category
    return string( "\\lcRawHtml{")
	+ handleHtmlIndexC( find_sort_key( param[0]), cc_string)
	+ '}';
}

string
cross_link( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);  // uses cc_string as crosslink key, param[0] is
                         // used as dummy to avoid parsing of spaces.
    return string( "\\lcRawHtml{") + handleHtmlCrossLink( cc_string) + '}';
}

string
pop_output( const string&, string [], size_t n, size_t opt) {
    NParamCheck( 0, 0);
    pop_current_output();
    return string();
}

string
push_output( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);
    push_current_output( param[0]);
    return string();
}

string
open_biblio( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);
    push_current_output();
    current_filename = tmp_path + macroX("\\lciBibFilename");
    current_ostream = open_file_for_write( current_filename);
    istream* in  = open_config_file( macroX( "\\lciBiblioHeader"));
    filter_config_file( *in, *current_ostream);
    delete in;
    return string();
}

string
close_biblio( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);
    *current_ostream << "</TD></TR>" << endl;
    istream* in  = open_config_file( macroX( "\\lciBiblioFooter"));
    filter_config_file( *in, *current_ostream);
    assert_file_write( *current_ostream, macroX( "\\lciBibFilename"));
    delete in;
    delete current_ostream;
    pop_current_output();
    return string();
}

string
html_class_file_end( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);
    handleClassFileEnd();
    return string();
}

// Initialize
// ======================================================================
void init_internal_macros() {
    insertInternalGlobalMacro( "{", "\\begingroup");
    insertInternalGlobalMacro( "}", "\\endgroup");

    insertInternalGlobalMacro( "\\newcommand", newcommand, 2);
    insertInternalGlobalMacro( "\\newcommand@mom", newcommand_opt);
    insertInternalGlobalMacro( "\\def", newcommand, 2);
    insertInternalGlobalMacro( "\\gdef", global_newcommand, 2);
    insertInternalGlobalMacro( "\\lciUndef", undef, 1);
    insertInternalGlobalMacro( "\\lciOptionalParameterAtEnd", 
			       optional_param_at_end, 2);
    insertInternalGlobalMacro( "\\lciBeginScope", begin_scope, 0);
    insertInternalGlobalMacro( "\\lciEndScope", end_scope, 0);
    insertInternalGlobalMacro( "\\lciIfDefined", if_defined, 1);
    insertInternalGlobalMacro( "\\lciIfEqual", if_equal, 2);
    insertInternalGlobalMacro( "\\lciIfLess", if_less, 2);
    insertInternalGlobalMacro( "\\lciIfLessOrEqual", if_less_or_equal, 2);
    insertInternalGlobalMacro( "\\lciAddTo", add_to, 2);
    insertInternalGlobalMacro( "\\lciMultTo", mult_to, 2);
    insertInternalGlobalMacro( "\\lciGlobalAddTo", global_add_to, 2);
    insertInternalGlobalMacro( "\\lciGlobalMultTo", global_mult_to, 2);
    insertInternalGlobalMacro( "\\lciToUpper", string_to_upper, 1);
    insertInternalGlobalMacro( "\\lciError", html_error, 1);
    insertInternalGlobalMacro( "\\lciMessage", html_message, 1);

    insertInternalGlobalMacro( "\\lciBibCite", bib_cite, 2);

    insertInternalGlobalMacro( "\\lciStyle", cc_style, 1);
    insertInternalGlobalMacro( "\\lciCCParameter", cc_parameter, 0);

    insertInternalGlobalMacro( "\\lciTwoColumnLayout",   
			       two_column_layout,   1);
    insertInternalGlobalMacro( "\\lciThreeColumnLayout", 
			       three_column_layout, 2);

    insertInternalGlobalMacro( "\\lciHtmlIndex", html_index, 2);
    insertInternalGlobalMacro( "\\lciHtmlIndexC", html_index_C, 1);
    insertInternalGlobalMacro( "\\lciHtmlCrossLink", cross_link, 1);
    insertInternalGlobalMacro( "\\lciHtmlClassFileEnd", html_class_file_end,0);

    insertInternalGlobalMacro( "\\lciPopOutput",  pop_output,  0);
    insertInternalGlobalMacro( "\\lciPushOutput", push_output, 1);
    insertInternalGlobalMacro( "\\lciOpenBibliography",  open_biblio, 0);
    insertInternalGlobalMacro( "\\lciCloseBibliography", close_biblio, 0);
}

// EOF //






