/**************************************************************************
 
  internal_macros.C
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Internal macro definitions.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
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

const string  sort_key_concept          = "<!sort B";
const string  sort_key_class            = "<!sort C";
const string  sort_key_struct           = "<!sort C";  // structs like classes
const string  sort_key_nested_type      = "<!sort D";
const string  sort_key_enum             = "<!sort E";
const string  sort_key_enum_tags        = "<!sort F";
const string  sort_key_typedef          = "<!sort G";
const string  sort_key_macro            = "<!sort H";
const string  sort_key_variable         = "<!sort I";
const string  sort_key_function         = "<!sort J";
const string  sort_key_member_function  = "<!sort K";

const string& find_sort_key( string txt) {
    if ( txt.size() > 0)
	txt[0] = tolower( txt[0]);
    if ( txt == "concept")
	return sort_key_concept;
    if ( txt == "functionObjectConcept")
	return sort_key_concept;
    if ( txt == "class")
	return sort_key_class;
    if ( txt == "functionObjectClass")
	return sort_key_class;
    if ( txt == "struct")
	return sort_key_struct;
    if ( txt == "nested_type")
	return sort_key_nested_type;
    if ( txt == "enum")
	return sort_key_enum;
    if ( txt == "enum_tags")
	return sort_key_enum_tags;
    if ( txt == "typedef")
	return sort_key_typedef;
    if ( txt == "macro")
	return sort_key_macro;
    if ( txt == "variable")
	return sort_key_variable;
    if ( txt == "constant")
	return sort_key_variable;
    if ( txt == "function")
	return sort_key_function;
    if ( txt == "member_function")
	return sort_key_member_function;
    printErrorMessage( UnknownIndexCategoryError);
    return sort_key_class;
}

void write_headers_to_index( ostream& out){
    out << sort_key_concept     << "0 !><P><LI><B>Concepts</B><P>" << endl;
    out << sort_key_class       << "0 !><P><LI><B>Classes</B><P>" << endl;
    //out << sort_key_struct      << "0 !><P><LI><B>Structs</B><P>" << endl;
    out << sort_key_nested_type << "0 !><P><LI><B>Nested Types</B><P>" << endl;
    out << sort_key_enum        << "0 !><P><LI><B>Enums</B><P>" << endl;
    out << sort_key_enum_tags   << "0 !><P><LI><B>Enum Tags</B><P>" << endl;
    out << sort_key_typedef     << "0 !><P><LI><B>Typedefs</B><P>" << endl;
    out << sort_key_macro       << "0 !><P><LI><B>Macros</B><P>" << endl;
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


// Chapter and Class Files
// =================================================

static int next_class_link_counter = 0;
static int next_class_link_last    = 0;
static int chapter_num = 0;

string chapter_title;
string part_title = "";

void handleChapter(  const Buffer_list& T) {
    next_class_link_last = 0;
    chapter_num++;
    //string new_main_filename = macroX( "\\lciInputPath")
    //                         + macroX( "\\lciChapterPrefix")
    //                         + macroX( "\\lciInputFilenameBase")
    //                         + macroX( "\\lciHtmlSuffix");
    string new_main_filename = macroX( "\\lciChapterPrefix")
                             + macroX( "\\lciInputFilenameBase")
	                     + macroX( "\\lciHtmlSuffix");
    if ( new_main_filename == main_filename) {
        printErrorMessage( ChapterStructureError);
	return;
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
        



    if (macroIsTrue("\\lciIfRefCrossLink"))
      *main_stream << "<TABLE WIDTH=100%> <TR> <TD ALIGN=LEFT VALIGN=TOP> <H1>" 
             << chapter_title << "</H1> </TD>" << endl;
    else
      *main_stream << "<H1>" << chapter_title  << "</H1>" << endl;


    // table of contents
   
    if (macroIsTrue("\\lciIfMultipleParts"))
       *contents_stream << " <TR> <TD>"<< chapter_num <<"    <A HREF=\"" 
                     << main_filename 
		     << "\">" << chapter_title << "</A> </TD> </TR>" << endl;
     else
       *contents_stream << "    <LI> <A HREF=\"" << main_filename
		     << "\">" << chapter_title << "</A>" << endl;     
}



void handlePart(  const Buffer_list& T) {
  

    if (!part_title.empty()) // end previous part 
    {
       *contents_stream << " </TABLE></TD><TD VALIGN=TOP> <TABLE>"  << endl;
       *contents_stream << "<!-- End of manual part -->"  << endl;
    }
    else
    {
       *contents_stream << "<TABLE WIDTH=100%><TD VALIGN=TOP><TABLE>" << endl;
    } 
    
    part_title = string( text_block_to_string( T));

    // add new part title to table of contents
    *contents_stream << "<!-- Start of new manual part -->"  << endl;
    *contents_stream << "<TR> <TH ALIGN=LEFT VALIGN=TOP><H3>" << part_title 
                     << "</H3></TH> </TR>" << endl; 
    if ( macroIsTrue( "\\lciIfNumberChaptersByPart") )
        chapter_num=0;
}



void handleBiblio(  const Buffer_list& T) {
    
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

// Opens a new classfile. Only a filename and a HTML formatted reference
// text are given.
void handleClassFile( const string& filename, 
		      const string& formatted_reference) {
    if ( next_class_link_last != 0) {
        *anchor_stream << "\"<!Next_class_link_" << next_class_link_last
		       << "!>\"    { fputs( \"<HR> <B>Next:</B> "
		       << convert_to_C_printable( formatted_reference) 
		       << "\\n\", stdout);}" << endl;
    }
    class_filename = filename;
    class_stream = open_file_for_write( tmp_path + class_filename);
    open_html( *class_stream);
    // Make a hyperlink in the chapter to the class file.
    if ( main_stream != &cout) {
         if (macroIsTrue("\\lciIfMultipleParts"))
           *main_stream  << "<TR> <TD> <UL><LI>" << formatted_reference
		      << ".</UL> </TD> </TR>\n" << endl;
         else
           *main_stream  << "<UL><LI>\n" << formatted_reference
		      << ".</UL>\n" << endl;
    }

    // table of contents
    if ( macroIsTrue( "\\lciIfHtmlClassToc"))
        if (macroIsTrue("\\lciIfMultipleParts"))
            *contents_stream << "<TR> <TD> <UL><LI> " << formatted_reference
			 << "</UL> </TD> </TR>" << endl;
        else
            *contents_stream << "          <UL><LI> " << formatted_reference
			 << "</UL> " << endl;
    current_ostream  = class_stream;
    current_filename = class_filename;
    insertInternalGlobalMacro( "\\lciOutputFilename", current_filename);
}

void handleClassFileEnd( void) {
    if ( class_stream != 0) {
	// implements the link from one class to the next class
	++next_class_link_counter;
	next_class_link_last = next_class_link_counter;
	*class_stream << "<!Next_class_link_" << next_class_link_counter 
		      << "!>";
        close_html( *class_stream);
	assert_file_write( *class_stream, class_filename);
	delete   class_stream;
	class_filename = string();
	class_stream = 0;
    }
    current_ostream  = main_stream;
    current_filename = main_filename;
    insertInternalGlobalMacro( "\\lciOutputFilename", current_filename);
}

void handleClassEnvironment() {
    string ref_scope_name = macroX("\\ccPureRefScope");
//    template_class_name = ref_scope_name + macroX("\\ccPureClassTemplateName");
    template_class_name = macroX("\\ccPureClassTemplateName");
    string formatted_template_class_name = 
	convert_C_to_html( template_class_name);
    class_name = macroX( "\\ccPureClassName");

    if ( macroIsTrue( "\\lciIfHtmlClassNotInline") &&
	 macroIsTrue( "\\lciIfHtmlClassFile")) {
	// Start a new class file.
	string filename = replace_template_braces_and_colons(
                            remove_font_commands( class_name))
	                + macroX( "\\lciHtmlSuffix");
        filename = replace_asterisks(filename);
	string contents( " Class declaration of ");
	contents += string("<A HREF=\"") + filename + "\">"
	            + formatted_template_class_name + "</A> ";
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

void handleHtmlClassFile( const string& filename, const Buffer_list& T) {
    char*  s = text_block_to_string( T);
    string fixed_filename = replace_template_braces_and_colons(filename);
    fixed_filename = replace_asterisks(fixed_filename);
    string p = string("<A HREF=\"") + fixed_filename + "\">" + s + "</A>";
    handleClassFile( fixed_filename, p);
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
    remove_separator( macroname);
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
    remove_separator( macroname);
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
    remove_separator( macroname);
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
    remove_separator( macroname);
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
    for ( size_t i = 0; i < s.size(); i++)
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
    else
	cerr << endl << "*** Error: " << convert_quoted_string(param[0]) <<'.';
    printErrorMessage( UserDefinedError);
    return string();
}

string
html_message( const string&, string param[], size_t n, size_t opt) {
    if ( n != 1 && opt != 0)
	printErrorMessage( MacroParamNumberError);
    else if ( ! quiet_switch)
	cerr << convert_quoted_string( param[0]);
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
    crop_string( param[0]);
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
    crop_string( param[0]);
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

string
line_number( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);
    if ( in_file)
	return int_to_string( in_file->line());
    return "0";
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

    insertInternalGlobalMacro( "\\lciStyle", cc_style, 1);
    insertInternalGlobalMacro( "\\lciCCParameter", cc_parameter, 0);

    insertInternalGlobalMacro( "\\lciTwoColumnLayout",  two_column_layout,  1);
    insertInternalGlobalMacro( "\\lciThreeColumnLayout",three_column_layout,2);

    insertInternalGlobalMacro( "\\lciHtmlIndex", html_index, 2);
    insertInternalGlobalMacro( "\\lciHtmlIndexC", html_index_C, 1);
    insertInternalGlobalMacro( "\\lciHtmlCrossLink", cross_link, 1);
    insertInternalGlobalMacro( "\\lciHtmlClassFileEnd", html_class_file_end,0);

    insertInternalGlobalMacro( "\\lciPopOutput",  pop_output,  0);
    insertInternalGlobalMacro( "\\lciPushOutput", push_output, 1);
    insertInternalGlobalMacro( "\\lciOpenBibliography",  open_biblio, 0);
    insertInternalGlobalMacro( "\\lciCloseBibliography", close_biblio, 0);

    insertInternalGlobalMacro( "\\lciLineNumber",  line_number,  0);
}

// EOF //






