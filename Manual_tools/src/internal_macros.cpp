/**************************************************************************

  internal_macros.cpp
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Internal macro definitions.
  System    : bison, flex, C++ (g++)  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$

**************************************************************************/


#include <internal_macros.h>

#include <mstring.h>
#include <ctype.h>
#include <string_conversion.h>
#include <macro_dictionary.h>
#include <cpp_formatting.h>
#include <input.h>
#include <lexer.h>
#include <error.h>
#include <config.h>
#include <output.h>
#include <sstream>
#include <cassert>

#include <regex_dictionary.h>

#include <map>


typedef std::map< string, string, Case_insensitive_string_greater_than > String_map;
hash_map< string, String_map >  sorted_map_table;

hash_map< string, hash_map< string, string > > hash_map_table;

// ======================================================================
// External visible functions
// ======================================================================


static int index_anchor_counter = 0;

static int cross_link_anchor_counter = 0;

string handleHtmlCrossLink( string key, bool tmpl_class) {
    crop_string( key);
    if ( key.empty()) {
      printErrorMessage( EmptyCrossLinkError);
      exit( 1);
    }

    char *tmp_namex = convert_fontified_ascii_to_html( key);
    string tmp_name( tmp_namex );
    delete[] tmp_namex;

    ostringstream replacement_text;
    replacement_text << REPLACE_WITH_CURRENT_PATH_TOKEN
                     << current_basename
                     << "#Cross_link_anchor_" << cross_link_anchor_counter;


    string anchor_key = "c ";
    do {
      *anchor_stream << anchor_key << tmp_name << '\t';
      wrap_anchor( replacement_text.str(), tmp_name, *anchor_stream );
      *anchor_stream << endl;

      string::size_type t = tmp_name.rfind( "&lt;" );
      string::size_type t2 = tmp_name.rfind( "::" );
      if( t != string::npos &&
          (t2 == string::npos || t > t2) )
      {
        string tmp_name2 = tmp_name.substr( 0, t );
        *anchor_stream << anchor_key << tmp_name2 << '\t';
        wrap_anchor( replacement_text.str(), tmp_name2, *anchor_stream );
        *anchor_stream << endl;
      }
      t2 = tmp_name.find( "::" );
      if( t2 != string::npos ) {
        const string prefix = tmp_name.substr( 0, t2 );
        if( prefix != "CGAL" )
          anchor_key = "localc ";
        tmp_name = tmp_name.substr( t2 + 2 );
      } else
        break;
    } while( tmp_name.length() > 2 );

    return string("\n<A NAME=\"Cross_link_anchor_")
      + int_to_string( cross_link_anchor_counter++) + "\"></A>\n";
}


// ========================================================================
// Index
// ========================================================================


string correct_name (string s) {
// corrects index entries for HTML (with <> and </> commands) so that
// the result is a correct command in HTML.
// Example: ``<I> index'' is changed into ``<I> index </I>''.

  string tab[20];
  string search_item;
  string correct_s=s;
  int tab_item=0;
  int k=1;
  bool found;
  size_t i = 0;

  while (i< s.size()) {
     if (s[i]=='<') {
        if (s[i+1]=='/' && (i+1<s.size())) {
             tab_item++;
             while (i< s.size() && s[i] != '>') {
               tab[tab_item]+=s[i];
               i++;
             }
             tab[tab_item].replace(0,2,"");
             search_item = tab[tab_item];
             tab[tab_item--]="";
             found = false;
             while ( k<=tab_item && !(found) ) {
                if (tab[k] == search_item) {
                   tab[k]="";
                   found = true;
                }
                k++;
             }
             if (!found)  correct_s = "<" + search_item + ">" + correct_s;
        }
        else {
          tab_item++;
          while (i< s.size() && s[i] != '>') {
             tab[tab_item]+=s[i];
             i++;
          }
          tab[tab_item].replace(0,1,"");
        }
     }
     i++;
   }
   for (int j=1;j<=tab_item;j++) {
      if (!(tab[j]==""))
          correct_s += "</" + tab[j] + ">";
   }
   return correct_s;

}


void mainTextParse(string s,
                   string& index_name,
                   string& modifier) {
// This function separates the normal entry and the modified entry from s.

   remove_leading_spaces(s);
   size_t i = 0;
   while (i< s.size() && s[i] != ',')  index_name+=s[i++];
   if (s[i] ==',') {
      s.replace(0,++i,"");
      remove_leading_spaces(s);
      i=0;
      while (i < s.size()){
        modifier+=s[i];
        i++;
      }
   }
   index_name = correct_name(index_name);
   modifier = correct_name(modifier);
}


ostream* temp_main_stream  = 0;
ostream* temp_current_stream = 0;


static string temp_main_filename;
static string temp_current_filename;

void OpenFileforIndex() {
  string new_main_filename;
  temp_main_stream = main_stream;
  temp_current_stream = current_ostream;
  temp_main_filename = main_filename;
  temp_current_filename = current_filename;

  string WhichItem = macroX("\\lciWhichItem");
  if (WhichItem=="MainItem" || WhichItem=="index") {
     new_main_filename = macroX("\\lciMainItemFile");
  }
  else {
     if (WhichItem=="SubItem")
        new_main_filename = macroX("\\lciSubItemFile");
     else {
        if (WhichItem=="SubSubItem") {
          new_main_filename = macroX("\\lciSubSubItemFile");
        }
        else {
           cerr << " internal_macros.cpp : warning: unknown variable WhichItem"
                << endl;
           exit(1);
        }
     }
  }
  main_filename = new_main_filename;
  main_stream = open_file_for_write( tmp_path + main_filename);
  current_ostream  = main_stream;
  current_filename = main_filename;
  insertInternalGlobalMacro( "\\lciOutputFilename", current_filename);
  insertInternalGlobalMacro( "\\lciMainFilename",   main_filename);
}

void CloseFileforIndex() {
  assert_file_write( *main_stream, main_filename);
  delete   main_stream;
  main_stream = temp_main_stream;
  current_ostream = temp_current_stream;
  main_filename = temp_main_filename;
  current_filename = temp_current_filename;
  insertInternalGlobalMacro( "\\lciOutputFilename", current_filename);
  insertInternalGlobalMacro( "\\lciMainFilename",   main_filename);
}



string name_for_ordering(string s) {
// index entry without HTML commands

   string ord_s="";
   size_t i = 0;
   while (i< s.size()) {
     if (s[i]=='<') {
        while (i< s.size() && s[i] != '>')
             i++;
     }
     else {
        ord_s+=s[i];
     }
     i++;
   }
   return ord_s;
}

// hash map of modified entries
typedef hash_map< string, int> Modifier_map;
typedef Modifier_map::iterator Modifier_map_iterator;
Modifier_map modifier_map;


void handleIndex() {
    string category;
    string index_name, sub_index_name, sub_sub_index_name;
    string modifier="", sub_modifier="", sub_sub_modifier="";
    string ord_modifier="", ord_sub_modifier="", ord_sub_sub_modifier="";
    string ord_index_name="", ord_sub_index_name="", ord_sub_sub_index_name="";
    string s;

    string main_item = "";
    string WhichItem = macroX("\\lciWhichItem");
    if (WhichItem=="MainItem" || WhichItem=="SubItem" ||
          WhichItem=="SubSubItem") {
       open_file_to_string(tmp_path + macroX("\\lciMainItemFile"),s);
       mainTextParse(s, index_name, modifier);
       s="";
       ord_index_name = name_for_ordering(index_name);
       ord_modifier = name_for_ordering(modifier);
    }
    if (WhichItem=="SubItem" || WhichItem=="SubSubItem") {
       open_file_to_string(tmp_path + macroX("\\lciSubItemFile"),s);
       mainTextParse(s, sub_index_name, sub_modifier);
       s="";
       ord_sub_index_name = name_for_ordering(sub_index_name);
       ord_sub_modifier = name_for_ordering(sub_modifier);
    }
   if (WhichItem=="SubSubItem") {
       open_file_to_string(tmp_path + macroX("\\lciSubSubItemFile"),s);
       mainTextParse(s, sub_sub_index_name, sub_sub_modifier);
       s="";
       ord_sub_sub_index_name = name_for_ordering(sub_sub_index_name);
       ord_sub_sub_modifier = name_for_ordering(sub_sub_modifier);
   }

   string sub_item="";
   string sub_sub_item="";

   if (sub_index_name!="") sub_index_name="??? "+sub_index_name;
   if (sub_sub_index_name!="") sub_sub_index_name="??? "+sub_sub_index_name;

   if (!(sub_index_name=="")) {
       if (sub_modifier=="")
          sub_item = ord_sub_index_name+"@ "+sub_index_name;
       else
         sub_item = ord_sub_index_name+" "+ord_sub_modifier +"@ "
                    + sub_index_name+", "+sub_modifier;
   }
   if (!(sub_sub_index_name=="")) {
       if (sub_sub_modifier=="")
          sub_sub_item = ord_sub_sub_index_name+"@ "+sub_sub_index_name;
       else
         sub_sub_item = ord_sub_sub_index_name +" " +ord_sub_sub_modifier
                        + "@ "+sub_sub_index_name+", " +sub_sub_modifier;
   }


   if (modifier=="") {
         handleIndex2( ord_index_name+"@ ??? "+ index_name, sub_item,
                       sub_sub_item, 0);
   }
   else {
       int number;
       Modifier_map_iterator mod_it = modifier_map.find( modifier+index_name);
       if ( mod_it != modifier_map.end()) {
            handleIndex2(ord_index_name+" "+ord_modifier+"@ ??? "+index_name+
                         ",  "+ modifier+"<A NAME=\""
                         + int_to_string( mod_it->second)
                         +"\"></A>",sub_item,sub_sub_item,0);
       } else {
            if (ord_modifier=="2D" || ord_modifier=="3D" ||
                 ord_modifier=="dD") {
                handleIndex2(ord_index_name+" "+ord_modifier+"@ ??? "
                        +index_name+ ",  "+ modifier,sub_item,sub_sub_item,0);
            } else {
                handleIndex2(ord_index_name+" "+ord_modifier+"@ ??? " +
                   index_name+ ",  " + modifier,sub_item,sub_sub_item,1);
                handleIndex2(ord_modifier+ " "+ord_index_name+" see "
                    +ord_index_name+", "+ord_modifier+"@"+modifier
                    + "  "+index_name+", <I> see </I> ??? "+index_name
                    + ", "+modifier ,"","",2);
                modifier_map[ modifier+index_name] = HREF_counter-4;
            }
        }
   }
}


void handleIndex2(string main_item, string sub_item, string sub_sub_item,
                  int modifier) {
    remove_separator( main_item);
    remove_separator( sub_item);
    remove_separator( sub_sub_item);
   if (modifier==1)
     *index_stream << "\\indexentry{"
		        << main_item<< "<A NAME=\""<<HREF_counter
                        <<"\"></A>";
   else {
     *index_stream << "\\indexentry{"
		        << main_item;
   }

   if (sub_item!=""){
        *index_stream <<"! " << sub_item;
        if (sub_sub_item!="") *index_stream <<"! " << sub_sub_item;
   }
   if ( !macroIsTrue("\\lciIfSeeAlso")) {
      switch (modifier) {
        case 0 :
             *HREF_stream << HREF_counter << " HREF=\"" << current_filename
                          << "#Index_anchor_" << index_anchor_counter << "\""
                          << endl;
             break;
        case 1 :
             *HREF_stream << HREF_counter <<" HREF=\"" << current_filename
                          << "#Index_anchor_" << index_anchor_counter << "\""
                          << endl;
             break;
        case 2 :
             *HREF_stream << HREF_counter << " HREF=\"" << "#"
                          << HREF_counter-2 << "\"" << endl;
             break;
      }
   }


   *index_stream << "}{"<< HREF_counter << "}" << endl;

   HREF_counter+=2;
   *current_ostream << "\n<A NAME=\"Index_anchor_"
	 << int_to_string( index_anchor_counter++)
	 << "\"></A> \n";
}


void TraitsClassTextParse(string s, string p) {

   string index_name, index_class_name;
   remove_leading_spaces(p);
   remove_trailing_spaces(p);
   index_class_name="see also"+p+"@<I> see also </I> ??? " + p;
   remove_leading_spaces(s);
   size_t i = 0;
   while (i< s.size()) {
      if (s[i] ==';') {
         s.replace(0,++i,"");
         remove_leading_spaces(s);
         i=0;
         remove_trailing_spaces(index_name);
         if ( ! macroIsTrue("\\lciIfPackage"))
              index_name=name_for_ordering(index_name)+"@ ??? <I>"+
                         name_for_ordering(index_name)+"</I>";
         else   index_name = name_for_ordering(index_name)+"@ ??? "+
                             name_for_ordering(index_name);
         handleIndex2(index_name,"traits class@ ??? traits class",
                      index_class_name,0);
         index_name="";
      }
      else {
         index_name+=s[i++];
      }
   }
   remove_leading_spaces(index_name);
   remove_trailing_spaces(index_name);
   if ( ! macroIsTrue("\\lciIfPackage"))
       index_name=name_for_ordering(index_name)+"@ ??? <I>"+
                  name_for_ordering(index_name)+"</I>";
   else index_name = name_for_ordering(index_name)+"@ ??? "+
                             name_for_ordering(index_name);
   handleIndex2( index_name, "traits class@ ??? traits class",
                 index_class_name, 0);
}



void handleIndexTraitsClass() {
   string s1, s2;
   open_file_to_string(tmp_path + macroX("\\lciMainItemFile"),s1);
   open_file_to_string(tmp_path + macroX("\\lciSubItemFile"),s2);
   TraitsClassTextParse(s1, s2);
}


void handleIndexRefName() {
   string s1, s2, s3;
   string index_name;
   open_file_to_string(tmp_path + macroX("\\lciMainItemFile"),s1);
   open_file_to_string(tmp_path + macroX("\\lciSubItemFile"),s2);
   open_file_to_string(tmp_path + macroX("\\lciSubSubItemFile"),s3);
   // s1 = local scope (may be empty)
   // s2 = classname (functioname etc.)
   // s3 = Ref category
   if ((s3=="Class") || (s3=="FunctionObjectClass")) {
      index_name = name_for_ordering(s2)+"@ ??? "+s2;
      handleIndex2(index_name,"","",0);
   }
   else {
     if ((s3=="Concept") || (s3=="FunctionObjectConcept")) {
        index_name = name_for_ordering(s2)+"@ ??? "+name_for_ordering(s2);
        handleIndex2(index_name,"","",0);
     }
     else {
       if (s1=="") {
          index_name = name_for_ordering(s2)+"@ ??? "+s2;
          handleIndex2(index_name,"","",0);

       }
       else {
         for ( size_t i = 0; i < s1.size(); ++i) {
	    if ( s1[i]==':' ) {
              if ( (i+1)<s1.size() && s1[i+1]==':') {
	        s1.replace(i,s1.size()-i,"");
                break;
              }
            }
	  }
          index_name = name_for_ordering(s1)+"@ ??? "+s1;
          string index_name2 =  name_for_ordering(s2)+"@ ??? "+s2;
          if ((s3=="Function"))
             handleIndex2(index_name2,index_name,"",0);
          else
             handleIndex2(index_name,index_name2,"",0);
       }
     }
   }
}


// Opens a new classfile of name 'filename'
void handleClassFile( string filename) {
    string filepath = macroX( "\\lciInputPath");
    if ( filepath[0] == '/' ) {
        printErrorMessage( ClassAbsolutePathError);
        filepath = main_filepath;
    }
    filename = filepath + filename;
    class_filename   = filename;
    current_filename = class_filename;
    current_basename = basename_string(class_filename);
    current_rootname = rootname_string(current_basename);
    current_filepath = path_string(class_filename);
    current_uppath   = uppath_string( current_filepath);
    insertInternalGlobalMacro( "\\lciOutputFilename", current_filename);
    insertInternalGlobalMacro( "\\lciOutputBasename", current_basename);
    insertInternalGlobalMacro( "\\lciOutputRootname", current_rootname);
    insertInternalGlobalMacro( "\\lciOutputPath",     current_filepath);
    insertInternalGlobalMacro( "\\lciOutputUppath",   current_uppath);

    if ( current_filepath == "") {
        anchor_stream = global_anchor_stream;
    } else if ( current_filepath == main_filepath) {
        anchor_stream = main_anchor_stream;
    } else {
        anchor_stream = open_file_for_append_with_path(
            tmp_path + current_filepath + macroX( "\\lciAnchorFilename"));
    }
    if(class_stream != 0) {
        printErrorMessage( RefPageNotClosedError);
        exit(1);
    }
    class_stream = open_file_for_write_with_path( tmp_path + class_filename);
    current_ostream  = class_stream;

}

void handleClassFileEnd( void) {
    assert_file_write( *class_stream, class_filename);
    delete   class_stream;
    class_filename = string();
    class_stream = 0;

    if ( anchor_stream != 0 && anchor_stream != main_anchor_stream
         && anchor_stream != global_anchor_stream) {
        assert_file_write( *anchor_stream, macroX( "\\lciAnchorFilename"));
        delete anchor_stream;
    }
    anchor_stream = main_anchor_stream;
    current_ostream  = main_stream;
    current_filename = main_filename;
    current_basename = main_basename;
    current_rootname = main_rootname;
    current_filepath = main_filepath;
    current_uppath   = main_uppath;
    insertInternalGlobalMacro( "\\lciOutputFilename", current_filename);
    insertInternalGlobalMacro( "\\lciOutputBasename", current_basename);
    insertInternalGlobalMacro( "\\lciOutputRootname", current_rootname);
    insertInternalGlobalMacro( "\\lciOutputPath",     current_filepath);
    insertInternalGlobalMacro( "\\lciOutputUppath",   current_uppath);
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
div_to( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    string macroname( param[0]);
    crop_string( macroname);
    remove_separator( macroname);
    int prod = atoi( expandFirstMacro( macroname).c_str())
             / atoi(param[1].c_str());
    insertMacro( macroname, in_string->name(), in_string->line(),
		 int_to_string( prod));
    return string();
}

string
mod_to( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    string macroname( param[0]);
    crop_string( macroname);
    remove_separator( macroname);
    int prod = atoi( expandFirstMacro( macroname).c_str())
             % atoi(param[1].c_str());
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

string
global_div_to( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    string macroname( param[0]);
    crop_string( macroname);
    remove_separator( macroname);
    int prod = atoi( expandFirstMacro( macroname).c_str())
             / atoi(param[1].c_str());
    insertGlobalMacro( macroname, in_string->name(), in_string->line(),
		       int_to_string( prod));
    return string();
}

string
global_mod_to( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    string macroname( param[0]);
    crop_string( macroname);
    remove_separator( macroname);
    int prod = atoi( expandFirstMacro( macroname).c_str())
             % atoi(param[1].c_str());
    insertGlobalMacro( macroname, in_string->name(), in_string->line(),
		       int_to_string( prod));
    return string();
}


string
add( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    int sum = atoi( param[0].c_str() )
            + atoi( param[1].c_str() );
    return int_to_string( sum );
}

string
expand_n( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    int N      = atoi( param[0].c_str() );
    string cmd = param[1];
    string retval = string();
    //std::cout << std::endl << "expand_n " << N << " " << cmd << " -> " << std::endl;
    for( ; N > 0; --N )
        retval += macroX( cmd );
    //std::cout << retval << std::endl;
    return retval;
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

string
int_to_roman( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);
    return int_to_roman_string( atoi( param[0].c_str()));
}

string
int_to_roman_upper( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);
    string s = int_to_roman_string( atoi( param[0].c_str()));
    for ( size_t i = 0; i < s.size(); i++)
	s[i] = toupper( s[i]);
    return s;
}

string
int_to_alpha( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);
    int i = atoi( param[0].c_str());
    if ( i < 1 || i > 27) {
	printErrorMessage( AlphaOutOfBoundsError);
        return string( "[Alpha digit out of bounds]");
    }
    char s[2] = "-";
    s[0] = 'a' + i - 1;
    return string( s);
}

string
int_to_alpha_upper( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);
    int i = atoi( param[0].c_str());
    if ( i < 1 || i > 27) {
	printErrorMessage( AlphaOutOfBoundsError);
        return string( "[Alpha digit out of bounds]");
    }
    char s[2] = "-";
    s[0] = 'A' + i - 1;
    return string( s);
}

string
substring( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 3, 0 );
    int begin = atoi( param[0].c_str() );
    int end   = atoi( param[1].c_str() );
    string s  = param[2];
    if( end <= 0 ) {
      end = s.length() + end;
      std::cout << "len: " << s.length() << " end: " << end << std::endl;
    } else
      end = begin + end;
    if( end < 0 )
      end = 0;
    return s.substr( begin, end );
}


string
remove_newlines( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0 );
    string s  = param[0];
    string::size_type pos = 0, nl;
    while( ( nl=s.find( '\n' ) ) != string::npos ) {
      s.replace( nl, 1, " "  );
      pos = nl;
    }
    return s;
}


// Some special characters hard to print otherwise
// ======================================================================

string
format_chapter_author( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);
    // count if none, one, or several \and's are used
    std::string::size_type idx = param[0].find("\\and");
    // make sure it was not a macro that just starts with '\and...'
    while ( idx != std::string::npos && idx + 4 < param[0].size() &&
            std::isalpha( param[0][idx+4]))
        idx = param[0].find("\\and", idx+4);
    if ( idx == std::string::npos) // no \and's used, do nothing
        return param[0];
    std::string::size_type idx2 = param[0].find("\\and", idx+1);
    // make sure it was not a macro that just starts with '\and...'
    while ( idx2 != std::string::npos && idx2 + 4 < param[0].size() &&
            std::isalpha( param[0][idx2+4]))
        idx2 = param[0].find("\\and", idx2+4);
    if ( idx2 == std::string::npos) { // one \and used, replace with " and "
        param[0].replace( idx, 4, " and ");
        return param[0];
    }
    // else replace each idx position with ", " and wait for last \and
    do {
        param[0].replace( idx, 4, ", ");
        // eliminate whitespaces left of the comma
        while ( idx > 0 && param[0][idx-1] <= ' ') {
            --idx;
            --idx2;
            param[0].replace( idx, 1, "");
        }
        idx = idx2 - 2;
        idx2 = param[0].find("\\and", idx+1);
        // make sure it was not a macro that just starts with '\and...'
        while ( idx2 != std::string::npos && idx2 + 4 < param[0].size() &&
                std::isalpha( param[0][idx2+4]))
            idx2 = param[0].find("\\and", idx2+4);
    } while (idx2 != std::string::npos);
    // finish with replacing last \and with ", and "
    param[0].replace( idx, 4, ", and ");
    // eliminate whitespaces left of the comma
    while ( idx > 0 && param[0][idx-1] <= ' ') {
        --idx;
        param[0].replace( idx, 1, "");
    }
    return param[0];
}


// Some special characters hard to print otherwise
// ======================================================================
string
backslash_char( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);
    return string( "\\lciRawOutputN{1}\\");
}

string
open_brace_char( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);
    return string( "\\lciRawOutputN{1}{");
}

string
close_brace_char( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);
    return string( "\\lciRawOutputN{1}}");
}



// Error and message output
// ======================================================================
string
html_error( const string&, string param[], size_t n, size_t opt) {
    cerr << endl;
    if ( n != 1 && opt != 0)
	printErrorMessage( MacroParamNumberError);
    else
	cerr << endl << "ERROR: " << convert_quoted_string(param[0]) <<'.';
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

string
html_dump( const string&, string param[], size_t n, size_t opt) {
    if ( n != 1 && opt != 0)
	printErrorMessage( MacroParamNumberError);
    else if ( ! quiet_switch)
	cerr << convert_quoted_string_seps( param[0]);
    return string();
}


// C++ Formatting
// ======================================================================
string
cc_style( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);  // uses cc_string as parameter, param[0] is a dummy
                         // to avoid scanning of spaces behind \ccc
    string html = convert_C_to_html( cc_string );
    std::stringstream s;
    s << "\\lciRawOutputN{" << html.length() << "}" << html;
    return s.str();
}

string
cc_parameter( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);  // uses cc_string as parameter
    return cc_string;
}

string
cache_class_name( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);
    template_class_name = macroX("\\ccPureClassTemplateName");
    class_name          = macroX( "\\ccPureClassName");
    return string();
}

string
store_file_name( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);
    cc_filename = replace_asterisks( replace_template_braces_and_colons(
                                         remove_font_commands( cc_string)));
    insertInternalGlobalMacro( "\\lciNewFilename", cc_filename);
    return string();
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
    NParamCheck( 1, 0);  // param[0] is index text
    crop_string( param[0]);
    handleIndex2( param[0] +"@ ??? "+ param[0],"","",0);
    return "";
}

string
html_index_C( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);  // uses cc_string as index item, param[0] is dummy
    handleIndex2( cc_string +"@ ??? "+ convert_C_to_html(cc_string),"","",0);
    return "";
}

string
cross_link( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);  // uses cc_string as crosslink key, param[0] is
                         // used as dummy to avoid parsing of spaces.
    return string( "\\lcRawHtml{") + handleHtmlCrossLink(cc_string,false)+ '}';
}

string
cross_link_template( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);  // uses cc_string as crosslink key, param[0] is
                         // used as dummy to avoid parsing of spaces.
    return string( "\\lcRawHtml{") + handleHtmlCrossLink(cc_string,true) + '}';
}

string
handle_include_only( const string& s, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0 );
    include_only( param[0] ); // only single filename supported currently
    return string();
}

string
handle_to_be_included( const string& s, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0 );
    if( is_to_be_included( param[0] ) )
      return string("\\lcTrue");
    else
      return string("\\lcFalse");
}

// another hack:
string
RCSdef( string input, bool date ) {
  input = input.substr( 1, input.length() -2 );
  string::size_type pos, pos2;
  pos = input.find_first_of( " " );
  pos = input.find_first_of( " ", pos +1 );
  if( date )  // date entry is next to revision entry
    pos = input.find_first_of( " ", pos +1 );
  pos2 = input.find_first_of( " ", pos +1 );
  return input.substr( pos + 1, pos2 - pos - 1);
}

string
handle_RCSdef( const string& s, string param[], size_t n, size_t opt ) {
  NParamCheck( 1, 0 );
  return RCSdef( param[0], false );
}

string
handle_RCSdefDate( const string& s, string param[], size_t n, size_t opt ) {
  NParamCheck( 1, 0 );
  return RCSdef( param[0], true );
}


// Chapter File Handling
// ======================================================================

string                // macro
handle_chapter( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);
    push_current_output();
    //string chapter_title( param[0]);
    //crop_string( chapter_title);
    //remove_separator( chapter_title);

    string new_main_filename;
    string new_main_filepath = macroX( "\\lciInputPath");
    if ( new_main_filepath[0] == '/' ) {
        printErrorMessage( ChapterAbsolutePathError);
        new_main_filename =  macroX( "\\lciChapterPrefix")
                           + macroX( "\\lciInputRootname")
	                   + macroX( "\\lciHtmlSuffix");
        new_main_filepath = string();
    } else {
        new_main_filename = new_main_filepath
                          + macroX( "\\lciChapterPrefix")
                          + macroX( "\\lciInputRootname")
                          + macroX( "\\lciHtmlSuffix");
    }
    if ( new_main_filename == main_filename) {
        printErrorMessage( ChapterStructureError);
	return string();
    }
    if ( main_stream != &cout && main_stream != pre_stream) {
        assert( main_stream != 0 );
	assert_file_write( *main_stream, main_filename);
	delete   main_stream;
	main_stream = 0;
    }
    if ( anchor_stream != 0 && anchor_stream != main_anchor_stream
         && anchor_stream != global_anchor_stream) {
        assert_file_write( *anchor_stream, macroX( "\\lciAnchorFilename"));
        delete anchor_stream;
        anchor_stream = 0;
    }
    if ( main_anchor_stream != 0
         && main_anchor_stream != global_anchor_stream) {
        assert_file_write( *main_anchor_stream,
                           macroX( "\\lciAnchorFilename"));
        delete main_anchor_stream;
        main_anchor_stream = 0;
    }
    main_filename = new_main_filename;
    main_basename = basename_string( main_filename);
    main_rootname = rootname_string( main_basename);
    main_filepath = new_main_filepath;
    main_uppath   = uppath_string( main_filepath);
    main_stream = open_file_for_write_with_path( tmp_path + main_filename);
    current_ostream  = main_stream;
    current_filename = main_filename;
    current_basename = main_basename;
    current_rootname = main_rootname;
    current_filepath = main_filepath;
    current_uppath   = main_uppath;
    insertInternalGlobalMacro( "\\lciOutputFilename", current_filename);
    insertInternalGlobalMacro( "\\lciOutputRootname", current_rootname);
    insertInternalGlobalMacro( "\\lciOutputBasename", current_basename);
    insertInternalGlobalMacro( "\\lciOutputPath",     current_filepath);
    insertInternalGlobalMacro( "\\lciOutputUppath",   current_uppath);
    insertInternalGlobalMacro( "\\lciMainFilename",   main_filename);
    insertInternalGlobalMacro( "\\lciMainBasename",   main_basename);
    insertInternalGlobalMacro( "\\lciMainRootname",   main_rootname);
    insertInternalGlobalMacro( "\\lciMainPath",       main_filepath);
    insertInternalGlobalMacro( "\\lciMainUppath",     main_uppath);
    if ( current_filepath == "") {
        anchor_stream = global_anchor_stream;
    } else {
        anchor_stream = open_file_for_append( tmp_path + current_filepath +
                                              macroX( "\\lciAnchorFilename"));
    }
    main_anchor_stream = anchor_stream;
    return string();
}

// Reference page and class file handling
// ======================================================================


// File management
// ======================================================================

string
if_file_exists( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);
    string name( param[0]);
    remove_separator( name);
    if ( find_filename_with_suffix_w_input_dirs( name) != string(""))
	return string("\\lcTrue");
    return string("\\lcFalse");
}

string
copy_file( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 2, 0);
    string source( param[0]);
    string target( param[1]);
    remove_separator( source);
    remove_separator( target);
    istream* in  = open_file_for_read_w_input_dirs( source);
    ostream* out = open_file_for_write_with_path( target);
    if ( in != 0 && out != 0) {
        char c; // yes, probably slow, but then not that often used anyway
        while( in->get(c)) {
            out->put(c);
        }
    }
    delete in;
    delete out;
    return string("");
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
open_tmp_file( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);
    remove_separator( param[0]);
    push_current_output_w_filename( tmp_path + expandFirstMacro(param[0]) );
    return string();
}

string
open_file( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);
    remove_separator( param[0]);
    push_current_output_w_filename( expandFirstMacro(param[0]) );
    return string();
}

string
close_file( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);
    delete current_ostream;
    pop_current_output();
    return string();
}

string
open_biblio( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);
    push_current_output();
    current_filename = tmp_path + macroX("\\lciBibFilename");
    current_ostream = open_file_for_write( current_filename);
    return string();
}

string
close_biblio( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);
    assert_file_write( *current_ostream, macroX( "\\lciBibFilename"));
    delete current_ostream;
    pop_current_output();
    return string();
}

string
open_reference_file( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 0, 0);
    handleClassFile( cc_filename);
    return string();
}

string
close_reference_file( const string&, string param[], size_t n, size_t opt) {
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

string
to_html_width( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);
    string width=param[0];
    string::size_type pos = width.find( "\\textwidth" );
    if( pos != string::npos ) {
        width = width.substr( 0, pos );
        return float_to_string( atof( width.c_str() ) * 100.0 ) + "\\%";
    }
    return "0";
}



string
sorted_map_add_to( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 3, 0);
  string name  = param[0];
  string key   = param[1];
  string value = param[2];

  if( key[0] == '\\' )
    key = macroX( key );
  if( value[0] == '\\' )
    value = macroX( value );

  sorted_map_table[name][ key ] = value;
  return string();
}

string
sorted_map_foreach( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 2, 0);
  string name  = param[0];
  string func  = param[1];

  if( func[0] != '\\' )
    return string("!! Error: argument to \\lciSortedMapForeach has to be a function");

  const String_map& m = sorted_map_table[name];
  string func_params[2];
  Macro_item item = fetchMacro( func );
  for( String_map::const_iterator it = m.begin(); it != m.end(); ++it ) {
    func_params[0] = it->first;
    func_params[1] = it->second;
    include_stack.push_string(
      func,
      expandMacro(func, item, func_params, 2, 0),
      0 );
  }

  return string();
}

string
sorted_map_clear( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 1, 0);
  string name  = param[0];
  sorted_map_table[name].clear();
  return string();
}

string
sorted_map_is_empty( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 1, 0);
  string name  = param[0];
  if( sorted_map_table[name].empty() )
    return string("\\lcTrue");
  else
    return string("\\lcFalse");
}

/* hash maps */

string
hash_map_insert( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 3, 0);
  string name  = param[0];
  string key   = param[1];
  string value = param[2];

  crop_string( key );
  crop_string( value );

  if( key[0] == '\\' )
    key = macroX( key );
  if( value[0] == '\\' )
    value = macroX( value );
  crop_string( name );
  crop_string( key );
  crop_string( value );

  hash_map_table[name][ key ] = value;
  return string();
}

string
hash_map_is_defined( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 2, 0);
  string name  = param[0];
  string key   = param[1];
  if( hash_map_table[name].find( key ) == hash_map_table[name].end() )
    return string("\\lcFalse");
  else
    return string("\\lcTrue");
}

string
hash_map_get( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 2, 0);
  string name  = param[0];
  string key   = param[1];
  return hash_map_table[name][key];
}

string
hash_map_clear( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 1, 0);
  hash_map_table[ param[0] ].clear();
  return string();
}

/* regular expressions */


string
handle_regex_register( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 2, 0);
  regex_register( param[0], param[1] );
  return string();
}

string
handle_regex_does_match( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 2, 0);
  if( regex_does_match( param[0], param[1] ) )
    return string("\\lcTrue");
  else
    return string("\\lcFalse");
}

string
handle_regex_store_result( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 2, 0);
  unsigned int index = atoi( param[1].c_str() );
  string macroname = param[0];
  string result = regex_get_submatch(index);
  //std::cerr << "!! Warning. store regex result [" << macroX(macroname) << "] [" << result << "]" << std::endl;
  insertInternalGlobalMacro( macroname, macroX(macroname) + result, 0);
  return string();
}


string
handle_store_savebox( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 1, 0);
  insertInternalGlobalMacro( param[0], savebox_stream->str(), 0);
  // clear the string stream:
  savebox_stream->str("");
  return string();
}

string
handle_savestream_open( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 1, 0);
  savestream_open( param[0] );
  return string();
}

string
handle_savestream_use( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 2, 0);
  return savestream_use( param[0], param[1] );
}

string
handle_savestream_close( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 1, 0);
  savestream_close( param[0] );
  return string();
}

string
handle_minitoc_open( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 0, 0);
  minitoc_open();
  return string();
}

string
handle_minitoc_close( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 0, 0);
  minitoc_close();
  return string();
}

string
handle_shell_command( const string&, string param[], size_t n, size_t opt) {
  NParamCheck( 1, 0);
  string cmd = param[0];
  std::stringstream out, err, exitcode;
  int retval = execute_shell_command( cmd, out, err );
  exitcode << retval;
  insertInternalGlobalMacro( "\\lciShellExitCode", exitcode.str() );
  insertInternalGlobalMacro( "\\lciShellStdout", out.str() );
  insertInternalGlobalMacro( "\\lciShellStderr", err.str() );
  if( retval != 0 ) {
    return string("\\lciError{error while executing \"" + cmd + "\". error log:" + err.str());
  } else {
    return string();
  }
}


// Initialize
// ======================================================================
void init_internal_macros() {
    insertInternalGlobalMacro( "{", "\\begingroup");
    insertInternalGlobalMacro( "}", "\\endgroup");

    insertInternalGlobalMacro( "$$", "\\lciTexToggleDisplayMathMode");

    insertInternalGlobalMacro( "\\newcommand", newcommand, 2);
    insertInternalGlobalMacro( "\\newcommand@mom", newcommand_opt);
    insertInternalGlobalMacro( "\\def", newcommand, 2);
    insertInternalGlobalMacro( "\\gdef", global_newcommand, 2);
    insertInternalGlobalMacro( "\\lciRCSdef", handle_RCSdef, 1);
    insertInternalGlobalMacro( "\\lciRCSdefDate", handle_RCSdefDate, 1);
    insertInternalGlobalMacro( "\\lciUndef", undef, 1);
    insertInternalGlobalMacro( "\\lciBeginScope", begin_scope, 0);
    insertInternalGlobalMacro( "\\lciEndScope", end_scope, 0);
    insertInternalGlobalMacro( "\\lciIfDefined", if_defined, 1);
    insertInternalGlobalMacro( "\\lciIfEqual", if_equal, 2);
    insertInternalGlobalMacro( "\\lciIfLess", if_less, 2);
    insertInternalGlobalMacro( "\\lciIfLessOrEqual", if_less_or_equal, 2);
    insertInternalGlobalMacro( "\\lciAddTo",  add_to, 2);
    insertInternalGlobalMacro( "\\lciAdd",  add,    2);
    insertInternalGlobalMacro( "\\lciMultTo", mult_to, 2);
    insertInternalGlobalMacro( "\\lciDivTo",  div_to, 2);
    insertInternalGlobalMacro( "\\lciModTo",  mod_to, 2);
    insertInternalGlobalMacro( "\\lciGlobalAddTo",  global_add_to, 2);
    insertInternalGlobalMacro( "\\lciGlobalMultTo", global_mult_to, 2);
    insertInternalGlobalMacro( "\\lciGlobalDivTo",  global_div_to, 2);
    insertInternalGlobalMacro( "\\lciGlobalModTo",  global_mod_to, 2);

    insertInternalGlobalMacro( "\\lciExpandN",      expand_n, 2);

    insertInternalGlobalMacro( "\\lciToUpper",      string_to_upper, 1);
    insertInternalGlobalMacro( "\\lciToRoman",      int_to_roman, 1);
    insertInternalGlobalMacro( "\\lciToRomanUpper", int_to_roman_upper, 1);
    insertInternalGlobalMacro( "\\lciToAlpha",      int_to_alpha, 1);
    insertInternalGlobalMacro( "\\lciToAlphaUpper", int_to_alpha_upper, 1);
    insertInternalGlobalMacro( "\\lciSubstring",    substring, 3);
    insertInternalGlobalMacro( "\\lciRemoveNewlinesX",remove_newlines, 1);

    insertInternalGlobalMacro( "\\lciFormatChapterAuthor",
                                                   format_chapter_author, 1);

    insertInternalGlobalMacro( "\\lciBackslash",   backslash_char, 0);
    insertInternalGlobalMacro( "\\lciOpenBrace",   open_brace_char, 0);
    insertInternalGlobalMacro( "\\lciCloseBrace",  close_brace_char, 0);

    insertInternalGlobalMacro( "\\lciError",   html_error, 1);
    insertInternalGlobalMacro( "\\lciMessage", html_message, 1);
    insertInternalGlobalMacro( "\\lciDump",    html_dump, 1);

    insertInternalGlobalMacro( "\\lciCCStyle", cc_style, 1);
    insertInternalGlobalMacro( "\\lciCCParameter", cc_parameter, 0);
    insertInternalGlobalMacro( "\\lciCacheClassName", cache_class_name, 0);
    insertInternalGlobalMacro( "\\lciStoreFileName", store_file_name, 0);

    insertInternalGlobalMacro( "\\lciTwoColumnLayout",  two_column_layout,  1);
    insertInternalGlobalMacro( "\\lciThreeColumnLayout",three_column_layout,2);

    insertInternalGlobalMacro( "\\lciHtmlIndex", html_index, 1);
    insertInternalGlobalMacro( "\\lciHtmlIndexC", html_index_C, 1);
    insertInternalGlobalMacro( "\\lciHtmlCrossLink", cross_link, 1);
    insertInternalGlobalMacro( "\\lciHtmlCrossLinkTemplate",
                                                       cross_link_template, 1);
    insertInternalGlobalMacro( "\\lciOpenReferenceFile",
                                                       open_reference_file,0);
    insertInternalGlobalMacro( "\\lciCloseReferenceFile",
                                                       close_reference_file,0);

    insertInternalGlobalMacro( "\\lciSortedMapAddTo", sorted_map_add_to,3 );
    insertInternalGlobalMacro( "\\lciSortedMapForeach", sorted_map_foreach,2 );
    insertInternalGlobalMacro( "\\lciSortedMapClear", sorted_map_clear,1 );
    insertInternalGlobalMacro( "\\lciSortedMapIsEmpty", sorted_map_is_empty,1 );

    insertInternalGlobalMacro( "\\lciHashMapInsert", hash_map_insert,3 );
    insertInternalGlobalMacro( "\\lciHashMapIsDefined", hash_map_is_defined,2 );
    insertInternalGlobalMacro( "\\lciHashMapGet", hash_map_get,2 );
    insertInternalGlobalMacro( "\\lciHashMapClear", hash_map_clear,1 );

    insertInternalGlobalMacro( "\\lciIfFileExists", if_file_exists, 1);
    insertInternalGlobalMacro( "\\lciCopyFile",     copy_file, 2);

    insertInternalGlobalMacro( "\\lciChapter",     handle_chapter, 0);

    insertInternalGlobalMacro( "\\lciPopOutput",   pop_output,  0);
    insertInternalGlobalMacro( "\\lciPushOutput",  push_output, 1);
    insertInternalGlobalMacro( "\\lciOpenTmpFile", open_tmp_file, 1);
    insertInternalGlobalMacro( "\\lciOpenFile",    open_file, 1);
    insertInternalGlobalMacro( "\\lciCloseFile",   close_file, 0);

    insertInternalGlobalMacro( "\\lciLineNumber",  line_number,  0);

    insertInternalGlobalMacro( "\\lciIncludeOnly",    handle_include_only,  1);
    insertInternalGlobalMacro( "\\lciIfToBeIncluded", handle_to_be_included,  1);
    insertInternalGlobalMacro( "\\lciToHtmlWidth",    to_html_width,  1);

    insertInternalGlobalMacro( "\\lciRegexRegister",   handle_regex_register,   2 );
    insertInternalGlobalMacro( "\\lcRegexDoesMatch",   handle_regex_does_match, 2 );
    insertInternalGlobalMacro( "\\lcRegexStoreResult", handle_regex_store_result, 2 );

    insertInternalGlobalMacro( "\\lciStoreSavebox",    handle_store_savebox, 1 );
    insertInternalGlobalMacro( "\\lciSavestreamOpen",  handle_savestream_open, 1 );
    insertInternalGlobalMacro( "\\lciSavestreamUseX",  handle_savestream_use, 2 );
    insertInternalGlobalMacro( "\\lciSavestreamClose", handle_savestream_close, 1 );

    insertInternalGlobalMacro( "\\lciMinitocOpen",   handle_minitoc_open, 0 );
    insertInternalGlobalMacro( "\\lciMinitocClose",  handle_minitoc_close, 0 );

    insertInternalGlobalMacro( "\\lciShellCommand",  handle_shell_command, 1 );
}

// EOF //






