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


static int index_anchor_counter = 0;

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
 //new_main_filename+"\">" +chapter_title +"</A>" ;    
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
        

    if (macroIsTrue("\\lciIfRef"))
      *main_stream << "<TABLE WIDTH=100%> <TR> <TD ALIGN=LEFT VALIGN=TOP> <H1>" 
             << chapter_title << "</H1> </TD>" << endl;
    else
      *main_stream << "<H1>" << chapter_title  << "</H1>" << endl;


    // table of contents
   
    *contents_stream << "    <LI> <A HREF=\"" << main_filename
		     << "\">" << chapter_title << "</A>" << endl;      
}



void handlePart(  const Buffer_list& T) {
  

    if (!part_title.empty()) // end previous part 
    {
       *contents_stream << "</OL>"  << endl;
       *contents_stream << "<!-- End of manual part -->"  << endl;
     //  chapter_num=0; 
    }
       
    part_title = string( text_block_to_string( T));

    // add new part title to table of contents
    *contents_stream << "<!-- Start of new manual part -->"  << endl;
    *contents_stream << "<H3>" << part_title << "</H3>" << endl; 
    if ( macroIsTrue( "\\lciIfNumberChaptersByPart") )
       *contents_stream << "<OL>" << endl;
    else
       *contents_stream << "<OL START=" << chapter_num+1 << ">" << endl;
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
           cerr << " internal_macros.C : warning: unknown variable WhichItem"
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

// list of modified entries

typedef struct Lines{
  string text;
  int number; 
  struct Lines *next;
} sList, *pList;

pList modifier_list=NULL;


pList AddToList(pList list , string l, int n) {
  if (list==NULL) {
     list = new Lines;
     list->text= l;
     list->number= n;
     list->next = NULL;
  } else {
     pList p = new Lines;
     p->text= l;
     p->number= n;
     p->next = list;
     list = p;
 }
 return list;
}

bool search(pList list, string l, int& number) {
     if (list==NULL) return false;
        while ((list!=NULL) ) {
           if (list->text==l) {
               number = list->number;
               return true;
           } else  list=list->next;
        }
     return false;
} 



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
         handleIndex2(ord_index_name+"@ ??? "+index_name,sub_item ,sub_sub_item,
                      0);
   }
   else { 
       int number;
       if (search(modifier_list, modifier+index_name, number)) { 
            handleIndex2(ord_index_name+" "+ord_modifier+"@ ??? "+index_name+ 
                         ",  "+ modifier+"<A NAME=\""+ int_to_string(number) 
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
                modifier_list = AddToList(modifier_list, 
                                          modifier+index_name,HREF_counter-4);
            }
        }
   }
} 


void handleIndex2(string main_item, string sub_item, string sub_sub_item, int modifier) {
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
   } else *current_ostream<<"NIE1"<<main_item<<endl;
  
     
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
   handleIndex2(index_name,"traits class@ ??? traits class",index_class_name,0);
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
           *main_stream  << "<UL><LI>\n" << formatted_reference
		      << ".</UL>\n" << endl;
    }

    // table of contents
    if ( macroIsTrue( "\\lciIfHtmlClassToc"))
           *contents_stream << "<UL><LI> " <<  formatted_reference
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
         macroIsTrue( "\\lciIfHtmlIndex")  &&  macroIsTrue( "\\ccIndex")
         &&  macroIsTrue( "\\ccAutoIndex"))    // Index.
            handleIndex2(class_name +"@ ??? "+
                        convert_C_to_html(class_name),"","",0);


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
    handleIndex2(param[1]+"@ ??? "+convert_C_to_html(param[1]),"","",0);
    return "";
}

string
html_index_C( const string&, string param[], size_t n, size_t opt) {
    NParamCheck( 1, 0);  // uses cc_string as index item, param[0] is
                         // index category
    crop_string( param[0]);
    handleIndex2(cc_string+"@ ??? "+convert_C_to_html(cc_string),"","",0); 
    return "";
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






