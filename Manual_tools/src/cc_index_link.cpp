          
/**************************************************************************
 
  cc_index_link.cpp
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual_index.sty.
  Function  : adds links in HTML-index (combines manual.ind with HREF)
  System    : C++ (g++)
  Author    : 2003 Renata Krysta
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <algorithm>
#include <list>

struct Lines{
  int number;
  std::string text;

  Lines(int n, const std::string& s) : number(n), text(s) {}
};

class Lines_eq : public std::unary_function<Lines, bool> {
  int ref_number;
public:
  explicit Lines_eq(int number) : ref_number(number) {}
  bool operator() (const Lines& l) const {return l.number == ref_number;}
};

/* >main: main function with standard unix parameter input */
/* ------------------------------------------------------- */
 
main( int argc, char **argv) {
// argv[1] is a file with sorted index entries and their numbers
// argv[2] is a file with links and their numbers
// argv[3] is an output file
  if (argc < 4) {
      std::cerr << "*** Error: program needs an additional parameter"<<'\n';
      exit(1);
  }

  char* file_name1 = argv[1];
  char* file_name2 = argv[2];
  char* out_file_name = argv[3]; 

  std::ofstream out_file(out_file_name);
  if (!out_file)
     std::cerr<<"*** Error: cannot open file "<< out_file_name <<'\n';

  std::ifstream in_file1(file_name1);
  if (!in_file1)
     std::cerr<<"*** Error: cannot open file "<< file_name1 <<'\n';

  std::list<Lines> list_of_lines;
  std::list<Lines>::const_iterator lli;

  int number;
  std::string link;
  while (in_file1) {
    in_file1 >> number >> link;
    list_of_lines.push_front(Lines(number, link));
  }

  std::ifstream in_file2(file_name2);
  if (!in_file2)
     std::cerr<<"*** Error: cannot open file "<< file_name2 <<'\n';
 
  char ch;
  std::string name_link;

  bool end_name;
  while (in_file2.get(ch)) {
    if (ch=='?') {
      if (in_file2.get(ch)) {
         if (ch=='?') {
             if (in_file2.get(ch)) {
                if (ch=='?') {    // if index entry
                  end_name=0; 
		  std::string tmp;
                  bool page = 0;
                  while (in_file2 && !end_name) {
                    in_file2>> name_link;

                    if (name_link.compare("???")!=0 && 
			name_link.compare("</TD></TR>")!=0 && 
			name_link.compare("!!!")!=0) {
		      tmp+=" ";
                      tmp+=name_link;
                    } else {
                       if (name_link.compare("???")==0) page = 1;
                       else {
                          if (name_link.compare("!!!")==0) {
                            in_file2 >> name_link; 
                            out_file << "<H2><A NAME=\"Index" << name_link
                                     << "\">";                               
                            name_link+="</A></H2>"; 
                          }
                          else out_file << tmp <<" ";
                       } 
                       end_name=1;
                    }
                  }
                  if (page) {    //if link
                    in_file2 >> number;
                    lli = find_if(list_of_lines.begin(), list_of_lines.end(), 
                                  Lines_eq(number));
		    bool number_found = (lli != list_of_lines.end());
		    if (number_found) {
		      out_file <<"<A "<< lli->text <<">" <<tmp << "</A> ";
		    } 
		    else 
		      out_file << tmp;
		    
		    while (in_file2.get(ch) && ch!='<') {
		      if (ch=='|') {
			in_file2 >> number;
			if (number_found) {
			  lli = find_if(list_of_lines.begin(), 
                                        list_of_lines.end(), 
                                        Lines_eq(number));
			  out_file << "<A " << lli->text << ">&rarr;</A> ";
			}  
		      }
		    } 
 
		    out_file << ch;
                  } else out_file << name_link; 

                 } else out_file << ch;
              }
          } else out_file << ch;
      } 
    }  else out_file << ch;
  }

in_file1.close();
in_file2.close();
out_file.close();

return 0;
}

// EOF //

