          
/**************************************************************************
 
  cc_index_link.C
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual_index.sty.
  Function  :add links in HTML-index
  System    : bison, flex, C++ (g++)
  Author    : 2003 Lutz Kettner, Renata Krysta
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/


#include <stdio.h>
#include <iostream.h>
#include <fstream.h>


typedef struct Lines{
  int number;
  char *text;
  struct Lines *next;
} *pList;


pList list=NULL;


pList AddToList(pList list , const int n, const char *l) {
  if (list==NULL) {
     list = new Lines;
     list->number = n;
     list->text= new char[strlen(l)];
     strcpy(list->text,l);
     list->next = NULL;
  } else {
     pList p = new Lines;
     p->number = n;
     p->text= new char[strlen(l)];
     strcpy(p->text,l);
     p->next = list;
     list = p;
 }
 return list;
}



char* search(pList list, int number) {
     if (list==NULL) return '\0';
        while ((list!=NULL) ) {
           if (list->number==number) {
               return list->text;
           } else  list=list->next;
        }
} 



/* >main: main function with standard unix parameter input */
/* ------------------------------------------------------- */
 
main( int argc, char **argv) {
  if (argc < 5) {
      cerr << "*** Error: program needs an additional parameter"<<'\n';
      exit(1);
  }
  char* file_name1 = argv[1];
  char* file_name2 = argv[2];
  char* out_file_name = argv[3]; 

  ofstream out_file(out_file_name);
  if (!out_file)
     cerr<<"*** Error: cannot open file "<< out_file_name <<'\n';

  ifstream in_file1(file_name1);
  if (!in_file1)
     cerr<<"*** Error: cannot open file "<< file_name1 <<'\n';

  int number;
  char link[400];
  while (in_file1) {
    in_file1 >> number >> link;
    list=AddToList(list, number, link);                                         
  }


  ifstream in_file2(file_name2);
  if (!in_file2)
     cerr<<"*** Error: cannot open file "<< file_name2 <<'\n';
 
 
 char ch;

  char name_link[300];

  bool end_name;
  while (in_file2.get(ch)) {
    if (ch=='?') {
      if (in_file2.get(ch)) {
         if (ch=='?') {
             if (in_file2.get(ch)) {
                if (ch=='?') {
                  end_name=0; 
                  char tmp[300]="";
                  bool page = 0;
                  while (in_file2 && !end_name) {
                    in_file2>> name_link;
                    if (strcmp(name_link,"???")!=0 && 
                                strcmp(name_link,"</TD></TR>")!=0 && 
                                strcmp(name_link,"!!!")!=0) {
                       strcat(tmp," ");
                       strcat(tmp,name_link);
                    } else {
                       if (strcmp(name_link,"???")==0) page = 1;
                       else {
                          if (strcmp(name_link,"!!!")==0) {
                            in_file2 >> name_link; 
                            out_file << "<H2><A NAME=\"Index" << name_link
                                     << "\">";                               
                            strcat(name_link,"</A></H2>"); 
                          }
                          else out_file << tmp <<" ";
                       } 
                       end_name=1;
                    }
                  }
                
                  if (page) {
                    in_file2 >> number;
                    out_file <<"<A "<<search(list,number)<<">" <<tmp 
                             << "</A> ";
                    while (in_file2.get(ch) && ch!='<') {
                      if (ch=='|') {
                        in_file2 >> number;
                        out_file << "<A "<< search(list,number)<<">" 
                                 << "<img SRC=\""<< argv[4] << "/index_arrow.gif\" ALT=\"reference\" WIDTH=\"14\" HEIGHT=\"12\" VALIGN=BOTTOM BORDER=0></A> ";  
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

