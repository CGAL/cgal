// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-75 $
// release_date  : $CGAL_Date: 2001/06/21 $
//
// file          : src/CGALWin/_file_panel.C
// package       : cgal_window (1.0)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.0
// revision_date : 20 June 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================


#include <CGAL/LEDA/file_panel.h>
#include <CGAL/LEDA/file.h>
#include <CGAL/LEDA/impl/x_basic.h>

#include <string>
#include <algorithm>

namespace CGAL {


enum { file_load, file_save, file_exit, file_all, file_cancel };

 
file_panel::~file_panel()
{ //dir_list.print("dir_list: ");
  //cout << endl;
  //file_list.print("file_list: ");
  //cout << endl;
 }


void file_panel::change_dir(char* dname)
{ if (strcmp(dname,"..") == 0)
  { window* call_win = window::get_call_window();
    file_panel* fp = (file_panel*)call_win->get_inf();
    fp->chdir();
   }
}

void file_panel::update_dir_menu(int)
{ window* call_win = window::get_call_window();
  file_panel* fp = (file_panel*)call_win->get_inf();
  fp->chdir();
}


void file_panel::update_file_menu(int)
{ window* call_win = window::get_call_window();
  file_panel* fp = (file_panel*)call_win->get_inf();
  fp->chdir();
}

void file_panel::init(window* wptr, string& fname, string& dname)
{ wp = wptr;
  load_handler = 0;
  save_handler = 0;
  cancel_handler= 0;
  
  load_ptr = 0;
  save_ptr = 0;
  cancel_ptr=0;
  
  if (dname == "") dname = ".";
  if (fname == "") fname = " ";
  
  file_name  = &fname;
  dir_name  = &dname;
  P.set_inf(this);
  load_string = "Load from File";
  save_string = "Save to File";
  panel_init = 0;

#if defined(__win32__) && (__BORLANDC__ != 0x550)
  mswin = true;
#else
  mswin = false;
#endif

  char* s;
  dir_list0.push_back("..");
  if ((s = getenv("HOME")) != 0)
  { home_dir = s;
    dir_list0.push_back(s);
   }
  if ((s = getenv("LEDAROOT")) != 0) 
  { dir_list0.push_back(s);
   }
}


file_panel::file_panel(string& fname, string& dname) : P("File Panel")
{ init(0,fname,dname); }

file_panel::file_panel(window& W, string& fname, string& dname): P("File Panel")
{ init(&W,fname,dname); }


void file_panel::set_load_object(const file_panel_handle_base& f)
{ load_ptr = &f; }

void file_panel::set_save_object(const file_panel_handle_base& f)
{ save_ptr = &f; }

void file_panel::set_cancel_object(const file_panel_handle_base& f)
{ cancel_ptr = &f; }



void file_panel::init_panel()
{
  if (panel_init) return;

  P.set_item_width(250);
  P.buttons_per_line(4);

  if ((load_handler || load_ptr) && (save_handler || save_ptr))
     P.text_item("\\bf File Panel");
  else  // was ~load_string and ~save_string
  { if (load_handler || load_ptr) P.text_item(string("\\bf\\blue ").append(load_string));
    if (save_handler || save_ptr) P.text_item(string("\\bf\\blue ").append(save_string));
   }

  P.text_item("");

  dir_list.push_back(*dir_name);
  file_list.push_back(*file_name);

  dir_item  = P.string_item("directory",*dir_name,dir_list,8,change_dir);
  //dir_item  = P.string_item("directory",*dir_name,dir_list,8);
  if (pat_string != "") 
              P.string_item("filter",pat_string);
  file_item = P.string_item("file name",*file_name,file_list,8);

  P.set_item_menu_func(dir_item,update_dir_menu);
  P.set_item_menu_func(file_item,update_file_menu);  

  chdir();

  if (load_handler || load_ptr) P.fbutton("load",file_load);
  if (save_handler || save_ptr) P.fbutton("save",file_save);

  P.button("cancel",file_cancel);

  panel_init++;
}


void file_panel::chdir()
{
  string old_dir = get_directory();

  if (old_dir != *dir_name && 
      std::find(dir_list0.begin(),dir_list0.end(),old_dir) == dir_list0.end() ) 
  { //cout << "append " << old_dir << endl;
    dir_list0.push_back(old_dir);
   }
  

  if ((*dir_name)[0] == '$')
  { //char* s = getenv(dir_name->del(0));
  
    string hlp = dir_name->substr(1);
  
    char* s = getenv(hlp.c_str());
    if (s) *dir_name = s;
   }

  if (is_directory(*dir_name))
  {
    set_directory(*dir_name);

    *dir_name = get_directory();

    if (pat_string.length() > 0)
       file_list = get_files(*dir_name,pat_string);
    else
       file_list = get_files(*dir_name);

    file_list.sort();

    // delete hidden files 
    std::list<string>::iterator it = file_list.begin();
    std::list<string>::iterator it2,it3;
    
    for(;it != file_list.end(); it++)
    { string s = *it;
      it3 = it; // maybe too complicated
      it2 = it++;
      it = it3;
      if (s[0] == '.') file_list.erase(it);
      it = it2;
     }

    dir_list = get_directories(*dir_name);
    dir_list.sort();
    while (!dir_list.empty() && (*(dir_list.begin()))[0] == '.') dir_list.pop_front();
    string s;
    it = dir_list0.begin();
    for(;it != dir_list0.end();it++) {
      s = *it;
      dir_list.push_front(s);
    }
  }
  else
  { dir_list.clear();
    string s;
    std::list<string>::iterator it = dir_list0.begin();
    
    for(;it != dir_list0.end();it++) {
      s = *it;
      dir_list.push_front(s);
    }
    
    string::size_type ST;
    ST = dir_name->find("\\");
    
    while (ST != string::npos){
     // replace ...
     dir_name->replace(ST,ST+2,"/");
     // was  s = dir_name->replace_all("\\","/");
     // new find ...
     ST = dir_name->find("\\");
    }
    
    s = *dir_name;
    
    panel W("File Panel");
    W.text_item("");
    W.text_item(string("\\bf\\blue ").append(s).append(":~~no such directory."));
    W.text_item("");
    W.fbutton("continue");
    W.open();
    file_list.clear();
  }
  P.set_menu(dir_item,dir_list,8);
  P.set_menu(file_item,file_list,8);
  if (P.is_open()) P.redraw_panel();
}



void file_panel::x_open()
{ int result;
  init_panel();
  if (wp)
    result = P.open(*wp,window::center,window::center);
  else
    result = P.open(window::center,window::center);

  switch (result) {
  case file_load: chdir();
                  if (load_ptr) { load_ptr->operator()(*file_name); }
                  else load_handler(*file_name);
                  break;
  case file_save: chdir();
                  if (save_ptr) { save_ptr->operator()(*file_name); }  
                  else save_handler(*file_name);
                  break;
  }
}



void file_panel::x_open(int xpos, int ypos)
{ int result;
  init_panel();
  result = P.open(xpos,ypos);

  switch (result) {
  case file_load: chdir();
                  if (load_ptr) { load_ptr->operator()(*file_name); }  
                  else load_handler(*file_name);
                  break;
  case file_save: chdir();
                  if (save_ptr) { save_ptr->operator()(*file_name); }    
                  else save_handler(*file_name);
                  break;
  }
}






void file_panel::ms_open()
{ 
/*
  window W(425,265);
  W.display(*wp,window::center,window::center);
*/

  window& W = *wp;

  char dname[128];
  strcpy(dname,(*dir_name).c_str());

  char fname[128];
  strcpy(fname,(*file_name).c_str());

  string filter = descr_string + "|" + pat_string;

  filter += "|All Files (*.*)|*.*|";

  int len = filter.length();

  for(int i=0; i< len; i++)
  { char& c = filter[i];
    if (c == '|') c = '\0';
   }

  if (load_handler || load_ptr)
  { if (x_choose_file(W.draw_win,0,load_string.c_str(),filter.c_str(),dname,fname)) 
    { *file_name = fname;
      *dir_name = dname;
      if (load_ptr) { load_ptr->operator()(fname); }  
      else load_handler(fname);
     }
   }

  if (save_handler || save_ptr)
  { if (x_choose_file(W.draw_win,1,save_string.c_str(),filter.c_str(),dname,fname)) 
    { *file_name = fname;
      *dir_name = dname;
      if (save_ptr) { save_ptr->operator()(fname); }    
      else save_handler(fname);
     }
   }
}



void file_panel::ms_open(int, int) { ms_open(); }

}
