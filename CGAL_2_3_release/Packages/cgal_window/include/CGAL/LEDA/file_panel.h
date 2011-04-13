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
// file          : include/CGAL/LEDA/file_panel.h
// package       : cgal_window (1.0.3)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.0.3
// revision_date : 25 June 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================


#ifndef CGAL_WINDOW_FILE_PANEL_H
#define CGAL_WINDOW_FILE_PANEL_H

#if defined(CGAL_USE_CGAL_HEADERS)
#include <CGAL/basic.h>
#endif

#include <CGAL/LEDA/window.h>

namespace CGAL {


typedef void (*file_panel_handler)(string);

// handler object ...
class __exportC file_panel_handle_base {
protected:
  file_panel_handler  FPH;
public:
  file_panel_handle_base() { FPH = NULL; }
  file_panel_handle_base(file_panel_handler f) { FPH = f; }
  
  virtual void operator()(string s) const
  { if (FPH != NULL) FPH(s); }
};


class __exportC file_panel {

panel P;

window* wp;

string*      dir_name;
std::list<string> dir_list0;
std::list<string> dir_list;
panel_item   dir_item;

string*      file_name;
std::list<string> file_list;
panel_item   file_item;

string       load_string;
string       save_string;
string       pat_string;
string       descr_string;

int          panel_init;

bool         mswin;

file_panel_handler load_handler;
file_panel_handler save_handler;
file_panel_handler cancel_handler;

//load/save/cancel objects ...
const file_panel_handle_base* load_ptr; 
const file_panel_handle_base* save_ptr; 
const file_panel_handle_base* cancel_ptr; 


void*  dummy_ptr;

string home_dir;

void chdir();

void init(window*,string& ,string&);

static void change_dir(char* dname);
static void update_dir_menu(int);
static void update_file_menu(int);

public:

 panel& get_panel() { return P; }

 void init_panel();

 file_panel(string& fname, string& dname);
 file_panel(window& W, string& fname, string& dname);
~file_panel();

 void set_load_handler(file_panel_handler f)  { load_handler  = f;  }
 void set_save_handler(file_panel_handler f)  { save_handler  = f;  }
 void set_cancel_handler(file_panel_handler f){ cancel_handler = f; }
 
 void set_load_object(const file_panel_handle_base& f);
 void set_save_object(const file_panel_handle_base& f);
 void set_cancel_object(const file_panel_handle_base& f);

 void set_load_string(string s) { load_string = s; }
 void set_save_string(string s) { save_string = s; }
 void set_pattern(string s)     { pat_string = s; }
 void set_pattern(string s1, string s2) { descr_string = s1, pat_string = s2; }

 void set_mswin(bool b) { mswin = b; }

 void set_frame_label(string s) { P.set_frame_label(s); }

 void x_open();
 void x_open(int x, int y);

 void ms_open();
 void ms_open(int x, int y);

 void open()             { if (mswin) ms_open();    else x_open();    }
 void open(int x, int y) { if (mswin) ms_open(x,y); else x_open(x,y); }

};

}

#endif
