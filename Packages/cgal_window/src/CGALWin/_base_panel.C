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
// file          : src/CGALWin/_base_panel.C
// package       : cgal_window (1.0)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 0.9.7
// revision_date : 23 May 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================

#include <CGAL/LEDA/base_window.h>
#include <CGAL/LEDA/impl/x_basic.h>

#if defined(_MSC_VER)
#pragma warning(disable:4305)
#pragma warning(disable:4309)
#endif

#include "bitmaps/pstyle.h"
#include "bitmaps/lstyle.h"
#include "bitmaps/lwidth.h"


#include <cassert>
#include <cstdio>

#if defined(__BORLANDC__)
using std::isprint;
#endif


namespace CGAL {


// static members


BASE_WINDOW* BASE_WINDOW::call_window = 0;
panel_item   BASE_WINDOW::call_item = 0;


BASE_WINDOW* BASE_WINDOW::get_call_window() 
{ return BASE_WINDOW::call_window; }

panel_item BASE_WINDOW::get_call_item()   
{ return BASE_WINDOW::call_item;   }

int BASE_WINDOW::get_call_button() 
{ panel_item it = BASE_WINDOW::call_item;  
  return (it && it->kind == Button_Item) ? it->dat1 : -1;
}





void BASE_WINDOW::scroll_down_action(int)
{ ((BASE_WINDOW*)call_window->get_inf())->scroll(-1); }

void BASE_WINDOW::scroll_up_action(int)
{ ((BASE_WINDOW*)call_window->get_inf())->scroll(+1); }


void BASE_WINDOW::scroll_drag_action(int i)
{ 
  BASE_WINDOW* wp = (BASE_WINDOW*)call_window->get_inf();
  char* scroll_buf;

  wp->clipping(0);

  int menu_h = wp->item_count * wp->yskip;
  int pan_h  = wp->panel_height;
  int pan_w  = wp->panel_width;

  panel_item p0 = wp->Item[0];

  if (i == -1) // start
  { // init buffer
   int y0 = p0->ycoord;
   int j;
   for(j=0; j<wp->item_count; j++) wp->Item[j]->ycoord = j*wp->yskip;
   wp->start_buffering(pan_w,menu_h+pan_h);
   wp->panel_height = menu_h;
   wp->clear();
   wp->draw_panel_items();
   wp->stop_buffering(scroll_buf);
   wp->panel_height  = pan_h;
   wp->window_height = pan_h;
   for(j=0; j<wp->item_count; j++) wp->Item[j]->ycoord = y0 + j*wp->yskip;
   wp->set_inf(scroll_buf);
   return;
  }

  scroll_buf= (char*)wp->get_inf();

  if (i == -2) // finish
  { // round to next button position and delete buffer
    int y0 = p0->ycoord;
    int dy = y0 % wp->yskip;
    int y1 = y0 - dy;
    int d  = 1;

    if (dy < -wp->yskip/2)
    { y1 = y0 - (wp->yskip + dy);
      d = -1;
     }

    while (y0 != y1+d)
    { x_insert_pixrect(wp->draw_win,0,menu_h+y0+pan_h-1,scroll_buf);
      y0 += d;
     }

    y0 -= d;

    x_delete_pixrect(scroll_buf);
    for(int j=0; j < wp->item_count; j++) 
       wp->Item[j]->ycoord = y0 + j*wp->yskip;
    return;
   }

  // drag

  double f = i/1000.0; 
  p0->ycoord = -int(f*(menu_h-pan_h));
  x_insert_pixrect(wp->draw_win,0,menu_h+p0->ycoord+pan_h-1,scroll_buf);
}


//------------------------------------------------------------------------------
// some auxiliary functions for string manipulation
//------------------------------------------------------------------------------

static char* string_dup(const char* p)
{ if (p == 0) return 0;
  char* q = new char[strlen(p)+1];
  if (q==0) 
  { std::cerr << "string_dup: out of memory";
    abort();
   }
  strcpy(q,p);
  return q;
}

static char* make_string(int x)
{ char str[256];
  CGAL_CLIB_STD::sprintf(str,"%d",x);
  return string_dup(str);
 }

static char* make_string(double x)
{ char str[256];
  CGAL_CLIB_STD::sprintf(str,"%f",x);
  return string_dup(str);
 }


//------------------------------------------------------------------------------
// construction and destruction of panel items
//------------------------------------------------------------------------------


Panel_Item::Panel_Item(int k, char* s, void* addr, int h, int i)
{ kind = k;
  enabled = true;
  label_str = s;
  ref = addr;
  height  = h;
  index = i;
  data_str = 0;
  help_str = 0;
  xcoord = 0;
  ycoord = 0;
  width = 0;
  dat1 = 0;
  dat2 = 0;
  offset = 0;
  argc = 0;
  argv = 0;
  
  action = 0;
  str_action = 0;
  
  action_ptr = 0;
  str_action_ptr = 0;
  
  menu_win = 0;
  menu_but = 0;
  shortcut = -1;
}


Panel_Item::~Panel_Item()
{ 
  if (label_str) delete[] label_str;
  if (help_str)  delete[] help_str;

  if (kind != String_Item && kind != String_Menu_Item && data_str) 
    delete[] data_str;

  if (argc > 0)
  { for(int i = 0; i < argc; i++)
    { if (kind == Button_Item || 
          kind == Bitmap_Choice_Item || kind == Bitmap_Choice_Mult_Item) 
         x_delete_bitmap(argv[i]);
      else
         delete[] argv[i];
     }
    delete[] argv;
   }

  if (kind == Button_Item && ref) 
  { ((BASE_WINDOW*)ref)->owner_item = 0;
    //if (menu_win) delete menu_win;
   }

  if (menu_win) delete menu_win;

}



//------------------------------------------------------------------------------
// panel member functions
//------------------------------------------------------------------------------

void BASE_WINDOW::panel_init()
{ 
  item_count = 0;
  but_count = 0;

  item_bg_color = white;
  bitmap_color0 = black;
  bitmap_color1 = black;

  if ( mono() )
   { panel_bg_color = white;
     shadow_color = black;   
     item_bg_color = white;   
     press_color   = white;
     disable_color  = white;
    }
  else
   { panel_bg_color = grey1;
     press_color    = grey2;
     shadow_color   = grey3;
     disable_color  = grey2;
    }

  panel_width = 0;
  panel_height = 0;


  x_set_button_font(0);
  tw = x_text_width(0,"H");
  th = x_text_height(0,"H");

  ytskip   = th + 2;
  yskip    = th + 11;
  yoff     = th/2 - 1;

  xoff     = tw + 2;

  slider_h = 12;
  color_h  = 12;

  button_h = th + 3;
  //button_h = th + 4;

  button_d = th/2;
  choice_h = button_h;

  label_w   = 10;          // minimal label length
  sl_num_w  = 4*tw+tw/2;   // slider num field width
  slider_w  = 16*color_h;  // minimal slider length
  string_w  = slider_w;    // minimal string item length
  number_w  = string_w;    // minimal int/float item length
  choice_w  = 35;          // minimal choice field width
  button_w  = 20;          // minimal button width

  buts_per_line = 0;  
  center_button_label = 1;
  //menu_button_style = 0;
  menu_button_style = 1;
  menu_size = 0;

  act_str_item = 0;
  last_sel_button = 0;
  panel_menu_mode = 0;
  has_menu_bar = 0;
  menu_bar_height = 0;
  focus_button = 0;

  x_set_fixed_font(0);
  tw = x_text_width(0,"H");
  th = x_text_height(0,"H");

  string_h = th + 2;
}


void BASE_WINDOW::set_item_height(int h)  { yskip = h; button_h = h-2; }

void BASE_WINDOW::set_item_width(int w)
{ slider_w = w;
  string_w = w;
  number_w = w;
 }

void BASE_WINDOW::set_item_space(int s)   { yoff = s; }

void BASE_WINDOW::set_button_space(int d) { button_d = d; }


void BASE_WINDOW::item_error()
{ std::cerr << "sorry, too many items."; 
  exit(1);
 }


//------------------------------------------------------------------------------
// adding items
//------------------------------------------------------------------------------

panel_item BASE_WINDOW::new_panel_item(int k, const char* s, void* addr, int h,
                                                             const char* hlp)
{ 
  if (item_count >= MAX_ITEM_NUM) item_error();

  char* s1 = string_dup(s);
  char* spos = s1;
  while (*spos != '\0' && *spos != '&') spos++;
  if (*spos) {
    for(char* p = spos; *p; p++) *p = *(p+1);
   }

  panel_item p = new Panel_Item(k,s1,addr,h,item_count);

  p->shortcut = (*spos) ? (spos - s1) : -1;

  p->help_str = string_dup(hlp);

  Item[item_count++] = p;

  if (k != Text_Item && k != Button_Item)
  { x_set_button_font(draw_win);
    int w = x_text_width(draw_win,s) + 1;
    if (k == Slider_Item) 
       w += sl_num_w;
    else 
       if (k == String_Menu_Item) 
         w += 4*tw;
       else
         w += 2*tw;
    if (w > label_w) label_w = w;
    x_set_text_font(draw_win);
   }
  return p;
}


void BASE_WINDOW::hspace(int d)
{ panel_item p = new_panel_item(Text_Item,"",0,0,0);
  p->width = d;
  p->height = 0;
}

void BASE_WINDOW::vspace(int d)
{ panel_item p = new_panel_item(Text_Item,"",0,d,0);
  p->width = 0;
  p->height = d;
}


void BASE_WINDOW::new_line()
{ new_panel_item(Newline_Item,"",0,0,0); 
  but_count++;
 }

void BASE_WINDOW::separator()
{ new_panel_item(Newline_Item,"",0,0,0); 
  but_count++;
 }


void BASE_WINDOW::fill_line()
{ new_panel_item(Fill_Item,"",0,0,0); }


  

panel_item BASE_WINDOW::text_item(const char* s)
{ panel_item it = (item_count > 0)  ? Item[item_count-1] : 0;

  if (s==0) s = "";

  ytskip = x_text_height(draw_win,"H");

  if (strlen(s) == 0 || it == 0 || it->kind != Text_Item)
    { if (strlen(s) == 0)
         ytskip = x_text_height(draw_win,"H") + 4;
      it = new_panel_item(Text_Item,"",0,ytskip,0); 
      it->argv = new char*[2];
      it->argc = 0;
     }

  // split into words

  int    argc1 = 0;
  char** argv1;

  split_text(s,argc1,argv1);

  int    argc0 = it->argc;
  char** argv0 = it->argv;

  it->argc = argc0 + argc1;

  if (it->argc > 0)
    it->argv = new char*[it->argc];
  else
    it->argv = 0;

  int i;
  int n = 0;
  for(i = 0; i < argc0; i++) it->argv[n++] = argv0[i]; 
  for(i = 0; i < argc1; i++) it->argv[n++] = argv1[i]; 

  delete[] argv0;
  delete[] argv1;

  return it;
}

void BASE_WINDOW::set_text(panel_item it, const char* s)
{
  if (it->kind != Text_Item)
  { std::cerr << "illegal item in window::set_text\n";
    return;
   }

  delete[] it->argv;

  split_text(s,it->argc,it->argv);
}


panel_item BASE_WINDOW::string_item0(double x, double y, void* s)
{ panel_item p = new Panel_Item(String_Item,(char*)"",s,string_h,-1);
  p->xcoord = xpix(x) - label_w;
  p->ycoord = ypix(y);
  p->help_str = 0;
  p->data_str = access_str(s);
  p->offset = strlen(p->data_str);
  p->dat1 = 0;
  p->dat2 = string_w/tw - 1;
  p->str_action = 0;
  return p;
 }
 



panel_item BASE_WINDOW::string_item(const char* s, void* x, 
                                    panel_str_action_func F, const char* hlp)
{ panel_item p = new_panel_item(String_Item,s,x,string_h,hlp);
  p->data_str = access_str(x);
  p->offset = strlen(p->data_str);
  p->dat1 = 0;
  p->dat2 = string_w/tw - 1;
  //if (p->offset > p->dat2) p->offset = p->dat2;
  p->str_action = F;
  return p;
 }




void str_menu_selection(int code)
{ 
  int sel = code & 0xFFFF;
  int i   = (code>>16) & 0xFFFF;

  // I assume that str_menu_selection is called after closing the 
  // sub-menu and that at this moment active_window points to
  // the window that contains the string item to be changed


  BASE_WINDOW* str_menu_win = BASE_WINDOW::active_window;
  panel_item   str_menu_it = str_menu_win->Item[i];

  if (str_menu_it && sel >= 0)
  { str_menu_win->assign_str(str_menu_it->ref,str_menu_it->argv[sel]);
    str_menu_it->data_str = str_menu_win->access_str(str_menu_it->ref);
    str_menu_it->offset = strlen(str_menu_it->data_str);
    str_menu_it->dat1 = 0;
    str_menu_it->dat2 = str_menu_win->string_w/str_menu_win->tw - 1;
    str_menu_win->clipping(1);

    int draw_w = str_menu_win->draw_win;

    int          save_lw = x_set_line_width(draw_w,1);
    line_style   save_ls = x_set_line_style(draw_w,solid);
    text_mode    save_tm = x_set_text_mode(draw_w,transparent);
    drawing_mode save_mo = x_set_mode(draw_w,src_mode);
    
    //str_menu_win->draw_string_item(str_menu_it);

    str_menu_win->activate_string_item(str_menu_it,0);

    x_set_line_width(draw_w,save_lw);
    x_set_line_style(draw_w,save_ls);
    x_set_text_mode(draw_w,save_tm);
    x_set_mode(draw_w,save_mo);

    if (str_menu_it->str_action || str_menu_it->str_action_ptr) 
    { //clipping(2);
      BASE_WINDOW::call_window = str_menu_win;
      BASE_WINDOW::call_item = str_menu_it;
      
      if (str_menu_it->str_action_ptr){
        str_menu_it->str_action_ptr->set_params((window*)str_menu_win,(char*)str_menu_it->argv[sel] );
        str_menu_it->str_action_ptr->operator()();
      }
      else str_menu_it->str_action((char*)str_menu_it->argv[sel]);
      //clipping(1);
     }

   }
}


void BASE_WINDOW::set_menu(panel_item p, int argc,const char** argv, int sz)
{ 
  if (p->kind != String_Menu_Item)
  { std::cerr << "illegal item " << p->label_str << " in window::set_menu\n";
    return;
   }

  const char* none_str = "none";

  if (argc == 0) 
  { argc = 1; 
    argv = &none_str;
   }

  if (p->argc > 0) 
  { for(int i = 0; i < p->argc; i++) delete[] p->argv[i];
    delete[] p->argv;
   }

/*
  if (argc == 0) 
  { p->argc = 0;
    p->argv = 0;
    p->menu_but->ref = 0;
    if (p->menu_win) delete p->menu_win;
    p->menu_win = 0;
    return;
   }
*/

  p->argc = argc;
  p->argv = new char*[argc];

  int max_choice_len = 0;
  int i;
  for(i = 0; i < argc; i++) 
  { p->argv[i] = string_dup(argv[i]);
    int l = strlen(argv[i]);
    if (l > max_choice_len) max_choice_len = l;
   }

  // build menu panel

  int cols = 1;
  int rows = argc;

  if (sz == 0)
  { if (argc < 20) cols = 1;
    else if (argc < 40) cols = 2;
    else if (argc < 80) cols = 3;
    else cols = 4;
    rows = (argc + cols -1)/cols;
   }



  if (p->menu_win) delete p->menu_win;


  p->menu_win = new BASE_WINDOW(-1,-1);
  p->menu_win->buttons_per_line(cols);

  //p->menu_win->menu_button_style = (sz > 0) ? 2 : 1;
  p->menu_win->menu_button_style = (sz > 0) ? 2 : 0;

  p->menu_win->menu_size = sz;

  for(int r = 0; r<rows; r++)
    for(int c = 0; c<cols; c++) 
    { int i = c*rows + r;
      if (i < argc)
         p->menu_win->button(p->argv[i],(p->index<<16) + i, 
                                        str_menu_selection);
      else
         p->menu_win->button("",0,(BASE_WINDOW*)0);
     }

  if (argv == &none_str)
  { panel_item it = p->menu_win->get_item(argv[0]);
    it->enabled = 0;
   }

  p->menu_but->ref = p->menu_win;
  p->menu_win->owner_item = p->menu_but;


}



panel_item BASE_WINDOW::string_menu_item(const char* s, void* x, const char*,
                                  int argc,const char** argv, int size,
                                  panel_str_action_func F, const char* hlp)
{ 
  if (argc == 0)
     return string_item(s,x,0,hlp);

  panel_item p = new_panel_item(String_Menu_Item,s,x,string_h,hlp);
  p->data_str = access_str(x);
  p->offset = strlen(p->data_str);
  p->dat1 = 0;
  p->dat2 = string_w/tw - 1;
  if (p->offset > p->dat2) p->offset = p->dat2;

  button("",0,(BASE_WINDOW*)0);
  but_count--;

  panel_item q = Item[item_count-1];
  p->menu_but = q;
  q->width = string_h;
  q->height = string_h;
  q->menu_but = q;

  //p->menu_win->button_w += 2*string_h-8; 

  p->str_action = F;

  p->argc = 0;
  p->argv = 0;
  set_menu(p,argc,argv,size);

  return p;
}




panel_item BASE_WINDOW::int_item(const char* s, int* x, const char* hlp)
{ panel_item p = new_panel_item(Int_Item,s,x,string_h,hlp);
  p->data_str = make_string(*x);
  p->offset = strlen(p->data_str);
  p->dat1 = 0;
  p->dat2 = number_w/tw - 1;
  if (p->offset > p->dat2) p->offset = p->dat2;
  return p;
 }

panel_item BASE_WINDOW::float_item(const char* s, double* x, const char* hlp)
{ panel_item p = new_panel_item(Float_Item,s,x,string_h,hlp);
  p->data_str = make_string(*x);
  p->offset = strlen(p->data_str);
  p->dat1 = 0;
  p->dat2 = number_w/tw - 1;
  if (p->offset > p->dat2) p->offset = p->dat2;
  return p;
 }


panel_item BASE_WINDOW::slider_item(const char* s, int* x, int low, int high, 
                                                      panel_action_func F,
                                                      const char* hlp)
{ panel_item p = new_panel_item(Slider_Item,s,x,slider_h,hlp);
  p->data_str = make_string(*x);
  p->dat1 = low;
  p->dat2 = high;
  p->action = F;
  return p;
 }

panel_item BASE_WINDOW::choice_item(const char* s, int* x, int argc,
                             const char** argv, int step, int off,
                             panel_action_func F, const char* hlp)
{ 
  panel_item p = new_panel_item(Choice_Item,s,x,choice_h,hlp);
  p->data_str = 0;
  p->dat2 = step;
  p->offset = off;
  p->argc = argc;
  p->argv = new char*[argc];
  for(int i=0; i<argc; i++) 
  { int w;
    p->argv[i] = string_dup(argv[i]);
    if ((w = x_text_width(draw_win,argv[i])+10) > choice_w) choice_w = w;
   }
  p->action = F;
  return p;
 }


panel_item BASE_WINDOW::choice_mult_item(const char* s, int* x, int argc,
                                   const char** argv, panel_action_func F,
                                   const char* hlp)
{ choice_item(s,x,argc,argv,0,0,F,hlp); 
  panel_item p = Item[item_count-1];
  p->kind = Choice_Mult_Item;
  return p;
 }


panel_item BASE_WINDOW::bitmap_choice_item(const char* label, int *x, 
                                           int argc, int width, int height, 
                                                      unsigned char **argv, 
                                                      panel_action_func F,
                                                      const char* hlp)
{ int h = choice_h;
  if (height > th) h = choice_h + (height-th);
  panel_item p = new_panel_item(Bitmap_Choice_Item,label,x,h,hlp);  

  p->dat1=width+3;	   // width in dat1
  p->dat2=height+4;	   // height in dat2
  p->argc=argc;		   // number of bitmaps
  p->argv=new char*[argc]; // bitmaps
  for(int i=0; i < argc; i++) 
     p->argv[i] = x_create_bitmap(draw_win,width,height,argv[i]);
  p->action = F;
  return p;
}


panel_item BASE_WINDOW::bitmap_choice_mult_item(const char* label, int *x,
                                          int argc, int width, int height, 
                                          unsigned char **argv, 
                                          panel_action_func F,
                                          const char* hlp)
{ 
  bitmap_choice_item(label,x,argc,width,height,argv,F,hlp); 
  panel_item p = Item[item_count-1];
  p->kind = Bitmap_Choice_Mult_Item;
  return p;
 }



panel_item BASE_WINDOW::bool_item(const char* s, bool* x, panel_action_func F, 
                                                          const char* hlp)
{ panel_item p = new_panel_item(Bool_Item,s,x,color_h,hlp);
  p->action = F;
  return p;
}



panel_item BASE_WINDOW::color_item(const char* s, int* x, panel_action_func F,
                                                          const char* hlp)
{ panel_item p = new_panel_item(Color_Item,s,x,color_h,hlp);
  int w = 18*(color_h+2) - 2;
  if (w > slider_w)
  { slider_w  = w;
    string_w  = w;
    number_w  = w;
  }

/*
  button("",0,&color_menu);
  panel_item q = Item[item_count-1];
  p->menu_but = q;
  q->menu_but = q;
  but_count--;
*/

  p->action = F;
  return p;
}

 
panel_item BASE_WINDOW::pstyle_item(const char* label, point_style* x, 
                                        panel_action_func F, const char* hlp)
{
  unsigned char* pstyle_bits[7]; 
  pstyle_bits[0] = pixel_point_bits; 
  pstyle_bits[1] = cross_point_bits; 
  pstyle_bits[2] = plus_point_bits; 
  pstyle_bits[3] = circle_point_bits; 
  pstyle_bits[4] = disc_point_bits; 
  pstyle_bits[5] = rect_point_bits; 
  pstyle_bits[6] = box_point_bits;


  return bitmap_choice_item(label,(int*)x,pstyle_num,pstyle_width,
                                                     pstyle_height,
                                                     pstyle_bits,
                                                     F,hlp);
}



panel_item BASE_WINDOW::lstyle_item(const char* label, line_style* x, 
                                        panel_action_func F, const char* hlp)
{
  unsigned char* lstyle_bits[4]; 
  lstyle_bits[0] = solid_bits; 
  lstyle_bits[1] = dashed_bits; 
  lstyle_bits[2] = dotted_bits; 
  lstyle_bits[3] = dashed_dotted_bits;

  return bitmap_choice_item(label,(int*)x,lstyle_num,lstyle_width,
                                                     lstyle_height,
                                                     lstyle_bits,
                                                     F,hlp);
}


panel_item BASE_WINDOW::lwidth_item(const char* label, int* x, 
                                    panel_action_func F, const char* hlp)
{
  unsigned char* lwidth_bits[7];
  lwidth_bits[0] = width_zero_bits; 
  lwidth_bits[1] = width_one_bits; 
  lwidth_bits[2] = width_two_bits; 
  lwidth_bits[3] = width_three_bits; 
  lwidth_bits[4] = width_four_bits; 
  lwidth_bits[5] = width_five_bits; 
  lwidth_bits[6] = width_six_bits;

  return bitmap_choice_item(label,x,lwidth_num,lwidth_width, lwidth_height,
                                                             lwidth_bits,
                                                             F,hlp);
}




// buttons

int BASE_WINDOW::button(const char* s, int val, panel_action_func F, 
                                                const char* hlp)
{ panel_item p = new_panel_item(Button_Item,s,0,yskip,hlp);
  p->dat1 = (val == -1) ? but_count : val;
  p->action = F;
  x_set_button_font(draw_win);
  int w = x_text_width(draw_win,s) + x_text_height(draw_win,s);
  x_set_text_font(draw_win);
  if (w  > button_w) button_w = w;
  p->width = w;
  but_count++;
  return p->dat1;
 }


void BASE_WINDOW::set_focus_button(int but)
{ int i;
  for(i = 0; i<item_count; i++)
  { panel_item p = Item[i];
    if (p->kind != Button_Item || p->menu_but == p) continue;
    if (p->dat1 == but) break;
   }
  focus_button =  (i < item_count) ? Item[i] : 0;
}


int BASE_WINDOW::button(const char* s, int val, BASE_WINDOW* wp,const char* hlp)
{ 
  panel_item p = new_panel_item(Button_Item,s,0,yskip,hlp);

  p->dat1 = val;
  if (val == -1)
  { p->dat1 = but_count;
    p->dat2 = 2;  // second bit = 1: return result of sub-window wp
   }

  p->ref = wp;
  if (wp) wp->owner_item = p;
  x_set_button_font(draw_win);
  int w = x_text_width(draw_win,s) + 28;
  x_set_text_font(draw_win);
  if (w  > button_w) button_w = w;
  p->width = w;
  but_count++;
  return p->dat1;
 }


int BASE_WINDOW::button(int w, int h, unsigned char* bits, const char* s, 
                                                           int val, 
                                                           panel_action_func F,
                                                           const char*hlp)
{ if (h+1 > button_h) button_h = h+1;
  if (w+2 > button_w) button_w = w+2;

  if (hlp == 0) hlp = s;
  panel_item p = new_panel_item(Button_Item,s,0,yskip,hlp);
  p->argv = new char*[1];
  p->argv[0] = x_create_bitmap(draw_win,w,h,bits);
  p->argc = 1; // bitmap
  p->dat1 = (val == -1) ? but_count : val;
  p->action = F;
  p->width = w-1;
  but_count++;
  return p->dat1;
 }


int BASE_WINDOW::button(int w, int h, unsigned char* bits, const char* s, 
                                                           int val, 
                                                           BASE_WINDOW* wp,
                                                             const char*hlp)
{ if (h+2 > button_h) button_h = h+2;
  if (w+2 > button_w) button_w = w+2;

  if (hlp == 0) hlp = s;
  panel_item p = new_panel_item(Button_Item,s,0,yskip,hlp);
  p->argv = new char*[1];
  p->argv[0] = x_create_bitmap(draw_win,w,h,bits);
  p->argc = 1; // bitmap

  p->dat1 = val;
  if (val == -1)
  { p->dat1 = but_count;
    p->dat2 = 2;  // second bit = 1: return result of sub-window wp
   }

  p->ref = wp;
  if (wp) wp->owner_item = p;
  p->width = w-1;
  but_count++;
  return p->dat1;
 }


int BASE_WINDOW::button(char* pmap0, char* pmap1, const char* s, int val, 
                                                            panel_action_func F,
                                                            const char*hlp)
{ int w,h;
  x_pixrect_dimensions(pmap0,&w,&h);
  if (h+2 > button_h) button_h = h+2;
  if (w+2 > button_w) button_w = w+2;
  if (hlp == 0) hlp = s;
  panel_item p = new_panel_item(Button_Item,s,0,yskip,hlp);
  p->argv = new char*[2];
  p->argv[0] = pmap0;
  p->argv[1] = pmap1;
  p->argc = -1; // pmap
  p->dat1 = (val == -1) ? but_count : val;
  p->action = F;
  p->width = w-1;
  but_count++;
  return p->dat1;
 }


int BASE_WINDOW::button(char* pmap0, char* pmap1, const char* s, int val, 
                                                               BASE_WINDOW* wp,
                                                               const char*hlp)
{ int w,h;
  x_pixrect_dimensions(pmap0,&w,&h);
  if (h+2 > button_h) button_h = h+2;
  if (w+2 > button_w) button_w = w+2;
  if (hlp == 0) hlp = s;
  panel_item p = new_panel_item(Button_Item,s,0,yskip,hlp);
  p->argv = new char*[2];
  p->argv[0] = pmap0;
  p->argv[1] = pmap1;
  p->argc = -1; // pmap

  p->dat1 = val;
  if (val == -1)
  { p->dat1 = but_count;
    p->dat2 = 2;  // second bit = 1: return result of sub-window wp
   }

  p->ref = wp;
  if (wp) wp->owner_item = p;
  p->width = w-1;
  but_count++;
  return p->dat1;
 }




int BASE_WINDOW::menu_button(const char* s, int val, BASE_WINDOW* wp, 
                                                     const char* hlp)
{ panel_item p = new_panel_item(Button_Item,s,0,yskip,hlp);

  p->dat1 = val;
  p->dat2 = 1;  // first bit = 1:  menu (bar) button

  if (val == -1)
  { p->dat1 = but_count;
    p->dat2 += 2;  // second bit = 1: return result of sub-window wp
   }

  p->ref = wp;
  if (wp) wp->owner_item = p;
  x_set_button_font(draw_win);
  int w = x_text_width(draw_win,s) + 10;
  x_set_text_font(draw_win);

//if (w  > button_w) button_w = w;

  p->width = w;
  has_menu_bar = 1;
  but_count++;
  return p->dat1;
 }



//------------------------------------------------------------------------
// drawing buttons & items
//------------------------------------------------------------------------


void BASE_WINDOW::draw_label(panel_item it)
{ char* str = it->label_str;
  if (strlen(str) == 0) return;
  if (str[0] == '!') str++;
  int x = it->xcoord;
  int y = it->ycoord;

  x_set_button_font(draw_win);
//int txtw = x_text_width(draw_win,str);
  int txth = x_text_height(draw_win,str);

  x_set_color(draw_win,panel_bg_color); 
  x_box(draw_win,x,y,x+label_w,y+it->height);

  x_set_color(draw_win,it->enabled ? black : disable_color);

  int dy = it->height -txth;
  int yt = y + dy/2 + dy%2;

  x_text(draw_win,x,yt,str);

  x_set_text_font(draw_win);
}


void BASE_WINDOW::draw_box_with_shadow(int x1,int y1,int x2,int y2,int c,int w)
{ // below and right of box(x1,y1,x2,y2)
  x_set_color(draw_win,shadow_color);
  x_box(draw_win,x1+3,y1+3,x2+w,y2+w);
  x_set_color(draw_win,c);
  x_box(draw_win,x1,y1,x2,y2);
  x_set_color(draw_win,black);
  x_rect(draw_win,x1,y1,x2,y2);
 }


void BASE_WINDOW::draw_d3_box(int x1,int y1,int x2,int y2, int pressed, 
                                                           int enabled)
{
  int save_color = press_color;

  if (press_color == invisible) press_color = panel_bg_color;

  if (pressed && enabled) 
     x_set_color(draw_win,press_color);
  else
     x_set_color(draw_win,panel_bg_color);

  x_box(draw_win,x1,y1,x2+1,y2+1);

  int c1 = (pressed) ? black : white;
  int c2 = (pressed) ? shadow_color : panel_bg_color;
  int c3 = (pressed) ? white : shadow_color;
  int c4 = (pressed) ? press_color : black;

  if (!enabled) 
  { //c1 = white;
    c1 = panel_bg_color;
    c2 = panel_bg_color;
    //c3 = white;
    c3 = disable_color;
    c4 = disable_color;
   }

  x_set_color(draw_win,c1);
  x_line(draw_win,x1,y1,x2+(pressed?2:1),y1);		// outer top corner
  x_line(draw_win,x1,y1+1,x1,y2+(pressed?2:1));		// outer left corner

  x_set_color(draw_win,c2);
  x_line(draw_win,x1+1,y1+1,x2+(pressed?1:0),y1+1);	// inner top
  x_line(draw_win,x1+1,y1+2,x1+1,y2+(pressed?1:0));	// inner left

  x_set_color(draw_win,c3);
  x_line(draw_win,x1+(pressed?1:0),y2+1,x2+1,y2+1);	// outer bottom 
  x_line(draw_win,x2+1,y1+(pressed?1:0),x2+1,y2+2);	// outer right

  x_set_color(draw_win,c4);
  x_line(draw_win,x1+(pressed?2:1),y2,x2,y2);		// inner bottom
  x_line(draw_win,x2,y1+(pressed?2:1),x2,y2+1); 	// inner right

  press_color = save_color;
}


int BASE_WINDOW::draw_text_item(panel_item it, int win)
{ 
  if (win)
  { x_set_color(draw_win,panel_bg_color);
    x_box(draw_win,xoff,it->ycoord,panel_width,it->ycoord+it->height+1);
   }

  int y= format_text(it->argc,it->argv,xoff,panel_width-xoff,
                     it->ycoord, ytskip,win) + 2;
  return y - it->ycoord + 2;
}


void BASE_WINDOW::draw_string_item(panel_item it, const char* s)
{ 
  if (it->kind == String_Item || it->kind == String_Menu_Item)  
      it->data_str = access_str(it->ref);

  int active = it == act_str_item;

  int x1 = it->xcoord + label_w;
  int y1 = it->ycoord;
  int x2 = x1 + string_w;
  int y2 = y1 + string_h;
  int yt = y1 + 1;

  int use_buf = x_window_opened(draw_win) && buf_level == 0 && it->index >= 0;

  //int use_buf = x_window_opened(draw_win) && buf_level == 0;

  if (it->index < 0)
  { x_set_color(draw_win,grey1);
    x_box(draw_win,x1-1,y1-1,x2+1,y2+1);
    x_set_color(draw_win,black);
    x_rect(draw_win,x1-1,y1-1,x2+1,y2+1);
   }

  if (use_buf)
  { x_start_buffering(draw_win);
    x_set_color(draw_win,panel_bg_color);
    x_box(draw_win,x1-label_w-1,y1-2,x2+1,y2+2);
   }

  draw_label(it);

  
  if (s == 0)  s = it->data_str;

  if (it->index >=0)
  { 
    if (it->kind != String_Menu_Item)
    { int xr = x1-8;
      int yr = it->ycoord+it->height/2;
  
      if (it->dat1 > 0)
        draw_left_menu_button(xr,yr,-1,it->enabled);
      else
       { x_set_color(draw_win,panel_bg_color);
         x_box(draw_win,xr-6,yr-8,xr+6,yr+8);
        }
  
      xr = x2+8;

      int slen = strlen(s);
  
      if (slen > it->dat2)
        draw_right_menu_button(xr,yr,-1,it->enabled);
      else
       { x_set_color(draw_win,panel_bg_color);
         x_box(draw_win,xr-6,yr-8,xr+6,yr+8);
        }
     }
  
    int save = press_color;
    //press_color = (active) ? blue : ivory;
    press_color = ivory;
    draw_d3_box(x1,y1,x2,y2,1,it->enabled);
    press_color = save;
  }


  // write text into box

  while (it->offset > it->dat2)
  { it->dat1++;
    it->dat2++;
   }

  s += it->dat1;

  //int xtoff = active ? 5 : 4;
  //int ytoff = active ? 1 : 0;

  int xtoff = 4;
  int ytoff = 1;

  x_set_fixed_font(draw_win);
  x_set_color(draw_win, black);

  if (!it->enabled)
  { x_set_color(draw_win,white);
    x_text(draw_win,x1+xtoff+1,yt+ytoff+1,s,it->dat2-it->dat1);
    x_set_color(draw_win, disable_color);
   }
  x_text(draw_win,x1+xtoff,yt+ytoff,s,it->dat2-it->dat1);

  x_set_text_font(draw_win);

  if (active)
  { // draw cursor
    int xc = x1 + xtoff + tw*(it->offset - it->dat1);
    x_line(draw_win,xc-3,y2-0,xc-3,y2+1);
    x_line(draw_win,xc-2,y2-1,xc-2,y2+1);
    x_line(draw_win,xc-1,y2-2,xc-1,y2+1);
    x_line(draw_win,xc+0,y2-3,xc+0,y2+1);
    x_line(draw_win,xc+1,y2-2,xc+1,y2+1);
    x_line(draw_win,xc+2,y2-1,xc+2,y2+1);
    x_line(draw_win,xc+3,y2-0,xc+3,y2+1);
/*
    x_line(draw_win,xc-4,y2-0,xc-4,y2+1);
    x_line(draw_win,xc-3,y2-1,xc-3,y2+1);
    x_line(draw_win,xc-2,y2-2,xc-2,y2+1);
    x_line(draw_win,xc-1,y2-3,xc-1,y2+1);
    x_line(draw_win,xc+0,y2-4,xc+0,y2+1);
    x_line(draw_win,xc+1,y2-3,xc+1,y2+1);
    x_line(draw_win,xc+2,y2-2,xc+2,y2+1);
    x_line(draw_win,xc+3,y2-1,xc+3,y2+1);
    x_line(draw_win,xc+4,y2-0,xc+4,y2+1);
*/
   }

  if (use_buf)
  { if (it->kind == String_Menu_Item)
    { x_set_button_font(draw_win);
      int xl1 = it->xcoord;
      int xl2 = xl1 + x_text_width(draw_win,it->label_str);
      x_flush_buffer(draw_win,xl1,y1-1,xl2,y2+1);
      x_flush_buffer(draw_win,x1,y1-1,x2+1,y2+1);
     }
    else
      if (it->index < 0)
         x_flush_buffer(draw_win,x1-1,y1-1,x2+1,y2+1);
      else
         x_flush_buffer(draw_win,x1-label_w-1,y1-1,x2+1,y2+1);
    x_stop_buffering(draw_win);
   }
}



void BASE_WINDOW::activate_string_item(panel_item it, int x)
{ 
  if (!it->enabled) return;

  int sl = strlen(it->data_str);
  int x0 = it->xcoord + label_w;
  int  j = it->dat1 + (x-x0)/tw;

  if (x > x0 && j <= it->dat2)
  { if (j > sl) j = sl;
    it->offset = j;
   }


  if (act_str_item && act_str_item != it)
  { panel_item it0 = act_str_item;
    act_str_item = it;
    if (it0->kind == Int_Item)
        { delete[] it0->data_str;
          it0->data_str = make_string(*(int*)it0->ref);
         }
    if (it0->kind == Float_Item)
        { delete[] it0->data_str;
          it0->data_str = make_string(*(double*)it0->ref);
         }
    draw_string_item(it0);
/*
    if (it0->str_action)
    { clipping(2);
      x_set_text_font(draw_win);
      call_window = this;
      call_item = it0;
      it0->str_action(it0->data_str);
      clipping(1);
     }
*/
   }

  draw_string_item(it);
}
   
   

bool BASE_WINDOW::panel_text_edit(panel_item it)
{ 
  int len = 2*strlen(it->data_str);
  if (len < 256) len = 256;

  char*  str = new char[len];

  strcpy(str,it->data_str);

  int k,val;
  int xc,yc;
  int w;
  unsigned long t;

  do k = x_get_next_event(&w,&val,&xc,&yc,&t);
  while (w != draw_win || (k != key_press_event && k != button_press_event));

  if (k == button_press_event) 
  { if (it->index < 0)
    { int x1 = it->xcoord + label_w;
      int y1 = it->ycoord;
      int x2 = x1 + string_w;
      int y2 = y1 + string_h;
      bool res = false;
      if (yc > y1 && yc < y2 && xc > x1 && xc < x2) 
      { activate_string_item(it,xc);
        res = true;
       }
      delete[] str;
      return res;
     }
    x_put_back_event();
    delete[] str;
    return false;
   }

  if (val == KEY_RETURN || val == KEY_UP || val == KEY_DOWN)
  { if (it->index < 0) return false;
    if (it->kind == String_Menu_Item && (it->str_action || it->str_action_ptr))
    { clipping(2);
      x_set_text_font(draw_win);
      call_window = this;
      call_item = it;
      
      if (it->str_action_ptr) {
        it->str_action_ptr->set_params((window*)this,str);
        it->str_action_ptr->operator()();
      }
      else it->str_action(str);
      
      clipping(1);
     }
    x_put_back_event();
    delete[] str;
    return false;
   }


  if (it->kind != String_Item && it->kind != String_Menu_Item) 
   delete[] it->data_str;


  char c    = (char)val;
  int  j    = it->offset;
  int  strl = strlen(str);

  if (isprint(c))
  { for(int i=strl; i>=j; i--) str[i+1] = str[i];
    str[j]=c;
    j++;
    if (j > it->dat2)
    { it->dat1++;
      it->dat2++;
     }
   }
  else
   switch (c) {

   case KEY_BACKSPACE:    
                   if (j > 0) 
                   { for(int i=j; i<=strl; i++) str[i-1] = str[i];
                     if (--j < it->dat1)
                     { it->dat1--;
                       it->dat2--;
                      }
                    }
                   break;

   case KEY_LEFT:  if (j == 0) break;
                   if (--j < it->dat1) 
                   { it->dat1--;
                     it->dat2--;
                   }
                   break;

   case KEY_RIGHT: if (j == strl) break;
                   if (++j > it->dat2)
                   { it->dat1++;
                     it->dat2++;
                   }
                   break;

   case KEY_HOME:  j = 0;
                   it->dat2 -= it->dat1;
                   it->dat1 = 0;
                   break;

   case KEY_END:   j = strl;
                   if (it->dat2 < j) 
                   { it->dat1 += (j-it->dat2);
                     it->dat2 = j;
                   }
                   break;
   }

 it->offset = j;
 draw_string_item(it,str);

 it->data_str = string_dup(str);


 if ((it->str_action || it->str_action_ptr) && /* it->kind == String_Item && */
     (isprint(c) || c == KEY_BACKSPACE))
 { clipping(2);
   x_set_text_font(draw_win);
   call_window = this;
   call_item = it;
   
   if (it->str_action_ptr) {
     it->str_action_ptr->set_params((window*)this,str);
     it->str_action_ptr->operator()();
   }
   else it->str_action(str);
   
   clipping(1);
  }

 delete[] str;
 return true;
}



void BASE_WINDOW::draw_choice_item(panel_item it, int val)
{ draw_label(it);

  int x = it->xcoord + label_w;
  int y = it->ycoord;
  int n = it->argc;
  int c;

  if (it->kind == Choice_Mult_Item) 
     c = val;
  else
    { c = (val - it->offset)/it->dat2;
      if (val != it->offset + c*it->dat2) c = -1;
     }


  for(int j=0; j<n; j++)
  { bool pressed = (j == c); 

    if (it->kind == Choice_Mult_Item) pressed = ((c&(1 << j)) != 0);

    draw_d3_box(x,y,x+choice_w-3,y+choice_h,pressed,it->enabled);
    int xx = x + choice_w/2 - (pressed ? 1 : 2);
    int yy = y + choice_h/2 - (pressed ? 0 : 1);

    yy++;

    x_set_color(draw_win,black);
    x_set_button_font(draw_win);

    if (!it->enabled)
    { x_set_color(draw_win,white);
      x_ctext(draw_win,xx+1,yy+1,it->argv[j]);
      x_set_color(draw_win,disable_color);
     }
    x_ctext(draw_win,xx,yy,it->argv[j]);

    x_set_text_font(draw_win);

    x += choice_w;
   }

}




void BASE_WINDOW::draw_bitmap_choice_item(panel_item it) 
{
  draw_label(it);

  int n = it->argc;
  int w = it->dat1;
  int h = it->dat2;
  int x = it->xcoord + label_w;
  int y = it->ycoord;
  int c = *(int*)it->ref;

  for(int j=0; j<n; j++) 
  { bool pressed = (j == c); 
    if (it->kind == Bitmap_Choice_Mult_Item) pressed = ((c&(1 << j)) != 0);
    draw_d3_box(x,y+1,x+w,y+h,pressed,it->enabled);
    if (it->enabled)
       x_set_color(draw_win,pressed ? bitmap_color1 : bitmap_color0);
    else
       x_set_color(draw_win, disable_color);
    int dx = pressed ? 3 : 2;
    int dy = pressed ? 1 : 2;
    x_insert_bitmap(draw_win,x+dx,y+h-dy,it->argv[j]);
    x += (w+3);
   }
}



void BASE_WINDOW::change_bitmap_choice_item(panel_item it, int j)
{

  int old = *(int*)it->ref;

  int n = it->argc;
  int w = it->dat1;
  int h = it->dat2;
  int y = it->ycoord;

  int c = old;

  if (it->kind == Bitmap_Choice_Mult_Item)
    c ^= (1 << j);
  else
    c = j;

  if (c == old) return;


  int x = it->xcoord + label_w + j*(w+3);

  if (it->kind == Bitmap_Choice_Mult_Item)
  { int pressed =  (c & (1<<j));
    draw_d3_box(x,y+1,x+w,y+h,pressed,it->enabled);
    if (it->enabled)
       x_set_color(draw_win,pressed ? bitmap_color1 : bitmap_color0);
    else
       x_set_color(draw_win, disable_color);
    if (pressed)
       x_insert_bitmap(draw_win,x+3,y+h-1,it->argv[j]);
    else
       x_insert_bitmap(draw_win,x+2,y+h-2,it->argv[j]);
   }
  else
  { draw_d3_box(x,y+1,x+w,y+h,1,it->enabled);
    x_set_color(draw_win,it->enabled ? bitmap_color1 : disable_color);
    x_insert_bitmap(draw_win,x+3,y+h-1,it->argv[j]);

    if (old >= 0 && old < n)
    { int x = it->xcoord + label_w + old * (w+3);
      draw_d3_box(x,y+1,x+w,y+h,0,it->enabled);
      x_set_color(draw_win,it->enabled ? bitmap_color0 : disable_color);
      x_insert_bitmap(draw_win,x+2,y+h-2,it->argv[old]);
     }
  }
}


void BASE_WINDOW::draw_color_button(int xt, int yt, int w, int col, int pressed,
                                                                    int enabled)
{ int x1  = xt - w;
  int y1  = yt - w;
  int x2  = xt + w;
  int y2  = yt + w;

  int save = press_color;
  press_color = col;
  draw_d3_box(x1,y1,x2,y2,pressed,enabled);
  press_color = save;

  if (col == invisible || (pressed && enabled)) return;

  if (col == panel_bg_color)
    { x_set_color(draw_win,press_color);
      x_rect(draw_win,xt-2,yt-2,xt+2,yt+2);
     }
   else
    { x_set_color(draw_win, enabled ? col : disable_color);
      x_box(draw_win,xt-2,yt-2,xt+2,yt+2);
     }

  x_set_color(draw_win,black);
}



void BASE_WINDOW::draw_color_item(panel_item it)
{ 
  draw_label(it);
  int dx = color_h+2;
  int yt = it->ycoord + it->height/2;
  int xt = it->xcoord + label_w +color_h/2;
  int c = *(int*)it->ref;
  for(int j=-1; j < 17; j++)
  { draw_color_button(xt,yt,color_h/2,j,j==c,it->enabled);
    xt += dx;
   }
}


void BASE_WINDOW::change_color_item(panel_item it, int j)
{  
  int yt = it->ycoord + it->height/2;
  int xt = it->xcoord + label_w + color_h/2;
  int c = *(int*)it->ref;

  if (j == c) return; 

  draw_color_button(xt + (j+1)*(color_h+2),yt,color_h/2,j,1,it->enabled);

  *(int*)it->ref = j;

  if (c >= -1 && c <= 16) 
    draw_color_button(xt + (c+1)*(color_h+2),yt,color_h/2,c,0,it->enabled);
}


void BASE_WINDOW::draw_bool_item(panel_item it)
{  
  draw_label(it);
  int yt = it->ycoord + it->height/2;
  int xt = it->xcoord + label_w + color_h/2;
  int c  = *(bool*)it->ref;
  int w  = color_h/2;
  draw_d3_box(xt-w,yt-w,xt+w,yt+w,c,it->enabled);
}



void BASE_WINDOW::draw_slider_item(panel_item it, int x)
{ 
  int x0 = it->xcoord + label_w;
  int y0 = it->ycoord;

  int x1 = x0 + slider_w;
  int y1 = y0 + slider_h;

  int first_time = 0;

  int mi = it->dat1;
  int ma = it->dat2;
  float d = (ma > mi) ? float(slider_w)/(ma-mi+1) : 1;


  int val0 = *(int*)it->ref;

  bool undef = false;

  if (x == 0) // first time
  { first_time = 1;
    draw_label(it);
    x = x0  + (int)(d * (val0 - mi + 0.5)) + 1;
    it->offset = x;
    if (val0 < mi || val0 > ma) undef = true;
   }

  if (x < x0) x = x0;
  if (x > x1) x = x1;

  int val  = mi + int((x-x0-d/2)/d);

  if (!undef) *(int*)it->ref = val;

  if (it->label_str[0] != '!')
  { char text[16];

    if (undef)
		CGAL_CLIB_STD::sprintf(text,"   ?");
    else
        CGAL_CLIB_STD::sprintf(text,"%4d",val);

    x_set_color(draw_win,panel_bg_color);
    x_box(draw_win,x0-sl_num_w,y0,x0-1,y1);

    x_set_fixed_font(draw_win);
    int dy = slider_h - x_text_height(draw_win,"0");
    int yt = y0 + dy/2 + dy%2;
    x_set_color(draw_win,it->enabled ? black : disable_color);
    x_text(draw_win,x0-sl_num_w,yt,text);
    x_set_text_font(draw_win);
  }

  // slider

  int slid_x = int(x+0.5) + 1;

  if (slid_x < x1-1) 
     draw_d3_box(slid_x,y0,x1+1,y1,1,it->enabled);

  if (slid_x > x0+2)
     draw_d3_box(x0,y0,slid_x,y1,0,it->enabled);


  it->offset = x;

  x_set_color(draw_win,black);

  if (first_time == 0)
    if ((it->action || it->action_ptr) && val != val0) 
    { clipping(2);
      call_window = this;
      call_item = it;
      
      if (it->action_ptr) {
        it->action_ptr->set_params((window*)this,NULL,val);
        it->action_ptr->operator()();
      }
      else it->action(val);
      
      clipping(1);
     }

}

void BASE_WINDOW::draw_separator(panel_item it)
{ x_set_color(draw_win,shadow_color);
  x_line(draw_win,1,it->ycoord,panel_width-1,it->ycoord);
  x_set_color(draw_win,white);
  x_line(draw_win,1,it->ycoord+1,panel_width-1,it->ycoord+1);
}


void BASE_WINDOW::draw_button(panel_item it, int pressed)
{ 
   //if (it == active_button) pressed = 1;

   int enab = it->enabled;
   if (!enab) pressed = 0;

   bool is_menu = (panel_height >= window_height && item_count == but_count);


   if (it->dat2 % 2)
   { draw_menu_button(it,pressed);
     return;
    }


   if (has_menu_bar && (it->argc != 0)) // before: == -1
   { int bar = has_menu_bar;
     has_menu_bar = 0;
     draw_menu_button(it,pressed);
     has_menu_bar = bar;
     return;
   }


   if (it->menu_but == 0 && menu_button_style == 2) 
   { draw_menu_item(it,pressed);
     return;
    }

   char* s = it->label_str;

   int x1 = it->xcoord;
   int y1 = it->ycoord;
   int x2 = x1 + it->width;
   int y2 = y1 + it->height;
   int xt = x1 + it->width/2;
   int yt = y1 + it->height/2; //+ it->height%2;

   int dt = 1;

   x_set_color(draw_win, panel_bg_color);

/*
   if (!it->enabled)
      x_box(draw_win,x1,y1,x2+1,y2+1);
   else
*/
   { if (is_menu && menu_button_style == 1)
       { y1++; y2++;
         if (pressed) 
            draw_d3_box(x1+1,y1,x2,y2-2,0,it->enabled);
         else
            if (buts_per_line == 1)
              x_box(draw_win,x1+1,y1,x2,y2-1);
            else
              x_box(draw_win,x1+1,y1,x2+1,y2-1);
         dt = -1;
        }
     else
       if (it->menu_but == it && strlen(it->label_str) == 0) 
         { x_box(draw_win,x1,y1,x2+2,y2+2);
           if (pressed || it->ref == 0) 
             //draw_d3_box(x1,y1,x2,y2,pressed,it->enabled);
             draw_d3_box(x1,y1,x2,y2,pressed,it->ref != 0);
          }
       else
         draw_d3_box(x1,y1,x2,y2,pressed);
  
     if (!pressed && it == focus_button)
     { int r1,g1,b1,r2,g2,b2;
       x_get_rgb(grey1,&r1,&g1,&b1);
       x_get_rgb(grey2,&r2,&g2,&b2);
       int fcol = x_new_color((r1+r2)/2,(g1+g2)/2,(b1+b2)/2);
       if (fcol == -1) fcol = press_color;
       x_set_color(draw_win,fcol);
       x_box(draw_win,x1+1,y1+1,x2-1,y2-1);
      }
   }

   x_set_color(draw_win, (enab) ? black : disable_color);

   x_set_button_font(draw_win);

   int scut = it->shortcut;

   if (it->argc == 0) // no bitmap/pixmap button
     { if (it->ref != this) // insert text
       { int ch = x_text_height(draw_win,"H");
         if (center_button_label)
           { if (pressed) { xt += dt; yt += dt; }
             x_ctext(draw_win,xt,yt,s);
             if (scut >= 0) x_ctext_underline(draw_win,xt,yt,s,scut,scut);
            }
         else
           { xt = x1 + ch/2;
             yt = yt - ch/2;
             if (pressed) { xt += dt; yt += dt; }
             x_text(draw_win,xt,yt,s);
             if (scut >= 0) x_text_underline(draw_win,xt,yt,s,scut,scut);
            }
        }

       if (it->ref && it->enabled) // sub window
       { int d = (y2-y1)/2;
         int x = (pressed) ? dt : 0;
         if (dt == 0) pressed = !pressed;
         if (it->ref == this) 
         { 
            if (strcmp(it->label_str,"DOWN") == 0)
               draw_down_menu_button(x1+d+x+1,y1+d+x-1,pressed,enab);
            else
            if (strcmp(it->label_str,"UP") == 0)
               draw_up_menu_button(x1+d+x+1,y1+d+x+1,pressed,enab);
          }
         else
         if (buts_per_line == 1 && item_count == but_count) // sub menu
            draw_right_menu_button(x2-d+x,y1+d+x,-1,enab);
         else
            draw_down_menu_button(x1+d+1,y1+d-1,pressed,enab);
        }

     }
   else 
   { 
    if (it->argc == 1)
     { // bitmap
       if (pressed)
         x_insert_bitmap(draw_win,x1+3,y2,it->argv[0]);
       else
         x_insert_bitmap(draw_win,x1+2,y2-1,it->argv[0]);
      }

      if (it->argc == -1)
      { // pixmap
        if (pressed)
          x_insert_pixrect(draw_win,x1+3,y2,it->argv[1]);
        else
          x_insert_pixrect(draw_win,x1+2,y2-1,it->argv[0]);
      }
    }

 x_set_text_font(draw_win);
 x_set_color(draw_win,black);
 x_flush_display();
}


void BASE_WINDOW::draw_menu_button(panel_item it, int pressed)
{ 
  if (it == active_button) pressed = 1;

  char* s = it->label_str;

  int x1 = it->xcoord;
  int y1 = it->ycoord;

  int x2 = x1 + it->width;
  int y2 = y1 + it->height+1;

  int xt = x1 + it->width/2;
  int yt = y1 + it->height/2; //+ it->height%2;

  if (!pressed) 
    { x_set_color(draw_win,panel_bg_color);
      x_box(draw_win,x1-1,y1-1,x2+2,y2+2);
     }
  else
    draw_d3_box(x1,y1,x2,y2,0,it->enabled);


  x_set_color(draw_win,it->enabled ? black : disable_color);

  if (it->argc == 1) // bitmap
  {
    if (it->enabled)
       x_set_color(draw_win,pressed ? bitmap_color1 : bitmap_color0);

    if (pressed)
       x_insert_bitmap(draw_win,x1+2,y2-1,it->argv[0]);
    else
       x_insert_bitmap(draw_win,x1+3,y2,it->argv[0]);
   }
  else
  if (it->argc == -1) // pixmap
  { if (pressed)
       x_insert_pixrect(draw_win,x1+2,y2-1,it->argv[0]);
    else
       x_insert_pixrect(draw_win,x1+3,y2,it->argv[1]);
   }
  else
  { if (!pressed) { xt++; yt++; }
    if (s[0] == '<' && s[1] == '<') 
      draw_left_menu_button(xt,yt+1,0,it->enabled);
    else
    if (s[0] == '>' && s[1] == '>')
      draw_right_menu_button(xt,yt+1,0,it->enabled);
    else
    { x_set_button_font(draw_win);
      int scut = it->shortcut;
      x_ctext(draw_win,xt,yt,s);
      if (scut >= 0) x_ctext_underline(draw_win,xt,yt,s,scut,scut);
      x_set_text_font(draw_win);
     }
   }

  x_flush_display();
}



void BASE_WINDOW::draw_menu_item(panel_item it, int pressed)
{ 
   char* s = it->label_str;

   int x1 = it->xcoord;
   int y1 = it->ycoord;
   int x2 = x1 + it->width;
   int y2 = y1 + it->height+1;
   int xt = x1 + 4;

   //if (but_count <= menu_size) xt += button_h+1;

   if (it == act_str_item) pressed = true;

   //int col  =  (pressed) ? blue : ivory; 

   int col  =  (pressed) ? blue : white; 

   x_set_color(draw_win,col);
   x_box(draw_win,x1,y1,x2,y2);

   if (it->enabled)
    { x_set_color(draw_win,black);
      x_rect(draw_win,x1,y1,x2,y2);
      x_set_color(draw_win,text_color(col));
     }
   else 
     x_set_color(draw_win,disable_color);

   x_set_fixed_font(draw_win);
   int dy = it->height - x_text_height(draw_win,"H");
   int yt = y1 + dy/2 + dy%2;
   x_text(draw_win,xt,yt,s);
   x_set_color(draw_win,black);
   x_set_text_font(draw_win);

   x_flush_display();
}



void BASE_WINDOW::draw_up_menu_button(int x, int y, int pressed, int enabled)
{ 
/*
  if (pressed >= 0) 
    draw_d3_box(x-9,y-9,x+9,y+9,pressed,enabled);
*/

  y++;

  int X[4];
  int Y[4];
  X[0] = x;
  Y[0] = y - 6;
  X[1] = x + 4;
  Y[1] = y + 2;
  X[2] = x - 4;
  Y[2] = y + 2;
  X[3] = X[0];
  Y[3] = Y[0];

  //x_set_color(draw_win, (pressed>0) ? press_color : panel_bg_color);
  x_set_color(draw_win, panel_bg_color);
  x_fill_polygon(draw_win,4,X,Y);

  color c1 = shadow_color;
  color c2 = white;

  if (pressed)
  { c1 = white;
    c2 = shadow_color;
   }

  if (!enabled) c1 = c2 = disable_color;

  x_set_color(draw_win,c1);
  x_line(draw_win,X[2],Y[2],X[3],Y[3]);

  x_set_color(draw_win,c2);
  x_polyline(draw_win,3,X,Y);

  x_set_color(draw_win,black);
  x_flush_display();
}



void BASE_WINDOW::draw_down_menu_button(int x, int y, int pressed, int enabled)
{ 

/*
  if (pressed >= 0) 
    draw_d3_box(x-9,y-9,x+9,y+9,pressed,enabled);
*/

  int X[4];
  int Y[4];
  X[0] = x;
  Y[0] = y + 6;
  X[1] = x - 4;
  Y[1] = y - 2;
  X[2] = x + 4;
  Y[2] = y - 2;
  X[3] = X[0];
  Y[3] = Y[0];

  //x_set_color(draw_win, (pressed>0) ? press_color : panel_bg_color);

  x_set_color(draw_win, panel_bg_color);
  x_fill_polygon(draw_win,4,X,Y);

  color c1 = shadow_color;
  color c2 = white;

  if (pressed)
  { c1 = white;
    c2 = shadow_color;
   }

  if (!enabled) c1 = c2 = disable_color;

  X[0]--; Y[0]--;
  x_set_color(draw_win,c1);
  x_polyline(draw_win,3,X,Y);

  X[2]--;
  x_set_color(draw_win,c2);
  x_line(draw_win,X[2],Y[2],X[3],Y[3]);

  x_set_color(draw_win,black);
  x_flush_display();
}


void BASE_WINDOW::draw_right_menu_button(int x, int y, int pressed, int enabled)
{ 
/*
  if (pressed >= 0) 
    draw_d3_box(x-9,y-9,x+9,y+9,pressed,enabled);
  else
*/
  x++;

  int X[4];
  int Y[4];
  X[0] = x - 4;
  Y[0] = y + 4;
  X[1] = x + 4;
  Y[1] = y;
  X[2] = x - 4;
  Y[2] = y - 4;
  X[3] = X[0];
  Y[3] = Y[0];

  color c0 = (pressed>0) ? press_color : panel_bg_color;
  color c1 = white;
  color c2 = shadow_color;

  if (!enabled) c0 = c1 = c2 = disable_color;

  x_set_color(draw_win,c0);
  x_fill_polygon(draw_win,4,X,Y);

  if (enabled)
  { x_set_color(draw_win,c1);
    x_polyline(draw_win,3,X+1,Y+1);
    x_set_color(draw_win,c2);
    x_line(draw_win,X[0]+1,Y[0]-1,X[1],Y[1]);
   }

  x_set_color(draw_win,black);
  x_flush_display();
}


void BASE_WINDOW::draw_left_menu_button(int x, int y, int pressed, int enabled)
{ 
/*
  if (pressed >= 0) 
    draw_d3_box(x-9,y-9,x+9,y+9,pressed,enabled);
  else
*/

  x--;

  int X[4];
  int Y[4];
  X[0] = x + 4;
  Y[0] = y + 4;
  X[1] = x - 4;
  Y[1] = y;
  X[2] = x + 4;
  Y[2] = y - 4;
  X[3] = X[0];
  Y[3] = Y[0];


  color c0 = (pressed>0) ? press_color : panel_bg_color;
  color c1 = white;
  color c2 = shadow_color;

  if (!enabled) c0 = c1 = c2 = disable_color;

  x_set_color(draw_win,c0);
  x_fill_polygon(draw_win,4,X,Y);

  if (enabled)
  { x_set_color(draw_win,c1);
    x_polyline(draw_win,3,X,Y);
    x_set_color(draw_win,c2);
    x_line(draw_win,X[2],Y[2]+1,X[3],Y[3]);
   }

  x_set_color(draw_win,black);
  x_flush_display();
}




void BASE_WINDOW::make_menu_bar(int kind)
{ // collects all buttons inserted so far into a menu bar

  if (has_menu_bar) return;

  for(int i=0; i<item_count; i++)
  { panel_item it = Item[i];
    if (it->kind == Button_Item && it->menu_but != it)
    { it->dat2 |= 1; // goes into menu bar
#if defined(__win32__)
      it->height = button_h;
#else
      it->height = button_h + 1;
#endif
      if (it->argc == 0)
      { x_set_button_font(draw_win);
        it->width = x_text_width(draw_win,it->label_str) + 10;
        x_set_text_font(draw_win);
       }
      else
        it->width += 4;

      if (kind > 0)
      { int bh = it->height + yoff + 1;
        if (menu_bar_height < bh) menu_bar_height = bh;
       }
     }
   }

  if (kind > 0)
  { has_menu_bar = kind;
    button_w = 20; 
    button_d = th/2;
    button_h = th + 3;
   }
}


void BASE_WINDOW::place_panel_items()
{
  bool is_menu = (p_win && item_count == but_count);

  // adjust panel width

  if (is_menu)  // menu
    { xoff  = 0;
      yoff  = 0;
      yskip = button_h+2;
      panel_width = buts_per_line*(button_w+1) + 1;
      center_button_label = 0;

      if (menu_button_style == 1) 
      { button_h = string_h + 4;
        yskip = button_h+2;
       }

      if (menu_button_style == 2) 
      { if (but_count > menu_size) panel_width += button_h+1;
        panel_width--;
        yskip = button_h+1;
        yoff = -2;
       }
     }
  else
  { if (yskip < button_h+2) yskip = button_h+2;
    int w = x_text_width(draw_win,default_frame_label);
    if (panel_width == -1 && item_count == but_count) // only buttons
    { if (w > 300) w = 300;
      if (w > panel_width) panel_width = w;
     }
    //else if (panel_width < 0) panel_width = 50;

    if (has_menu_bar) xoff = 3;
   }

  if (buts_per_line > but_count) buts_per_line = but_count;


  // set width for all items and adjust panel width

  int max_w = 0;
  int max_choice_width = 0;
  panel_item max_choice_it = 0;

  int non_button_items = 0;

  int i;

  for(i = 0; i<item_count; i++)
  { panel_item it = Item[i];

    if (has_menu_bar && it->kind == Button_Item 
                     && (it->dat2 % 2)) continue;

    switch (it->kind)
    {
      case  Button_Item: 
             if (it->menu_but != it) 
                it->width  = button_w;
             break;

      case  Text_Item:
             { it->width = 0;
               if (it->argc > 0)
               { if ((it->argv[0][0]&0xf0) == 0x10) x_set_bold_font(draw_win);
                 if ((it->argv[0][0]&0xf0) == 0x20) x_set_fixed_font(draw_win);
                }
               for(int j=0; j<it->argc; j++) 
                   it->width += x_text_width(draw_win,it->argv[j]);
               x_set_text_font(draw_win);
               it->width += xoff;
               if (it->width > 400) it->width =  400;
               non_button_items = 1;
               break;
              }

      case Bool_Item:
             it->width = label_w + color_h;
             non_button_items = 1;
             break;

      case Choice_Item:
      case Choice_Mult_Item: 
           { int w = choice_w*it->argc - 3;
             if (w > max_choice_width) 
             { max_choice_it = it;
               max_choice_width = w; 
              }
             it->width = label_w + w;
             non_button_items = 1;
             break;
            }

      case Color_Item:
             it->width = label_w + 18*(color_h+2)-2;
             non_button_items = 1;
             break;

      case Bitmap_Choice_Item:
      case Bitmap_Choice_Mult_Item:
             it->width = label_w + it->argc * (it->dat1+3) - 3;
             non_button_items = 1;
             break;

      case Slider_Item:
             it->width = label_w + slider_w;
             non_button_items = 1;
             break;

      case String_Item:
      case String_Menu_Item:
             it->width = label_w + string_w;
             non_button_items = 1;
             break;

      case Float_Item:
      case Int_Item:
             it->width = label_w + number_w;
             non_button_items = 1;
             break;
  
     }


    int w = it->width;
    if (it->kind == String_Item) w += xoff;

    if (w > max_w) max_w = w;

  }


  int w = max_w + 2*xoff;
  if (w > panel_width) panel_width = w;

  // adjust width of string and slider items

  if (slider_w < max_choice_width && slider_w + choice_w > max_choice_width)
  { slider_w = max_choice_width;
    string_w = max_choice_width;
    number_w = max_choice_width;
   }

  for(i = 0; i<item_count; i++)
  { panel_item it = Item[i];
    switch (it->kind) {
      case String_Item:
      case String_Menu_Item:
              it->dat2 = string_w/tw - 1;
      case Slider_Item:
      case Float_Item:
      case Int_Item:
              it->width = label_w + string_w;
              break;
      }
   }


  if (slider_w > max_choice_width && slider_w - choice_w < max_choice_width)
  { int dw = (slider_w - max_choice_width)/max_choice_it->argc;
    choice_w += dw;
    for(i = 0; i<item_count; i++)
    { panel_item it = Item[i];
      switch (it->kind) {
        case Choice_Item:
        case Choice_Mult_Item:
                it->width = label_w + choice_w*it->argc - 3;
                break;
      
       }
     }
   }


  //  button layout

  int bxoff = xoff;      // left and right button boundary space
  int bskip = button_d;  // horizontal space between buttons

  if (is_menu)
    bskip = 1;
  else
    if (has_menu_bar)
      bxoff = xoff;
    else
    { int bw = button_w + bskip;

      int w = but_count*bw+2*xoff-bskip;
      if (buts_per_line == 0 /* && non_button_items */&& w < panel_width) 
          buts_per_line = but_count;

      if (buts_per_line) // increase panel width if necessary 
      { w = buts_per_line*bw+2*xoff-bskip;
        if (w > panel_width) panel_width = w;
        // center buttons
        bxoff = (panel_width - (buts_per_line * bw - bskip))/2;
       }
     }



  // place items

  // menu bar


  int bar_y = yoff-1;

  if (button_h < string_h+1) button_h = string_h+1;

  // int bspace = (menu_bar_height-button_h-2)/2;

  if (has_menu_bar == 1) 
  { // use user defined ordering of buttons
    int x = xoff + 4;
    int i;
    for(i=0; i<item_count; i++)
    { panel_item it = Item[i];
      if (it->kind == Button_Item && (it->dat2 % 2))
      { it->xcoord = x;
        it->ycoord = bar_y;
        it->ycoord += (menu_bar_height-it->height)/2;
        x += it->width + bskip;
      }
     }
   }


  if (has_menu_bar == 2) 
  { // place menu buttons left of normal buttons
    // place non-menu buttons at right border
    int x1 = panel_width;
    for(i=item_count-1; i>=0; i--)
    { panel_item it = Item[i];
      if (it->kind == Button_Item && (it->dat2 % 2) && it->ref==0)
      { x1 -= (it->width + bskip);
        it->xcoord = x1;
        it->ycoord = bar_y;
        it->ycoord += (menu_bar_height-it->height)/2;
       }
     }

    // menu buttons (place at left border)

    int x = xoff+4;

    for(i=0; i<item_count; i++)
    { panel_item it = Item[i];
      if (it->kind == Button_Item && (it->dat2 % 2) && it->ref)
      { it->xcoord = x;
        it->ycoord = bar_y;
        it->ycoord += (menu_bar_height-it->height)/2;
        x += it->width;
        if (x > x1) it->xcoord = panel_width+1; // overlap
        x += bskip;
      }
     }

  }


  int cur_xcoord = panel_width+1;
  int cur_ycoord = 2*yoff;
  int cur_skip   = 0;

  if (item_count > 0 && Item[0]->kind == Text_Item) cur_ycoord = yoff;

  if (cur_ycoord < 0) cur_ycoord = -1;


  int non_bar_items = 0;
  int non_bar_buttons = 0;

  for(i = 0; i<item_count; i++)
  { panel_item it = Item[i];
    if (it->kind != Button_Item || (it->dat2 % 2) == 0 || !has_menu_bar) 
    {  non_bar_items++;
       if (it->kind == Button_Item) non_bar_buttons++;
     }
   }

  if (has_menu_bar) 
  { cur_ycoord = menu_bar_height + yoff/2;
    if (non_bar_buttons) cur_ycoord += 3;
    if (item_count > but_count) cur_ycoord += 2*yoff;
    if (non_bar_items || but_count == 0) cur_ycoord += yoff;
   }


  for(i = 0; i<item_count; i++)
  { panel_item it = Item[i];

    if (it->kind == Button_Item)  continue;
    if (it->kind == Space_Item)   continue;
    if (it->kind == Fill_Item)    continue;
    if (it->kind == Newline_Item) continue;


    if ( it ->kind == Text_Item ||
        (cur_xcoord > xoff && cur_xcoord + it->width > panel_width) )
    { 
      cur_xcoord = xoff;
      cur_ycoord += cur_skip;

      switch (it->kind) { 
        case Bitmap_Choice_Item:
        case Bitmap_Choice_Mult_Item:
              cur_ycoord += (it->dat2 - yskip)/2;
              cur_skip =  (yskip + it->dat2)/2;
              if (cur_skip < yskip) cur_skip = yskip;
              break;
        case Text_Item: 
              if (cur_skip == yskip) cur_ycoord += th/2;
              cur_skip = draw_text_item(it,0);
              break;
        default: 
              cur_skip = yskip;
              break;
        }
     }

    /*if (it->xcoord == 0)*/ it->xcoord = cur_xcoord;
    /*if (it->ycoord == 0)*/ it->ycoord = cur_ycoord + (yskip - it->height)/2;

    if (it->kind == Bool_Item)
      cur_xcoord += it->width + 2*xoff;
    else
      cur_xcoord += max_w+4*xoff;

    if (it->menu_but)
    { panel_item p = it->menu_but;
      p->xcoord = it->xcoord + label_w - p->width - 3;
      p->ycoord = it->ycoord;
     }
   }

  if (non_bar_items)
     cur_ycoord += cur_skip;


  // place buttons

  int x = panel_width + 1;
  int y = cur_ycoord - yskip;

  if (is_menu)
  { x = 0;
    y = 0;
   }

  if (non_button_items) y += (yoff-1);


  for(i=0; i<item_count; i++)
  { panel_item it = Item[i];

    switch (it->kind) {

     case Newline_Item: x = bxoff;
                        y += yskip;
                        it->ycoord = y;
                        break;

     case Space_Item:   if (it->width > 0) x += it->width;
                        break;


     case Fill_Item:  { if (x > bxoff && x+it->width > panel_width-bxoff) 
                        { x = bxoff;
                          y += yskip;
                         }
                        int x1 = panel_width - button_w - bxoff;
                        if (x1 > x) x = x1;
                        break;
                       }

     case Button_Item: { if (it->menu_but == it || 
                             ((it->dat2 % 2) && has_menu_bar)) continue;
                         if (x > bxoff && x+it->width > panel_width-bxoff) 
                         { x = bxoff;
                           y += yskip;
                          }
                         it->height = button_h;
                         it->xcoord = x;
                         it->ycoord = y + (yskip - it->height)/2;
                         x += (it->width + bskip);
                         break;
                       }
     }
  }

  y += yskip; 

  panel_height = y + 1;

  if (is_menu == false) panel_height += 4;

  if (menu_button_style == 2)  
    { // use scrollbar (if but_count > menu_size)
      if (but_count > menu_size) 
      { panel_height = menu_size*yskip + 1;
        if (act_str_item)
        { // center active item
          int j = (menu_size/2 - act_str_item->index);
          if (menu_size%2 == 0) j--;
          if (j > 0) j = 0;
          if (j < menu_size-item_count) j = menu_size-item_count;
          for(i=0; i< item_count; i++) Item[i]->ycoord += j*yskip;
         }
       }
     }
}



void BASE_WINDOW::draw_panel_items()
{ drawing_mode save_mode  = x_set_mode(draw_win,src_mode);
  line_style   save_style = x_set_line_style(draw_win,solid);
  int          save_lw    = x_set_line_width(draw_win,1);

  bool is_menu = (panel_height >= window_height && item_count == but_count);

  if (is_menu && menu_button_style == 1)
   { //draw_d3_box(0,0,window_width-2,window_height-2,0,1);
     x_set_bg_color(draw_win,shadow_color);
     x_clear_window(draw_win,0,0,window_width-1,window_height-1);
     x_set_color(draw_win,panel_bg_color);
     x_box(draw_win,1,1,window_width-2,window_height-2);
     x_set_color(draw_win,white);
     x_line(draw_win,0,0,window_width,0);
     x_line(draw_win,0,0,0,window_height);
    }
  else
   { x_set_color(draw_win,panel_bg_color);
     x_box(draw_win,0,0,window_width-1,panel_height-1);
    }


  if (window_height > panel_height)
  { 
    //int y = yoff-1;
    int y = yoff-2;

    x_set_color(draw_win,shadow_color);
    x_line(draw_win,0,panel_height-1,window_width-1,panel_height-1);
    x_line(draw_win,window_width-2,y,window_width-2,panel_height);
    //x_line(draw_win,0,y-1,window_width-2,y-1);

    x_set_color(draw_win,white);
    x_line(draw_win,0,y,0,panel_height);
    x_line(draw_win,0,y,window_width-1,y);
    //x_line(draw_win,window_width-1,y,window_width-1,panel_height);


    if (has_menu_bar)
    { int non_bar_items = 0;
      for(int i = 0; i<item_count; i++)
      { panel_item it = Item[i];
        if (it->kind != Button_Item || (it->dat2 % 2) == 0) non_bar_items++;
       }


      if (non_bar_items)
      { int y = yoff + menu_bar_height + 1;
        x_set_color(draw_win,shadow_color);
        x_line(draw_win,0,y,window_width,y);
        x_set_color(draw_win,white);
        x_line(draw_win,0,y+1,window_width-1,y+1);
       }
/*
      else
      { int x0 = xoff;
        int y0 = yoff-1;
        int y1 = y0 + menu_bar_height+1;
        int x1 = window_width-xoff;

        x_set_color(draw_win,white);
        x_line(draw_win,x0,y0,x1,y0);
        x_line(draw_win,x0,y0,x0,y1);
        x_set_color(draw_win,shadow_color);
        x_line(draw_win,x0+1,y1,x1+1,y1);
        x_line(draw_win,x1,y1,x1,y0);
       }
*/
    }
  }

  redraw_panel(0);
  clipping(0); 

  x_set_mode(draw_win,save_mode);
  x_set_line_style(draw_win,save_style);
  x_set_line_width(draw_win,save_lw);
}


void BASE_WINDOW::redraw_panel(panel_item x)
{
  // redraw x or (if x==0) all items

  int          save_lw = x_set_line_width(draw_win,1);
  line_style   save_ls = x_set_line_style(draw_win,solid);
  text_mode    save_tm = x_set_text_mode(draw_win,transparent);
  drawing_mode save_mo = x_set_mode(draw_win,src_mode);

  clipping(0);

  x_set_color(draw_win,black);

  int i;
  for(i=0; i<item_count; i++)
  { panel_item it=Item[i];

    if (x != 0 && it != x) continue;

    switch (it->kind) {

    case Newline_Item:
        { if (menu_button_style != 1) break;
          draw_separator(it);
          break;
         }

    case Button_Item:
        { draw_button(it,it == last_sel_button);
          break;
         }

    case Text_Item:
        { draw_text_item(it,draw_win);
          break;
         }

    case Choice_Item:
    case Choice_Mult_Item:
        { draw_choice_item(it,*(int*)it->ref);
          break;
         }

    case Bool_Item:
        { draw_bool_item(it);
          break;
         }

    case Color_Item:
        { draw_color_item(it);
          break;
         }

    case Bitmap_Choice_Item:
    case Bitmap_Choice_Mult_Item:
        { draw_bitmap_choice_item(it);
          break;
         }

    case Slider_Item:
        { draw_slider_item(it,0);
          break;
         }

    case Int_Item:
        { if (act_str_item == 0) act_str_item = it;
          delete[] it->data_str;
          it->data_str = make_string(*(int*)it->ref);
          draw_string_item(it);
          break;
         }
 
    case Float_Item:
        { if (act_str_item == 0) act_str_item = it;
          delete[] it->data_str;
          it->data_str = make_string(*(double*)it->ref);
          draw_string_item(it);
          break;
         }
  
    case String_Item:
    case String_Menu_Item:
        { it->data_str = access_str(it->ref);
          if (act_str_item == 0) act_str_item = it;
          it->offset = strlen(it->data_str);
          it->dat1 = 0;
          it->dat2 = string_w/tw - 1;
          draw_string_item(it);
          break;
         }
   }

 }

 clipping(2);
 x_set_line_width(draw_win,save_lw);
 x_set_line_style(draw_win,save_ls);
 x_set_mode(draw_win,save_mo);
 x_set_text_mode(draw_win,save_tm);

 x_flush_display();
}


    
void BASE_WINDOW::open_sub_panel(panel_item it, bool anim)
{ 
  BASE_WINDOW* wp = (BASE_WINDOW*)it->ref;

  wp->last_sel_button = 0;

  if (wp->is_open()) return;

  bool is_menu = (buts_per_line == 1 && item_count == but_count);

  // find "owner" p of (menu)button
  int j = 0;
  while (j<item_count && Item[j]->menu_but != it) j++;
  panel_item owner = (j < item_count) ? Item[j] : 0;


  int x,y;
  BASE_WINDOW* wpp = p_root;

  wp->p_win  = this;
  wp->p_root = p_root;

  //if (wp->menu_button_style == 0) 
  { if (wp->owner_item == 0 || wp->owner_item->menu_but != wp->owner_item)
      wp->menu_button_style = p_root->menu_button_style;
   }

  //if (this != p_root) 
  if (is_menu)
    { x = window_width + 1;
      y = it->ycoord + 1;
      if (wp->menu_button_style == 1) { x--; y++; }
      wp->hotx1 = -window_width;
      wp->hotx2 = 0;
      wp->hoty1 = 0;
      wp->hoty2 = yskip;
      x_set_border_width(wp->draw_win,0);
    }
  else
   { x = it->xcoord; 
     y = it->ycoord + it->height + 2;
   //y = it->ycoord + (yskip+it->height)/2 + 2;

     wp->hotx1 = 0;
     wp->hotx2 = it->width;
     wp->hoty1 = it->ycoord - y;
     wp->hoty2 = 0;

     if (wp->menu_button_style == 2) // scroll box menu
       { wp->yskip = 0;
         wp->button_w = string_w - 1;
         if (wp->but_count > wp->menu_size) wp->button_w -= wp->button_h+1;

         // find active item
         wp->act_str_item = 0;
         for(j = 0; j < wp->item_count; j++) 
           if (strcmp(wp->Item[j]->label_str,owner->data_str) == 0)
           { wp->act_str_item = wp->Item[j];
             break;
            }
         
         x = owner->xcoord + label_w;
         //y += 2;
        }
     else x_set_border_width(wp->draw_win,0);
    }
  
    if (owner && owner->kind==String_Menu_Item)
      activate_string_item(owner,panel_width);

    x_window_to_screen(draw_win,&x,&y);
    x_screen_to_window(p_root->draw_win,&x,&y);

    wp->place_panel_items();

    if (p_root->window_height < y + wp->panel_height + 2 ||
        p_root->window_width  < x + wp->panel_width + 2)
    { 
      wpp = wp;
      x_window_to_screen(p_root->draw_win,&x,&y);

      int dx = x_display_width()  - (x + wp->panel_width);
      int dy = x_display_height() - (y + wp->panel_height);
  
      if (dx < 0) 
      { x += dx;
        wp->hotx1 -= dx;
        wp->hotx2 -= dx;
       }
  
      if (dy < 0)
      { y += dy;
        wp->hoty1 -= dy;
        wp->hoty2 -= dy;
       }
    }

  active_window = wp;


  if (!anim)
     wp->display(x,y,wpp);
  else
   { int w = wp->panel_width;
     int h = wp->panel_height;
     wp->window_height = 1;
     wp->display(x,y,wpp);
     wp->window_width  = w;
     wp->window_height = h;
     wp->start_buffering(w,h);
     wp->draw_panel_items();
     char* pr;
     wp->stop_buffering(pr);
     wp->set_bg_pixrect(pr);
     for(int i=1; i<h; i++) wp->resize(x,y,w,i); 
     wp->draw_frame();
    }
  

  if (wp->menu_button_style == 2)
    { double sz = double(wp->menu_size)/wp->item_count;
      if (sz < 1)
      { int  h = wp->item_count * wp->yskip;
        double pos = double(wp->Item[0]->ycoord)/(wp->panel_height - h);
        wp->open_scrollbar(scroll_up_action,scroll_down_action,
                           scroll_drag_action, sz,pos);
       }
     }


  x_grab_pointer(wp->draw_win);
  x_set_focus(wp->draw_win);


/*
  { // skip pending events from wpp
    int e,w,v,x,y;
    unsigned long t;
    do { e = x_check_next_event(&w, &v, &x, &y, &t);
    } while (e != no_event && wp->draw_win != w);
    if (e != no_event) x_put_back_event();
  }
*/

  wp->last_sel_button = 0;
  wp->panel_menu_mode = 1;
}

    
    
void BASE_WINDOW::close_sub_panel(BASE_WINDOW* wp)
{ wp->last_sel_button = 0;
  wp->panel_menu_mode = 0;

  if (wp->is_open())
  { 
    // read all remaining events 
    int e,w,v,x,y;
    unsigned long t;
    while (x_check_next_event(&w, &v, &x, &y, &t) != no_event);

    wp->close();

    wp = wp->p_win;

    if (wp == 0) return;

    x_set_focus(wp->draw_win);

    if (!x_display_bits_saved() && win_parent != this)
       { while (x_get_next_event(&w, &v, &x, &y, &t) != exposure_event);
         x_put_back_event();
         event_handler(wp,1);
        }
    else
       while ((e=x_check_next_event(&w, &v, &x, &y, &t)) != no_event)
       { if (e == exposure_event)
         { x_put_back_event();
           event_handler(wp,1);
          }
        }

  }
}


void BASE_WINDOW::scroll(int steps)
{ 
  panel_item p1 = Item[0];
  panel_item p2 = Item[item_count-1];

  int dy = (steps > 0) ? yskip : -yskip;

  if (steps < 0) steps = -steps;

  if (p1->ycoord + steps*dy > button_h/2) return;
  if (p2->ycoord + steps*dy < panel_height-3*button_h/2) return;

  clipping(0);
  x_start_buffering(draw_win);
  clipping(0);

  if (steps == 0)
    { draw_panel_items();
      x_flush_buffer(draw_win,0,0,panel_width,panel_height);
     }
  else
    while (steps--)
    { for(int i=0; i<item_count; i++) Item[i]->ycoord += dy;
      draw_panel_items();
      x_flush_buffer(draw_win,0,0,panel_width,panel_height);
     }
  x_stop_buffering(draw_win);

  double f = double(p1->ycoord)/(panel_height - item_count*yskip);
  set_scrollbar_pos(f);
}



int BASE_WINDOW::panel_event_handler(int w, int k, int b, int x, int y, 
                                                         unsigned long t)
{
  //printf("panel event: w = %d  k = %d  b = %d  x = %d  y = %d\n",w,k,b,x,y);

  int          save_lw = x_set_line_width(draw_win,1);
  line_style   save_ls = x_set_line_style(draw_win,solid);
  text_mode    save_tm = x_set_text_mode(draw_win,transparent);
  drawing_mode save_mo = x_set_mode(draw_win,src_mode);
  clipping(1);

  if (w != draw_win) 
  { x = -1;
    y = -1;
   }

  int but = -1;

  panel_item it = 0;

  if (menu_button_style == 2)
  { int i = 0;
    int max_i = item_count-1;

    if (act_str_item) 
       i = act_str_item->index;
    else
      { if (b == KEY_DOWN) i = -1; 
        if (b == KEY_UP)   i = max_i + 1;  
       }

    if (k == key_press_event)
    {
      x_start_buffering(draw_win);

      switch (b) {

      case KEY_UP   : if (i == 0) break;
                      act_str_item = Item[i-1];
                      if (act_str_item->ycoord < 0) 
                         scroll(+1);
                      else
                         scroll(0);
                      break;

      case KEY_DOWN : if (i == max_i) break;
                      act_str_item = Item[i+1];
                      if (act_str_item->ycoord+button_h > panel_height) 
                         scroll(-1);
                      else
                         scroll(0);
                      break;

      case KEY_HOME : scroll(-item_count);
                      break;

      case KEY_END  : scroll(+item_count);
                      break;
      }

      x_stop_buffering(draw_win);

      if ( b == KEY_RETURN && act_str_item)
      { // simulate button release on last selected button
        x = act_str_item->xcoord + button_w/2;
        y = act_str_item->ycoord + button_h/2;
        last_sel_button = act_str_item;
        k = button_release_event;
       }
      else return -1;
     }

   }


  if (k == key_press_event && b != KEY_RETURN && act_str_item != 0) 
     { // input for active string item

        if (/*b == KEY_RETURN ||*/ b == KEY_DOWN || b == KEY_UP)
           { int i = act_str_item->index;
             int ki;
             do { if (b == KEY_UP)
                    { if (--i < 0) i = item_count-1; }
                  else
                    { if (++i >= item_count) i = 0; }
                  ki =  Item[i]->kind; 
                 }
             while (!Item[i]->enabled ||
                    (ki != String_Item && ki != String_Menu_Item && 
                     ki != Int_Item    && ki != Float_Item));
             activate_string_item(Item[i],panel_width);
            }
         else // character
            { it = act_str_item;
              x = it->xcoord + string_w + label_w;
              y = it->ycoord;
              x_put_back_event(); // event will be handled by panel_text_edit
             }
       }
    else
      { 
        if (focus_button && k == key_press_event && b == KEY_RETURN) 
        { panel_item it = focus_button;
          draw_button(it,1);
          { // wait for key release event
            int w,v,x,y;
            unsigned long t;
            while (x_get_next_event(&w,&v,&x,&y,&t) != key_release_event);
           }
          // simulate button release on focus button
          x = it->xcoord + button_w/2;
          y = it->ycoord + button_h/2;
          last_sel_button = it;
          k = button_release_event;
        }

        int i;
        for(i=0; i<item_count; i++) // search for selected item
        { it = Item[i];
          if (!it->enabled) continue;
          int x1 = it->xcoord;
          int x2 = x1 + it->width; 
          int y1 = it->ycoord; 
          int y2 = y1 + it->height;
          if (it->kind != Button_Item) x1 += label_w;
          if (it->kind == Color_Item)  x1 -= (button_h+2);
          if (panel_menu_mode) 
          { y1 -= 2;
            y2 += 2;
           }
          if (x < x1 || x > x2 || y < y1 || y > y2) continue;
          break;
         }

        if (i == item_count)  it = 0;  // no item selected

        if (it && k == motion_event && it->kind == Button_Item
               && ((it->dat2 & 1) || it->help_str || 
                   (it->menu_but==it && strlen(it->label_str) == 0)))
          { active_button_win = this;
            active_button = it;
           }
        else
           active_button = 0;

        if (k != button_press_event && !panel_menu_mode) it = 0;

        if (k == button_press_event && last_sel_button
                                    && (it==0 || it->kind!=Button_Item)) 
        { // buttons are handled later
          // close open sub window (if not selected again)
          panel_item lsb = last_sel_button;
          last_sel_button = 0;
          draw_button(lsb,0);
          if (lsb->ref)
          { close_sub_panel((BASE_WINDOW*)lsb->ref);
            active_window = this;
           }
         }


       // repeaters, e.g. scrolling (it->offset > 0)

       if (it && it->kind == Button_Item && it->offset > 0 && (it->action || it->action_ptr)
                                         && k == button_press_event)
       { int v;
         unsigned long t;
         draw_button(it,1);
         last_sel_button = it;
         call_window = this;
         call_item = it;
         clipping(2);
	 
	 if (it->action_ptr) {
	   it->action_ptr->set_params((window*)this,NULL,it->dat1);
	   it->action_ptr->operator()();
	 }
         else it->action(it->dat1);
	 
         clipping(1);
         int t1 = 500;
         int t2 = it->offset;
         // return if button released after t1 msec
         if (x_get_next_event(&w,&v,&x,&y,&t,t1) != button_release_event)
         { // otherwise execute action every t2 msec until button released
           draw_button(it,1);
           while (x_get_next_event(&w,&v,&x,&y,&t,t2)!=button_release_event)
           { if (it->ref == 0) draw_button(it,0);
             clipping(2);
	     
	     if (it->action_ptr) {
	       it->action_ptr->set_params((window*)this,NULL,it->dat1); 
	       it->action_ptr->operator()();
	     }
             else it->action(it->dat1);
	     
             clipping(1);
             if (it->ref == 0) draw_button(it,1);
            }
          }
         last_sel_button = 0;
         clipping(1);
         draw_button(it,0);
         clipping(2);
         return -1;
       }

       // check special buttons (scrollbar slider: it->ref == this)

       if (it && it->kind == Button_Item && it->ref == this
                                         && k == button_press_event)
       { int e,v;
         unsigned long t;
         draw_button(it,1);
         last_sel_button = it;

         if (it->action || it->action_ptr) 
         { call_window = this;
           call_item = it;
           //int but_dy = it->ycoord + (yskip - it->height)/2 - y;
           int but_dy = it->ycoord - y;
           int yoff1 = button_h;
           //int yoff2 = it->height + button_h - 1;
           int yoff2 = it->height + button_h + 1;
           int y0 = 0;

           BASE_WINDOW* owner = (BASE_WINDOW*)get_inf();

           if (owner->p_win == 0)
               x_grab_pointer(draw_win);
	       
           if (it->action_ptr) {
	     it->action_ptr->set_params((window*)this,NULL,-1);
	     it->action_ptr->operator()();
	   }
           else it->action(-1); // start

           while ((e=x_get_next_event(&w,&v,&x,&y,&t))!=button_release_event)
           { if (e != motion_event || y == y0) continue;
             y0 = y;
             y += but_dy;
             if (y < yoff1) y = yoff1;
             if (y > panel_height - yoff2) y = panel_height - yoff2;
             //it->ycoord  = y + (it->height - yskip)/2;
             it->ycoord  = y;
             int D = (1000*(y - yoff1))/(panel_height-yoff1-yoff2);
	     
	     if (it->action_ptr) {
	       it->action_ptr->set_params((window*)this,NULL,D);
	       it->action_ptr->operator()();
	     }
             else it->action(D);
	     
             //clipping(0);
             x_start_buffering(draw_win);
             clipping(0);
             draw_panel_items();
             x_flush_buffer(draw_win,0,0,panel_width,panel_height);
             x_stop_buffering(draw_win);
            }

           if (owner->p_win == 0) x_ungrab_pointer();

           last_sel_button = 0;
           draw_button(it,0);

           if (it->action_ptr) {
	     it->action_ptr->set_params((window*)this,NULL,-2);
	     it->action_ptr->operator()();
	   }
           else it->action(-2); // finish

           return -1;
          }
       }


       if (k == button_release_event || 
          (k == button_press_event && panel_menu_mode))
       { 
         if (last_sel_button == 0)
           { active_window = this;
             panel_menu_mode = 0; }
         else
           { it = last_sel_button;
             if (it->ref == 0)
               { draw_button(it,0);
                 but = it->dat1;
                }
             else
               { BASE_WINDOW* wp =  ((BASE_WINDOW*)it->ref);
                 active_window = wp;
                 wp->panel_menu_mode = 0;
                 it = 0;
                }
            }
         }

      }


    if (but == -1 && it) // item selected
    {
      switch (it->kind) {
  
      case Text_Item: break;
  
      case Slider_Item:
      { int xoff1 = it->xcoord + label_w;
        int w = draw_win;
        int xlast = 0;
        if (x < xoff1 || x > xoff1+slider_w) break;
        //while (w == draw_win)
        for(;;)
        { if (xlast != x)
          { xlast = x;
            draw_slider_item(it,xlast);
           }
          if (x_get_next_event(&w,&b,&x,&b,&t) == button_release_event) break;
         }
        break;
       }
  
      case Int_Item:
      { activate_string_item(it,x);
        panel_text_edit(it);
        *(int*)it->ref = atoi(it->data_str);
        break;
       }
  
      case Float_Item:
      { activate_string_item(it,x);
        panel_text_edit(it);
        *(double*)it->ref = atof(it->data_str);
        break;
       }
  
      case String_Item:
      case String_Menu_Item:
      { activate_string_item(it,x);
        panel_text_edit(it);
        assign_str(it->ref,it->data_str);
        it->data_str = access_str(it->ref);
        break;
       }
  
     case Choice_Item:
     { int d = (x-it->xcoord-label_w)/choice_w;
       if (d < it->argc)
       { int x = it->offset + d * it->dat2;
         draw_choice_item(it,x);
         if (it->action || it->action_ptr)
         { clipping(2);
           call_window = this;
           call_item = it;
	   
	   if (it->action_ptr) {
	     it->action_ptr->set_params((window*)this,NULL,x);
	     it->action_ptr->operator()();
	   }
           else it->action(x);
           clipping(1);
          }
         *(int*)it->ref = x;
        }
       break;
      }

     case Choice_Mult_Item:
     { int d = (x-it->xcoord-label_w)/choice_w;
       if (d < it->argc)
       { int x = *(int*)it->ref;
         x ^= (1 << d);
         draw_choice_item(it,x);
         if (it->action || it->action_ptr)
         { clipping(2);
           call_window = this;
           call_item = it;
	   
	   if (it->action_ptr) {
	     it->action_ptr->set_params((window*)this,NULL,x);
	     it->action_ptr->operator()();
	   }
           else it->action(x);
	   
           clipping(1);
          }
         *(int*)it->ref = x;
        }
       break;
      }
  
  
     case Bool_Item:
     { int b = (*(bool*)it->ref) ? 0 : 1;
       *(bool*)it->ref = (b != 0);
       draw_bool_item(it);
       if (it->action || it->action_ptr) 
       { clipping(2);
         call_window = this;
         call_item = it;
	 
	 if (it->action_ptr) {
	   it->action_ptr->set_params((window*)this,NULL,b);
	   it->action_ptr->operator()();
	 }
         else it->action(b);
         clipping(1);
       }
       break;
      }
  
     case Color_Item:
     { int dx = color_h + 2;
       int x0 = it->xcoord + label_w;
       int d = (x - x0)/dx - 1;
       if (d < 17) 
       { if (b == 3 && d >= 0)
         { x_choose_color(0,d);
           draw_color_item(it);
           break;
          }
         if (it->action || it->action_ptr) 
         { clipping(2);
           call_window = this;
           call_item = it;
	   
	   if (it->action_ptr) {
	     it->action_ptr->set_params((window*)this,NULL,d);
	     it->action_ptr->operator()();
	   }
           else it->action(d);
           clipping(1);
         }
         change_color_item(it,d);
       }
       break;
      }
  
     case Bitmap_Choice_Mult_Item:
     case Bitmap_Choice_Item:
     { int d = (x - it->xcoord - label_w)/(it->dat1 + 3);
       if (d < it->argc) 
       { int c = d;
         if (it->kind == Bitmap_Choice_Mult_Item) 
            c = (*(int*)it->ref) ^ (1 << d);
         if (c != *(int*)it->ref) change_bitmap_choice_item(it,d);
         if (it->action || it->action_ptr) 
         { clipping(2);
           call_window = this;
           call_item = it;
	   
	   if (it->action_ptr) {
	     it->action_ptr->set_params((window*)this,NULL,c);
	     it->action_ptr->operator()();
	   }
           else it->action(c);
           clipping(1);
         }
         *(int*)it->ref = c;
        }
       break;
     }

     case Button_Item:
     { 
       panel_item lsb = last_sel_button;

       if (menu_button_style == 2 && act_str_item)
       { panel_item tmp = act_str_item;
         act_str_item =0;
         draw_button(tmp,0);
        }

       if (it != lsb) 
       { draw_button(it,1);
         last_sel_button = it;
        }
       else  // second click on active button
         if (k == button_press_event) last_sel_button = 0;
    
       if (lsb && (it != lsb || k==button_press_event))
       { draw_button(lsb,0);
         if (lsb->ref)
         { close_sub_panel((BASE_WINDOW*)lsb->ref);
           active_window = this;
          }
        }

       if (it->menu_but == it)
       { if (it->action || it->action_ptr)
         { call_window = this;
           call_item = it;
	   
	   if (it->action_ptr) {
	     it->action_ptr->set_params((window*)this,NULL,0);
	     it->action_ptr->operator()();
	   }
           else it->action(0);
          }
         if (it->ref == 0)
         { draw_button(it,0);
           last_sel_button = 0;
          }
        }

       if (it != lsb && it->ref && it->ref != this) open_sub_panel(it,b==2);

       it = 0;
       break;
     }

    }
  }

  int result = but;

  if (it  && it->kind == Button_Item  && it->ref == 0 ) 
  { BASE_WINDOW* wp = this;
    while(wp->p_win) 
    { BASE_WINDOW* wpp = wp->p_win;
      close_sub_panel(wp);
      wp->p_win = 0;
      wp = wpp;
     }
    active_window = wp;

    panel_item lsb = wp->last_sel_button;

    if (lsb)
    { drawing_mode save_mo = x_set_mode(wp->draw_win,src_mode);
      wp->clipping(1);
      wp->draw_button(lsb,0);
      x_set_mode(wp->draw_win,save_mo);
      wp->clipping(2);
      wp->last_sel_button = 0;
      if (lsb->dat2 < 2)  // second bit = 0: return number of lsb 
           result = lsb->dat1;
     }

     if (lsb == 0 || lsb->menu_but == lsb) result = -1;

   }

  x_set_line_width(draw_win,save_lw);
  x_set_line_style(draw_win,save_ls);
  x_set_mode(draw_win,save_mo);
  x_set_text_mode(draw_win,save_tm);
  clipping(2);

  if (it && it->kind == Button_Item && (it->action || it->action_ptr) && it->menu_but == 0)
  { call_window = this;
    call_item = it;
    
    if (it->action_ptr) {
      it->action_ptr->set_params((window*)this,NULL,but);
      it->action_ptr->operator()();
    }
    else it->action(but);
   }

 return result;
}




panel_item BASE_WINDOW::get_item(const char* s)
{ for(int i = 0; i<item_count; i++)
  { panel_item p = Item[i];
    if (strcmp(p->label_str,s) == 0) return p;
   }
  return NULL;
}

int BASE_WINDOW::get_button(const char* s)
{ panel_item p = get_item(s);
  if (p && p->kind == Button_Item)
     return p->dat1;
  else
     return -1;
}


panel_item BASE_WINDOW::get_button_item(int but)  
{ for(int i = 0; i<item_count; i++)
  { panel_item p = Item[i];
    if (p->kind == Button_Item && p->dat1 == but)  return p;
   }
  return 0;
 }




BASE_WINDOW* BASE_WINDOW::get_window(panel_item it)
{ return (it->kind == Button_Item) ? (BASE_WINDOW*)it->ref : 0; }

BASE_WINDOW* BASE_WINDOW::get_window(int but)
{ panel_item p = get_button_item(but);
  return p ? (BASE_WINDOW*)p->ref : 0;
}



BASE_WINDOW* BASE_WINDOW::set_window(int but, BASE_WINDOW* wp)
{ BASE_WINDOW* ret = 0;
  panel_item p = get_button_item(but);
  if (p)
  { ret =  (BASE_WINDOW*)p->ref;
    p->ref = wp;
    if (wp) wp->owner_item = p;
    if (ret) ret->owner_item = 0;
   }
  return ret;
}

panel_action_func BASE_WINDOW::set_action(int but, 
                                                 panel_action_func func)
{ panel_action_func ret = 0;
  panel_item p = get_button_item(but);
  if (p)
  { ret =  p->action;
    p->action = func;
   }
  return ret;
}

const window_handler* BASE_WINDOW::set_action(int but, const window_handler& obj)
{ const window_handler* ret = 0;
  panel_item p = get_button_item(but);
  if (p)
  { ret =  p->action_ptr;
    p->action_ptr = & (window_handler&) obj;
   }
  return ret;
}


int BASE_WINDOW::set_value(int but, int val)
{ int ret = -1;
  panel_item p = get_button_item(but);
  if (p)
  { ret =  p->dat1;
    p->dat1 = val;
   }
  return ret;
}



void BASE_WINDOW::enable_item(panel_item it) 
{ if (it && !it->enabled)
  { it->enabled = 1; 
    panel_enabled = 1;
    if (is_open()) redraw_panel(it);
   }
}


void BASE_WINDOW::disable_item(panel_item it) 
{ if (it && it->enabled)
  { it->enabled = 0; 
    if (is_open()) redraw_panel(it);
   }
}


void BASE_WINDOW::disable_panel(bool disable_items)  
{ panel_enabled = 0; 
  if (!disable_items) return;
  for (panel_item it = first_item(); it; it = next_item(it))
     it->enabled = 0;
  if (is_open()) redraw_panel(0);
}


void BASE_WINDOW::enable_panel()  
{ for (panel_item it = first_item(); it; it = next_item(it))
     it->enabled = 1;
  if (is_open()) redraw_panel(0);
  panel_enabled = 1; 
}


const char* BASE_WINDOW::get_item_label(panel_item p) { return p->label_str; }


void BASE_WINDOW::set_item_label(panel_item p, const char* s)
{ delete[] p->label_str;
  p->label_str = string_dup(s);
  if (is_open()) redraw_panel(p);
 }



const char* BASE_WINDOW::get_button_label(int but)
{ panel_item p = get_button_item(but);
  return  (p) ? p->label_str : 0;
}

void BASE_WINDOW::set_button_label(int but, const char* s)
{ panel_item p = get_button_item(but);
  if (p) 
  { delete[] p->label_str;
    p->label_str = string_dup(s);
    if (is_open()) redraw_panel(p);
   }
}


void BASE_WINDOW::set_button_help_str(int but, const char* hlp)
{ panel_item p = get_button_item(but);
  if (p)
  { delete[] p->help_str;
    p->help_str = string_dup(hlp);
   }
}


void BASE_WINDOW::set_button_repeat(int but, int msec)
{ panel_item p = get_button_item(but);
  if (p) p->offset = msec;
}


void BASE_WINDOW::set_button_pixrects(int but,char*pmap0,char*pmap1)
{ panel_item p=get_button_item(but);
  if(p) 
  { delete[] p->argv; 
    p->argv = new char*[2];
    p->argv[0] = pmap0;
    p->argv[1] = pmap1;
    p->argc = -1; // pmap
    if (is_open()) redraw_panel(p);
  }
}


// Iteration

panel_item BASE_WINDOW::first_item() const { return Item[0]; } 

panel_item BASE_WINDOW::next_item(panel_item it)  const
{ if (it == 0) return 0;
  int i = it->index + 1;
  return (i < item_count) ? Item[i] : 0;
 }




void BASE_WINDOW::open_scrollbar(void (*scroll_up)(int),
                                 void (*scroll_down)(int),
                                 void (*scroll_drag)(int), 
                                 double sz, double pos)
{ if (sz > 1) sz = 1;
  if (sz < 0) sz = 0;

  int w = button_h - 2;

  if (scroll_bar) // adjust slider size of existing scroll bar
  { panel_item slid = scroll_bar->Item[2];
    slid->height = int(sz*(scroll_bar->panel_height - 2*w));
    set_scrollbar_pos(pos);
    return;
   }

  BASE_WINDOW* wp = this;
  BASE_WINDOW* sp = new BASE_WINDOW(-1,-1);

  sp->menu_button_style = 0;

  sp->set_inf(wp);

  sp->panel_width  = w+2; 

  if (window_height == panel_height)
     sp->panel_height = window_height-2;
  else 
     sp->panel_height = window_height - panel_height-1;

  sp->window_width  = sp->panel_width;
  sp->window_height = sp->panel_height;
  sp->yskip         = w;

  sp->button("UP",-1,scroll_up);
  sp->button("DOWN",-1,scroll_down);
  sp->button("SLIDER",-1,scroll_drag);

  panel_item p;

  p = sp->Item[0];
  p->ref    = sp;
  p->width  = w;
  p->height = w;
  p->xcoord = 0;
  p->ycoord = 0;
  p->menu_but = p;
  p->offset = 50;  // repeat every 50 msec
    
  p = sp->Item[1];
  p->ref    = sp;
  p->width  = w;
  p->height = w;
  p->xcoord = 0;
  p->ycoord = sp->panel_height - w - 1;
  p->menu_but = p;
  p->offset = 50;  // repeat every 50 msec

  p = sp->Item[2];
  p->ref    = sp;
  p->width  = w;
  p->height = int(sz*(sp->panel_height - 2*w));
  p->xcoord = 0;
  p->ycoord = int(w + pos *(sp->panel_height - 2*w - p->height - 5)) + 2;
  p->menu_but = p;

  sp->window_xpos = window_width-w-4;

  if (window_height == panel_height)
    sp->window_ypos = 0;
  else
    sp->window_ypos = panel_height-1;


  x_open_window(sp->draw_win,sp->window_xpos,sp->window_ypos,
                             sp->window_width,sp->window_height,draw_win);
  sp->configure();
  sp->draw_frame();
  sp->scroll(0);

  if (scroll_bar) delete scroll_bar;
  scroll_bar = sp;
  
  draw_frame();  
}


void BASE_WINDOW::close_scrollbar()
{ if (scroll_bar) delete scroll_bar;
  scroll_bar = 0;
 }



void BASE_WINDOW::set_scrollbar_pos(double p)
{
  if (scroll_bar == 0) return;

  if (p > 1) p = 1;
  if (p < 0) p = 0;

  BASE_WINDOW* sp = scroll_bar;

  panel_item it1 = sp->Item[0];
  panel_item it2 = sp->Item[1];
  panel_item it3 = sp->Item[2];

  int scroll_h = sp->panel_height-it1->height-it2->height-it3->height - 5;

  it3->ycoord = int(it1->height + p*scroll_h) + 2;

  if (sp->is_open())
  { sp->clipping(0);
    x_start_buffering(sp->draw_win);
    sp->clipping(0);
    sp->draw_panel_items();
    x_flush_buffer(sp->draw_win,0,0,sp->panel_width,sp->panel_height);
    x_stop_buffering(sp->draw_win);
  }
}
  

}


