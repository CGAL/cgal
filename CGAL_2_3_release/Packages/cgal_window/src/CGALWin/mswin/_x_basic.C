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
// release       : 
// release_date  : 2001, May 23
//
// file          : src/CGALWin/mswin/_x_basic.C
// package       : cgal_window (0.9.7)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 0.9.7
// revision_date : 23 May 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================


//-----------------------------------------------------------------------------
// Basic Graphics for Win32 (nt, win95, win32s)
//
// (c) Algorithmic Solutions 1996-1999
//-----------------------------------------------------------------------------

#include <CGAL/LEDA/basic.h>
#include <CGAL/LEDA/system.h>

#if !defined(__BORLANDC__)
#include <cctype>
#endif

#include <CGAL/LEDA/impl/x_basic.h>

#if defined(_MSC_VER)
#pragma warning(disable:4305)
#pragma warning(disable:4309)
#endif



#include <CGAL/LEDA/pixmaps/win_icon1.xpm>
#include <CGAL/LEDA/pixmaps/win_small_icon1.xpm>

#include <cstdio>
#include <cstdarg>
#include <string>
#include <cctype>
#include <ctime>
#include <cassert>

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <commdlg.h>


#include <shellapi.h>

#define LEDA_TASKBAR_ICON

#define MAX_WIN    256
#define MAX_COLORS 1024


#define XWM_TASKBAR_ICON (WM_USER+1)
#define XWM_FONT_DIALOG  (WM_USER+2)
#define XWM_COLOR_DIALOG (WM_USER+3)


namespace CGAL {


typedef int Window;

typedef void (*redraw_func)(void*,int,int,int,int);
typedef void (*drop_func)(void*,const char*,int,int);

struct win_struct 
{
  HWND         hwnd;
  HDC          hdc;
  HDC          hdcMem;
  HDC          hdc1;
  HDC          hdc2;
  HENHMETAFILE hmf;
  HBITMAP      hbm;
  LOGPEN       pen_data;
  HBITMAP      stip_bm;
  int          stip_bgcol;
  HFONT        font;
  HFONT        tmp_font;
  UINT         timer_id;


  int          width;
  int          height;
  int          border_w;

  int          buf_width;
  int          buf_height;

  int          COLOR;
  int          B_COLOR;
  char*        B_PIXMAP;
  int          LWIDTH;
  int          JSTYLE;
  line_style   LSTYLE;
  text_mode    TMODE;
  drawing_mode MODE;

  int          save_co;
  int          save_lw;
  line_style   save_ls;
  text_mode    save_tm;
  drawing_mode save_mo;

  char         font_name[8];

  char*        header;
  void*        inf;

  int          cursor;
  
  int wm_paint_count;

  redraw_func repaint;
  drop_func   drop_handler;


win_struct(void* ptr, int w, int h, int bg_col, const char* label, 
                                                redraw_func redraw)
{
  ZeroMemory(this,sizeof(win_struct));

  header = new char[strlen(label) + 1];
  strcpy(header,label);

  width   = w;
  height  = h;
  B_COLOR = bg_col;
  inf     = ptr;
  repaint = redraw;
  drop_handler = NULL;


  hwnd   = NULL;
  hdc    = NULL;
  hdcMem = NULL;
  hdc1   = NULL;
  hdc2   = NULL;
  hmf    = NULL;

  pen_data.lopnStyle   = PS_SOLID;
  pen_data.lopnWidth.x = 1;
  pen_data.lopnWidth.y = 1;
  pen_data.lopnColor   = 0;

  stip_bm    = NULL;
  stip_bgcol = white;

  B_PIXMAP = NULL;
  COLOR    = black;
  LSTYLE   = solid;
  JSTYLE   = 1;
  LWIDTH   = 1;
  MODE     = src_mode;
  TMODE    = transparent;

  font     = NULL;
  tmp_font = NULL;
  font_name[0] = '\0';

  timer_id = 0;
  border_w = 1;
  cursor   = -1;

  save_co = COLOR;
  save_ls = LSTYLE;
  save_lw = LWIDTH;
  save_mo = MODE;
  save_tm = TMODE;
}

~win_struct() { delete[] header; }

};


struct x_pixrect {
 int w;
 int h;
 HBITMAP map;
};


struct event {
  int win;
  int kind;
  int val;
  int x;
  int y;
  unsigned long t;
};



static event cur_event;
static event last_event;
static int putback;
static int trace_events;

static win_struct* wlist[MAX_WIN];
static int wcount = 0;


static COLORREF rgb_custom[16];

static COLORREF rgb_table[MAX_COLORS];
static int      rgb_count = 0;

static HFONT root_font;
static HFONT text_font;
static HFONT italic_font;
static HFONT bold_font;
static HFONT fixed_font;
static HFONT button_font;

static LOGFONT  root_lf;
static LOGFONT  text_lf;
static LOGFONT  italic_lf;
static LOGFONT  bold_lf;
static LOGFONT  fixed_lf;
static LOGFONT  button_lf;

static HCURSOR  dot_box_cursor;
static HCURSOR  hand_cursor;
static HICON    leda_icon;
static HICON    leda_small_icon;

static HWND     taskbar_icon_win;
static HMENU    taskbar_icon_menu;

static char szAppName1[] = "LEDA Window";
static char szAppName2[] = "LEDA Subwindow";

static Window grab_win;

static int sys_key_down = 0;

static char dot_box_bits[] = 
  { 0xFF,0xF8,0x00,0x00,
    0x80,0x08,0x00,0x00,
    0x80,0x08,0x00,0x00,
    0x80,0x08,0x00,0x00,
    0x80,0x08,0x00,0x00,
    0x80,0x08,0x00,0x00,
    0x82,0x08,0x00,0x00,
    0x80,0x08,0x00,0x00,
    0x80,0x08,0x00,0x00,
    0x80,0x08,0x00,0x00,
    0x80,0x08,0x00,0x00,
    0x80,0x08,0x00,0x00,
    0xFF,0xF8,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00};


static char hand_bits[] = 
  { 0x7F,0xF0,0x00,0x00,
    0x80,0x08,0x00,0x00,
    0x7F,0x84,0x00,0x00,
    0x08,0x02,0x00,0x00,
    0x07,0x82,0x00,0x00,
    0x08,0x02,0x00,0x00,
    0x07,0x85,0x00,0x00,
    0x08,0x08,0x80,0x00,
    0x07,0x90,0x40,0x00,
    0x00,0xE0,0x80,0x00,
    0x00,0x49,0x00,0x00,
    0x00,0x22,0x00,0x00,
    0x00,0x14,0x00,0x00,
    0x00,0x08,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00};



void x_add_taskbar_icon(HWND hwnd, HICON hicon, char* tip)
{
#if defined(LEDA_TASKBAR_ICON)
  NOTIFYICONDATA icon_data;
  icon_data.cbSize = sizeof(NOTIFYICONDATA);
  icon_data.hWnd = hwnd;
  icon_data.uID = 9999;
  icon_data.uFlags = NIF_MESSAGE | NIF_ICON | NIF_TIP;
  icon_data.uCallbackMessage = XWM_TASKBAR_ICON;
  icon_data.hIcon = hicon;
  strcpy(icon_data.szTip,tip);
  Shell_NotifyIcon(NIM_ADD,&icon_data);
  HMENU hm = CreatePopupMenu();
  char win_label[128];
  GetWindowText(hwnd,win_label,128);
  AppendMenu(hm,MF_STRING,0,"LEDA Setup");
  AppendMenu(hm,MF_SEPARATOR,0,NULL);
  AppendMenu(hm,MF_STRING,XWM_FONT_DIALOG,"Fonts");
  AppendMenu(hm,MF_STRING,XWM_COLOR_DIALOG,"Colors");
  taskbar_icon_menu = hm;
#endif
}
  
 
void x_change_taskbar_icon(HWND hwnd,HICON hicon,char* tip)
{
#if defined(LEDA_TASKBAR_ICON)
  NOTIFYICONDATA icon_data;
  icon_data.cbSize = sizeof(NOTIFYICONDATA);
  icon_data.hWnd = hwnd;
  icon_data.uID = 9999;
  icon_data.uFlags = NIF_MESSAGE | NIF_ICON | NIF_TIP;
  icon_data.uCallbackMessage = XWM_TASKBAR_ICON;
  icon_data.hIcon = hicon;
  strcpy(icon_data.szTip,tip);
  Shell_NotifyIcon(NIM_MODIFY,&icon_data);
#endif
}
  
   
 
void x_del_taskbar_icon(HWND hwnd)
{
#if defined(LEDA_TASKBAR_ICON)
  NOTIFYICONDATA icon_data;
  icon_data.cbSize = sizeof(NOTIFYICONDATA);
  icon_data.hWnd = hwnd;
  icon_data.uID = 9999;
  Shell_NotifyIcon(NIM_DELETE,&icon_data);
#endif
}




char* XpmBitmapBits(const char** xpm) 
{ 
  int width;
  int height;
  int colors;
  int chars;
  CGAL_CLIB_STD::sscanf(*xpm,"%d %d %d %d",&width,&height,&colors, &chars);

  if (chars > 2)
  { std::cerr << "xpm: sorry, chars_per_pixel > 2 not implemented.\n";
    exit(1);
   }

  char black_c1 = 0;
  char black_c2 = 0;

  int i;
  for(i=0; i<colors; i++)
  { xpm++;
    char c1=0;
    char c2=0;
    char rgb_str[16];

    if (chars == 1) CGAL_CLIB_STD::sscanf(*xpm,"%c c %s",&c1,rgb_str);
    if (chars == 2) CGAL_CLIB_STD::sscanf(*xpm,"%c%c c %s",&c1,&c2,rgb_str);

    if (strcmp(rgb_str,"None") == 0)
    { black_c1 = c1;
      black_c2 = c2;
     }
   }

  int   bytes = height*width/8;
  char* bits = new char[bytes];

  for(i=0; i<bytes; i++) bits[i] = 0;

  int pos = 0;
  for(int y=0; y<height; y++)
  { const char* line = *++xpm;
    for(int x=0; x<width; x++)
    { char c1 = line[x*chars];
      char c2 = 0;
      if (chars == 2) c2 = line[x*chars+1];
      if (c1 == black_c1 && c2 == black_c2)  bits[pos/8] |= (1 << (7-(pos%8)));
      pos++;
     }
   }

  return bits;
}



HICON CreateXpmIcon(int w, int h, const char** icon_xpm)
{

  char* icon_pr  = x_create_pixrect(0,icon_xpm); 

  HBITMAP icon_hbm = ((x_pixrect*)icon_pr)->map;

  BITMAP icon_bm;
  GetObject(icon_hbm,sizeof(BITMAP),&icon_bm);
  long icon_bytes = icon_bm.bmWidthBytes * icon_bm.bmHeight * icon_bm.bmPlanes; 
  char* icon_bits = new char[icon_bytes];

  GetBitmapBits(icon_hbm,icon_bytes,icon_bits);


  char* mask_bits = XpmBitmapBits(icon_xpm);


  HICON icon = CreateIcon(NULL,w,h,(char)icon_bm.bmPlanes,
                                   (char)icon_bm.bmBitsPixel,
                                   (const unsigned char*)mask_bits,
                                   (const unsigned char*)icon_bits);


  x_delete_pixrect(icon_pr);

  delete[] mask_bits;
  delete[] icon_bits;

  return icon;
}

  
   



inline void SWAP(int& x1, int& x2)
{ int t = x1; x1 = x2; x2 = t; }




static long WINAPI WndProc (HWND hwnd, UINT message, UINT wParam, LONG lParam);

void x_do_not_open_display(int) {}

void x_open_display(void)  
{ 
  if (text_font) return;

  // HDC hdc = CreateIC("DISPLAY",NULL,NULL,NULL); 
  // int bits_per_pixel = GetDeviceCaps(hdc,BITSPIXEL);
  // printf("bits per pixel: %d\n",bits_per_pixel);
  // DeleteDC(hdc);


  // create fonts

  ZeroMemory(&root_lf,  sizeof(LOGFONT));
  root_lf.lfHeight = -4;
  root_lf.lfPitchAndFamily = FF_SWISS | VARIABLE_PITCH;
  root_lf.lfWeight = FW_NORMAL;
  root_lf.lfItalic = 0;
  root_font = CreateFontIndirect(&root_lf);

  int disp_w = x_display_width();

  int sz = -13;
  if (disp_w <= 800)  sz = -12;


  ZeroMemory(&text_lf,  sizeof(LOGFONT));
  text_lf.lfHeight = sz;
  text_lf.lfWeight = FW_NORMAL;
  text_lf.lfItalic = 0;
  strcpy(text_lf.lfFaceName,"MS Sans Serif");
  text_font = CreateFontIndirect(&text_lf);

  ZeroMemory(&bold_lf,  sizeof(LOGFONT));
  bold_lf.lfHeight = sz; 
  bold_lf.lfWeight = FW_BOLD;
  bold_lf.lfItalic = 0;
  strcpy(bold_lf.lfFaceName,"MS Sans Serif");
  bold_font = CreateFontIndirect(&bold_lf);

  ZeroMemory(&italic_lf,sizeof(LOGFONT));
  italic_lf.lfHeight = sz;
  italic_lf.lfWeight = FW_NORMAL;
  italic_lf.lfItalic = 1;
  strcpy(italic_lf.lfFaceName,"MS Sans Serif");
  italic_font = CreateFontIndirect(&italic_lf);


  ZeroMemory(&fixed_lf, sizeof(LOGFONT));
  fixed_lf.lfHeight = sz+1;
  fixed_lf.lfWeight = FW_NORMAL;
  fixed_lf.lfItalic = 0;
  strcpy(fixed_lf.lfFaceName,"Lucida Console");
  fixed_font = CreateFontIndirect(&fixed_lf);


  ZeroMemory(&button_lf,sizeof(LOGFONT));
  button_lf.lfHeight = sz+1;
  button_lf.lfWeight = FW_NORMAL;
  button_lf.lfItalic = 0;
  strcpy(button_lf.lfFaceName,"Arial");
  button_font = CreateFontIndirect(&button_lf);

  //button_font= (HFONT)GetStockObject(DEFAULT_GUI_FONT);
  //button_font= (HFONT)GetStockObject(SYSTEM_FONT);


  // initialize rgb table

  rgb_table[white]  = RGB(255,255,255);
  rgb_table[black]  = RGB(  0,  0,  0);
  rgb_table[red]    = RGB(255,  0,  0);
  rgb_table[green]  = RGB(  0,255,  0);
//rgb_table[blue]   = RGB(  0,  0,128);
  rgb_table[blue]   = RGB(  0,  0,160);
  rgb_table[yellow] = RGB(255,255,  0);
  rgb_table[violet] = RGB(128,  0,255);
  rgb_table[orange] = RGB(255,128,  0);
  rgb_table[cyan]   = RGB(  0,255,255);
  rgb_table[brown]  = RGB(128,  0,  0);
  rgb_table[pink]   = RGB(255,  0,255);
  rgb_table[green2] = RGB(  0,128,128);
//rgb_table[blue2]  = RGB(  0,  0,255);
  rgb_table[blue2]  = RGB(100,100,255);

  //rgb_table[grey1]  = RGB(205,205,205);
  rgb_table[grey1]  = RGB(220,220,220);
  rgb_table[grey2]  = RGB(175,175,175);
  rgb_table[grey3]  = RGB(128,128,128);
  rgb_table[ivory]  = RGB(255,255,230);

  rgb_count = 17;


  int i;
  for(i=0; i<16; i++) rgb_custom[i] = RGB(255,255,255);


  // wlist[0]
  win_struct* wp = new win_struct(0,0,0,0,"root window",0);
  wp->font = fixed_font;
  wp->hwnd = 0;
  wlist[0] = wp;
  wcount = 0;

  for(int j = 1; j < MAX_WIN; j++) wlist[j] = 0;


  // cursor
  // int cw = GetSystemMetrics(SM_CXCURSOR);
  // int ch = GetSystemMetrics(SM_CYCURSOR);

  for(i=0;i<128;i++) dot_box_bits[i] = ~dot_box_bits[i];
  for(i=0;i<128;i++) hand_bits[i] = ~hand_bits[i];

  char* xor_bits = new char[512];
  for(i=0;i<512;i++) xor_bits[i] = 0x00;

  dot_box_cursor = CreateCursor(NULL,6,6,32,32,dot_box_bits,xor_bits);
  hand_cursor = CreateCursor(NULL,0,1,32,32,hand_bits,xor_bits);

  delete[] xor_bits;

  leda_icon = LoadIcon (NULL,MAKEINTRESOURCE(101));

  if (leda_icon == NULL)
    leda_icon = CreateXpmIcon(32,32,win_icon1_xpm);

  leda_small_icon = CreateXpmIcon(16,16,win_small_icon1_xpm);



  // register window classes

  WNDCLASSEX    wndclass;

  wndclass.cbSize = sizeof(WNDCLASSEX);


  wndclass.lpszClassName = szAppName1;
  wndclass.style         = CS_OWNDC |CS_HREDRAW | CS_VREDRAW;
  wndclass.lpfnWndProc   = (WNDPROC)WndProc;
  wndclass.cbClsExtra    = 0;
  wndclass.cbWndExtra    = 4;
  wndclass.hInstance     = NULL;
  wndclass.hIcon         = leda_icon;
  wndclass.hIconSm       = leda_small_icon;
  wndclass.hCursor       = LoadCursor (NULL, IDC_ARROW);
  wndclass.hbrBackground = (HBRUSH)GetStockObject (WHITE_BRUSH);
  wndclass.lpszMenuName  = NULL;
  RegisterClassEx(&wndclass);

  wndclass.lpszClassName = szAppName2;
  wndclass.style         = CS_SAVEBITS | CS_OWNDC | CS_HREDRAW | CS_VREDRAW;
  wndclass.lpfnWndProc   = (WNDPROC)WndProc;
  wndclass.cbClsExtra    = 0;
  wndclass.cbWndExtra    = 4;
  wndclass.hInstance     = NULL;
//wndclass.hIcon         = LoadIcon (NULL, IDI_APPLICATION);
  wndclass.hIcon         = NULL;
  wndclass.hIconSm       = NULL;
  wndclass.hCursor       = LoadCursor (NULL, IDC_ARROW);
  wndclass.hbrBackground = (HBRUSH)GetStockObject (WHITE_BRUSH);
  wndclass.lpszMenuName  = NULL;
  RegisterClassEx(&wndclass);

}



int   x_root_color() { return white; }

char* x_root_pixrect(int x1 ,int y1, int x2, int y2)
{ HDC hdc = CreateDC("DISPLAY",NULL,NULL,NULL); 
  HDC hdcMem  = CreateCompatibleDC (hdc) ;
  x_pixrect* im = new x_pixrect;
  if (x1 > x2) SWAP(x1,x2);
  if (y1 > y2) SWAP(y1,y2);
  im->w = x2-x1+1;
  im->h = y2-y1+1;
  im->map = CreateCompatibleBitmap (hdc,im->w,im->h) ;
  SelectObject (hdcMem,im->map);
  BitBlt (hdcMem,0,0,im->w,im->h,hdc,x1,y1,SRCCOPY) ;
  DeleteDC(hdcMem) ;
  DeleteDC(hdc) ;
  return (char*)im;
}


void x_close_display(void) 
{ 
  if (taskbar_icon_win) 
  { x_del_taskbar_icon(taskbar_icon_win);
    taskbar_icon_win = 0;
   }

  DeleteObject(text_font);  
  DeleteObject(italic_font);  
  DeleteObject(bold_font);  
  DeleteObject(fixed_font); 
  DeleteObject(button_font); 
  DestroyCursor(dot_box_cursor);
  DestroyCursor(hand_cursor);
  DestroyIcon(leda_icon);
  DestroyIcon(leda_small_icon);

  text_font   = NULL;
  italic_font = NULL;
  bold_font   = NULL;
  fixed_font  = NULL;
  button_font = NULL;

  if (wlist[0]) delete wlist[0];

/*
  for(int i=0; i <= wcount; i++) 
    if (wlist[i]) delete wlist[i];
*/

}


static void display_info(void)
{ HDC hdc = CreateIC("DISPLAY",NULL,NULL,NULL); 
  int d = GetDeviceCaps(hdc,BITSPIXEL);
  int c = GetDeviceCaps(hdc,NUMCOLORS);
  DeleteDC(hdc);
  CGAL_CLIB_STD::printf("display: bitspixel = %d  numcolors = %d\n",c,d);
}


int x_display_width(void)  { return GetSystemMetrics(SM_CXSCREEN); }
int x_display_height(void) { return GetSystemMetrics(SM_CYSCREEN); }

int x_display_depth(void)
{ // color display: return something greater than 1
  return 4;
}


int x_display_bits_saved(void) { return 1; }


int x_create_buffer(Window win, int w, int h) 
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC  hdc  = wp->hdc;
  x_delete_buffer(win);

  HBITMAP hbm = CreateCompatibleBitmap(hdc,w+2,h+2);

  if (hbm == NULL) 
  { std::cerr << "\nCannot create bitmap w =" << w << " h =" << h << "\n";
    return 0;
   }

  SIZE sz;
  SetBitmapDimensionEx(hbm,w,h,&sz);

  HDC hdcMem = CreateCompatibleDC(hdc) ;
  wp->hbm = (HBITMAP)SelectObject(hdcMem,hbm);
  wp->hdcMem = hdcMem;

  HBRUSH hBrush = CreateSolidBrush(rgb_table[wp->B_COLOR]);

  RECT r;
  SetRect(&r,0,0,w+1,h+1);
  FillRect(hdcMem,&r,hBrush);
  DeleteObject(hBrush);
  ReleaseDC(hwnd,hdc);
  wp->buf_width = w;
  wp->buf_height= h;
  return 1;
}


int x_create_buffer(Window win) 
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  RECT r;
  GetClientRect(hwnd,&r);
  return x_create_buffer(win,r.right,r.bottom);
}


void x_delete_buffer(Window win)
{ win_struct* wp = wlist[win];
  HDC  hdcMem  = wp->hdcMem;
  if (hdcMem != NULL)
  { HWND    hwnd = wp->hwnd;
    HBITMAP  hbm = wp->hbm;
    HDC      hdc = GetDC(hwnd);
    DeleteObject(SelectObject(hdcMem,hbm));
    DeleteDC(hdcMem);
    wp->hdc = hdc;
    wp->hdcMem = NULL;
    wp->hbm = NULL;
    ReleaseDC(hwnd,hdc);
  }
}


int x_start_buffering(Window win) 
{ win_struct* wp = wlist[win];
  HDC  hdcMem = wp->hdcMem;
  int new_buf = 0;
  if (hdcMem == NULL) 
  { x_create_buffer(win);
    hdcMem = wp->hdcMem;
    new_buf = 1;
   }
  DeleteObject(SelectObject(hdcMem,CreatePenIndirect(&(wp->pen_data))));
  wp->hdc = hdcMem;
  return new_buf;
 }

int x_start_buffering(Window win, int w, int h) 
{ win_struct* wp = wlist[win];
  HDC  hdcMem = wp->hdcMem;
  x_stop_buffering(win);
  x_delete_buffer(win);
  x_create_buffer(win,w,h);
  hdcMem = wp->hdcMem;
  DeleteObject(SelectObject(hdcMem,CreatePenIndirect(&(wp->pen_data))));
  wp->hdc = hdcMem;
  return 1;
 }


void x_set_buffer(Window win, char* pr) 
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC  hdc = GetDC(hwnd);
  if (pr == 0) 
  { if (wp->hdc1)
    { DeleteDC(wp->hdc);
      wp->hdc = wp->hdc1;
      wp->hdc1 = 0;
     }
   }
  else
  { if (wp->hdc1)
  { std::cerr << "x_set_buffer: nested call.";
      FatalExit(0);
     }
    x_pixrect* im = (x_pixrect*)pr;
    HDC hdc1 = CreateCompatibleDC(hdc) ;
    SelectObject(hdc1,im->map);
    wp->hdc1 = wp->hdc;
    DeleteObject(SelectObject(hdc1,CreatePenIndirect(&(wp->pen_data))));
    wp->hdc = hdc1;
   }
  ReleaseDC(hwnd,hdc);
} 



int x_test_buffer(int win)
{ win_struct* wp = wlist[win];
  return wp->hwnd && (wp->hdcMem == wp->hdc); 
}


void x_flush_buffer(Window win,int x1,int y1,int x2,int y2,int xoff,int yoff)
{ win_struct* wp = wlist[win];
  
  HDC  hdcMem  = wp->hdcMem;

  if (hdcMem)
  { HWND hwnd    = wp->hwnd;
    HDC  hdc     = GetDC(hwnd);
    if (x1 > x2) SWAP(x1,x2);
    if (y1 > y2) SWAP(y1,y2);
    if (x1 < 0) x1 = 0;
    if (y1 < 0) y1 = 0;
    int w = x2-x1+1;
    int h = y2-y1+1;
    BitBlt(hdc,x1,y1,w,h,hdcMem,x1+xoff,y1+yoff,SRCCOPY) ;
    ReleaseDC(hwnd,hdc);
  }
}


void x_flush_buffer(int w, int x1, int y1, int x2, int y2) 
{ x_flush_buffer(w,x1,y1,x2,y2,0,0); }



void x_stop_buffering(Window win)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = GetDC(hwnd); 
  wp->hdc = hdc;
  DeleteObject(SelectObject(hdc,CreatePenIndirect(&(wp->pen_data))));
  ReleaseDC(hwnd,hdc);
}



void x_stop_buffering(Window win, char** pr)
{ win_struct* wp = wlist[win];
  x_stop_buffering(win);
  x_pixrect* im = NULL;
  HDC  hdcMem  = wp->hdcMem;
  if (hdcMem != NULL)
  { HWND    hwnd = wp->hwnd;
    HBITMAP  hbm = wp->hbm;
    HDC      hdc = GetDC(hwnd);
    im = new x_pixrect;
    im->map = (HBITMAP)SelectObject(hdcMem,hbm);
    SIZE sz;
    GetBitmapDimensionEx(im->map,&sz);
    im->w = sz.cx;
    im->h = sz.cy;
    DeleteDC(hdcMem);
    wp->hdc = hdc;
    wp->hdcMem = NULL;
    wp->hbm = NULL;
    ReleaseDC(hwnd,hdc);
  }
  *pr = (char*)im;
 }


/* windows */

void x_grab_pointer(Window win) 
{ win_struct* wp = wlist[win];
  grab_win = win; 
  SetCapture(wp->hwnd); 
}

void x_ungrab_pointer() 
{ grab_win = 0; 
  ReleaseCapture(); 
}

void x_set_focus(Window win) 
{ win_struct* wp = wlist[win];
  SetFocus(wp->hwnd); 
}

void x_move_pointer(Window win, int x, int y) 
{ x_window_to_screen(win,&x,&y);
  SetCursorPos(x,y);
}


void x_window_to_screen(Window win, int* x, int* y)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  POINT pt;
  pt.x = *x;
  pt.y = *y;
  ClientToScreen(hwnd,&pt);
  *x = pt.x;
  *y = pt.y;
}


void x_screen_to_window(Window win, int* x, int* y)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  POINT pt;
  pt.x = *x;
  pt.y = *y;
  ScreenToClient(hwnd,&pt);
  *x = pt.x;
  *y = pt.y;
}




static HWND center_hwnd(Window win, int w, int h)
{
  win_struct* wp = wlist[win];

  // open an invisible window centered over w

  HWND hwnd;

  if (win > 0) 
     hwnd = wp->hwnd;
  else
     hwnd = GetDesktopWindow();


  RECT r;
  GetClientRect(hwnd,&r);
  POINT pt;
  pt.x = (r.right - w)/2;
  pt.y = (r.bottom - h)/2;
  ClientToScreen(hwnd,&pt);

  return CreateWindow (szAppName2,"",WS_POPUPWINDOW,pt.x,pt.y,w,h,hwnd,
                                                           NULL,NULL,NULL);
}





int x_choose_color(Window win, int c)
{ //win_struct* wp = wlist[win];
  //HWND hwnd = wp->hwnd;
  HWND hwndc = center_hwnd(win,400,400); 

  CHOOSECOLOR cc;
  cc.lStructSize = sizeof(CHOOSECOLOR);
  cc.hwndOwner = hwndc;
  cc.hInstance = 0;
  cc.lpCustColors = rgb_custom;
  cc.rgbResult = rgb_table[c];
  cc.Flags = CC_RGBINIT | CC_FULLOPEN;

  int result = 0;

  for(int i=0; i<16; i++) rgb_custom[i] = rgb_table[i];

  if (ChooseColor(&cc)) 
  { rgb_table[c] = cc.rgbResult;
    for(int i=0; i<16; i++) rgb_table[i] = rgb_custom[i];
    //InvalidateRect(hwnd,NULL,TRUE);
    result = 1;
   }


  DestroyWindow(hwndc);
  return result;
}





int x_choose_file(Window win, int mode, const char* title, const char* filt, 
                                                           char* dname,
                                                           char* fname)
{ //win_struct* wp = wlist[win];
  //HWND hwnd = wp->hwnd;
  HWND hwndc = center_hwnd(win,430,290); 


  char file_path[256];

  CGAL_CLIB_STD::sprintf(file_path,"%s\\%s",dname,fname);

  // static char custfilt[128];
  // if (custfilt[0] == 0) strcpy(custfilt,filt);

  OPENFILENAME ofn;       // common dialog box structure

  ZeroMemory(&ofn, sizeof(OPENFILENAME));

  ofn.lStructSize = sizeof(OPENFILENAME);
//ofn.hwndOwner = GetDesktopWindow();
  ofn.hwndOwner = hwndc;
  ofn.lpstrFile = file_path;
  ofn.nMaxFile = 256;
  ofn.lpstrFilter = filt;
  ofn.nFilterIndex = 1;

  ofn.lpstrInitialDir = dname;

//ofn.nFilterIndex = 0;
//ofn.lpstrCustomFilter = custfilt;
//ofn.nMaxCustFilter = 128;

  ofn.lpstrFileTitle = NULL;
  ofn.nMaxFileTitle = 0;
  ofn.lpstrInitialDir = NULL;
  ofn.lpstrTitle = title;
  ofn.Flags = OFN_PATHMUSTEXIST | OFN_HIDEREADONLY | OFN_FILEMUSTEXIST 
                                                   | OFN_OVERWRITEPROMPT;
  
  int result = 0;
  if (mode == 0)
    { if (GetOpenFileName(&ofn))
      { int off = ofn.nFileOffset;
        strcpy(fname,file_path+off);
        file_path[off-1] = 0;
        strcpy(dname,file_path);
        result = 1;
       }
     }
  else
    { if (GetSaveFileName(&ofn))
      { int off = ofn.nFileOffset;
        strcpy(fname,file_path+off);
        file_path[off-1] = 0;
        strcpy(dname,file_path);
        result = 1;
       }
     }

  DestroyWindow(hwndc);
  return result;
}




static void choose_font(Window win, HFONT* hf, LOGFONT* lf)
{ win_struct* wp = wlist[win];
  HWND hwnd  = wp->hwnd;
  HWND hwndc = center_hwnd(win,440,550); 

  CHOOSEFONT cf;
  cf.lStructSize = sizeof(CHOOSEFONT);
  //cf.lpstrTitle = "Choose Font";
  cf.hwndOwner = hwndc;
  cf.lpLogFont = lf;
  cf.Flags = CF_SCREENFONTS | CF_INITTOLOGFONTSTRUCT | CF_NOSCRIPTSEL;
  if (ChooseFont(&cf))
  { DeleteObject(*hf);
    *hf = CreateFontIndirect(lf);
    InvalidateRect(hwnd,NULL,TRUE);
    CGAL_CLIB_STD::printf("\n");
    CGAL_CLIB_STD::printf("FaceName       = %s\n", lf->lfFaceName);
    CGAL_CLIB_STD::printf("Height         = %ld\n", lf->lfHeight);
    CGAL_CLIB_STD::printf("Width          = %ld\n", lf->lfWidth);
    CGAL_CLIB_STD::printf("Escapement     = %ld\n", lf->lfEscapement);
    CGAL_CLIB_STD::printf("Orientation    = %ld\n", lf->lfOrientation);
    CGAL_CLIB_STD::printf("Weight         = %ld\n", lf->lfWeight);
    CGAL_CLIB_STD::printf("Italic         = %d\n", lf->lfItalic);
    CGAL_CLIB_STD::printf("Underline      = %d\n", lf->lfUnderline);
    CGAL_CLIB_STD::printf("StrikeOut      = %d\n", lf->lfStrikeOut);
    CGAL_CLIB_STD::printf("CharSet        = %d\n", lf->lfCharSet);
    CGAL_CLIB_STD::printf("OutPrecision   = %d\n", lf->lfOutPrecision);
    CGAL_CLIB_STD::printf("ClipPrecision  = %d\n", lf->lfClipPrecision);
    CGAL_CLIB_STD::printf("Quality        = %d\n", lf->lfQuality);
    CGAL_CLIB_STD::printf("PitchAndFamily = %d\n", lf->lfPitchAndFamily);
    CGAL_CLIB_STD::printf("\n");
   }
  DestroyWindow(hwndc);
}



 

static void set_cursor(int i) 
{ 
  LPCSTR c;

  switch (i) {

  case XC_sb_v_double_arrow : 
           c = IDC_SIZENS;
           break;

  case XC_sb_h_double_arrow : 
           c = IDC_SIZEWE;
           break;

  case XC_fleur: 
           c = IDC_SIZEALL;
           break;

  case XC_crosshair: 
           c = IDC_CROSS;
           break;

  case XC_watch: 
           c = IDC_WAIT;
           break;

  case XC_hand2: 
           SetCursor(hand_cursor);
           return;

  case XC_dotbox: 
           SetCursor(dot_box_cursor);
           return;

  default: c = IDC_ARROW;
           break;
  }

  SetCursor(LoadCursor(NULL,c));
}




static void x_create_tmp_file(char* dname, char* fname)
{ char path[256];
  CGAL_CLIB_STD::sprintf(path,"%s%s",dname,fname);
  CreateFile(path,GENERIC_READ|GENERIC_WRITE,FILE_SHARE_READ,0,OPEN_ALWAYS,
                                                     FILE_ATTRIBUTE_NORMAL,0);
}


static void setup_fonts(Window win)
{
  char tmp_dir[256];
  GetTempPath(256,tmp_dir);

  CGAL_CLIB_STD::sprintf(tmp_dir+strlen(tmp_dir),"leda_fonts");

  CreateDirectory(tmp_dir,NULL);

  x_create_tmp_file(tmp_dir,"\\text.font");
  x_create_tmp_file(tmp_dir,"\\italic.font");
  x_create_tmp_file(tmp_dir,"\\bold.font");
  x_create_tmp_file(tmp_dir,"\\fixed.font");
  x_create_tmp_file(tmp_dir,"\\button.font");


  char fname[64];
  strcpy(fname,"text.font");
  if (x_choose_file(win,0,"LEDA Font Setup","Leda Fonts \0*.font\0",tmp_dir,
                                                                    fname))
  {
    if (strcmp(fname,"text.font")  == 0) 
      choose_font(win,&text_font,&text_lf);

    if (strcmp(fname,"italic.font")== 0) 
      choose_font(win,&italic_font,&italic_lf);

    if (strcmp(fname,"bold.font")  == 0) 
      choose_font(win,&bold_font,&bold_lf);

    if (strcmp(fname,"fixed.font") == 0) 
      choose_font(win,&fixed_font,&fixed_lf);

    if (strcmp(fname,"button.font")== 0) 
      choose_font(win,&button_font,&button_lf);
   }
}


static void setup_colors(Window win) 
{ 
  int c = white;
  x_choose_color(win,c); 
}




static long WINAPI WndProc (HWND hwnd, UINT message, UINT wParam, LONG lParam)
{
  if (trace_events)
     CGAL_CLIB_STD::printf("hwnd = %d mess = %d  wParam = %d lParam = %ld\n",
             (int)hwnd,message,wParam,lParam);

  int win = GetWindowLong(hwnd,0);

  cur_event.win = win;

  win_struct* wp = wlist[win];



  switch (message) { 

      case XWM_TASKBAR_ICON:
         {
           switch(lParam) {

           case WM_MOUSEMOVE:   break;
           case WM_RBUTTONDOWN: 
           case WM_LBUTTONDOWN: { POINT pos;
                                  GetCursorPos(&pos);
                                  TrackPopupMenu(taskbar_icon_menu,0,
                                                 pos.x,pos.y,0,hwnd,NULL);
                                  PostMessage(hwnd,WM_USER,0,0);
                                  break;
                                 }
           }
           break;
          }


      case WM_COMMAND:
        { 
          switch (LOWORD(wParam)) {

          case XWM_FONT_DIALOG:  setup_fonts(win);
                                 break;
          case XWM_COLOR_DIALOG: setup_colors(win);
                                 break;

          }
          break;
         }


      case WM_DROPFILES:
         { HDROP hDrop = (HDROP)wParam;
           int count = DragQueryFile(hDrop,0xFFFFFFFF,NULL,0);
           POINT pt;
           DragQueryPoint(hDrop,&pt);
           if (wp->drop_handler)
           { for(int i=0; i<count; i++)
             { char fname[80];
               DragQueryFile(hDrop,i,fname,sizeof(fname));
               wp->drop_handler(wp->inf,fname,pt.x,pt.y);
              }
            }
           else
           { char label[32];
             CGAL_CLIB_STD::sprintf(label,"Files dropped at (%ld,%ld)",pt.x,pt.y);
             char txt[1024];
             CGAL_CLIB_STD::sprintf(txt,"\n");
             for(int i=0; i<count; i++)
             { char fname[80];
               DragQueryFile(hDrop,i,fname,sizeof(fname));
               CGAL_CLIB_STD::sprintf(txt+strlen(txt),"%s\n",fname);
              }
             CGAL_CLIB_STD::sprintf(txt+strlen(txt),"\n");
             MessageBox(NULL,txt,label,MB_OK);
            }
           DragFinish(hDrop);
           break;
          }
      


      case WM_CREATE :
         { if (trace_events) 
             CGAL_CLIB_STD::printf("CREATE win = %d  wm_paint_count = %d\n", win,
                                      wp->wm_paint_count);
           break;
          }

      case WM_SETCURSOR:
         { if (trace_events) 
             CGAL_CLIB_STD::printf("SETCURSOR win = %d   hit code = %d\n",win,LOWORD(lParam));
           if (LOWORD(lParam) == HTCLIENT) 
           { set_cursor(wp->cursor);
             return 0;
            }
           break;
          }

      case WM_TIMER: 
         { if (trace_events) CGAL_CLIB_STD::printf("TIMER win = %d  id = %d\n", win, wParam);
           cur_event.kind = timer_event;
           cur_event.x = 0;
           cur_event.y = 0;
           cur_event.val = wParam;
           cur_event.t = 0; 
           return 0;
          }

      case WM_ERASEBKGND:
         { return 1;
          }

      case WM_PAINT :
         { //if (!GetUpdateRect(hwnd,NULL,FALSE)) return 0;

           RECT rect;
           GetUpdateRect(hwnd,&rect,TRUE);

            if (trace_events) CGAL_CLIB_STD::printf("PAINT win = %d  %ld %ld %ld %ld\n", win,
                                     rect.left,rect.top,rect.right,rect.bottom);


           if (wp->wm_paint_count) // skip first paint event
           { 
             int x = rect.left; 
             int y = rect.top;
             int w = rect.right - rect.left + 1;
             int h = rect.bottom - rect.top + 1;

             if (wp->repaint)
                (wp->repaint)(wp->inf,x,y,w,h);
             else
              { cur_event.kind = exposure_event;
                cur_event.x   = x;
                cur_event.y   = y;
                cur_event.val = w;
                cur_event.t   = (unsigned long)h;
               }
            }

           ValidateRect(hwnd,NULL);
           wp->wm_paint_count++;
           return 0;
          }


/*
      case WM_SIZE:
         { if (trace_events) CGAL_CLIB_STD::printf("SIZE  win = %d\n", w);
           if (w > wcount) 
           { //fprintf(stderr,"SIZE: win out of range.\n");
             return 0;
            }
           if (wp->hdcMem)
           { int buffering = (wp->hdc == wp->hdcMem);
             x_delete_buffer(w);
             x_create_buffer(w);
             if (buffering) wp->hdc = wp->hdcMem;
            }
           return 0;
          }
*/


      case WM_MOVE:
         { cur_event.kind = configure_event;
           cur_event.x = LOWORD (lParam);
           cur_event.y = HIWORD (lParam);
           cur_event.t = 0;

           if (trace_events) 
               CGAL_CLIB_STD::printf("MOVE win = %d  x = %d  y = %d\n", 
                       win, cur_event.x, cur_event.y);
           return 0;
          }


      case WM_MOUSEMOVE:
         { cur_event.kind = motion_event;
            short x = LOWORD (lParam);
            short y = HIWORD (lParam);
            cur_event.x = x;
            cur_event.y = y;
            cur_event.t = 0;

            if (trace_events) 
              CGAL_CLIB_STD::printf("MOUSEMOVE win = %d  x = %d  y = %d\n", win, cur_event.x, cur_event.y);
            return 0;
          }


      case WM_KILLFOCUS:
         { if (trace_events) CGAL_CLIB_STD::printf("KILLFOCUS\n");
           //if (grab_win == w) x_ungrab_pointer();
           return 0;
          }


      case WM_LBUTTONDBLCLK : // CGAL_CLIB_STD::printf("left double click\n");
      case WM_LBUTTONDOWN : {
           cur_event.kind = button_press_event;
           short x = LOWORD (lParam);
           short y = HIWORD (lParam);
           cur_event.x = x;
           cur_event.y = y;
           cur_event.t = 0;
           cur_event.val = 1;
           if (sys_key_down) cur_event.val = 2;
           if (wParam & MK_SHIFT) cur_event.val |= 256;
           if (wParam & MK_CONTROL) cur_event.val |= 512;

           //x_grab_pointer(win);

           if (trace_events)
            CGAL_CLIB_STD::printf("LBUTTONDOWN win = %d  x = %d  y = %d\n", win, cur_event.x, 
                                                                  cur_event.y);
           return 0;
         }

      case WM_MBUTTONDBLCLK : // printf("middle double click\n");
      case WM_MBUTTONDOWN : {
           cur_event.kind = button_press_event;
           short x = LOWORD (lParam);
           short y = HIWORD (lParam);
           cur_event.x = x;
           cur_event.y = y;
           cur_event.val = 2;
           cur_event.t = 0;
           if (wParam & MK_SHIFT) cur_event.val |= 256;
           if (wParam & MK_CONTROL) cur_event.val |= 512;

           if (trace_events)
            CGAL_CLIB_STD::printf("MBUTTONDOWN win = %d  x = %d  y = %d\n", win, cur_event.x, 
                                                                  cur_event.y);
           return 0;
         }


      case WM_RBUTTONDBLCLK : // CGAL_CLIB_STD::printf("right double click\n");
      case WM_RBUTTONDOWN : {
           cur_event.kind = button_press_event;
           short x = LOWORD (lParam);
           short y = HIWORD (lParam);
           cur_event.x = x;
           cur_event.y = y;
           cur_event.val = 3;
           cur_event.t = 0;
           if (sys_key_down) cur_event.val = 2;
           if (wParam & MK_SHIFT) cur_event.val |= 256;
           if (wParam & MK_CONTROL) cur_event.val |= 512;

           if (trace_events)
            CGAL_CLIB_STD::printf("RBUTTONDOWN win = %d  x = %d  y = %d\n", win, cur_event.x, 
                                                                  cur_event.y);
           return 0;
         }

      case WM_LBUTTONUP : {
           cur_event.kind = button_release_event;
           short x = LOWORD (lParam);
           short y = HIWORD (lParam);
           cur_event.x = x;
           cur_event.y = y;
           cur_event.val = 1;
           cur_event.t = 0;
           if (sys_key_down) cur_event.val = 2;
           if (wParam & MK_SHIFT) cur_event.val |= 256;
           if (wParam & MK_CONTROL) cur_event.val |= 512;

           //x_ungrab_pointer();

           if (trace_events)
             CGAL_CLIB_STD::printf("LBUTTONUP win = %d  x = %d  y = %d\n", win, cur_event.x, 
                                                                 cur_event.y);
           return 0;
         }

      case WM_MBUTTONUP : {
           cur_event.kind = button_release_event;
           short x = LOWORD (lParam);
           short y = HIWORD (lParam);
           cur_event.x = x;
           cur_event.y = y;
           cur_event.val = 2;
           cur_event.t = 0;
           if (wParam & MK_SHIFT) cur_event.val |= 256;
           if (wParam & MK_CONTROL) cur_event.val |= 512;
           if (trace_events)
             CGAL_CLIB_STD::printf("MBUTTONUP win = %d  x = %d  y = %d\n", win, cur_event.x, 
                                                                 cur_event.y);
           return 0;
         }

      case WM_RBUTTONUP : {
           cur_event.kind = button_release_event;
           short x = LOWORD (lParam);
           short y = HIWORD (lParam);
           cur_event.x = x;
           cur_event.y = y;
           cur_event.val = 3;
           cur_event.t = 0;
           if (sys_key_down) cur_event.val = 2;
           if (wParam & MK_SHIFT) cur_event.val |= 256;
           if (wParam & MK_CONTROL) cur_event.val |= 512;
           if (trace_events)
             CGAL_CLIB_STD::printf("RBUTTONUP win = %d  x = %d  y = %d\n", win, cur_event.x, 
                                                                 cur_event.y);
           return 0;
         }

      case WM_SYSKEYDOWN: {
              if (trace_events) 
                  CGAL_CLIB_STD::printf("SYSKEYDOWN win = %d wparam = %d\n",win,wParam);
              sys_key_down = 1;
              return 0;
      }

      case WM_SYSKEYUP: {
             if (trace_events) 
                  CGAL_CLIB_STD::printf("SYSKEYUP win = %d wparam = %d\n",win,wParam);
             cur_event.kind = key_release_event;
             cur_event.val = sys_key_down | 1024;
             cur_event.t = 0;
             sys_key_down = 0;
             return 0;
      }

      case WM_SYSCHAR: {
              if (trace_events) 
                  CGAL_CLIB_STD::printf("SYSCHAR win = %d code = %d\n",win,(TCHAR)wParam);
             sys_key_down = (TCHAR)wParam;
             cur_event.kind = key_press_event;
             cur_event.val = sys_key_down | 1024;
             cur_event.t = 0;
             return 0;
      }


      case WM_KEYUP:
      case WM_KEYDOWN: {

             if (message == WM_KEYDOWN)
             { cur_event.kind = key_press_event;
               if (trace_events) 
                  CGAL_CLIB_STD::printf("KEYDOWN win = %d  wParam = %x \n",win,wParam);
              }
             else
             { cur_event.kind = key_release_event;
               if (trace_events) 
                  CGAL_CLIB_STD::printf("KEYUP   win = %d  wParam = %x \n",win,wParam);
              }

             cur_event.t = 0;

             switch (wParam) {

              case VK_F12:   setup_fonts(win);
                             cur_event.kind = no_event;
                             break;

              case VK_F1:     cur_event.val = KEY_F1;
                              break;
              case VK_F2:     cur_event.val = KEY_F2;
                              break;
              case VK_F3:     cur_event.val = KEY_F3;
                              break;
              case VK_F4:     cur_event.val = KEY_F4;
                              break;
              case VK_F5:     cur_event.val = KEY_F5;
                              break;
              case VK_F6:     cur_event.val = KEY_F6;
                              break;
              case VK_F7:     cur_event.val = KEY_F7;
                              break;
              case VK_F8:     cur_event.val = KEY_F8;
                              break;

              case VK_F9:     //cur_event.val = KEY_F9;
                              if (message == WM_KEYUP) break;
                              trace_events = !trace_events;
                              CGAL_CLIB_STD::printf("TRACE EVENTS = %d\n",trace_events);
                              display_info();
                              break;


              case VK_HOME:   cur_event.val = KEY_HOME;
                              break;
              case VK_END:    cur_event.val = KEY_END;
                              break;
              case VK_UP:     cur_event.val = KEY_UP;
                              break;
              case VK_DOWN:   cur_event.val = KEY_DOWN;
                              break;
              case VK_LEFT:   cur_event.val = KEY_LEFT;
                              break;
              case VK_RIGHT:  cur_event.val = KEY_RIGHT;
                              break;
              case VK_ESCAPE: cur_event.val = KEY_ESCAPE;
                              //exit(0);
                              break;
              case VK_RETURN: cur_event.val = KEY_RETURN;
                              break;
              case VK_BACK:   cur_event.val = KEY_BACKSPACE;
                              break;

              case VK_SCROLL: 
              case VK_PRINT:  cur_event.val = KEY_PRINT;
                              break;

              default:        cur_event.kind = no_event;
                              break;

             }
             return 0;
           }



      case WM_CHAR: {
           char c = (char)wParam;
           if (trace_events) 
             CGAL_CLIB_STD::printf("CHAR win = %d  c = %d\n",win,c);
           if (isprint(c))
           { cur_event.kind = key_press_event;
             cur_event.val = c;
             cur_event.t = 0;
            }
           return 0;
         }
         


      case WM_CLOSE:
           if (trace_events) CGAL_CLIB_STD::printf("WM_CLOSE\n");
           cur_event.kind = destroy_event;
           return 0;


      case WM_DESTROY :
           if (trace_events) CGAL_CLIB_STD::printf("DESTROY: w = %d\n",win);
           //if (wlist[win] > 0 && GetParent(hwnd) == NULL) exit(0);
           return 0;

     }

  return DefWindowProc (hwnd, message, wParam, lParam);
}




Window x_create_window(void* inf, int width,int height,int bg_col, 
                       const char* label, const char* icon_label, 
                       int pwin, redraw_func redraw)
{
  int i = 1;
  while (i <= wcount && wlist[i] != 0) i++;

  if (i > wcount) wcount = i;

  if (wcount >= MAX_WIN) 
  { std::cerr << "\n Too many Windows (" << wcount << ").\n";
    FatalExit(0);
   }


   win_struct* wp = new win_struct(inf,width,height,bg_col,label,redraw);
   wlist[i] = wp;
   return i;
}


void x_set_clip_rectangle(Window win, int x0, int y0, int w, int h) 
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  if (hwnd && wp->hdc)
  { HDC hdc = wp->hdc;
    HRGN hrgn = CreateRectRgn(x0,y0,x0+w,y0+h);
    SelectClipRgn(hdc,hrgn);
    DeleteObject(hrgn);
    ReleaseDC(hwnd,hdc);
   }
 }


void x_bezier_ellipse(Window win, int x, int y, int r1, int r2)
{
  win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;

  int f1 = int(1.33*r1);
  int f2 = int(0.93*r2);

  POINT p[7];
  p[0].x = x;    p[0].y = y-r2;
  p[1].x = x+f1; p[1].y = y-f2;
  p[2].x = x+f1; p[2].y = y+f2;
  p[3].x = x;    p[3].y = y+r2;
  p[4].x = x-f1; p[4].y = y+f2;
  p[5].x = x-f1; p[5].y = y-f2;
  p[6].x = x;    p[6].y = y-r2;

  PolyBezier(hdc,p,7);
  ReleaseDC(hwnd,hdc);
}




void x_clip_mask_ellipse(Window win, int x, int y, int r1, int r2, int mode) 
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;
  if (hwnd == NULL) return;

  OSVERSIONINFO osinf;
  osinf.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
  GetVersionEx(&osinf);

  BeginPath(hdc);
  if (osinf.dwMajorVersion >= 4 && osinf.dwMinorVersion == 0)
    Ellipse(hdc,x-r1,y-r2,x+r1+1,y+r2+1);
  else 
    x_bezier_ellipse(win,x,y,r1,r2);
  EndPath(hdc);

  switch (mode) {
  case 0: SelectClipPath(hdc,RGN_DIFF);
          break;
  case 1: SelectClipPath(hdc,RGN_OR);
          break;
  case 2: SelectClipPath(hdc,RGN_COPY);
          break;
  }

  ReleaseDC(hwnd,hdc);
}



void x_clip_mask_polygon(Window win, int n, int* xcoord, int* ycoord, int mode) 
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;

  if (hwnd == NULL) return;

  if (n == 0)
  { SelectClipRgn(hdc,NULL);
    return;
   }

  POINT* p = new POINT[n];
  for(int i=0; i < n; i++) 
  { p[i].x = xcoord[i];
    p[i].y = ycoord[i];
   }

  BeginPath(hdc);
  Polygon(hdc,p,n);
  EndPath(hdc);

  switch (mode) {
  case 0: SelectClipPath(hdc,RGN_DIFF);
          break;
  case 1: SelectClipPath(hdc,RGN_OR);
          break;
  case 2: SelectClipPath(hdc,RGN_COPY);
          break;
  }

  ReleaseDC(hwnd,hdc);
  delete[] p;
 }



void x_set_clip_polygon(Window win, int n, int* xcoord, int* ycoord) 
{ x_clip_mask_polygon(win,n,xcoord,ycoord,2); }


void x_set_clip_ellipse(Window win, int x, int y, int r1, int r2) 
{ x_clip_mask_ellipse(win,x,y,r1,r2,2); }


int  x_window_opened(Window win) 
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  return hwnd && IsWindowVisible(hwnd); 
 }


void x_resize_window(Window win, int xpos, int ypos, int width,int height,int)  
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
//HWND hwnd_parent = GetParent(hwnd);

  char class_name[32];
  GetClassName(hwnd,class_name,32);


  if (strcmp(class_name,szAppName1) == 0)
  { width  += 2*GetSystemMetrics(SM_CXFRAME);
    height += 2*GetSystemMetrics(SM_CYFRAME);
    height += GetSystemMetrics(SM_CYCAPTION);
   }
  else
  { width  += 2*GetSystemMetrics(SM_CXBORDER);
    height += 2*GetSystemMetrics(SM_CYBORDER);
   }

  SetWindowPos(hwnd,HWND_TOP,xpos,ypos,width,height,
                                          SWP_DRAWFRAME | SWP_NOACTIVATE);

  x_clear_window(win,0,0,width,height);

  if (wp->repaint)
     (wp->repaint)(wp->inf,0,0,width,height);

  wp->wm_paint_count = 0;
}


void x_open_window(Window win, int x, int y, int width,int height,int pwin) 
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;

  HWND hwndparent = 0;
  DWORD win_style = 0;
  DWORD win_ex_style = 0;
  char* class_name = 0;

  int w0 = width;
  int h0 = height;

  if (width  < 0) width  = -width;
  if (height < 0) height = -height;
  
  if (pwin == 0)
  { width  += 2*GetSystemMetrics(SM_CXFRAME);
    height += 2*GetSystemMetrics(SM_CYFRAME);
    height += GetSystemMetrics(SM_CYCAPTION);
    hwndparent = NULL;
    win_style = WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN;
    win_ex_style = WS_EX_APPWINDOW;
    class_name = szAppName1;
   }
  else
  { if (pwin < 0)
    { width  += 2*GetSystemMetrics(SM_CXBORDER);
      height += 2*GetSystemMetrics(SM_CYBORDER);
      hwndparent = NULL;
      win_style = WS_POPUPWINDOW | WS_CLIPCHILDREN;
      win_ex_style = WS_EX_TOOLWINDOW;
      class_name = szAppName2;
     }

    if (pwin > 0)
    { width  += 2*GetSystemMetrics(SM_CXBORDER);
      height += 2*GetSystemMetrics(SM_CYBORDER);
      hwndparent = wlist[pwin]->hwnd;
      win_style = WS_CHILD | WS_CLIPSIBLINGS | WS_CLIPCHILDREN; 
      win_ex_style = WS_EX_TOOLWINDOW;
      if (wp->border_w > 0)
         win_style |= WS_BORDER; 
      else
        { width  -= 2;
          height -= 2;
         }
      class_name = szAppName2;
     }
   }

  if (w0 < 0) x = x - (width-1);
  if (h0 < 0) y = y - (height-1);


  // restrict window size to display size

  int disp_width  = x_display_width();
  int disp_height = x_display_height();

  int dy = disp_height - height;
  if (dy < 0) 
  { float f = float(disp_height)/height; 
    height = int(f*height); 
    width  = int(f*width);
    dy = 0; 
   }

  int dx = disp_width - width;
  if (dx < 0) 
  { float f = float(disp_width)/width; 
    width  = int(f*width);
    height = int(f*height); 
    dx = 0; 
   }

  if (pwin == 0)
  { if ( x < 0 ) x = 0;
    if ( y < 0 ) y = 0;
    if ( x > dx ) x = dx;
    if ( y > dy ) y = dy;
   }


  if (!hwnd)
  { 
    //win_style |= WS_VISIBLE;

    hwnd = CreateWindowEx(win_ex_style,         // extended style
                          class_name,           // window class name
                          wp->header,           // frame label
                          win_style,            // window style
                          x,y,                  // window position
                          width,height,         // window size
                          hwndparent,           // handle parent window 
                          NULL,                 // handle menu
                          NULL,                 // handle program copy
                          NULL);                // special parameter

    wp->hwnd = hwnd;
    SetWindowLong(hwnd,0,win);
  }

  wp->wm_paint_count = 0;


  if (pwin == 0)
   { MoveWindow(hwnd,x,y,width,height,FALSE);
     ShowWindow (hwnd, SW_SHOWNORMAL);
    }
  else
   { //HWND save_focus = GetFocus();
     SetParent(hwnd,hwndparent);
     if (pwin < 0) 
     { //SetFocus(save_focus);
       SetWindowPos(hwnd,HWND_TOP,x,y,width,height,
                                            SWP_NOACTIVATE|SWP_SHOWWINDOW);
       //SetFocus(hwnd);
      }
     if (pwin > 0) 
     { //MoveWindow(hwnd,x,y,width,height,FALSE);
       //ShowWindow(hwnd,SW_SHOWNOACTIVATE);
       SetWindowPos(hwnd,HWND_TOP,x,y,width,height,
                                            SWP_NOACTIVATE|SWP_SHOWWINDOW);
       //SetFocus(hwnd);
      }
    }


  while (wp->wm_paint_count==0)
  { MSG msg;
    GetMessage(&msg,NULL,0,0);
    TranslateMessage(&msg);
    DispatchMessage(&msg);
   }


  HDC hdc = GetDC(hwnd);

  wp->font = text_font;
  //SelectObject(hdc, text_font);

  wp->hdc = hdc;
  ReleaseDC(hwnd,hdc);

  x_set_color(win,wp->COLOR);
  x_set_mode(win,wp->MODE);
  x_set_text_mode(win,wp->TMODE);
  x_set_line_width(win,wp->LWIDTH);
  x_set_line_style(win,wp->LSTYLE);

/*
  if (pwin == 0 && taskbar_icon_win == 0) 
  { x_add_taskbar_icon(hwnd,leda_small_icon,wp->header);
    taskbar_icon_win = hwnd;
   }
*/

  DragAcceptFiles(hwnd,TRUE);
}


void x_close_window(Window win) 
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;

  if (taskbar_icon_win == hwnd) 
  { x_del_taskbar_icon(hwnd);
    taskbar_icon_win = 0;
   }

  ReleaseDC(hwnd,hdc);
  if (hwnd) ShowWindow(hwnd,SW_HIDE);
  if (grab_win == win) x_ungrab_pointer();
  wp->hdc = NULL;
 }


void x_destroy_window(Window win) 
{ 
  win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;

  x_delete_buffer(win);
  x_stop_timer(win);

  wlist[win] = 0;

  if (wp->hmf) DeleteEnhMetaFile(wp->hmf);

  if (hwnd) DestroyWindow(hwnd);

  delete wp;

  if (wcount == win) wcount--;
 }


int x_set_cursor(Window win, int i) 
{ win_struct* wp = wlist[win];
  if (wp->cursor == i) return i;
  int i0 = wp->cursor;
  wp->cursor = i;
  set_cursor(i);
  return i0;
 }

void x_set_label(Window win, const char *s)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  SetWindowText(hwnd,s);
 }



void x_start_timer(Window win, int msec)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  x_stop_timer(win);
  wp->timer_id = SetTimer(hwnd,1,msec,NULL);
 }

void x_stop_timer(Window win)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  UINT id = wp->timer_id;
  if (id) KillTimer(hwnd,id);
  wp->timer_id = 0;
 }



void* x_window_inf(Window win) 
{ win_struct* wp = wlist[win];
  return wp->inf; 
}

int x_window_width(Window win)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  if (hwnd)
     { RECT r;
       GetClientRect(wp->hwnd,&r);
       return r.right; 
      }
   else
     return wp->width;
 }

int x_window_height(Window win)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  if (hwnd)
     { RECT r;
       GetClientRect(wp->hwnd,&r);
       return r.bottom; 
      }
   else
     return wp->height;
 }


void x_window_frame(Window win, int* x0, int* y0, int* x1, int* y1)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  if (hwnd)
    { RECT r;
      GetWindowRect(hwnd,&r);
      *x0 = r.left;
      *y0 = r.top;
      *x1 = r.right-1;
      *y1 = r.bottom-1;
     }
  else
    { *x0 = 0;
      *y0 = 0;
      *x1 = 0;
      *y1 = 0;
     }
}


void x_clear_window(Window win, int x0, int y0, int x1, int y1, int xorig, 
                                                                int yorig)
{ 
  if (!x_window_opened(win)) return;

  win_struct* wp = wlist[win];
 
  char* pm = wp->B_PIXMAP;

  if (x0 > x1)  SWAP(x0,x1);
  if (y0 > y1)  SWAP(y0,y1);

  if (pm == 0)
  { int save_col = x_set_color(win,wp->B_COLOR);
    x_box(win,x0,y0,x1,y1);
    x_set_color(win,save_col);
    return;
   }

  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;
  HDC hdcMem  = CreateCompatibleDC (hdc) ;

  x_pixrect* im = (x_pixrect*)pm;
  int w = im->w;
  int h = im->h;

  SelectObject (hdcMem,im->map);


  if (xorig > 0)
     while (xorig > 0) xorig -= w;
  else
     while (xorig+w < 0) xorig += w;

  if (yorig > 0)
     while (yorig > 0) yorig -= h;
  else
     while (yorig+h < 0) yorig += h;

  int xmax,ymax;
  if (x_test_buffer(win))
  { xmax = wp->buf_width;
    ymax = wp->buf_height;
   }
  else
  { RECT r;
    GetClientRect(hwnd,&r);
    xmax = r.right;
    ymax = r.bottom;
   }

  for(int y = yorig;  y < ymax; y += h)
    for(int x = xorig; x < xmax; x += w)
      if (x < x1 && x+w > x0 && y < y1 && y+h >y0)
        BitBlt(hdc,x,y,w,h,hdcMem,0,0,SRCCOPY);

  DeleteDC (hdcMem) ;
  ReleaseDC(hwnd,hdc);
}

void x_clear_window(Window win, int x0, int y0, int x1, int y1)
{ x_clear_window(win,x0,y0,x1,y1,0,0); }


//------------------------------------------------------------------------------
// drawing
//------------------------------------------------------------------------------

inline void h_line(HDC hdc, int x1, int x2, int y) 
{ MoveToEx(hdc,x1,y,NULL); 
  LineTo(hdc,x2+1,y); 
 }


static void adjust_line(int s, int& x1, int& y1, int& x2, int& y2)
{ int dx = x2 - x1;
  int dy = y2 - y1;
  if (dx == 0 && dy == 0) return;

  int xoff = s;
  int yoff = s;

  if (dx < 0) { dx = -dx; xoff = -s; }
  if (dy < 0) { dy = -dy; yoff = -s; }

  if ( dx >= dy) x2 += xoff;
  if ( dx <= dy) y2 += yoff;
 }



void x_line(Window win, int x1, int y1, int x2, int y2)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;

  int jstyle = wp->JSTYLE;
  int lw     = wp->LWIDTH;
  int jleft  = (jstyle & 1);
  int jright = (jstyle & 2);

  if (x1 > x2 || (x1 == x2 && y1 > y2)) 
  { SWAP(x1,x2);
    SWAP(y1,y2);
    SWAP(jleft,jright);
   }

  if (!jleft ) adjust_line(-lw/2-1,x2,y2,x1,y1);
  if (!jright) adjust_line(-lw/2,  x1,y1,x2,y2);
  if ( jright) adjust_line( lw/2,  x1,y1,x2,y2);

  MoveToEx(hdc,x1,y1,NULL);
  LineTo(hdc,x2,y2);
  ReleaseDC(hwnd,hdc);
 }


void x_lines(int w, int n, int *x1, int *y1, int* x2, int* y2)
{ while (n--) x_line(w,*x1++,*y1++,*x2++,*y2++); }



void x_rect(Window win, int x1, int y1, int x2, int y2)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC  hdc = wp->hdc;
  if (x1 > x2) SWAP(x1,x2);
  if (y1 > y2) SWAP(y1,y2);
  HBRUSH save_brush = (HBRUSH)SelectObject(hdc,GetStockObject(NULL_BRUSH));
  Rectangle(hdc,x1,y1,x2+1,y2+1);
  SelectObject(hdc,save_brush);
  ReleaseDC(hwnd,hdc);
 }


void x_box(Window win, int x1, int y1, int x2, int y2)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC  hdc = wp->hdc;
  if (x1 > x2) SWAP(x1,x2);
  if (y1 > y2) SWAP(y1,y2);
  COLORREF col = wp->pen_data.lopnColor;

  RECT r;
  SetRect(&r,x1,y1,x2+1,y2+1);

  HBRUSH hBrush;
  if (wp->stip_bm)
   { hBrush = CreatePatternBrush(wp->stip_bm);
     SetBkColor(hdc,rgb_table[wp->stip_bgcol]);
    }
  else
    hBrush = CreateSolidBrush(col);

 
  FillRect(hdc,&r,hBrush);
  DeleteObject(hBrush);

  ReleaseDC(hwnd,hdc);
}


void x_circle0(Window win, int x0,int y0,int r0)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;
  int col = wp->pen_data.lopnColor;
  int lw  = wp->LWIDTH;

  int rin  = r0-lw/2;
  int rout = r0+lw/2;

  if ( (lw % 2) == 0) rout--;

  for (int r = rin; r <= rout; r++)
  { int y = r;
    int x = 0;
    int e = 3-2*y;

    SetPixel(hdc,x0,y0+r,col);
    SetPixel(hdc,x0,y0-r,col);
    SetPixel(hdc,x0+r,y0,col);
    SetPixel(hdc,x0-r,y0,col);


    for (x=1;x<y;)
      { SetPixel(hdc,x0+x,y0+y,col);
        SetPixel(hdc,x0+x,y0-y,col);
        SetPixel(hdc,x0-x,y0+y,col);
        SetPixel(hdc,x0-x,y0-y,col);
        SetPixel(hdc,x0+y,y0+x,col);
        SetPixel(hdc,x0+y,y0-x,col);
        SetPixel(hdc,x0-y,y0+x,col);
        SetPixel(hdc,x0-y,y0-x,col);
        x++;
        if (e>=0) { y--; e = e - 4*y; }
        e = e + 4*x + 2;
       }

    if (x == y)
    { SetPixel(hdc,x0+x,y0+y,col);
      SetPixel(hdc,x0+x,y0-y,col);
      SetPixel(hdc,x0-x,y0+y,col);
      SetPixel(hdc,x0-x,y0-y,col);
     }
  }
 ReleaseDC(hwnd,hdc);
}


void x_ellipse(Window win, int x, int y, int r1, int r2)
{ win_struct* wp = wlist[win];
  
  if (r1 == r2 && wp->LWIDTH == 1 && wp->hdc2 == NULL)
  { x_circle0(win,x,y,r1);
    return;
   }
   
  HWND hwnd = wp->hwnd;
  HDC  hdc = wp->hdc;
  HBRUSH save_brush = (HBRUSH)SelectObject(hdc,GetStockObject(NULL_BRUSH));
  Ellipse(hdc,x-r1,y-r2,x+r1+1,y+r2+1);
  SelectObject(hdc,save_brush);
  ReleaseDC(hwnd,hdc);
}



void x_fill_ellipse(Window win, int x, int y, int r1, int r2)
{ win_struct* wp = wlist[win];
  int lw = x_set_line_width(win,1); 
  HWND hwnd = wp->hwnd;
  HDC  hdc = wp->hdc;
  COLORREF col = wp->pen_data.lopnColor;

  HBRUSH hBrush;
  if (wp->stip_bm)
   { hBrush = CreatePatternBrush(wp->stip_bm);
     SetBkColor(hdc,rgb_table[wp->stip_bgcol]);
    }
  else
    hBrush = CreateSolidBrush(col);

  HBRUSH save_brush = (HBRUSH)SelectObject(hdc,hBrush);
  
  Ellipse(hdc,x-r1,y-r2,x+r1+1,y+r2+1);

  DeleteObject(SelectObject(hdc,save_brush));
  ReleaseDC(hwnd,hdc);
  x_set_line_width(win,lw);
}



void x_circle(Window win, int x, int y, int r)
{ x_ellipse(win,x,y,r,r); }


void x_fill_circle(Window win, int x, int y, int r)
{ x_fill_ellipse(win,x,y,r,r); }


void x_polyline(Window win, int n, int* xcoord, int* ycoord, int adjust)
{ if (adjust == 1) return;
  win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;
  POINT* p = new POINT[n];
  for(int i=0; i < n; i++) 
  { p[i].x = xcoord[i];
    p[i].y = ycoord[i];
   }
  Polyline(hdc,p,n);
  ReleaseDC(hwnd,hdc);
  delete[] p;
}


void x_polygon(Window win, int n, int* xcoord, int* ycoord)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC  hdc = wp->hdc;
  POINT* p = new POINT[n];
  for(int i=0; i < n; i++) 
  { p[i].x = xcoord[i];
    p[i].y = ycoord[i];
   }

  HBRUSH save_brush = (HBRUSH)SelectObject(hdc,GetStockObject(NULL_BRUSH));

  Polygon(hdc,p,n);

  DeleteObject(SelectObject(hdc,save_brush));

  ReleaseDC(hwnd,hdc);
  delete[] p;
}




void x_fill_polygon(Window win, int n, int* xcoord, int* ycoord)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC  hdc = wp->hdc;

  COLORREF col = wp->pen_data.lopnColor;

  POINT* p = new POINT[n];
  for(int i=0; i < n; i++) 
  { p[i].x = xcoord[i];
    p[i].y = ycoord[i];
   }

  LOGPEN lp = wp->pen_data;
  lp.lopnStyle = PS_NULL;
  DeleteObject(SelectObject(hdc,CreatePenIndirect(&lp)));

  HBRUSH hBrush;
  if (wp->stip_bm)
    { hBrush = CreatePatternBrush(wp->stip_bm);
      SetBkColor(hdc,rgb_table[wp->stip_bgcol]);
     }
  else
     hBrush = CreateSolidBrush(col);

  HBRUSH save_brush = (HBRUSH)SelectObject(hdc,hBrush);

  SetPolyFillMode(hdc,ALTERNATE);
  Polygon(hdc,p,n);

  DeleteObject(SelectObject(hdc,save_brush));
  DeleteObject(SelectObject(hdc,CreatePenIndirect(&wp->pen_data)));

  ReleaseDC(hwnd,hdc);
  delete[] p;
}




void x_pixel(Window win, int x, int y)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC  hdc = wp->hdc;
  SetPixel(hdc,x,y,wp->pen_data.lopnColor);
  ReleaseDC(hwnd,hdc);
}


int x_get_pixel(Window win, int x, int y)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC  hdc = wp->hdc;
  int col = GetPixel(hdc,x,y); 
  ReleaseDC(hwnd,hdc);
  return col;
}



void x_pixels(Window win, int n, int* x, int* y)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC  hdc = wp->hdc;
  while (n--) SetPixel(hdc,x[n],y[n],wp->pen_data.lopnColor);
  ReleaseDC(hwnd,hdc);
}


void x_point(Window win, int x, int y)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC  hdc = wp->hdc;
  COLORREF col = wp->pen_data.lopnColor;
  for(int i = -2; i <= 2; i++)
  { SetPixel(hdc,x+i,y+i,col);
    SetPixel(hdc,x+i,y-i,col);
   }
  ReleaseDC(hwnd,hdc);
}


void x_plus(Window win, int x, int y)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC  hdc = wp->hdc;
  COLORREF col = wp->pen_data.lopnColor;
  for(int i = -3; i <= 3; i++)
  { SetPixel(hdc,x,y+i,col);
    SetPixel(hdc,x+i,y,col);
   }
  ReleaseDC(hwnd,hdc);
}




void x_arc(Window win,int mx,int my,int r1,int r2,double start,double angle)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;

  if (angle > 0)
  { start += angle;
    angle = - angle;
   }

  int R  = r1+r2;
  int x0 = int(mx + R*cos(start+angle));
  int y0 = int(my - R*sin(start+angle));
  int x1 = int(mx + R*cos(start));
  int y1 = int(my - R*sin(start));
  Arc(hdc,mx-r1,my-r2,mx+r1+1,my+r2+1,x0,y0,x1,y1);
  ReleaseDC(hwnd,hdc);
}


void x_text_underline(Window win, int x, int y, const char* s, int l, int r)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC  hdc = wp->hdc;
  y += (x_text_height(win,s)-1);
  int x1 = x+x_text_width(win,s,l);
  int x2 = x+x_text_width(win,s,r+1);
  MoveToEx(hdc,x1,y,NULL);
  LineTo(hdc,x2,y);
  ReleaseDC(hwnd,hdc);
 }


void x_text(Window win, int x, int y, const char* s0, int l)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC  hdc = wp->hdc;
  COLORREF col = wp->pen_data.lopnColor;
  COLORREF bcol = rgb_table[wp->B_COLOR];

  if (wp->font == fixed_font) y++;

  if (l > (int)strlen(s0)) l = strlen(s0);

  const char* s = s0;

  HFONT save_font = (HFONT)SelectObject(hdc,wp->font);

  if (wp->TMODE == opaque)
  { SetBkColor(hdc,bcol);
    TextOut(hdc,x,y,s,l);
    SelectObject(hdc,save_font);
    ReleaseDC(hwnd,hdc);
    return;
   }

  if (wp->MODE == src_mode)
  { SetBkMode(hdc,TRANSPARENT);
    TextOut(hdc,x,y,s,l);
    SetBkMode(hdc,OPAQUE);
    SelectObject(hdc,save_font);
    ReleaseDC(hwnd,hdc);
    return;
   }

  SIZE sz;
  GetTextExtentPoint(hdc,s,l,&sz);
  SelectObject(hdc,save_font);

  int w = sz.cx;
  int h = sz.cy;

  HDC     hdcMem    = CreateCompatibleDC(hdc) ;
  HBITMAP save_hbm  = (HBITMAP)SelectObject(hdcMem,CreateCompatibleBitmap(hdc,w,h));

  save_font = (HFONT)SelectObject(hdcMem,wp->font);

  SetBkMode(hdcMem,OPAQUE);
  SetBkColor(hdcMem,RGB(0,0,0));
  SetTextColor(hdcMem,RGB(255,255,255));
  TextOut(hdcMem,0,0,s,l);


  HBRUSH save_brush = (HBRUSH)SelectObject(hdc,CreateSolidBrush(col));

  if (wp->MODE == xor_mode)
     BitBlt(hdc,x,y,w,h,hdcMem,0,0,0x00A60706);
  else
     BitBlt(hdc,x,y,w,h,hdcMem,0,0,0x00E20746);

  DeleteObject(SelectObject(hdcMem,save_hbm));
  SelectObject(hdcMem,save_font);
  DeleteDC(hdcMem);

  DeleteObject(SelectObject(hdc,save_brush));
  ReleaseDC(hwnd,hdc);
}

  
void x_text(Window win, int x, int y, const char* s)
{ x_text(win, x, y, s, strlen(s)); }



void x_ctext(Window win, int x, int y, const char* s)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;
  HFONT save_font = (HFONT)SelectObject(hdc,wp->font);
  SIZE sz;
  GetTextExtentPoint(hdc,s,strlen(s),&sz);
  SelectObject(hdc,save_font);
  ReleaseDC(hwnd,hdc);
  x_text(win,x-sz.cx/2,y-sz.cy/2/*+sz.cy%2*/,s);
}


void x_ctext_underline(Window win, int x, int y, const char* s, int l, int r)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;
  HFONT save_font = (HFONT)SelectObject(hdc,wp->font);
  SIZE sz;
  GetTextExtentPoint(hdc,s,strlen(s),&sz);
  SelectObject(hdc,save_font);
  ReleaseDC(hwnd,hdc);
  x_text_underline(win,x-sz.cx/2,y-sz.cy/2/*+sz.cy%2*/,s,l,r);
}



//------------------------------------------------------------------------------
// pixrects and bitmaps 
//------------------------------------------------------------------------------

static char rev_byte(char c)
{ char c1 = 0x00;
  if (c & 0x01) c1 |= 0x80;
  if (c & 0x02) c1 |= 0x40;
  if (c & 0x04) c1 |= 0x20;
  if (c & 0x08) c1 |= 0x10;
  if (c & 0x10) c1 |= 0x08;
  if (c & 0x20) c1 |= 0x04;
  if (c & 0x40) c1 |= 0x02;
  if (c & 0x80) c1 |= 0x01;
  return c1;
 }


char* x_create_bitmap(Window win, int w, int h, unsigned char* data)
{ 
//win_struct* wp = wlist[win];
//HWND hwnd = wp->hwnd;
//HDC   hdc = wp->hdc;

  int   bw0 = (w+7)/8;
  int   bw1 = 2*((bw0+1)/2);
  char* buf = new char[bw1*h];
  char* p   = buf;
  char* q   = (char*)data;

  for(int i=0; i<h; i++)
    for(int j=0; j<bw1; j++) 
       *p++ = (j >= bw0) ? 0 : rev_byte(*q++);

  BITMAP bm;
  bm.bmType = 0;
  bm.bmWidth = w;
  bm.bmHeight = h;
  bm.bmWidthBytes = bw1;
  bm.bmPlanes = 1;
  bm.bmBitsPixel = 1;
  bm.bmBits = buf;

  x_pixrect* im = new x_pixrect;
  im->w = w;
  im->h = h;
  im->map = CreateBitmapIndirect(&bm);

  delete[] buf;

//ReleaseDC(hwnd,hdc);
  return (char*)im;
}


void x_pixrect_to_ps(Window win, char* prect, char* buf, int full_color,
                                                         void (*progress)(int))
{ 
  HDC hdc = CreateDC("DISPLAY",NULL,NULL,NULL); 
  HDC hdcMem  = CreateCompatibleDC (hdc) ;
  x_pixrect* im = (x_pixrect*)prect;
  SelectObject (hdcMem,im->map);

  for(int y=0; y<im->h; y++)
  { if (progress) progress(y);
    for(int x=0; x<im->w; x++)
    { COLORREF col = GetPixel(hdcMem,x,y); 
      int r = GetRValue(col);
      int g = GetGValue(col);
      int b = GetBValue(col);
      if (full_color)
       { CGAL_CLIB_STD::sprintf(buf,"%02x%02x%02x",r,g,b);
         buf += 6;
        }
      else
       { CGAL_CLIB_STD::sprintf(buf,"%02x",(r+g+b)/3);
         buf += 2;
        }
     }
   }

  DeleteDC(hdcMem);
  DeleteDC(hdc);
}


void x_pixrect_to_matrix(Window, char* prect, int* matrix)
{ 
  HDC hdc = CreateDC("DISPLAY",NULL,NULL,NULL); 
  HDC hdcMem = CreateCompatibleDC (hdc) ;
  x_pixrect* im = (x_pixrect*)prect;
  SelectObject(hdcMem,im->map);

  for(int y=0; y<im->h; y++)
    for(int x=0; x<im->w; x++)
    { COLORREF col = GetPixel(hdcMem,x,y);
      int i = 0;
      while (i < rgb_count && rgb_table[i] != col) i++;
      if (i == rgb_count)
      { rgb_table[i] = col;
        rgb_count++;
       }
      *matrix++ = i;
     }

  DeleteDC(hdcMem);
  DeleteDC(hdc);
}


void x_matrix_to_pixrect(Window,int* matrix, int width, int height, char** pr)
{ 
  HDC hdc = CreateDC("DISPLAY",NULL,NULL,NULL); 

  x_pixrect* im = new x_pixrect;
  im->w = width;
  im->h = height;
  im->map = CreateCompatibleBitmap(hdc,im->w,im->h) ;

  HDC hdcMem = CreateCompatibleDC(hdc) ;
  SelectObject(hdcMem,im->map);

  for(int y=0; y<height; y++)
    for(int x=0; x<width; x++)
      SetPixel(hdcMem,x,y,*matrix++); 

  DeleteDC(hdcMem);
  DeleteDC(hdc);

  *pr = (char*)im;
}


void x_pixrect_to_clipboard(Window win, char* rect)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;
  x_pixrect* im = (x_pixrect*)rect;
  OpenClipboard(hwnd);
  EmptyClipboard();
  SetClipboardData(CF_BITMAP,im->map);
  im->map = 0;
  CloseClipboard();
  ReleaseDC(hwnd,hdc);
 }


char* x_pixrect_from_clipboard(Window win)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;
  OpenClipboard(hwnd);

  //UINT fmt = EnumClipboardFormats(0);
  //while (fmt != 0 && fmt != CF_BITMAP) EnumClipboardFormats(fmt);
  //if (fmt == 0) return 0;

  HBITMAP hbm = (HBITMAP)GetClipboardData(CF_BITMAP);

  if (hbm == NULL) 
  { CloseClipboard();
    return 0;
   }

  //BITMAP bm;
  //GetObject(hbm,sizeof(BITMAP),(LPSTR)&bm);
  //hbm = CreateBitmapIndirect(&bm);

  SIZE sz;
  GetBitmapDimensionEx(hbm,&sz);
  x_pixrect* im = new x_pixrect;
  im->map = hbm;
  im->w = sz.cx;
  im->h = sz.cy;

  CloseClipboard();
  ReleaseDC(hwnd,hdc);
  return (char*)im;
 }


char* x_create_pixrect(Window win, int w, int h, unsigned char* data, int fc, 
                                                                      int bc)
{ return x_create_bitmap(win,w,h,data); }

char* x_create_pixrect(Window win, int x1, int y1, int x2, int y2)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;
  HDC hdcMem  = CreateCompatibleDC (hdc) ;
  x_pixrect* im = new x_pixrect;
  if (x1 > x2) SWAP(x1,x2);
  if (y1 > y2) SWAP(y1,y2);
  im->w = x2-x1+1;
  im->h = y2-y1+1;
  im->map = CreateCompatibleBitmap (hdc,im->w,im->h) ;
  SelectObject (hdcMem,im->map);
  BitBlt (hdcMem,0,0,im->w,im->h,hdc,x1,y1,SRCCOPY) ;
  DeleteDC (hdcMem) ;
  ReleaseDC(hwnd,hdc);
  return (char*)im;
}


void x_insert_pixrect(Window win, int x, int y, char* rect)
{ // (x,y) lower left corner !
  win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;
  HDC hdcMem  = CreateCompatibleDC (hdc) ;
  x_pixrect* im = (x_pixrect*)rect;
  SelectObject (hdcMem,im->map);
  BitBlt (hdc,x,y-im->h+1,im->w,im->h,hdcMem,0,0,SRCCOPY) ;
  DeleteDC (hdcMem) ;
  ReleaseDC(hwnd,hdc);
}


void x_insert_pixrect(Window win, int x, int y, char* rect, int x0, int y0, 
                                                         int w, int h)
{ // (x,y) lower left corner !
  win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;
  HDC hdcMem  = CreateCompatibleDC (hdc) ;
  x_pixrect* im = (x_pixrect*)rect;
  SelectObject (hdcMem,im->map);
  BitBlt (hdc,x,y-h+1,w,h,hdcMem,x0,y0,SRCCOPY) ;
  DeleteDC (hdcMem) ;
  ReleaseDC(hwnd,hdc);
}



void x_insert_pixrect(Window win, char* rect)
{ // insert at (0,0)
  win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;
  HDC hdcMem  = CreateCompatibleDC (hdc) ;
  x_pixrect* im = (x_pixrect*)rect;
  SelectObject (hdcMem,im->map);
  BitBlt (hdc,0,0,im->w,im->h,hdcMem,0,0,SRCCOPY) ;
  DeleteDC (hdcMem) ;
  ReleaseDC(hwnd,hdc);
}



void x_insert_bitmap(Window win, int x, int y, char* rect)
{ // (x,y) lower left corner !
  win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;
  HDC hdcMem  = CreateCompatibleDC (hdc) ;
  x_pixrect* im = (x_pixrect*)rect;
  SelectObject (hdcMem,im->map);

  COLORREF col = wp->pen_data.lopnColor;
  HBRUSH  save_brush = (HBRUSH)SelectObject(hdc,CreateSolidBrush(col));

  if (wp->MODE == xor_mode)
     BitBlt (hdc,x,y-im->h+1,im->w,im->h,hdcMem,0,0,0x00A60706);
  else
     BitBlt (hdc,x,y-im->h+1,im->w,im->h,hdcMem,0,0,0x00E20746);

  DeleteDC (hdcMem);

  DeleteObject(SelectObject(hdc,save_brush));
  ReleaseDC(hwnd,hdc);
}



void x_delete_bitmap(char* rect)
{ if (rect)
  { x_pixrect* im = (x_pixrect*)rect;
    DeleteObject(im->map);
    delete im;
   }
 }

void x_delete_pixrect(char* rect) { x_delete_bitmap(rect); }



void x_copy_pixrect(Window win, int x1, int y1, int x2, int y2, int x, int y)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;
  HDC hdcMem  = CreateCompatibleDC (hdc) ;
  int w = x2-x1+1;
  int h = y2-y1+1;
  HBITMAP map = CreateCompatibleBitmap (hdc,w,h) ;
  SelectObject (hdcMem,map);
  BitBlt (hdcMem,0,0,w,h,hdc,x1,y1,SRCCOPY) ;
  BitBlt (hdc,x,y,w,h,hdcMem,0,0,SRCCOPY) ;
  DeleteDC (hdcMem) ;
  ReleaseDC(hwnd,hdc);
}




char* x_create_pixrect(Window win, const char** xpm) 
{ 
  HDC hdc = CreateDC("DISPLAY",NULL,NULL,NULL); 

  int width;
  int height;
  int colors;
  int chars;

  CGAL_CLIB_STD::sscanf(*xpm,"%d %d %d %d",&width,&height,&colors, &chars);

  if (chars > 2)
  { std::cerr << "xpm: sorry, chars_per_pixel > 2 not implemented.\n";
    exit(1);
   }

  x_pixrect* im = new x_pixrect;

  im->w = width;
  im->h = height;
  im->map = CreateCompatibleBitmap(hdc,im->w,im->h) ;

  HDC hdcMem = CreateCompatibleDC(hdc) ;
  SelectObject(hdcMem,im->map);

  COLORREF color_table[1<<16]; 
  int i = 1<<16; 
  while (i--) color_table[i] = rgb_table[white];


  for(i=0; i<colors; i++)
  { xpm++;
    char c1=0;
    char c2=0;
    char rgb_str[16];

    if (chars == 1) CGAL_CLIB_STD::sscanf(*xpm,"%c c %s",&c1,rgb_str);
    if (chars == 2) CGAL_CLIB_STD::sscanf(*xpm,"%c%c c %s",&c1,&c2,rgb_str);

    int len = strlen(rgb_str);

    if (strcmp(rgb_str,"None") == 0)
      if (win == 0)
         color_table[c1+256*c2] = RGB(0,0,0); // icon
      else
         color_table[c1+256*c2] = RGB(220,220,220); //grey1
    else
      { if (len == 13)
        { rgb_str[3] = rgb_str[5];
          rgb_str[4] = rgb_str[6];
          rgb_str[5] = rgb_str[9];
          rgb_str[6] = rgb_str[10];
          rgb_str[7] = 0;
         }
    
        int x;
        CGAL_CLIB_STD::sscanf(rgb_str+1,"%x",&x);
     
        int b = x % 256; x = x/256;
        int g = x % 256; x = x/256;
        int r = x % 256;
        color_table[c1+256*c2] = RGB(r,g,b);
       }
   }


  for(int y=0; y<height; y++)
  { const char* line = *++xpm;
    for(int x=0; x<width; x++)
    { int index = line[x*chars];
      if (chars == 2) index += 256*line[x*chars+1];
      SetPixel(hdcMem,x,y,color_table[index]); 
     }
   }

  DeleteDC(hdcMem);
  DeleteDC(hdc);

  return (char*)im;
}


void x_pixrect_dimensions(char* rect, int* w, int* h)
{ if (rect)
  { x_pixrect* im = (x_pixrect*)rect;
    *w = im->w;
    *h = im->h;
   }
}



// fonts and text

static HFONT x_create_font(const char* fname)
{ 
  LOGFONT lf;

  ZeroMemory(&lf,sizeof(LOGFONT));

  switch (fname[0]) {

    case 'T': lf.lfPitchAndFamily = FF_SWISS | VARIABLE_PITCH;
              lf.lfWeight = FW_NORMAL;
              break;
  
    case 'I': lf.lfPitchAndFamily = FF_SWISS | VARIABLE_PITCH;
              lf.lfWeight = FW_NORMAL;
              lf.lfItalic = 1;
              break;
  
    case 'B': lf.lfPitchAndFamily = FF_SWISS | VARIABLE_PITCH;
              lf.lfWeight = FW_BOLD;
              break;
  
    case 'F': lf.lfPitchAndFamily = FF_MODERN | FIXED_PITCH;
              lf.lfWeight = FW_NORMAL;
              break;

    case 'A':
              lf.lfPitchAndFamily = 82;
              lf.lfWeight = 400;
              lf.lfCharSet = 1;
              lf.lfOutPrecision = 3;
              lf.lfClipPrecision = 2;
              lf.lfQuality = 1;
              strcpy(lf.lfFaceName,"Empire BT");
              break;

    default : return NULL;
  }

  int h = 0;
  int i = 1;
  while (fname[i]) h = 10*h + fname[i++] - '0';

  if (i > 4 || h > 500) return NULL;

  lf.lfHeight = h+3;

  return CreateFontIndirect(&lf);
}



int x_set_font(Window win, const char* fname)
{ 
  win_struct* wp = wlist[win];

  if (strcmp(fname,wp->font_name) == 0) return 1;

  HFONT hf_new = x_create_font(fname);
  if (hf_new == NULL) return 0;

  strcpy(wp->font_name,fname);

  if (wp->tmp_font) DeleteObject(wp->tmp_font);

  wp->tmp_font = wp->font = hf_new;

  return 1;
}


void x_set_text_font(Window win)  
{ win_struct* wp = wlist[win];
  wp->font_name[0] = '\0'; 
  wp->font = text_font;
 }

void x_set_italic_font(Window win)
{ win_struct* wp = wlist[win];
  wp->font_name[0] = '\0'; 
  wp->font = italic_font;
 }

void x_set_bold_font(Window win)  
{ win_struct* wp = wlist[win];
  wp->font_name[0] = '\0'; 
  wp->font = bold_font;
 }

void x_set_fixed_font(Window win) 
{ win_struct* wp = wlist[win];
  wp->font_name[0] = '\0'; 
  wp->font = fixed_font;
 }


void x_set_button_font(Window win)
{ win_struct* wp = wlist[win];
  wp->font_name[0] = '\0'; 
  wp->font = button_font;
 }




static void x_text_dimensions(Window win,const char* s, int len, SIZE& sz) 
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC  hdc  = wp->hdc;

  if (hdc == 0)
  { // win still closed
    x_open_display();
    hdc = CreateDC("DISPLAY",NULL,NULL,NULL); 
   }


  HFONT save_font = (HFONT)SelectObject(hdc,wp->font);
  GetTextExtentPoint(hdc,s,len,&sz);
  SelectObject(hdc,save_font);

  if (hdc == wp->hdc)
    ReleaseDC(hwnd,hdc);
  else
    DeleteDC(hdc);
}


int x_text_width(Window win,const char* s, int len)
{ if (len > (int)strlen(s)) len = strlen(s);
  if (len == 0) return 0;
  SIZE sz;
  x_text_dimensions(win,s,len,sz);
  return sz.cx;
 }

int x_text_height(Window win,const char* s) 
{ SIZE sz;
  if (strlen(s) == 0)
    x_text_dimensions(win,"H",1,sz);
  else
    x_text_dimensions(win,s,strlen(s),sz);

  if (win == 0)
     return sz.cy - 2;
  else
     return sz.cy - 1;
 }

int x_text_width(Window win,const char* s)
{ return x_text_width(win,s,strlen(s)); }



/* drawing parameters */

int x_set_bg_color(Window win, int c) 
{ win_struct* wp = wlist[win];
  int save = wp->B_COLOR;
  wp->B_COLOR = c; 
  return save;
}

char* x_set_bg_pixmap(Window win, char* pm) 
{ win_struct* wp = wlist[win];
  char* save = wp->B_PIXMAP;
  wp->B_PIXMAP = pm; 
  return save; 
}


void x_set_stipple(Window win, char* bits, int c) 
{ win_struct* wp = wlist[win];

  if (wp->stip_bm) DeleteObject(wp->stip_bm);
  wp->stip_bm = NULL;

  if (bits)
  { BITMAP bm;
    bm.bmType = 0;
    bm.bmWidth = 16;
    bm.bmHeight = 16;
    bm.bmWidthBytes = 2;
    bm.bmPlanes = 1;
    bm.bmBitsPixel = 1;
    bm.bmBits = bits;
    wp->stip_bm = CreateBitmapIndirect(&bm);
    wp->stip_bgcol = c;
   }
}

int x_set_color(Window win, int col)
{ win_struct* wp = wlist[win];
  int old_col = wp->COLOR;

  //if (col == old_col) return col;

  wp->COLOR = col;

  HWND hwnd = wp->hwnd;

  if (hwnd) 
  { HDC hdc = wp->hdc;
    if (col > 12)
       wp->pen_data.lopnColor = GetNearestColor(hdc,rgb_table[col]);
    else
       wp->pen_data.lopnColor = rgb_table[col];

/*
   if (hdc == 0) 
      printf("win = %d   hwnd = %d  hdc = %d\n",win,hwnd,hdc);
*/

    DeleteObject(SelectObject(hdc,CreatePenIndirect(&wp->pen_data)));
    SetTextColor(hdc,wp->pen_data.lopnColor);
    ReleaseDC(hwnd,hdc);
   }
  return old_col;
 }


drawing_mode x_set_mode(Window win, drawing_mode mod) 
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  drawing_mode old_mod = wp->MODE;
  if (mod == old_mod) return old_mod;
  wp->MODE = mod;
  if (hwnd)
  { HDC hdc = wp->hdc;
    switch (mod) {
    case src_mode: SetROP2(hdc, R2_COPYPEN);
                   break;
    case and_mode: SetROP2(hdc, R2_MASKPEN);
                   break;
    case or_mode:  SetROP2(hdc, R2_MERGEPEN);
                   break;
    case xor_mode: SetROP2(hdc, R2_NOTXORPEN);
                   break;
    case diff_mode:
                   break;
    }
    ReleaseDC(hwnd,hdc);
   }

  return old_mod;
 }


int x_set_line_width(Window win, int lw) 
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  int old_lw = wp->LWIDTH;
  if (lw == old_lw) return old_lw;
  wp->LWIDTH = lw;
  if (hwnd)
  { HDC hdc = wp->hdc;
    wp->pen_data.lopnWidth.x = lw;
    wp->pen_data.lopnWidth.y = lw;
    DeleteObject(SelectObject(hdc,CreatePenIndirect(&wp->pen_data)));
    ReleaseDC(hwnd,hdc);
   }
  return old_lw;
}


line_style x_set_line_style(Window win, line_style ls)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  line_style old_ls = wp->LSTYLE;
  if (ls == old_ls) return old_ls;
  wp->LSTYLE = ls;
  if (hwnd)
  { HDC hdc = wp->hdc;
    switch (ls) {
    case solid:  wp->pen_data.lopnStyle = PS_SOLID;
                 break;
    case dashed: wp->pen_data.lopnStyle = PS_DASH;
                 break;
    case dotted: wp->pen_data.lopnStyle = PS_DOT;
                 break;
    case dashed_dotted: 
                 wp->pen_data.lopnStyle = PS_DASHDOT;
                 break;
    }
    DeleteObject(SelectObject(hdc,CreatePenIndirect(&wp->pen_data)));
    ReleaseDC(hwnd,hdc);
   }
  return old_ls;
 }


text_mode x_set_text_mode(Window win, text_mode tm) 
{ win_struct* wp = wlist[win];
  text_mode save = wp->TMODE;
  wp->TMODE = tm;  
  return save;
 }


int x_set_join_style(Window win, int js) 
{ win_struct* wp = wlist[win];
  int save = wp->JSTYLE;
  wp->JSTYLE = js;  
  return save;
 }



int x_get_color(Window win)      
{ win_struct* wp = wlist[win];
  return wp->COLOR;
}

drawing_mode x_get_mode(Window win)       
{ win_struct* wp = wlist[win];
  return wp->MODE; 
}

int x_get_line_width(Window win) 
{ win_struct* wp = wlist[win];
  return wp->LWIDTH;
}

line_style x_get_line_style(Window win) 
{ win_struct* wp = wlist[win];
  return wp->LSTYLE;
}

text_mode x_get_text_mode(Window win)  
{ win_struct* wp = wlist[win];
  return wp->TMODE;
}

int x_get_cursor(Window win)     
{ win_struct* wp = wlist[win];
  return wp->cursor;
}


int x_get_border_color(Window)
{ // not implemented
  return black; 
}

int x_get_border_width(Window) 
{ // not implemented
  return 1; 
}





void x_set_read_gc(Window win)
{ win_struct* wp = wlist[win];
  wp->save_mo = x_set_mode(win,xor_mode);
  wp->save_ls = x_set_line_style(win,solid);
  wp->save_lw = x_set_line_width(win,1);
  x_set_color(win,black);
 }


void x_reset_gc(Window win)
{ win_struct* wp = wlist[win];
  x_set_mode(win,wp->save_mo);
  x_set_line_style(win,wp->save_ls);
  x_set_line_width(win,wp->save_lw);
 }


/* colors */

int x_set_rgb(int i, int r, int g, int b)
{ x_open_display();
  rgb_table[i]  = RGB(r,g,b); 
  return 1;
}

void x_get_rgb(int i, int* red, int* green, int* blue)
{ x_open_display();

  if (i < 0 || i >= rgb_count)
  { *red = 0;
    *green = 0;
    *blue = 0;
    return;
   }

  *red   = GetRValue(rgb_table[i]);
  *green = GetGValue(rgb_table[i]);
  *blue  = GetBValue(rgb_table[i]);
 }


int x_new_color(int r, int g, int b) 
{ x_open_display();
  COLORREF col = RGB(r,g,b);
  int i = 0;
  while (i < rgb_count && rgb_table[i] != col) i++;
  if (i < rgb_count) return i;
  if (rgb_count < MAX_COLORS)
  { rgb_table[rgb_count] = col;
    return rgb_count++;
   }
  return -1;
}


int x_new_color(const char* name)
{ 
  if (name[0] != '#') return 0;

  char rgb_str[16];
  strcpy(rgb_str,name);

  if (strlen(rgb_str) == 13)
  { rgb_str[3] = rgb_str[5];
    rgb_str[4] = rgb_str[6];
    rgb_str[5] = rgb_str[9];
    rgb_str[6] = rgb_str[10];
    rgb_str[7] = 0;
  }
    
  int x;
  CGAL_CLIB_STD::sscanf(rgb_str+1,"%x",&x);
     
  int b = x % 256; x = x/256;
  int g = x % 256; x = x/256;
  int r = x % 256;

  return x_new_color(r,g,b);
}




/* events */

int x_get_next_event(Window* win, int* val, int* x, int* y,unsigned long *t)
{ 

  if (putback)
  { int w = last_event.win;
    if (wlist[w] == 0) putback = false;
   }


  if (putback) 
   { cur_event = last_event;
     putback = 0;
    }
  else
  { MSG msg;
    cur_event.kind = no_event;
    if (GetMessage (&msg,NULL,0,0))
      { TranslateMessage(&msg);
        DispatchMessage(&msg);
       }
    else 
       exit(0);

    if (cur_event.kind != exposure_event)
        cur_event.t = (unsigned long)msg.time;

    if (cur_event.kind != no_event) 
       last_event = cur_event;
   }

  *win = cur_event.win;
  *val = cur_event.val;
  *x   = cur_event.x;
  *y   = cur_event.y;
  *t   = cur_event.t;

  if ((wlist[*win] != 0))
    return cur_event.kind;
  else
    return no_event;
 }


int x_get_next_event(Window* win, int* val, int* x, int* y, unsigned long *t,
                                                            int msec)
{ // get next event with timeout

  MSG msg;
  if (putback || PeekMessage(&msg,NULL,0,0,PM_NOREMOVE))
  { int k = x_get_next_event(win,val,x,y,t); 
    if (k != no_event) return k;
   }

  HWND hwnd = wlist[*win]->hwnd;

  SetTimer(hwnd,2,msec,NULL);

  int k = x_get_next_event(win,val,x,y,t);

  KillTimer(hwnd,2);

  if (k == timer_event && *val == 2) 
   { *win = 0;
     k = no_event;
    }

  return k;
}


int x_check_next_event(Window* win, int* val, int* x, int* y, unsigned long *t)
{ MSG msg;
  if (putback || PeekMessage(&msg,NULL,0,0,PM_NOREMOVE))
     return x_get_next_event(win,val,x,y,t); 
  else
     return no_event;
 }

void x_put_back_event(void) { putback = 1; }

void x_set_border_width(Window win, int b) 
{ win_struct* wp = wlist[win];
  wp->border_w = b;
 }



//------------------------------------------------------------------------------
// not implemented
//------------------------------------------------------------------------------

void x_set_border_color(Window, int) 
{ /* not implemented */ }

void x_set_bg_origin(Window, int, int)
{ /* not implemented */ }

void x_fill_arc(Window, int, int, int, int, double, double)
{ /* not implemented */ }


void x_set_icon_pixmap(Window,char*)
{ /* not implemented */ }

void x_set_icon_label(Window, const char*) 
{ /* not implemented */ }

void x_set_icon_window(Window, Window)     
{ /* not implemented */ }

void x_iconify_window(Window) 
{ /* not implemented */ } 

void x_flush_display(void)
{ /* not implemented */ }

char* x_create_bitmap(Window win, int x1, int y1, int x2, int y2)
{ /* not implemented */
  return 0; 
 }


char* x_pixrect_to_bitmap(Window win, char* rect)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC  hdc = wp->hdc;

  x_pixrect* im = (x_pixrect*)rect;

  int w = im->w;
  int h = im->h;

  HDC hdcMem_src = CreateCompatibleDC(hdc);
  SelectObject(hdcMem_src,im->map);

  unsigned char* data = new unsigned char[h*(w/8+1)];
  char* bitmap = x_create_bitmap(win,w,h,data);

  HDC hdcMem_tgt = CreateCompatibleDC(hdc);
  SelectObject(hdcMem_tgt,((x_pixrect*)bitmap)->map);

  BitBlt(hdcMem_tgt,0,0,w,h,hdcMem_src,0,0,SRCCOPY) ;

  DeleteDC(hdcMem_src);
  DeleteDC(hdcMem_tgt);
  ReleaseDC(hwnd,hdc);

  return bitmap;
}


static int x_load_font(Window win, const char* fname, HFONT& hfont)   
{ win_struct* wp = wlist[win];
  HFONT hf = x_create_font(fname);
  if (hf == NULL) return 0;
  if (wp->font == hfont) wp->font = hf;
  DeleteObject(hfont);
  hfont = hf;
  return 1;
}

int x_load_text_font(Window win, const char* fn)   
{ return x_load_font(win,fn,text_font);  }

int x_load_italic_font(Window win, const char* fn) 
{ return x_load_font(win,fn,italic_font);}

int x_load_bold_font(Window win, const char* fn)   
{ return x_load_font(win,fn,bold_font);  }

int x_load_fixed_font(Window win, const char* fn)  
{ return x_load_font(win,fn,fixed_font); }

int x_load_button_font(Window win, const char* fn) 
{ return x_load_font(win,fn,button_font);}



void x_open_metafile(Window win, const char* fname)
{ win_struct* wp = wlist[win];
  x_close_metafile(win);
  wp->hdc2 = wp->hdc;
  if (wp->hmf) 
  { DeleteEnhMetaFile(wp->hmf);
    wp->hmf = NULL;
   }
  //RECT rect;
  //GetClientRect(wp->hwnd,&rect);
  //SetRect(&rect,0,0,100,100);
  //wp->hdc = CreateEnhMetaFile(NULL,fname,&rect,"LEDA\0window\0\0");
  wp->hdc = CreateEnhMetaFile(NULL,fname,NULL,"LEDA\0window\0\0");
}

void x_close_metafile(Window win)
{ win_struct* wp = wlist[win];
  if (wp->hdc2 == NULL) return;
  wp->hmf = CloseEnhMetaFile(wp->hdc);
  wp->hdc = wp->hdc2;
  wp->hdc2= NULL;
 }

void x_metafile_to_clipboard(Window win)
{ win_struct* wp = wlist[win];
  if (wp->hmf == NULL) return;
  OpenClipboard(wp->hwnd);
  EmptyClipboard();
  SetClipboardData(CF_ENHMETAFILE,wp->hmf);
  CloseClipboard();
 }

  

void x_load_metafile(Window win, int x0, int y0, int x1, int y1, 
                                                         const char* fname)
{ win_struct* wp = wlist[win];
  HWND hwnd = wp->hwnd;
  HDC hdc = wp->hdc;
  RECT rect;
  if (x0 > x1) SWAP(x0,x1);
  if (y0 > y1) SWAP(y0,y1);
  SetRect(&rect,x0,y0,x1,y1);
  HENHMETAFILE hmf = GetEnhMetaFile(fname);
  if (hmf)
  { PlayEnhMetaFile(hdc,hmf,&rect);
    DeleteEnhMetaFile(hmf);
   }
  ReleaseDC(hwnd,hdc);
}


void x_set_drop_handler(Window win, drop_func fun) 
{ wlist[win]->drop_handler = fun; }


int x_get_open_cmd(const char* suffix, char* buf, unsigned long buf_sz)
{ DWORD sz = buf_sz;
  HKEY hKey;
  RegOpenKeyEx(HKEY_CLASSES_ROOT, suffix, 0, KEY_QUERY_VALUE, &hKey);
  RegQueryValueEx(hKey,NULL, NULL,NULL, (unsigned char*)buf,&sz);
  RegCloseKey(hKey);
  RegOpenKeyEx(HKEY_CLASSES_ROOT, buf,    0, KEY_QUERY_VALUE, &hKey);
  RegOpenKeyEx(hKey,             "shell", 0, KEY_QUERY_VALUE, &hKey);
  RegOpenKeyEx(hKey,              "open", 0, KEY_QUERY_VALUE, &hKey);
  RegOpenKeyEx(hKey,           "command", 0, KEY_QUERY_VALUE, &hKey);
  sz = buf_sz;
  RegQueryValueEx(hKey,NULL, NULL,NULL, (unsigned char*)buf,&sz);
  RegCloseKey(hKey);
  return 1;
}


} // end namespace
