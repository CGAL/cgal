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
// release_date  : 
//
// file          : src/CGALWin/x11/_x_basic.C
// package       : cgal_window (1.0)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.0
// revision_date : 20 June 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================


//------------------------------------------------------------------------------
// basic graphics x_... functions for libW (declared in <LEDA/impl/x_basic.h>)
// implemented by X11 library functions
//
// S. N"aher (1994-99)
//------------------------------------------------------------------------------


#include <CGAL/LEDA/impl/x_basic.h>
#include <CGAL/LEDA/pixmaps/leda_icon.xpm>

#include <cstdio>
#include <cstdarg>
#include <string>
#include <cctype>
#include <ctime>
#include <cassert>

//#if defined(__KCC) || (defined(__sgi) && defined(_COMPILER_VERSION))
#include <cmath>
//#endif




#include <unistd.h>
#include <sys/stat.h>

#if defined(mips)
#include <bstring.h>
#endif


#include <sys/time.h>
#include <sys/times.h>


#include <X11/X.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/cursorfont.h>


namespace CGAL {


inline void SWAP(int& x1, int& x2)
{ int t = x1; x1 = x2; x2 = t; }


#define NUMBER_OF_COLORS (1<<12)


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int do_not_open_display = 0;

static Atom wm_protocols;
static Atom wm_delete_window;

static char dot_mask[2] = { 2,6 };
static char dash_mask[2] = { 6,4 };
static char dash_dot_mask[4] = { 6,2,2,2 };


#if defined(SMALL_FONTS)

static const char* text_font_name1  = "helvetica-12";
static const char* text_font_name2 = "-adobe-helvetica-medium-r-*-*-12-*-*-*-*-*-*-*";

static const char* bold_font_name1  = "helvetica-bold-12";
static const char* bold_font_name2  = "-adobe-helvetica-bold-r-*-*-12-*-*-*-*-*-*-*";

static const char* fixed_font_name1 = "lucidasanstypewriter-12";
static const char* fixed_font_name2 = "7x13";

static const char* italic_font_name1 = "-adobe-helvetica-medium-o-*-*-12-*-*-*-*-*-*-*";
static const char* italic_font_name2 = "7x12";

static const char* button_font_name1  = "lucidasans-12";
static const char* button_font_name2 = "-adobe-helvetica-medium-r-*-*-12-*-*-*-*-*-*-*";


#else


static const char* text_font_name1  = "helvetica-14";
static const char* text_font_name2 = "-adobe-helvetica-medium-r-*-*-14-*-*-*-*-*-*-*";

static const char* bold_font_name1  = "helvetica-bold-14";
static const char* bold_font_name2  = "-adobe-helvetica-bold-r-*-*-14-*-*-*-*-*-*-*";

static const char* fixed_font_name1 = "lucidasanstypewriter-12";
static const char* fixed_font_name2 = "7x14";

static const char* italic_font_name1 = "-adobe-helvetica-medium-o-*-*-14-*-*-*-*-*-*-*";
static const char* italic_font_name2 = "7x14";

static const char* button_font_name1 = "lucidasans-12";
static const char* button_font_name2 = "7x14";

#endif



static int trace_events = 0;

static Display* display = NULL;
static int      display_rel;
static int      screen;

static XEvent   event;

static Colormap color_map;

static XFontStruct* text_font;
static XFontStruct* italic_font;
static XFontStruct* bold_font;
static XFontStruct* button_font;
static XFontStruct* fixed_font;

static XFontStruct* tmp_font;
static char tmp_font_name[64];

static int Read_Only_Colors = 1;

static unsigned long color_pix[NUMBER_OF_COLORS];
static int           color_id[NUMBER_OF_COLORS];
static int           color_count= 0;

static int wm_frame_width =  7;
static int wm_title_width = 24;

struct  x11_win
{
  Drawable     win;
  Drawable     win_save;
  Drawable     win_save2;
  Drawable     buf;
  Pixmap       clip_mask;

  int          buf_w;
  int          buf_h;

  XGCValues    gc_val;
  GC           gc;

  XFontStruct* font;
  int          font_sz;

  int          LINEWIDTH;
  line_style   LINESTYLE;
  int          JOINSTYLE;
  drawing_mode MODE;
  text_mode    TEXTMODE;
  int          COLOR;
  int          B_COLOR;
  char*        B_PIXMAP;

  int          mapped;

  void* inf;

  void (*redraw)(void*,int,int,int,int);

  Window       icon_win;
  Window       frame_win;
  Cursor       cursor;
  int          cursor_id;

  int          BORDER_COLOR;
  int          BORDER_COLOR_SAVE;
};


struct x_pixrect {
 int w;
 int h;
 Pixmap P;
 Pixmap mask;
 char* buf;
 x_pixrect(int width, int height) : w(width), h(height) { mask = None; }
};



const int MAXWIN = 256;

static x11_win* wlist[MAXWIN];
static int wcount = 0;



//------------------------------------------------------------------------------
// auxiliary functions
//------------------------------------------------------------------------------

/*
static char* duplicate_string(const char* p)
{ char* q = new char[strlen(p)+1];
  strcpy(q,p);
  return q;
}

static int get_char_width(XFontStruct* fp)
{ XCharStruct char_struct;
  int asc, desc, dir;
  XQueryTextExtents(display,fp->fid,"s",1, &dir,&asc,&desc,&char_struct);
  return char_struct.width;
}
*/



//------------------------------------------------------------------------------
// display functions
//------------------------------------------------------------------------------

#if (__SUNPRO_CC < 0x500)
static
#else
extern "C" 
#endif
int x_error_handler(Display* disp, XErrorEvent* err)
{ char msg[80];
  XGetErrorText(disp,err->error_code,msg,80);
  fprintf(stderr, "\n X Error: %s\n\n",msg);
  abort();
  return 0;
}

static XFontStruct* x_load_font(const char *fname)
{ 
  if (strcmp(fname,tmp_font_name) == 0) 
    return tmp_font;

  strcpy(tmp_font_name,fname);

  char font_name[64];
  int  sz = 0;

  if (isdigit(fname[1]))
  { for(int i=1; fname[i]!=0; i++) 
       sz = 10*sz + fname[i] - '0';
   }

  switch (fname[0]) {

     case 'T': 
       sprintf(font_name,
              "-bitstream-charter-medium-r-normal--0-%d0-0-0-p-0-iso8859-1",sz);
       break;

     case 'I': 
       sprintf(font_name,
              "-bitstream-charter-medium-i-normal--0-%d0-0-0-p-0-iso8859-1",sz);
       break;

     case 'B':
       sprintf(font_name,
              "-bitstream-charter-bold-r-normal--0-%d0-0-0-p-0-iso8859-1",sz);
       break;

     case 'F':
       //sprintf(font_name,"lucidasanstypewriter-%d",sz);
       sprintf(font_name,
              "-bitstream-courier-medium-r-normal--0-%d0-0-0-m-0-iso8859-1",sz);
       break;

     default :
       strcpy(font_name,fname);
       break;
   }

  tmp_font = XLoadQueryFont(display,font_name);

  return tmp_font;
}


static XFontStruct* x_load_font(const char* name, const char* name2)
{ XFontStruct* ft;
  const char* name3 = "7x14";
  if ((ft = x_load_font(name)) == NULL)
  { if ((ft = x_load_font(name2)) == NULL)
     { fprintf(stderr,"Cannot load text font %s\n",name2);
       fprintf(stderr,"Trying to load %s instead.\n",name3);
       if ((ft = x_load_font(name3)) == NULL)
       { fprintf(stderr,"Cannot load font %s\n",name3);
         abort();
        }
     }
  }
  return ft;
}

void x_do_not_open_display(int x) { do_not_open_display = x; }

void x_open_display(void)
{ 
  if (display != NULL || do_not_open_display) return;

  if ((display = XOpenDisplay(0)) == NULL)	
  { fprintf(stderr, "Cannot open display.\n");
    exit(1);
   }

  display_rel  = XVendorRelease(display);

  if (display_rel < 1000) 
  { // XFree86 3.3.4 gives 334 instead of 3340
    display_rel = 10*display_rel;
   }
  

  screen    = DefaultScreen(display);
  color_map = DefaultColormap(display,screen);

  text_font   = x_load_font(text_font_name1,   text_font_name2);
  italic_font = x_load_font(italic_font_name1, italic_font_name2);
  bold_font   = x_load_font(bold_font_name1,   bold_font_name2);
  fixed_font  = x_load_font(fixed_font_name1,  fixed_font_name2);
  button_font = x_load_font(button_font_name1, button_font_name2);


/*
  if (DefaultDepth(display,screen) > 1)
  { unsigned long plane_masks;
    //while (NUMER_OF_COLORS > 0 && 
    //      !XAllocColorCells(display,color_map,False,&plane_masks,
    //                        0,color_pix, NUMER_OF_COLORS)) NUMBER_OF_COLORS--;

    if (!Read_Only_Colors)
    {
     if (!XAllocColorCells(display,color_map,False,&plane_masks,0,color_pix,64))
        Read_Only_Colors = 1;
     else
        NUMBER_OF_COLORS = 64;
     }
   }
*/


  // 17 predefined colors from /usr/lib/X11/rgb.txt

  color_count = 0;
  x_new_color("white");           // 0: white
  x_new_color("black");           // 1: black
  x_new_color("red");             // 2: red
  if (x_new_color("green2") == -1)// 3: green
  { color_count--;
    x_new_color("green");
   }
  if (x_new_color("blue3") == -1) // 4: blue
  { color_count--;
    x_new_color("blue");
   }
  x_new_color("yellow");          // 5: yellow
  x_new_color("purple");          // 6: violet
  x_new_color("darkorange");      // 7: orange
  x_new_color("cyan");            // 8: cyan
  x_new_color("sienna");          // 9: brown
  x_new_color("magenta");         //10: pink 
  x_new_color("#0cb3a0");         //11: green2
  x_new_color("cornflowerblue");  //12: blue2

  int failed = 0;

  //13: grey1  
  if (x_new_color("grey85") == -1)
  { color_count--;
    if (x_new_color("grey76") == -1) failed++;
   }

  //14: grey2  
  if (x_new_color("grey70") == -1)
  { color_count--;
    if (x_new_color("grey76") == -1) failed++;
   }

  //15: grey3  
  if (x_new_color("grey45") == -1)
  { color_count--;
    if (x_new_color("grey76") == -1) failed++;
   }


  if (failed) {
    fprintf(stderr,"\nCould not allocate all 16 basic colors.\
                      Probably another X-client is\
                      using too many colors.\n\n");
  }

  //16: ivory
  if (x_new_color("#ffffe4") == -1)
  { color_count--;
    x_new_color("LightYellow");
   }


  XSetErrorHandler(x_error_handler);

  wm_protocols     = XInternAtom(display,"WM_PROTOCOLS",False);
  wm_delete_window = XInternAtom(display,"WM_DELETE_WINDOW",False);


  // NULL window

  x11_win* wp = new x11_win;

  wp->win = RootWindow(display,screen);

  wp->font      = fixed_font;
  wp->LINEWIDTH = 1;
  wp->LINESTYLE = solid;
  wp->JOINSTYLE = 1;
  wp->MODE      = src_mode;
  wp->TEXTMODE  = transparent;
  wp->COLOR     = black;
  wp->redraw    = 0;
  wp->inf       = 0;
  wp->mapped    = 0;
  wp->cursor_id = -1;

  wp->BORDER_COLOR_SAVE = black;
  wp->BORDER_COLOR = black;

  wlist[0] = wp;
  wcount = 0;
}




char* x_root_pixrect(int x0, int y0, int x1, int y1)
{ 
  if (display == NULL) 
  { if (do_not_open_display) return 0;
    x_open_display();
   }

  Window win = XRootWindow(display,screen);

  XGCValues  gc_val;
  gc_val.background = color_pix[0];
  gc_val.foreground = color_pix[1];
  gc_val.function  = GXcopy; 
  GC gc = XCreateGC(display,win, GCBackground | GCForeground | GCFunction,
                    &gc_val);


  if (x0 > x1) SWAP(x0,x1);
  if (y0 > y1) SWAP(y0,y1);

  int w = x1-x0+1;
  int h = y1-y0+1;

  XImage* I = XGetImage(display,win,x0,y0,w,h,AllPlanes,ZPixmap);

  x_pixrect* im = new x_pixrect(w,h);

  im->P = XCreatePixmap(display,win,w,h,DefaultDepth(display,screen));

  XPutImage(display,im->P,gc,I,0,0,0,0,w,h);

  XDestroyImage(I);
  XFreeGC(display,gc);

  return (char*)im;

}



void x_close_display()
{ int i;
  if (display) 
  { /*
    XFreeFont(display,text_font);
    XFreeFont(display,italic_font);
    XFreeFont(display,bold_font);
    XFreeFont(display,fixed_font);
    XFreeFont(display,button_font);
    */
    XCloseDisplay(display); 
    for(i=0; i <= wcount; i++) 
        if (wlist[i]) delete wlist[i];
    display = 0; 
   } 
}


int x_display_width(void)
{ if (display == NULL)
  { if (do_not_open_display) return 1024;
    x_open_display();
   }
  return DisplayWidth(display,screen);  
 }

int x_display_height(void)
{ if (display == NULL)
  { if (do_not_open_display) return 768;
    x_open_display();
   }
  return DisplayHeight(display,screen); 
 }

int x_display_depth(void)
{ if (display == NULL)
  { if (do_not_open_display) return 8;
    x_open_display();
   }
  return DefaultDepth(display,screen); 
 }

int x_display_bits_saved(void)
{ return XDoesBackingStore(XScreenOfDisplay(display,screen)); }
 

void x_flush_display(void)
{ XFlush(display); }



//------------------------------------------------------------------------------
// colors
//------------------------------------------------------------------------------

int x_set_rgb(int i, int r, int g, int b)   
{ 
  if (r > 255) r = 255;
  if (g > 255) g = 255;
  if (b > 255) b = 255;

  int col_id = r;
  col_id = (col_id << 8) + g;
  col_id = (col_id << 8) + b;
  color_id[i] = col_id;

  if (display == NULL)
  { if (do_not_open_display) return 0;
    x_open_display();
   }

  if (Read_Only_Colors)
  { XColor xcolor;
    xcolor.red   = (unsigned short)r*256;
    xcolor.green = (unsigned short)g*256;
    xcolor.blue  = (unsigned short)b*256;
    XAllocColor(display,color_map,&xcolor);
    if (XAllocColor(display,color_map,&xcolor))
    { color_pix[i] = xcolor.pixel;
      return 1;
     }
    else return 0;
   }

  XColor color;
  color.pixel = color_pix[i];
  XQueryColor(display,color_map,&color);
  color.red   = (unsigned short)r*256;
  color.green = (unsigned short)g*256;
  color.blue  = (unsigned short)b*256;
  XStoreColor(display,color_map,&color);
  return 1;
}


void x_get_rgb(int i, int* r, int* g, int* b)   
{ 
  if (i < 0 || i >= color_count)
  { *r = 0;
    *g = 0;
    *b = 0;
    return;
   }

  int col_id = color_id[i];
  *b = col_id & 0xff; col_id >>= 8;
  *g = col_id & 0xff; col_id >>= 8;
  *r = col_id & 0xff;

/*
  XColor xcolor;
  xcolor.pixel = color_pix[i];
  XQueryColor(display,color_map,&xcolor);
  *r = xcolor.red/256;
  *g = xcolor.green/256;
  *b = xcolor.blue/256; 
*/
 }


int x_new_color(int r, int g, int b)
{ 
  if (display == NULL && !do_not_open_display) x_open_display();

  if (r < 0) r = 0;
  if (g < 0) g = 0;
  if (b < 0) b = 0;

  if (r > 255) r = 255;
  if (g > 255) g = 255;
  if (b > 255) b = 255;

  int col_id = r;

  col_id = (col_id << 8) + g;
  col_id = (col_id << 8) + b;

  for (int i = 0; i < color_count; i++)
     if (color_id[i] == col_id) return i;  // has been allocated before

  if (color_count >= NUMBER_OF_COLORS) 
  { fprintf(stderr,"Too many colors (only %d available)\n",NUMBER_OF_COLORS);
    abort();
   }

  int failed = 0;

  if (display)
  { if (DefaultDepth(display,screen) == 1) // monochrome display
      if (r == 255 && g == 255 &&  b == 255)
          color_pix[color_count] = WhitePixel(display,screen);
       else
          color_pix[color_count] = BlackPixel(display,screen);
    else
       { XColor xcolor;
         xcolor.red   = (unsigned short)r*256;
         xcolor.green = (unsigned short)g*256;
         xcolor.blue  = (unsigned short)b*256;
         if (Read_Only_Colors)
           { if (XAllocColor(display,color_map,&xcolor))
                color_pix[color_count] = xcolor.pixel;
             else
              { color_pix[color_count] = color_pix[(r+g+b>512) ? white : black];
                failed = 1;
               }
            }
         else
           { xcolor.pixel = color_pix[color_count];
             xcolor.flags = DoRed | DoGreen | DoBlue;
             XStoreColor(display,color_map,&xcolor);
            }
   
        }
  }
  
  color_id[color_count] = col_id;
  color_count++;

  return (failed) ? -1 : color_count-1;
}



int x_new_color(const char* name)
{ if (display == NULL)
  { if (do_not_open_display) return 0;
    x_open_display();
   }
  XColor xcolor;
  XParseColor(display, color_map, name, &xcolor);
  int r = xcolor.red/256;
  int g = xcolor.green/256;
  int b = xcolor.blue/256;
  return x_new_color(r,g,b);
}


static Pixmap xpm_bitmask(Window win, const char** xpm) 
{ 
  int width;
  int height;
  int colors;
  int chars;
  sscanf(*xpm,"%d %d %d %d",&width,&height,&colors, &chars);

  if (chars > 2)
  { fprintf(stderr,"xpm: sorry, chars_per_pixel > 2 not implemented.\n");
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

    if (chars == 1) sscanf(*xpm,"%c c %s",&c1,rgb_str);
    if (chars == 2) sscanf(*xpm,"%c%c c %s",&c1,&c2,rgb_str);

    if (strcmp(rgb_str,"None") == 0 || strcmp(rgb_str,"none") == 0)
    { black_c1 = c1;
      black_c2 = c2;
     }
   }

  if (black_c1 == 0 && black_c2 == 0) 
     return None;


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
      if (c1 == black_c1 && c2 == black_c2) bits[pos/8] |= (1 << (pos%8));
      pos++;
     }
   }
  for(i=0; i<bytes; i++) bits[i] = ~bits[i];

  Pixmap pm = XCreateBitmapFromData(display,win,bits,width,height);
  delete[] bits;
  return pm;
}


static Pixmap xpm_pixmap(Window win, const char** xpm, int& width, int& height, 
                                                       int& transp) 
{
  int failed = 0;
  transp = 0;

  int colors;
  int chars;
  sscanf(*xpm,"%d %d %d %d",&width,&height,&colors, &chars);

  //unsigned long color_table[1<<16];

  int cn = 1 << (8*chars);

  unsigned long* color_table = new unsigned long[cn];

  int i;
  for(i=0; i<cn; i++) color_table[i] = color_pix[0];

  for(i=0; i<colors; i++)
  { xpm++;
    char c1=0;
    char c2=0;
    char cstr[32];

    if (chars == 1) sscanf(*xpm,"%c c %s",&c1,cstr);
    if (chars == 2) sscanf(*xpm,"%c%c c %s",&c1,&c2,cstr);

    int col_index = c1 + 256*c2;

    if (strcmp(cstr,"None") == 0 || strcmp(cstr,"none") == 0)
      { transp = 1;
        color_table[col_index] = color_pix[white];
       }
    else
    { XColor xcolor;
      XParseColor(display, color_map, cstr, &xcolor);

      if (XAllocColor(display,color_map,&xcolor))
          color_table[col_index] = xcolor.pixel;
      else
        { if (xcolor.red+xcolor.green+xcolor.blue > 384*256) 
            color_table[col_index] = color_pix[white];
          else
            color_table[col_index] = color_pix[black];
          failed ++;
         }
     }
   }


  if (failed > 0)
  { fprintf(stderr,"\n");
    fprintf(stderr,"xpm: could not allocate %d of %d colors .\n",failed,colors);
   }

  XGCValues  gc_val;
  gc_val.background = color_pix[0];
  gc_val.foreground = color_pix[1];
  gc_val.function  = GXcopy; 
  GC gc = XCreateGC(display,RootWindow(display,screen),
                    GCBackground | GCForeground | GCFunction,
                    &gc_val);


  Pixmap pm = XCreatePixmap(display,win,
                            width, height, DefaultDepth(display,screen));
 
  //XFillRectangle(display,pm,gc,0,0,width,height);

  XImage* im = XGetImage(display,pm,0,0,width,height,AllPlanes,ZPixmap);


  for(int y=0; y<height; y++)
  { const char* line = *++xpm;
    for(int x=0; x<width; x++)
    { int index = line[x*chars];
      if (chars == 2) index += 256*line[x*chars+1];
      XPutPixel(im,x,y, color_table[index]);
     }
   }

  XPutImage(display,pm,gc,im,0,0,0,0,width,height);

  XDestroyImage(im);
  XFreeGC(display,gc);

  delete[] color_table;

  return pm;
}





static Window xpm_icon_window(const char** xpm) 
{ 
  Window win = XCreateSimpleWindow(display,RootWindow(display,screen),
                                           0, 0, 1, 1, 0,
                                           BlackPixel(display,screen),
                                           BlackPixel(display,screen));
  int width, height;
  int transp = 0;
  Pixmap pm = xpm_pixmap(win,xpm, width, height, transp);
  XResizeWindow(display,win,width,height);
  XSetWindowBackgroundPixmap(display,win,pm);
  XFreePixmap(display,pm);
  return win;
}


 
  

//------------------------------------------------------------------------------
// windows
//------------------------------------------------------------------------------

int x_create_window(void* inf, int width, int height, int bg_col,
                    const char* header, const char* label, int parent,
                    void (*func)(void*,int,int,int,int))
{
  if (display == NULL) return 0;

  x11_win* wp = new x11_win;

  XSetWindowAttributes attrib;
  attrib.override_redirect = True;
  attrib.backing_store = Always;

  if (width <= 0) width = 1;
  if (height <= 0) height = 1;

  wp->win = XCreateWindow(display, wlist[parent]->win, 0, 0,
                          width,height,1,DefaultDepth(display,screen),
                          InputOutput,DefaultVisual(display,screen), 
                          CWBackingStore, &attrib);


  XSelectInput(display,wp->win, EnterWindowMask | LeaveWindowMask    |
                                KeyPressMask    | KeyReleaseMask     |
                                ButtonPressMask | ButtonReleaseMask  |
                                PointerMotionMask  | 
                                ExposureMask    | StructureNotifyMask);
  
  if (DefaultDepth(display,screen) == 1) 
    XSetWindowBackground(display,wp->win, WhitePixel(display, screen)); 
  else
    XSetWindowBackground(display,wp->win,color_pix[bg_col]);
 
  XStoreName(display,wp->win,header);
  XSetIconName(display,wp->win,label);

  XWMHints wm_hints;

  wm_hints.flags = StateHint | InputHint;
  wm_hints.initial_state = NormalState;
  wm_hints.input = True;
  XSetWMProperties(display,wp->win,0,0,0,0,0,&wm_hints,0);

  XSetWMProtocols(display,wp->win,&wm_delete_window,1);
  

  wp->gc_val.background = color_pix[0];
  wp->gc_val.foreground = color_pix[1];
  wp->gc_val.function  = GXcopy; 
  wp->gc_val.line_style = LineSolid;
  wp->gc_val.cap_style = CapButt;
  //wp->gc_val.join_style = JoinMiter;
  wp->gc_val.join_style = JoinRound;
  wp->gc_val.line_width = 1;
  wp->gc_val.font = text_font->fid;
  wp->gc_val.graphics_exposures = False;

  wp->gc = XCreateGC(display,RootWindow(display,screen),
                    GCBackground | GCForeground |GCFunction | GCJoinStyle |
                    GCLineStyle  | GCLineWidth  |GCFont | GCGraphicsExposures,
                    &(wp->gc_val));

  wp->font = text_font;
  wp->LINEWIDTH = 1;
  wp->JOINSTYLE = 1;
  wp->LINESTYLE = solid;
  wp->MODE      = src_mode;
  wp->COLOR     = black;
  wp->B_COLOR   = white;
  wp->B_PIXMAP  = 0;
  wp->TEXTMODE  = transparent;
  wp->redraw    = func;
  wp->mapped    = 0;
  wp->inf       = inf;
  wp->buf       = 0;
  wp->buf_w     = 0;
  wp->buf_h     = 0;
  wp->win_save  = 0;
  wp->win_save2 = 0;
  wp->clip_mask = None;
  wp->icon_win  = 0;
  wp->frame_win  = 0;
  wp->cursor_id = -1;

  int i = 1;
  while (i <= wcount)
  { if (wlist[i] == 0) break;
    i++;
   }
  if (i > wcount) wcount = i;

  assert(wcount < MAXWIN);

  wlist[i] = wp;

  return i;
}


void x_resize_window(int w,int xpos,int ypos,int width,int height,int parent)
{ x11_win* wp = wlist[w];

  if (parent == 0)
  { xpos += wm_frame_width;
    ypos += wm_title_width;
   }

  if (wp->win_save != 0)
     XMoveResizeWindow(display,wp->win_save,xpos,ypos,width,height);
  else
     XMoveResizeWindow(display,wp->win,xpos,ypos,width,height);
  XFlush(display);
}


int x_window_opened(int w) 
{ return wlist[w] && wlist[w]->mapped; }


void x_open_window(int w, int xpos, int ypos, int width, int height, int pw)
{
  x11_win* wp = wlist[w];

  if (wp->mapped) return;

  if (width < 0)  
  { width = -width;
    xpos -= (width-1);
    if (pw == 0) xpos -= (2 * wm_frame_width);
   }

  if (height < 0)  
  { height = -height;
    ypos -= (height-1);
    if (pw == 0) ypos -= (wm_frame_width + wm_title_width);
   }


  if (pw == 0 && wp->icon_win == 0)
  { // assign icon window
    wp->icon_win = xpm_icon_window(leda_icon);
    XWMHints wm_hints;
    wm_hints.icon_window = wp->icon_win;
    wm_hints.flags = IconWindowHint;
    XSetWMProperties(display,wp->win,0,0,0,0,0,&wm_hints,0);
   }


  if (pw > 0)
  { x11_win* pwp = wlist[pw];
    if (pwp->win_save != 0)
       XReparentWindow(display,wp->win,pwp->win_save,xpos,ypos);
    else
       XReparentWindow(display,wp->win,pwp->win,xpos,ypos);
   }

  
  if (pw < 0)
  { XSetWindowAttributes attrib;
    attrib.override_redirect = True;
    XChangeWindowAttributes(display,wp->win, CWOverrideRedirect,&attrib);
    x_set_cursor(w,XC_arrow);
   }
  else
  { XSetWindowAttributes attrib;
    attrib.override_redirect = False;
    XChangeWindowAttributes(display,wp->win, CWOverrideRedirect,&attrib);
    XUndefineCursor(display,wp->win);
   }


  XMoveResizeWindow(display,wp->win,xpos,ypos,width,height);

  XSizeHints size_hints;
  size_hints.flags = PPosition;
  size_hints.x = xpos;
  size_hints.y = ypos;
  
  XSetWMProperties(display,wp->win,0,0,0,0,&size_hints,0,0);

  if (wp->B_PIXMAP)
  { Pixmap pm = ((x_pixrect*)wp->B_PIXMAP)->P;
    XSetWindowBackgroundPixmap(display,wp->win,pm);
   }
  else
    XSetWindowBackground(display,wp->win,color_pix[wp->B_COLOR]);

  XMapRaised(display,wp->win);

  wp->mapped = 1;

  XEvent e;
  do XMaskEvent(display,StructureNotifyMask,&e);
  while (e.type != MapNotify);

  while (XCheckWindowEvent(display, wp->win, ExposureMask, &e));

/*
  if (pw < 0)
    XSetInputFocus(display,wp->win,RevertToParent,CurrentTime);
*/


  Window wm_frame = wp->win;

//if (wp == 0)

 if (pw >= 0)
  { for(;;) 
    { Window root,parent;
      Window *childlist = 0;
      unsigned int u;
      Status status = XQueryTree(display,wm_frame,&root,&parent,&childlist,&u);
      if (childlist) XFree((char*)childlist);
      if (parent == root || !parent || !status) break;
      wm_frame = parent;
     }
   }
  
  wlist[w]->frame_win = wm_frame;

  if (pw == 0)
  { int wx0,wy0,wx1,wy1;
    x_window_frame(w,&wx0,&wy0,&wx1,&wy1);
    wx1 = wy1 = 0;
    x_window_to_screen(w,&wx1,&wy1);
    wm_frame_width = wx1 - wx0;
    wm_title_width = wy1 - wy0;
   }
}
  


void x_close_window(int w)
{ x11_win* wp = wlist[w];
  if (wp->mapped)
  { wp->mapped = 0;
    XUnmapWindow(display,wp->win);
  }
}

void x_iconify_window(int w)
{ x11_win* wp = wlist[w];
  if (wp->mapped) 
    XIconifyWindow(display,wp->win,screen);
}


void x_destroy_window(int w)
{ x11_win* wp = wlist[w];

  if (display)
  { x_stop_buffering(w);
    x_close_window(w);
    if (wp->icon_win) XDestroyWindow(display,wp->icon_win);
  
    // disconnect children
    Window r_win,p_win;
    Window* child=0;
    unsigned int n;
    XQueryTree(display,wp->win,&r_win,&p_win,&child,&n);
    for(unsigned int i=0; i < n; i++)
       XReparentWindow(display,child[i],RootWindow(display,screen),0,0);
    if (child) XFree((char*)child);
  
    XFreeGC(display,wp->gc);
    XDestroyWindow(display,wp->win); 
  }

  delete wlist[w];
  wlist[w] = 0;
  if (w == wcount) wcount--;
 }


static void ClearPixmap(int w, Pixmap pmap, int x0, int y0, int x1, 
                                                            int y1, 
                                                            int xorig,
                                                            int yorig)
{ x11_win* wp = wlist[w];

 if (wp->B_PIXMAP)
 { x_pixrect* im = (x_pixrect*)wp->B_PIXMAP;
   Pixmap pm = im->P;
   int wi = im->w;
   int he = im->h;

   if (xorig > 0)
      while (xorig > 0) xorig -= wi;
   else
      while (xorig+wi < 0) xorig += wi;

   if (yorig > 0)
      while (yorig > 0) yorig -= he;
   else
      while (yorig+he < 0) yorig += he;

   Window win1;
   int xpos,ypos;
   unsigned width,height,bw,dep;
   XGetGeometry(display,pmap,&win1,&xpos,&ypos,&width,&height,&bw,&dep);

   int xmax = width;
   int ymax = height;

   for(int y = yorig;  y < ymax; y += he)
     for(int x = xorig; x < xmax; x += wi)
      if (x < x1 && x+wi > x0 && y < y1 && y+he >y0)
        XCopyArea(display,pm,pmap,wp->gc,0,0,wi,he,x,y);
  }
 else
  { x_set_color(w,wp->B_COLOR);
    XSetClipMask(display,wp->gc,None);
    XFillRectangle(display,pmap,wp->gc,x0,y0,x1-x0+1,y1-y0+1);
   }
}

 

void x_clear_window(int w, int x0, int y0, int x1, int y1, int xorig, int yorig)
{ x11_win* wp = wlist[w];
  if (x0 > x1)  SWAP(x0,x1);
  if (y0 > y1)  SWAP(y0,y1);

  wp->gc_val.ts_x_origin = xorig;
  wp->gc_val.ts_y_origin = yorig;
  XChangeGC(display,wp->gc,GCTileStipXOrigin|GCTileStipYOrigin,&(wp->gc_val));

  if (wp->win_save == 0 && ((xorig == 0 && yorig == 0) || wp->B_PIXMAP == 0 ))
  { XClearArea(display,wp->win,x0,y0,x1-x0+1,y1-y0+1,False);
    //int c = x_set_color(w,grey1);
    //x_text(w,x1-100,y0,"LEDA Research 1998");
    //x_set_color(w,c);
   }
  else
    ClearPixmap(w,wp->win,x0,y0,x1,y1,xorig,yorig);
}

void x_clear_window(int w, int x0, int y0, int x1, int y1)
{ x_clear_window(w,x0,y0,x1,y1,0,0); }



void* x_window_inf(int w) { return wlist[w]->inf; }

int x_window_height(int w)
{ x11_win* wp = wlist[w];
  Window win1;
  int xpos,ypos;
  unsigned width,height,bw,dep;
  XGetGeometry(display,wp->win,&win1,&xpos,&ypos,&width,&height,&bw,&dep);
  return height;
 }

int x_window_width(int w)
{ x11_win* wp = wlist[w];
  Window win1;
  int xpos,ypos;
  unsigned width,height,bw,dep;
  XGetGeometry(display,wp->win,&win1,&xpos,&ypos,&width,&height,&bw,&dep);
  return width;
 }



void x_window_frame(int w, int* x0, int* y0, int* x1, int* y1)
{ x11_win* wp = wlist[w];
  Window win;
  int xpos,ypos;
  unsigned width,height,bw,dep;
  XGetGeometry(display,wp->frame_win,&win,&xpos,&ypos,&width,&height,&bw,&dep);

  *x0 = xpos;
  *y0 = ypos;
  *x1 = xpos + width - 1;
  *y1 = ypos + height - 1;
 }




//------------------------------------------------------------------------------
// drawing functions
//------------------------------------------------------------------------------

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



void x_pixel(int w, int x, int y)
{ x11_win* wp = wlist[w];
  XDrawPoint(display,wp->win,wp->gc,x,y); 
 }


int x_get_pixel(int, int, int) { return 0; }



void x_point(int w, int x, int y)
{ x11_win* wp = wlist[w];
  for(int i = -2; i <= 2; i++)
  { XDrawPoint(display,wp->win,wp->gc,x+i,y+i); 
    XDrawPoint(display,wp->win,wp->gc,x+i,y-i); 
   }
 }


void x_plus(int w, int x, int y)
{ x11_win* wp = wlist[w];
  for(int i = -3; i <= 3; i++)
  { XDrawPoint(display,wp->win,wp->gc,x+i,y); 
    XDrawPoint(display,wp->win,wp->gc,x,y+i); 
   }
 }



void x_pixels(int w, int n, int *x, int *y)
{ x11_win* wp = wlist[w];
  XPoint* points = new XPoint[n];
  int i;
  for(i=0; i<n; i++)
  { points[i].x = (short)x[i];
    points[i].y = (short)y[i];
   }
  XDrawPoints(display,wp->win,wp->gc,points,n,CoordModeOrigin);
  delete[] points;
 }


inline void bound_coord(int& x, int max_c)
{ if (x >  max_c) x =  max_c;
  else if (x < -max_c) x = -max_c;
}


void x_line(int w, int x1, int y1, int x2, int y2)
{ x11_win* wp = wlist[w];

  int jstyle = wp->JOINSTYLE;
  //int lwidth = (wp->LINEWIDTH+1)/2;
  int lwidth = 1;

  int max_c = (1 << 15) - 1;

  bound_coord(x1,max_c);
  bound_coord(x2,max_c);
  bound_coord(y1,max_c);
  bound_coord(y2,max_c);


// Server <= X11R4 ???

  if ((display_rel <  2100 && (x1 < x2 || (x1 == x2 && y1 < y2))) ||
      (display_rel >= 2100 && (x1 > x2 || (x1 == x2 && y1 > y2))))
  { SWAP(x1,x2);
    SWAP(y1,y2);
    if (jstyle == 1) 
       jstyle = 2; 
    else 
       if (jstyle == 2) jstyle = 1; 
   }

  if ((jstyle & 1) == 0) adjust_line(-lwidth,x2,y2,x1,y1);
  if ((jstyle & 2) == 1) adjust_line(+lwidth,x1,y1,x2,y2);

  XDrawLine(display,wp->win,wp->gc,x1,y1,x2,y2); 
 }


void x_lines(int w, int n, int *x1, int *y1, int* x2, int* y2)
{ 
  x11_win* wp = wlist[w];
  XSegment* segs = new XSegment[n];

  for(int i=0; i<n; i++)
  { segs[i].x1 = x1[i];
    segs[i].y1 = y1[i];
    segs[i].x2 = x2[i];
    segs[i].y2 = y2[i];
  }

  XDrawSegments(display,wp->win,wp->gc,segs,n);
  delete[] segs;
}
 


void x_rect(int w, int x1, int y1, int x2, int y2)
{ x11_win* wp = wlist[w];
  if (x1 > x2) SWAP(x1,x2);
  if (y1 > y2) SWAP(y1,y2);
  //XDrawRectangle(display,wp->win,wp->gc,x1,y1,x2-x1+1,y2-y1+1);
  XDrawRectangle(display,wp->win,wp->gc,x1,y1,x2-x1,y2-y1);
}
 
 
void x_box(int w, int x1, int y1, int x2, int y2)
{ x11_win* wp = wlist[w];
  if (x1 > x2) SWAP(x1,x2);
  if (y1 > y2) SWAP(y1,y2);
  XFillRectangle(display,wp->win,wp->gc,x1,y1,x2-x1+1,y2-y1+1);
}



#define setarcpix(win,x,y)\
{ if (orient >= 0)\
   { if (ax*y >= ay*x && x*by >= y*bx)\
        XDrawPoint(display,wp->win,wp->gc,x0+x,y0+y); }\
  else  if (ax*y >= ay*x || x*by >= y*bx)\
          XDrawPoint(display,wp->win,wp->gc,x0+x,y0+y);\
 }
 
using namespace std; 

void x_arc0(int w, int x0, int y0, int r0, int, double start,double angle)
{ x11_win* wp = wlist[w];

  if (angle > 0)
  { start += angle;
    angle = -angle;
   }

  double ax =  cos(start);
  double ay = -sin(start);
  double bx =  cos(start+angle);
  double by = -sin(start+angle);
  double orient = ax*by - ay*bx;

  for (int r = r0-wp->LINEWIDTH/2; r <= r0+wp->LINEWIDTH/2; r++)
  { int y = r;
    int x = 1;
    int e = 3 - 2*y;

    setarcpix(win, 0, r);
    setarcpix(win, 0,-r);
    setarcpix(win, r, 0);
    setarcpix(win,-r, 0);

    while (x < y)
    { setarcpix(win, x, y);
      setarcpix(win, x,-y);
      setarcpix(win,-x, y);
      setarcpix(win,-x,-y);
      setarcpix(win, y, x);
      setarcpix(win, y,-x);
      setarcpix(win,-y, x);
      setarcpix(win,-y,-x);
      x++;
      if (e>=0) { y--; e = e - 4*y; }
      e = e + 4*x + 2;
     }

    if (x == y)
    { setarcpix(win, x, y);
      setarcpix(win, x,-y);
      setarcpix(win,-x, y);
      setarcpix(win,-x,-y);
     }
  }
}


void x_arc(int w, int x0, int y0, int r1, int r2, double start, double angle)
{ double len = angle*r1;
  if (len > 2000 || len < -2000) 
  { x_arc0(w,x0,y0,r1,r2,start,angle);
    return;
   }
  x11_win* wp = wlist[w];
  int s = (int)(360*32*start/M_PI);
  int a = (int)(360*32*angle/M_PI);
  XDrawArc(display,wp->win,wp->gc,x0-r1,y0-r2,2*r1+1,2*r2+1,s,a);
}



void x_ellipse(int w, int x0, int y0, int r1, int r2)
{ x11_win* wp = wlist[w];
  XDrawArc(display,wp->win,wp->gc,x0-r1,y0-r2,2*r1+1,2*r2+1,0,360*64); }


void x_fill_arc(int w, int x0, int y0, int r1, int r2, double start, double angle)
{ x11_win* wp = wlist[w];
  int s = (int)(360*32*start/M_PI);
  int a = (int)(360*32*angle/M_PI);
  XFillArc(display,wp->win,wp->gc,x0-r1,y0-r2,2*r1+1,2*r2+1,s,a);
}


void x_circle(int w, int x0, int y0, int r)
{ x11_win* wp = wlist[w];
  XDrawArc(display,wp->win,wp->gc,x0-r,y0-r,2*r,2*r,0,360*64); }

void x_fill_circle(int w, int x0, int y0, int r)
{ x11_win* wp = wlist[w];
  //XFillArc(display,wp->win,wp->gc,x0-r,y0-r,2*r+1,2*r+1,0,360*64); 
  XFillArc(display,wp->win,wp->gc,x0-r,y0-r,2*r,2*r,0,360*64); 
}


/*
#define SETPIX(x,y)\
 XDrawPoint(display,wp->win,wp->gc,x,y); 

void x_circle(int w, int x0,int y0,int r0)
{ x11_win* wp = wlist[w];

  for (int r = r0-wp->LINEWIDTH/2; r <= r0+wp->LINEWIDTH/2; r++)
  { int y = r;
    int x = 0;
    int e = 3-2*y;

    XDrawPoint(display,wp->win,wp->gc,x,y); 

    SETPIX(x0,y0+r);
    SETPIX(x0,y0-r);
    SETPIX(x0+r,y0);
    SETPIX(x0-r,y0);

    for (x=1;x<y;)
      { SETPIX(x0+x,y0+y);
        SETPIX(x0+x,y0-y);
        SETPIX(x0-x,y0+y);
        SETPIX(x0-x,y0-y);
        SETPIX(x0+y,y0+x);
        SETPIX(x0+y,y0-x);
        SETPIX(x0-y,y0+x);
        SETPIX(x0-y,y0-x);
        x++;
        if (e>=0) { y--; e = e - 4*y; }
        e = e + 4*x + 2;
       }

    if (x == y)
    { SETPIX(x0+x,y0+y);
      SETPIX(x0+x,y0-y);
      SETPIX(x0-x,y0+y);
      SETPIX(x0-x,y0-y);
     }
  }
}



#define HLINE(x,y)\
XDrawLine(display,wp->win,wp->gc,x0-(x),y0+(y),x0+(x)+1,y0+(y))

void x_fill_circle(int w, int x0, int y0, int r)
{ x11_win* wp = wlist[w];

  int y = 1;
  int x = r;
  int e = 3-2*r;

  HLINE(x,0);

  while (y <= x)
  { HLINE(x,+y);
    HLINE(x,-y);
    if (y < x && e >= 0)
    { HLINE(y,+x);
      HLINE(y,-x);
      x--;
      e = e - 4*x;
     }
    y++;
    e = e + 4*y + 2;
   }
}
*/



void x_fill_ellipse(int w, int x0, int y0, int r1, int r2)
{ x11_win* wp = wlist[w];
  //if (wp->gc_val.fill_style == FillOpaqueStippled) { r1++; r2++; }
  XFillArc(display,wp->win,wp->gc,x0-r1,y0-r2,2*r1+1,2*r2+1,0,360*64); 
 }


void x_fill_polygon(int w, int n, int *xcoord, int *ycoord)
{ x11_win* wp = wlist[w];
  XPoint* edges = new XPoint[n];
  for(int i=0;i<n;i++) 
  { edges[i].x = (short) xcoord[i];
    edges[i].y = (short) ycoord[i];
   }

  XFillPolygon(display,wp->win,wp->gc,edges,n,Nonconvex,CoordModeOrigin);

  delete[] edges;
}


void x_polyline(int w, int n, int *xcoord, int *ycoord, int adjust)
{ x11_win* wp = wlist[w];
  int jstyle = wp->JOINSTYLE;

  if (xcoord[0] == xcoord[n-1] && ycoord[0] == ycoord[n-1]) jstyle = -1;

  if (adjust && (jstyle == 0 || jstyle == 3))
  { int d1,d2;
    if (jstyle == 0) 
    { d1 = -1;
      if (xcoord[0] > xcoord[1] || 
         (xcoord[0] == xcoord[1] && ycoord[0] > ycoord[1])) d1++;
      d2 = -1;
      if (xcoord[n-1] > xcoord[n-2] || 
         (xcoord[n-1] == xcoord[n-2] && ycoord[n-1] > ycoord[n-2])) d2++;
     }
  
    if (jstyle == 3) 
    { d1 = 0;
      if (xcoord[0] > xcoord[1] || 
         (xcoord[0] == xcoord[1] && ycoord[0] > ycoord[1])) d1++;
      d2 = 0;
      if (xcoord[n-1] > xcoord[n-2] || 
         (xcoord[n-1] == xcoord[n-2] && ycoord[n-1] > ycoord[n-2])) d2++;
    }
    adjust_line(d1,xcoord[1],ycoord[1],xcoord[0],ycoord[0]);
    adjust_line(d2,xcoord[n-2],ycoord[n-2],xcoord[n-1],ycoord[n-1]);
  }
  
  if (adjust == 1) return;

  XPoint* P = new XPoint[n];

  for(int i=0;i<n;i++) 
  { P[i].x = (short) xcoord[i];
    P[i].y = (short) ycoord[i];
   }

  XDrawLines(display,wp->win,wp->gc,P,n,CoordModeOrigin);
  delete[] P;
}




//------------------------------------------------------------------------------
// text
//------------------------------------------------------------------------------

void x_text_underline(int w, int x, int y, const char* s, int l,int r)
{ x11_win* wp = wlist[w];
  y += wp->font->ascent + 1;
  int x1 = x + x_text_width(w,s,l);
  int x2 = x + x_text_width(w,s,r+1);
  XDrawLine(display,wp->win,wp->gc,x1,y,x2,y);
}

void x_text(int w, int x, int y, const char* s, int l)
{ x11_win* wp = wlist[w];
  y += wp->font->ascent;
  if (unsigned(l) > strlen(s)) l = strlen(s);
  if (wp->TEXTMODE == transparent)
     XDrawString(display,wp->win,wp->gc,x,y,s,l);
  else
     XDrawImageString(display,wp->win,wp->gc,x,y,s,l);


  if (s[0] != '-') return;

  const char* ptr = s;
  while (*ptr == '-') ptr++;

  if (*ptr == '\0')
  { x += x_text_width(w,"-")/2;
    if (wp->TEXTMODE == transparent)
       XDrawString(display,wp->win,wp->gc,x,y,s,l-1);
    else
       XDrawImageString(display,wp->win,wp->gc,x,y,s,l-1);
   }
}

void x_text(int w, int x, int y, const char* s)
{ x_text(w,x,y,s,strlen(s)); }


void x_ctext(int w, int x, int y, const char* s)
{ int tw = x_text_width(w,s);
  int th = x_text_height(w,s);
  x_text(w, x-tw/2, y-th/2+th%2, s);
}

void x_ctext_underline(int w, int x, int y, const char* s, int l, int r)
{ int tw = x_text_width(w,s);
  int th = x_text_height(w,s);
  x_text_underline(w, x-tw/2, y-th/2+th%2, s, l, r);
}

 
int x_text_width(int w,const char* s)
{ return (display) ? XTextWidth(wlist[w]->font,s,strlen(s)) : 1; }

int x_text_width(int w,const char* s, int l)
{ return (display) ? XTextWidth(wlist[w]->font,s,l) : 1; }
 
 
int x_text_height(int w, const char*)
{ if (display == NULL) return 1;
  XFontStruct* fp = wlist[w]->font;
  return fp->ascent+(2*fp->descent)/3; 
}




//------------------------------------------------------------------------------
// pixrects
//------------------------------------------------------------------------------



char* x_create_pixrect(int w, const char** xpm) 
{ x11_win* wp = wlist[w];
  x_pixrect* im = new x_pixrect(0,0);
  int transp = 0;
  im->P = xpm_pixmap(wp->win,xpm,im->w,im->h,transp);
  if (transp)
     im->mask = xpm_bitmask(wp->win,xpm);
  else
     im->mask = None;
  return (char*)im;
}



char* x_create_pixrect(int w, int width, int height, unsigned char* data,
                                                     int fg_col, int bg_col) 
{ 
  x11_win* wp = wlist[w];
  x_pixrect* im=new x_pixrect(width,height);
  im->P=XCreatePixmapFromBitmapData(display,wp->win,(char*)data, width,height, 
                                    color_pix[fg_col], color_pix[bg_col],
                                    DefaultDepth(display,screen));
  return (char*)im;
 }



char* x_create_pixrect(int w, int x1, int y1, int x2, int y2)
{ x11_win* wp = wlist[w];
  drawing_mode save = x_set_mode(w,src_mode);
  if (x1 > x2) SWAP(x1,x2);
  if (y1 > y2) SWAP(y1,y2);

  x_pixrect* im = new x_pixrect(x2-x1+1,y2-y1+1);
  im->P = XCreatePixmap(display,wp->win,im->w,im->h,
                         DefaultDepth(display,screen));

  XCopyArea(display,wp->win,im->P,wp->gc,x1,y1,im->w,im->h,0,0);

  x_set_mode(w,save);
  return (char*)im;
 }

void x_pixrect_dimensions(char* pr, int* w, int* h)
{ if (pr)
  { x_pixrect* im = (x_pixrect*)pr;
    *w = im->w;
    *h = im->h;
   }
 }


void x_insert_pixrect(int w, int x, int y, char* prect)
{ // (x,y) lower left corner !

  assert(prect !=0);

  x11_win* wp = wlist[w];
  drawing_mode save = x_set_mode(w,src_mode);
  x_pixrect* p_im = (x_pixrect*)prect;

  y -= (p_im->h - 1);

  if (p_im->mask != None)
  { XSetClipOrigin(display,wp->gc,x,y);
    XSetClipMask(display,wp->gc,p_im->mask);
   }

  XCopyArea(display,p_im->P,wp->win,wp->gc,0,0,p_im->w,p_im->h,x,y);
  XFlush(display);
  x_set_mode(w,save);

  if (p_im->mask != None) XSetClipMask(display,wp->gc,None);
 }


void x_insert_pixrect(int w, int x, int y, char* prect, int x0, int y0, int wi, int he)
{ x11_win* wp = wlist[w];
  drawing_mode save = x_set_mode(w,src_mode);
  x_pixrect* im = (x_pixrect*)prect;
  XCopyArea(display,im->P,wp->win,wp->gc,x0,y0,wi,he,x,y-he+1);
  XFlush(display);
  x_set_mode(w,save);
 }


void x_insert_pixrect(int w, char* prect)
{ // inssert at (0,0)
  x11_win* wp = wlist[w];
  drawing_mode save = x_set_mode(w,src_mode);
  x_pixrect* im = (x_pixrect*)prect;
  XCopyArea(display,im->P,wp->win,wp->gc,0,0,im->w,im->h,0,0);
  XFlush(display);
  x_set_mode(w,save);
 }


void x_pixrect_to_matrix(int, char* prect, int* matrix)
{ 
  x_pixrect* im = (x_pixrect*)prect;
  XImage* I = XGetImage(display,im->P,0,0,im->w,im->h,AllPlanes,ZPixmap);

  int* p = matrix;
  for(int y=0; y<im->h; y++)
    for(int x=0; x<im->w; x++) 
    { unsigned long pix = XGetPixel(I,x,y);
      int i = 0;
      while (i < color_count && color_pix[i] != pix) i++;
      if (i == color_count && color_count < NUMBER_OF_COLORS)
      { color_pix[i] = pix;
        color_count++;
       }
      *p++ = i;
     }

  XDestroyImage(I);
}

void x_matrix_to_pixrect(int w, int* matrix, int width, int height, char** pr)
{ x11_win* wp = wlist[w];
  XGCValues  gc_val;
  gc_val.background = color_pix[0];
  gc_val.foreground = color_pix[1];
  gc_val.function  = GXcopy; 
  GC gc = XCreateGC(display,RootWindow(display,screen),
                    GCBackground | GCForeground | GCFunction,
                    &gc_val);

  x_pixrect* prect = new x_pixrect(width,height);

  prect->P = XCreatePixmap(display,wp->win,
                            width, height, DefaultDepth(display,screen));
 
  XImage* im = XGetImage(display,prect->P,0,0,width,height,AllPlanes,ZPixmap);

  for(int y=0; y<height; y++)
    for(int x=0; x<width; x++)
      XPutPixel(im,x,y,color_pix[*matrix++]);

  XPutImage(display,prect->P,gc,im,0,0,0,0,width,height);

  XDestroyImage(im);
  XFreeGC(display,gc);

  *pr = (char*)prect;
}





void x_pixrect_to_ps(int, char* prect, char* buf, int full_color, 
                                                    void (*progress)(int))
{ 

  x_pixrect* im = (x_pixrect*)prect;

  XImage* I = XGetImage(display,im->P,0,0,im->w,im->h,AllPlanes,ZPixmap);

  XColor* xcol = new XColor[im->w];
  XColor* xcstop = xcol + im->w;
  XColor* xcp;

  for(int y=0; y<im->h; y++)
  { if (progress) progress(y);
    int x = 0;
    for(xcp=xcol; xcp < xcstop; xcp++) xcp->pixel = XGetPixel(I,x++,y);
    XQueryColors(display,color_map,xcol,im->w);

    for(xcp=xcol; xcp < xcstop; xcp++)
    { int r = xcp->red/256;
      int g = xcp->green/256;
      int b = xcp->blue/256; 
      if (full_color)
       { sprintf(buf,"%02x%02x%02x",r,g,b);
         buf += 6;
        }
      else
       { sprintf(buf,"%02x",(r+g+b)/3);
         buf += 2;
        }
     }
   }
  XDestroyImage(I);
  delete[] xcol;
}



void x_delete_pixrect(char* prect)
{ x_pixrect* im = (x_pixrect*)prect;
  if (im)
  { XFreePixmap(display,im->P);
    if (im->mask != None) XFreePixmap(display,im->mask);
    delete im;
   }
 }


void x_copy_pixrect(int w, int x1, int y1, int x2, int y2, int x, int y)
{ char* im = x_create_pixrect(w,x1,y1,x2,y2); 
  x_insert_pixrect(w,x,y,im);
  x_delete_pixrect(im);
 }


char* x_create_bitmap(int w, int width, int height, unsigned char* bits)
{ x11_win* wp = wlist[w];
  x_pixrect* im = new x_pixrect(width,height);

  if (bits)
  { im->P = XCreateBitmapFromData(display,wp->win,(char*)bits,width,height);
    return (char*)im;
   }

  im->P = XCreatePixmap(display,wp->win,width,height,1);

  XGCValues  gc_val;
  gc_val.foreground = 0;
  gc_val.function  = GXcopy; 
  GC gc = XCreateGC(display,im->P, GCForeground | GCFunction, &gc_val);

  XFillRectangle(display,im->P,gc,0,0,width,height);

  XFreeGC(display,gc);
  return (char*)im;
}


char* x_pixrect_to_bitmap(int w, char* pr)
{ x11_win* wp = wlist[w];
  drawing_mode save = x_set_mode(w,src_mode);

  x_pixrect* im = (x_pixrect*)pr;
  int wi = im->w;
  int he = im->h;

  Pixmap bm = XCreatePixmap(display,wp->win,wi,he,1);

  XGCValues  gc_val;
  gc_val.background = 0;
  gc_val.foreground = 1;
  gc_val.function  = GXcopy; 

  GC gc0 = XCreateGC(display,bm,
                     GCBackground | GCForeground | GCFunction, &gc_val);

  XCopyPlane(display,im->P,bm,gc0,0,0,wi,he,0,0,1);

  gc_val.function = GXor;
  XChangeGC(display,gc0,GCFunction,&gc_val);
  XCopyPlane(display,im->P,bm,gc0,0,0,wi,he,0,0,2);
  XCopyPlane(display,im->P,bm,gc0,0,0,wi,he,0,0,4);
  XCopyPlane(display,im->P,bm,gc0,0,0,wi,he,0,0,8);

  XFreeGC(display,gc0);

  x_set_mode(w,save);

  x_pixrect* im1 = new x_pixrect(wi,he);
  im1->P = bm;


  return (char*)im1;
}



char* x_create_bitmap(int w, int x1, int y1, int x2, int y2)
{ x11_win* wp = wlist[w];
  drawing_mode save = x_set_mode(w,src_mode);

  if (x1 > x2) SWAP(x1,x2);
  if (y1 > y2) SWAP(y1,y2);

  x_pixrect* im = new x_pixrect(x2-x1+1,y2-y1+1);
  im->P = XCreatePixmap(display,wp->win,im->w,im->h,1);

  XGCValues  gc_val;
  gc_val.background = WhitePixel(display,screen);
  gc_val.foreground = BlackPixel(display,screen);
  gc_val.function  = GXcopy; 

  GC gc0 = XCreateGC(display,im->P,
                     GCBackground | GCForeground | GCFunction, &gc_val);

  XCopyPlane(display,wp->win,im->P,gc0,x1,y1,im->w,im->h,0,0,1);
  gc_val.function = GXor;
  XChangeGC(display,gc0,GCFunction,&gc_val);
  XCopyPlane(display,wp->win,im->P,gc0,x1,y1,im->w,im->h,0,0,2);
  XCopyPlane(display,wp->win,im->P,gc0,x1,y1,im->w,im->h,0,0,4);
  XCopyPlane(display,wp->win,im->P,gc0,x1,y1,im->w,im->h,0,0,8);

  XFreeGC(display,gc0);

  x_set_mode(w,save);

  return (char*)im;
}




void x_insert_bitmap(int w, int x, int y, char* bm)
{ x11_win* wp = wlist[w];
  x_pixrect* im = (x_pixrect*)bm;
  y -= (im->h - 1);
  XSetClipOrigin(display,wp->gc,x,y);
  XSetClipMask(display,wp->gc,im->P);
  x_box(w,x,y,x+im->w-1,y+im->h-1);
  XSetClipMask(display,wp->gc,None);
}


void x_delete_bitmap(char* bmap) 
{ x_pixrect* im = (x_pixrect*)bmap;
  if (im)
  { XFreePixmap(display,im->P);
    delete im;
   }
 }



//------------------------------------------------------------------------------
// fonts
//------------------------------------------------------------------------------


static void x_set_font(int w, XFontStruct* font)
{ if (display==NULL) return;
  x11_win* wp = wlist[w];
  if (wp->font != font)
  { wp->font = font;
    if (w != 0)
    { wp->gc_val.font = font->fid;
      XChangeGC(display,wp->gc,GCFont,&(wp->gc_val));
     }
   }
 }

static int x_load_font(int w, const char* font_name, XFontStruct*& xfont)
{ if (display == NULL) return 0;
  x11_win* wp = wlist[w];
  XFontStruct* fp = x_load_font(font_name);
  if (fp && fp != xfont)  
  { if (wp->font == xfont)  x_set_font(w,fp);
    XFreeFont(display,xfont);
    xfont = fp;
   }
  return (fp != NULL);
 }

int x_load_text_font(int w, const char* font_name)
{ return x_load_font(w,font_name,text_font); }

int x_load_italic_font(int w, const char* font_name)
{ return x_load_font(w,font_name,italic_font); }

int x_load_bold_font(int w, const char* font_name)
{ return x_load_font(w,font_name,bold_font); }

int x_load_fixed_font(int w, const char* font_name)
{ return x_load_font(w,font_name,fixed_font); }

int x_load_button_font(int w, const char* font_name)
{ return x_load_font(w,font_name,button_font); }


void x_set_text_font(int w)  { x_set_font(w,text_font);  }
void x_set_italic_font(int w){ x_set_font(w,italic_font);  }
void x_set_bold_font(int w)  { x_set_font(w,bold_font);  }
void x_set_fixed_font(int w) { x_set_font(w,fixed_font); }
void x_set_button_font(int w){ x_set_font(w,button_font); }


int x_set_font(int w, const char *fname)
{ if (display==NULL) return 0;
  XFontStruct* fp = x_load_font(fname);
  if (fp)  x_set_font(w,fp);
  return fp != NULL;
}



//------------------------------------------------------------------------------
// setting parameters
//------------------------------------------------------------------------------

int x_set_cursor(int w, int id)
{ x11_win* wp = wlist[w];
  Window win = (wp->win_save) ? wp->win_save : wp->win;
  if (wp->cursor_id == id) return id;
  if (wp->cursor_id >= 0) XFreeCursor(display,wp->cursor);
  if (id >= 0)
  { wp->cursor = XCreateFontCursor(display,id);
    XDefineCursor(display,win,wp->cursor);
   }
  else
    XUndefineCursor(display,win);

  XFlush(display);

  int old_id = wp->cursor_id;
  wp->cursor_id = id;
  return old_id;
}


void x_set_border_width(int w, int width)
{ x11_win* wp = wlist[w];
  XWindowChanges changes;
  changes.border_width = width;
  if (wp->win_save)
    XConfigureWindow(display,wp->win_save,CWBorderWidth,&changes);
  else
    XConfigureWindow(display,wp->win,CWBorderWidth,&changes);
}

void x_set_border_color(int w, int col) 
{ x11_win* wp = wlist[w];
  if (wp->win_save)
    { XSetWindowBorder(display,wp->win_save,color_pix[col]);
      wp->BORDER_COLOR_SAVE = col;
     }
  else
    { XSetWindowBorder(display,wp->win,color_pix[col]);
      wp->BORDER_COLOR = col;
     }
}
 
 
void x_set_label(int w, const char* label)
{ x11_win* wp = wlist[w];
  if (wp->win_save)
    XStoreName(display,wp->win_save,label); 
  else
    XStoreName(display,wp->win,label); 
  XFlush(display);
}

 
void x_set_icon_label(int w, const char* label)
{ x11_win* wp = wlist[w];
  if (wp->win_save)
    XSetIconName(display,wp->win_save,label);
  else
    XSetIconName(display,wp->win,label);
}




int x_set_bg_color(int w, int col)
{ x11_win* wp = wlist[w];
  int save = wp->B_COLOR;
  if (col != 0 && DefaultDepth(display,screen) == 1) col = 0; 
  Window win = (wp->win_save) ? wp->win_save : wp->win;
  XSetWindowBackground(display,win,color_pix[col]);
  wp->B_COLOR = col;
  return save;
 }


char* x_set_bg_pixmap(int w, char* prect)
{ x11_win* wp = wlist[w];
  char* save = wp->B_PIXMAP;
  wp->B_PIXMAP = prect;
  Window win = (wp->win_save) ? wp->win_save : wp->win;
  if (prect)
  { Pixmap pm = ((x_pixrect*)prect)->P;
    XSetWindowBackgroundPixmap(display,win,pm);
   }
  else
    XSetWindowBackground(display,win,color_pix[wp->B_COLOR]);
  return save;
 }

void x_set_bg_origin(int w, int xorig, int yorig)
{ x11_win* wp = wlist[w];
  wp->gc_val.ts_x_origin = xorig;
  wp->gc_val.ts_y_origin = yorig;
  XChangeGC(display,wp->gc,GCTileStipXOrigin|GCTileStipYOrigin,&(wp->gc_val));
}



int x_set_color(int w, int col)
{ x11_win* wp = wlist[w];
  int save = wp->COLOR;
  if (col != 0 && DefaultDepth(display,screen) == 1) col = 1; 
  wp->COLOR = col;
  wp->gc_val.foreground = color_pix[col];
  if (wp->MODE == xor_mode) wp->gc_val.foreground  ^= color_pix[white];
  XChangeGC(display,wp->gc,GCForeground,&(wp->gc_val));
  return save;
 }


void x_set_stipple(int w, char* bits, int c)
{ x11_win* wp = wlist[w];
  if (bits)
   { wp->gc_val.background = color_pix[c];
     Pixmap stip_pm = XCreateBitmapFromData(display,wp->win,bits,16,16);
     wp->gc_val.stipple = stip_pm;
     //wp->gc_val.fill_style = FillOpaqueStippled;
     wp->gc_val.fill_style = FillStippled;
     wp->gc_val.background = color_pix[c];
     XChangeGC(display,wp->gc,GCBackground|GCStipple|GCFillStyle,&(wp->gc_val));
    }
  else
   { wp->gc_val.background = color_pix[wp->B_COLOR];
     wp->gc_val.fill_style = FillSolid;
     XChangeGC(display,wp->gc,GCBackground|GCFillStyle,&(wp->gc_val));
    }
}



drawing_mode x_set_mode(int w, drawing_mode m)
{ x11_win* wp = wlist[w];

  if (wp->MODE == m) return m;

  drawing_mode save = wp->MODE;

  wp->MODE = m;

  wp->gc_val.foreground = color_pix[wp->COLOR];

  switch (m)  {

   case 0 : wp->gc_val.function = GXcopy;
            break;
   case 1 : wp->gc_val.function = GXxor;
            break;
   case 2 : wp->gc_val.function = GXor;
            wp->gc_val.foreground ^= color_pix[white];
            break;
   default: break;
  }

  XChangeGC(display,wp->gc,GCFunction,&(wp->gc_val));

  return save;
}


text_mode x_set_text_mode(int w, text_mode tm) 
{ x11_win* wp = wlist[w];
  text_mode save = wp->TEXTMODE;
  wp->TEXTMODE = tm;  
  return save;
 }

int x_set_join_style(int w, int js) 
{ x11_win* wp = wlist[w];
  int save = wp->JOINSTYLE;
  wp->JOINSTYLE = js;  
  return save;
 }



int x_set_line_width(int w, int lw)
{ x11_win* wp = wlist[w];
  if (wp->LINEWIDTH == lw) return lw;
  int save = wp->LINEWIDTH;
  wp->LINEWIDTH = lw;
  wp->gc_val.line_width = lw;
  XChangeGC(display,wp->gc,GCLineWidth,&(wp->gc_val));
  return save;
}


line_style x_set_line_style(int w, line_style s)
{ x11_win* wp = wlist[w];
  if (wp->LINESTYLE == s) return s;
  line_style save = wp->LINESTYLE;

  wp->LINESTYLE = s;

  switch (s)  {

   case solid  : wp->gc_val.line_style = LineSolid;
                 break;
   case dashed : wp->gc_val.line_style = LineOnOffDash;
                 XSetDashes(display,wp->gc,0,dash_mask,2);
                 break;
   case dotted : wp->gc_val.line_style = LineOnOffDash;
                 XSetDashes(display,wp->gc,0,dot_mask,2);
                 break;
   case dashed_dotted : 
                 wp->gc_val.line_style = LineOnOffDash;
                 XSetDashes(display,wp->gc,0,dash_dot_mask,4);
                 break;
   }

  XChangeGC(display,wp->gc,GCLineStyle,&(wp->gc_val));
  return save;
}


int           x_get_color(int w)      { return wlist[w]->COLOR;     }
drawing_mode  x_get_mode(int w)       { return wlist[w]->MODE;      }
int           x_get_line_width(int w) { return wlist[w]->LINEWIDTH; }
line_style    x_get_line_style(int w) { return wlist[w]->LINESTYLE; }
text_mode     x_get_text_mode(int w)  { return wlist[w]->TEXTMODE;  }
int           x_get_cursor(int w)     { return wlist[w]->cursor_id;      }



int x_get_border_width(int w)
{ x11_win* wp = wlist[w];
  XWindowAttributes attributes;
  if(wp->win_save) 
     XGetWindowAttributes(display,wp->win_save,&attributes);
  else
     XGetWindowAttributes(display,wp->win,&attributes);
  return attributes.border_width;
}

int x_get_border_color(int w)
{ x11_win* wp = wlist[w];
  if(wp->win_save)
     return wp->BORDER_COLOR_SAVE;
  else
     return wp->BORDER_COLOR;
}



//------------------------------------------------------------------------------
// event handling
//------------------------------------------------------------------------------

static int handle_event(int *w,int *val,int *x,int *y,unsigned long *t)
{
  KeySym keysym;
  XComposeStatus status;

  int  kind = no_event;

  Window win = event.xany.window;


  wlist[0]->win = win; //stopper

  int i = wcount;
  while (wlist[i] == 0 || 
         (wlist[i]->win != win && wlist[i]->win_save != win)) i--;

  *w = i;
  *val = 0;

  wlist[0]->win = RootWindow(display,screen);


  switch (event.type) {

  case ClientMessage:
                { Atom type_at = event.xclient.message_type;
                  Atom data_at = event.xclient.data.l[0];
                  if (type_at == wm_protocols && data_at == wm_delete_window) 
                  { //printf("delete window %d\n",*w);
                    kind = destroy_event;
                   }
                  break;
                 }

  case Expose:  //printf("expose: count = %d\n",event.xexpose.count);
                kind = exposure_event;
                *x = event.xexpose.x;
                *y = event.xexpose.y;
                *val = event.xexpose.width;
                *t = event.xexpose.height;
                break;


  case ConfigureNotify: kind = configure_event;
                        *x = event.xconfigure.x;
                        *y = event.xconfigure.y;
                        if (*x == 0 && *y == 0) x_window_to_screen(*w,x,y);
                        break;


  case DestroyNotify: //kind = destroy_event;
                      //printf("destroy notify %d\n",i);
                      break;


  case ButtonPress: *val = event.xbutton.button;
                    *x = event.xbutton.x;
                    *y = event.xbutton.y;
                    *t = event.xbutton.time;
                    kind = button_press_event;
                    if (event.xbutton.state & ShiftMask)   *val |= 256; 
                    if (event.xbutton.state & ControlMask) *val |= 512;
                    if (event.xkey.state & Mod1Mask)       *val |= 1024; // alt
                    if (event.xkey.state & Mod4Mask)       *val |= 1024; // sun
                    //XUngrabPointer(display,CurrentTime);
                    //x_ungrab_pointer();
                    break;

  case ButtonRelease: 
                    *val = event.xbutton.button;
                    *x = event.xbutton.x;
                    *y = event.xbutton.y;
                    *t = event.xbutton.time;
                    if (event.xbutton.state & ShiftMask)   *val |= 256; 
                    if (event.xbutton.state & ControlMask) *val |= 512;
                    if (event.xkey.state & Mod1Mask)       *val |= 1024; // alt
                    if (event.xkey.state & Mod4Mask)       *val |= 1024; // sun
                    kind = button_release_event;
                    break;

  case LeaveNotify:
  case EnterNotify:
  case MotionNotify: *x = event.xmotion.x;
                     *y = event.xmotion.y;
                     *t = event.xbutton.time;
                     kind = motion_event;
                     break;

  case KeyRelease: 
  case KeyPress: { *x = event.xkey.x;
                   *y = event.xkey.y;
                   *t = event.xkey.time;

                   char c = 0;
                   XLookupString(&event.xkey,&c,1, &keysym, &status);

                   if (c == 0) c = char(-1);
  
                   switch (keysym) {
  
                     case XK_R2:        // sun keyboard
                     case XK_Print:     c = KEY_PRINT;
                                        break;
                     case XK_R3:        c = KEY_PRINT1;
                                        break;
                     case XK_BackSpace: c = KEY_BACKSPACE;
                                        break;
                     case XK_Return:    c = KEY_RETURN;
                                        break;
                     case XK_Escape:    c = KEY_ESCAPE;
                                        break;
                     case XK_Left:      c = KEY_LEFT;
                                        break;
                     case XK_Right:     c = KEY_RIGHT;
                                        break;
                     case XK_Up:        c = KEY_UP;
                                        break;
                     case XK_Down:      c = KEY_DOWN;
                                        break;
                     case XK_Home:      c = KEY_HOME;
                                        break;
                     case XK_End:       c = KEY_END;
                                        break;
                     case XK_F1:        c = KEY_F1;
                                        break;
                     case XK_F2:        c = KEY_F2;
                                        break;
                     case XK_F3:        c = KEY_F3;
                                        break;
                     case XK_F4:        c = KEY_F4;
                                        break;
                     case XK_F5:        c = KEY_F5;
                                        break;
                     case XK_F6:        c = KEY_F6;
                                        break;
                     case XK_F7:        c = KEY_F7;
                                        break;
                     case XK_F8:        c = KEY_F8;
                                        break;
                     case XK_F9:        c = KEY_F9;
                                        break;
  
                     case XK_F10: if (event.type != KeyPress) break;
                                  if (!trace_events)
                                    { printf("START TRACING EVENTS\n");
                                      printf("DISPLAY: ");
                                      printf("X11 RELEASE = %d ",display_rel);
                                      printf("DEPTH = %d\n\n", 
                                              DefaultDepth(display,screen)); 
                                      trace_events = 1;
                                     }
                                  else
                                    { printf("STOP TRACING EVENTS\n");
                                      trace_events = 0;
                                     }
                                  break;
  
  
                     default:  if (c == 3) kind = destroy_event;
                               break;
                          
                    }

                   if (c == char(-1) || kind == destroy_event) break;

                   *val = c;
  
                   if (event.type == KeyPress) 
                     kind = key_press_event;
                   else
                     kind = key_release_event;
  
                   if (event.xkey.state & ShiftMask)   *val |= 256; 
                   if (event.xkey.state & ControlMask) *val |= 512;
                   if (event.xkey.state & Mod1Mask)    *val |= 1024; // alt
                   if (event.xkey.state & Mod4Mask)    *val |= 1024; // sun
                   break;
                 }
  
  case MappingNotify:
                 XRefreshKeyboardMapping((XMappingEvent*)&event);
                 break;

  }

  //if (trace_events) {
  // printf("%22s: w = %d  v = %d  x = %3d  y = %3d  t = %lu\n", 
  //                   event_name[kind],*w,*val,*x,*y,*t);
  //}

  return kind;
}


int x_check_next_event(int *w, int *val, int *x, int *y, unsigned long *t) 
{ // non-blocking 

  if (XCheckMaskEvent(display, 
                      EnterWindowMask | LeaveWindowMask    |
                      KeyPressMask    | KeyReleaseMask     |
                      ButtonPressMask | ButtonReleaseMask  |
                      PointerMotionMask  | 
                      ExposureMask    | StructureNotifyMask, &event) == 0) 
  { *w = 0;
    return no_event; 
   }

  return handle_event(w,val,x,y,t);
}


int x_get_next_event(int *w, int *val, int *x, int *y, unsigned long *t)
{ // blocking

  XNextEvent(display, &event);
  return handle_event(w,val,x,y,t);
}


// dummy select
static int select(size_t,void*,void*,void*, const timeval*) { return 0; }


int x_get_next_event(int *w, int *val, int *x, int *y, unsigned long *t, int msec)
{
  int fd = ConnectionNumber(display);

  timeval polltime;
  fd_set rdset, wrset, xset;

  if (XPending(display) == 0) 
  { 
    polltime.tv_sec  = msec / 1000; 
    polltime.tv_usec = 1000 * (msec % 1000);

    FD_ZERO(&rdset);
    FD_SET(fd,&rdset);
    FD_ZERO(&wrset);
    FD_ZERO(&xset);
    FD_SET(fd,&xset);

    int sel = select(fd+1,&rdset,&wrset,&xset,&polltime)
            + select(fd+1,(int*)&rdset,(int*)&wrset,(int*)&xset,&polltime);

    if (sel <= 0)
    { *w = 0;
      return no_event;
     }

  }

  XNextEvent(display, &event);

  int e = handle_event(w,val,x,y,t);
  if (e == no_event && *w == 0) *w = 1;
  return e;
}





void x_put_back_event(void)
{ XPutBackEvent(display,&event); }


int x_create_buffer(int w, int wi, int he)
{ x11_win* wp = wlist[w];
  if (wp->buf == 0 || wi != wp->buf_w || he != wp->buf_h)
  { if (wp->buf) XFreePixmap(display,wp->buf);
    wp->buf_w = wi;
    wp->buf_h = he;
    wp->buf = XCreatePixmap(display,wp->win,wi,he,DefaultDepth(display,screen));
    ClearPixmap(w,wp->buf,0,0,wi,he,0,0);
    if (wp->win_save) wp->win = wp->buf;
   }
  return 1;
}


int x_create_buffer(int w)
{ int wi = x_window_width(w);
  int he = x_window_height(w);
  x_create_buffer(w,wi,he);
  return 1;
}

int x_start_buffering(int w, int wi, int he)
{ x11_win* wp = wlist[w];
  x_stop_buffering(w);
  x_delete_buffer(w);
  x_create_buffer(w,wi,he);
  wp->win_save = wp->win;
  wp->win = wp->buf;
  return 1;
}



int x_start_buffering(int w) 
{ x11_win* wp = wlist[w];
  if (wp->win_save) return 0;
  int new_buf = wp->buf == 0;
  if (new_buf) x_create_buffer(w);
  wp->win_save = wp->win;
  wp->win = wp->buf;
  return new_buf;
} 


void x_set_buffer(int w, char* pr) 
{ x11_win* wp = wlist[w];

  if (pr == 0) 
  { if (wp->win_save2)
    { wp->win = wp->win_save2;
      wp->win_save2 = 0;
     }
   }
  else
  { if (wp->win_save2)
    { fprintf(stderr,"x_set_buffer: nested call.");
      abort();
     }
    x_pixrect* im = (x_pixrect*)pr;
    wp->win_save2 = wp->win;
    wp->win =  im->P;
   }
} 


int x_test_buffer(int w)
{ x11_win* wp = wlist[w];
  return wp->win_save != 0;
 } 

void x_flush_buffer(int w, int x1, int y1, int x2, int y2, int xoff, int yoff) 
{ x11_win* wp = wlist[w];
  if (wp->buf == 0) return;

  Window win = wp->win_save;
  if (win == 0) win = wp->win;

  drawing_mode save = x_set_mode(w,src_mode);
  if (x1 > x2) SWAP(x1,x2);
  if (y1 > y2) SWAP(y1,y2);
  int wi = x2-x1+1;
  int he = y2-y1+1;
  XCopyArea(display,wp->buf,win,wp->gc,x1+xoff,y1+yoff,wi,he,x1,y1);
  XFlush(display);
  x_set_mode(w,save);
}

void x_flush_buffer(int w, int x1, int y1, int x2, int y2) 
{ x_flush_buffer(w,x1,y1,x2,y2,0,0); }


void x_stop_buffering(int w)
{ x11_win* wp = wlist[w];
  if (wp->win_save) wp->win  = wp->win_save;
  wp->win_save = 0;
  x_set_bg_pixmap(w,wp->B_PIXMAP);
} 

void x_stop_buffering(int w, char** pr)
{ x11_win* wp = wlist[w];
  if (!wp->win_save) 
  { // not in buffering mode
    *pr = 0;
    return;
   }
  int wi = x_window_width(w);
  int he = x_window_height(w);
  x_pixrect* im = new x_pixrect(wi,he);
  im->P = wp->buf;  
  x_stop_buffering(w);
  wp->buf = 0;
  *pr = (char*)im;
} 

void x_delete_buffer(int w)  
{ x11_win* wp = wlist[w];
  x_stop_buffering(w);
  if (wp->buf) XFreePixmap(display,wp->buf);
  wp->buf = 0;
}



void x_start_timer(int, int ) {}
void x_stop_timer(int) {}



//------------------------------------------------------------------------------
// other functions
//------------------------------------------------------------------------------


void x_set_read_gc(int w)
{ XGCValues gc_val;
  gc_val.function = GXxor; 
  gc_val.foreground = BlackPixel(display,screen); 
  gc_val.line_style = LineSolid;
  gc_val.line_width = 1;
  XChangeGC(display,wlist[w]->gc,
            GCForeground|GCFunction|GCLineStyle|GCLineWidth,&gc_val);
  x_flush_display();
}

void x_reset_gc(int w)
{ x11_win* wp = wlist[w];
  XChangeGC(display,wp->gc,
            GCForeground|GCFunction|GCLineStyle|GCLineWidth,&(wp->gc_val));
  x_flush_display();
}



void x_grab_pointer(int w)
{ x11_win* wp = wlist[w];

  XGrabPointer(display,wp->win,False,ButtonPressMask|ButtonReleaseMask|
                                                     PointerMotionMask, 
               GrabModeAsync,GrabModeAsync,None,None,CurrentTime); 

//XGrabKeyboard(display,wp->win,False,GrabModeAsync,GrabModeAsync,CurrentTime); 
//XSetInputFocus(display,wp->win,RevertToParent,CurrentTime);
 }


void x_ungrab_pointer()
{ XUngrabPointer(display,CurrentTime);
//  XUngrabKeyboard(display,CurrentTime); 
}

void x_set_focus(int w)
{ x11_win* wp = wlist[w];
  if (wp->win_save == 0)
    XSetInputFocus(display,wp->win,RevertToParent,CurrentTime);
  else
    XSetInputFocus(display,wp->win_save,RevertToParent,CurrentTime);
 }


void x_set_icon_pixmap(int w, char* prect)
{ x11_win* wp = wlist[w];
  Window win = wp->icon_win;
  if (win == 0)
  { win = XCreateSimpleWindow(display,RootWindow(display,screen),
                                           0, 0, 1, 1, 0,
                                           BlackPixel(display,screen),
                                           BlackPixel(display,screen));
    wp->icon_win = win;
    XWMHints wm_hints;
    wm_hints.icon_window = win;
    wm_hints.flags = IconWindowHint;
    XSetWMProperties(display,wp->win,0,0,0,0,0,&wm_hints,0);
   }
  x_pixrect* im = (x_pixrect*)prect;
  XSetWindowBackgroundPixmap(display,win,im->P);
  XResizeWindow(display,win,im->w,im->h);
}


void x_set_icon_window(int w, int icon_w)
{ x11_win* wp = wlist[w];
  if (wp->icon_win) XDestroyWindow(display,wp->icon_win);
  wp->icon_win = wlist[icon_w]->win;
  XWMHints wm_hints;
  wm_hints.icon_window = wp->icon_win;
  wm_hints.flags = IconWindowHint;
  XSetWMProperties(display,wp->win,0,0,0,0,0,&wm_hints,0);
  wlist[icon_w]->mapped = 1;
}




void x_window_to_screen(int w, int* x, int* y)
{ x11_win* wp = wlist[w];
  Window src_win  = wp->win_save ? wp->win_save : wp->win;
  Window dest_win = RootWindow(display,screen);
  Window child_win;
  int x1,y1;
  XTranslateCoordinates(display,src_win,dest_win,*x,*y,&x1,&y1,&child_win);
  *x = x1;
  *y = y1;
}


void x_screen_to_window(int w, int* x, int* y)
{ x11_win* wp = wlist[w];
  Window src_win  = RootWindow(display,screen);
  Window dest_win = wp->win_save ? wp->win_save : wp->win;
  Window child_win;
  int x1,y1;
  XTranslateCoordinates(display,src_win,dest_win,*x,*y,&x1,&y1,&child_win);
  *x = x1;
  *y = y1;
}


void x_set_clip_rectangle(int w, int x0, int y0, int width, int height)
{ x11_win* wp = wlist[w];
  XRectangle rect;
  rect.x = 0;
  rect.y = 0;
  rect.width = width+1;
  rect.height = height+1;
  XSetClipMask(display,wp->gc,None);
  XSetClipRectangles(display,wp->gc,x0,y0,&rect,1,0);
}


void x_set_clip_ellipse(int win, int x, int y, int r1, int r2)
{ x11_win* wp = wlist[win];

  int width  = 2*r1+1;
  int height = 2*r2+1;

  Pixmap pm = XCreatePixmap(display,wp->win,width,height,1);

  XGCValues  gc_val;
  gc_val.background = 0;
  gc_val.foreground = 0;
  gc_val.function  = GXcopy; 
  GC gc = XCreateGC(display,pm,
                     GCBackground | GCForeground | GCFunction, &gc_val);

  XFillRectangle(display,pm,gc,0,0,width,height);

  gc_val.foreground = 1;
  XChangeGC(display,gc,GCForeground,&gc_val);

  XFillArc(display,pm,gc,0,0,width,height,0,360*64);

  XSetClipOrigin(display,wp->gc,x-r1,y-r2);
  XSetClipMask(display,wp->gc,pm);
  XFreePixmap(display,pm);
  XFreeGC(display,gc);
}



void x_set_clip_polygon(int w, int n, int *xcoord, int *ycoord)
{ x11_win* wp = wlist[w];
  XPoint* edges = new XPoint[n];
  int x0 = xcoord[0];
  int x1 = xcoord[0];
  int y0 = ycoord[0];
  int y1 = ycoord[0];

  int i;
  for(i=0;i<n;i++) 
  { int x = xcoord[i];
    int y = ycoord[i];
    if (x < x0) x0 = x;
    if (y < y0) y0 = y;
    if (x > x1) x1 = x;
    if (y > y1) y1 = y;
    edges[i].x = x;
    edges[i].y = y;
   }

  for(i=0;i<n;i++) 
  { edges[i].x -= x0;
    edges[i].y -= y0;
   }

  int width  = x1-x0 + 1;
  int height = y1-y0 + 1;

  Pixmap pm = XCreatePixmap(display,wp->win,width,height,1);


  XGCValues  gc_val;
  gc_val.background = 0;//BlackPixel(display,screen);
  gc_val.foreground = 0;//BlackPixel(display,screen);
  gc_val.function  = GXcopy; 
  GC gc = XCreateGC(display,pm,
                     GCBackground | GCForeground | GCFunction, &gc_val);


  XFillRectangle(display,pm,gc,0,0,width,height);

  gc_val.foreground = 1;//WhitePixel(display,screen);
  XChangeGC(display,gc,GCForeground,&gc_val);

  XFillPolygon(display,pm,gc,edges,n,Nonconvex,CoordModeOrigin);
  XSetClipOrigin(display,wp->gc,x0,y0);
  XSetClipMask(display,wp->gc,pm);
  XFreePixmap(display,pm);
  XFreeGC(display,gc);
  delete[] edges;
}



void x_clip_mask_polygon(int w, int n, int* xcoord, int* ycoord, int mode)
{ x11_win* wp = wlist[w];

  // mode = 0:  add to clip mask
  // mode = 1:  subtract

  if (n == 0)
  { if (wp->clip_mask)
    { XFreePixmap(display,wp->clip_mask);
      wp->clip_mask = None;
     }
    XSetClipMask(display,wp->gc,None);
    return;
   }
      

  XPoint* edges = new XPoint[n];

  if (wp->clip_mask == None)
  { int wi = x_window_width(w);
    int he = x_window_height(w);
    Pixmap pm = XCreatePixmap(display,wp->win,wi,he,1);
    wp->clip_mask = pm;
    XGCValues  gc_val;
    gc_val.foreground = 1;
    gc_val.function  = GXcopy; 
    GC gc = XCreateGC(display,pm, GCForeground | GCFunction, &gc_val);
    XFillRectangle(display,pm,gc,0,0,wi,he);
    XFreeGC(display,gc);
  }

  int i;
  for(i=0;i<n;i++) 
  { edges[i].x = xcoord[i];
    edges[i].y = ycoord[i];
   }

  Pixmap pm = wp->clip_mask;


  XGCValues  gc_val;
  gc_val.foreground = mode;
  gc_val.function  = GXcopy; 
  GC gc = XCreateGC(display,pm,GCForeground | GCFunction, &gc_val);
  XFillPolygon(display,pm,gc,edges,n,Nonconvex,CoordModeOrigin);
  XFreeGC(display,gc);

  XSetClipOrigin(display,wp->gc,0,0);
  XSetClipMask(display,wp->gc,pm);
  delete[] edges;
}


void x_clip_mask_ellipse(int w, int x, int y, int r1, int r2, int mode)
{ x11_win* wp = wlist[w];

  // mode = 0:  add to clip mask
  // mode = 1:  subtract

  if (wp->clip_mask == None)
  { int wi = x_window_width(w);
    int he = x_window_height(w);
    Pixmap pm = XCreatePixmap(display,wp->win,wi,he,1);
    wp->clip_mask = pm;
    XGCValues  gc_val;
    gc_val.foreground = 1;
    gc_val.function  = GXcopy; 
    GC gc = XCreateGC(display,pm, GCForeground | GCFunction, &gc_val);
    XFillRectangle(display,pm,gc,0,0,wi,he);
    XFreeGC(display,gc);
  }

  Pixmap pm = wp->clip_mask;

  XGCValues  gc_val;
  gc_val.foreground = mode;
  gc_val.function  = GXcopy; 
  GC gc = XCreateGC(display,pm,GCForeground | GCFunction, &gc_val);
  XFillArc(display,pm,gc,x-r1,y-r2,2*r1+1,2*r2+1,0,360*64); 
  XFreeGC(display,gc);

  XSetClipOrigin(display,wp->gc,0,0);
  XSetClipMask(display,wp->gc,pm);
}




// moving the pointer

void x_move_pointer(int win, int x, int y)
{ x11_win* wp = wlist[win];
  Drawable w = wp->win;
  if (wp->win_save != 0) w = wp->win_save;
  XEvent e;
  while (XCheckWindowEvent(display, wp->win, PointerMotionMask, &e));
  XWarpPointer(display,None,w,0,0,0,0,x,y);
 }



// not implemented

void  x_pixrect_to_clipboard(int, char*) {}
char* x_pixrect_from_clipboard(int)      { return 0; }

void x_open_metafile(int,const char*) {}
void x_close_metafile(int)            {}
void x_load_metafile(int,int,int,int,int,const char*) {}
void x_metafile_to_clipboard(int)     {}

int x_choose_file(int,int,const char*,const char*,char*,char*)  { return 0; }
int x_choose_color(int,int) { return 0; }


void x_set_drop_handler(int, void (*)(void*,const char*,int,int)) {}


int x_get_open_cmd(const char* suffix, char* buf, unsigned long /*buf_sz*/)
{ 

  if (strcmp(suffix,".htm") && strcmp(suffix,".html")) return 0;

  strcpy(buf,"netscape %1");

  char* home = getenv("HOME");
  if (home)
  { char lock_file[256];
    sprintf(lock_file,"%s/.netscape/lock",home);
    struct stat stat_buf;
    if (lstat(lock_file,&stat_buf) == 0)
        strcpy(buf,"netscape -remote 'openURL(%1)'");
   }
  return 1;
}



} // end namespace
