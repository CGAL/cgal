#include <CGAL/basic.h>

#if defined(CGAL_USE_CGAL_WINDOW)

#include <CGAL/LEDA/window.h>
#include <CGAL/LEDA/file.h>
#include <CGAL/LEDA/bitmaps/button32.h>

static unsigned char* bm_bits [] = {
  triang_bits,
  voro_bits,
  empty_circle_bits,
  encl_circle_bits,
  grid_bits,
  hull_bits,
};

using namespace std;


int main(int argc, char *argv[])
{ 
  if (argc >= 2) { return 0; }  

  bool   B     = false;
  CGAL::color  col   = CGAL::blue2;
  double R     = 3.1415;
  int    c     = 0;
  int    c1    = 0;
  int    c2    = 0;
  int    cm    = 5;
  int    N     = 100;
  string s0    = "string0";
  string s     = "dummy";
  string s1    = "menu";

  int but_num = 5;

  std::list<string> M1 = CGAL::get_files(".");
  std::list<string> M2 = CGAL::get_files(".");


  std::list<string> CML;
  CML.push_back("0");
  CML.push_back("1");
  CML.push_back("2");
  CML.push_back("3");
  CML.push_back("4");


  for(;;)
  {
    CGAL::panel P("PANEL DEMO");

    string text;
    text += " The panel section of a window is used for displaying text";
    text += " and for updating the values of variables. It consists";
    text += " of a list of panel items and a list of buttons.";
    text += " All operations adding panel items or buttons to the panel";
    text += " section of a window have to be called before";
    text += " the window is displayed for the first time.";
  
    P.text_item("");
    P.text_item("\\bf\\blue A Text Item");
    P.text_item("");
    P.text_item(text);
  
    P.text_item("");
    P.text_item("\\bf\\blue A Bool Item");
    P.bool_item("bool item",B);

    P.text_item("\\bf\\blue A Color Item");
    P.color_item("color item",col);

    P.text_item("");
    P.text_item("\\bf\\blue A Slider Item");
    P.int_item("slider item(1,20)",but_num,-1,20);

    P.text_item("");
    P.text_item("\\bf\\blue Simple Items");
    P.string_item("string item",s);
    P.int_item("int item",N);
    P.double_item("double item",R);

    P.text_item("");
    P.text_item("\\bf\\blue String Menu Items");
    P.string_item("string menu",s0,M1,8);

    P.text_item("");
    P.text_item("\\bf\\blue Choice Items");
    P.choice_item("simple choice",c,"one","two","three","four","five");
    P.choice_mult_item("multiple choice",cm,CML);
    P.int_item("integer choice",c1,0,80,20);
    P.choice_item("bitmap choice",c2,6,32,32,bm_bits);

    P.text_item("");
    P.text_item("\\bf\\blue Buttons");
  
    int i;
    for(i=0; i < but_num; i++) P.button(string("button"));
  
    P.display();

    if (P.open() == 0) break;
  }
 
  return 0;
}
#else
#include <iostream>

int main()
{
 std::cout << "CGAL::window is not available !\n";

 return 0;
}

#endif
  
