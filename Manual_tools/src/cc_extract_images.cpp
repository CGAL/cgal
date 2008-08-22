/**************************************************************************
 
  cc_extract_images.cpp
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Extracts inline filenames from a HTML stream given on stdin,
              matches '["/][^"/<>]*.(gif|jpg|png)"'. The option '-cc' selects
              between names starting with 'cc_' prefix, and those without.
              Implemented as a finite state automaton.
  System    : C++ (g++)
  Author    : (c) 1998,2004 Lutz Kettner
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

#include <stdlib.h>
#include <iostream>
#include <string>

extern int yylex();

bool      cc_option = false; // if true, select only images with 'cc_' prefix


int main( int argc, char* argv[]) {
    if ( argc == 2 && argv[1] == std::string("-cc") )
        cc_option = true;
    if ( argc > 1 && ! cc_option) {
        std::cerr << "Usage: cc_extract_images\n";
        std::cerr << "Filter to extract inline image names from a HTML stream "
                "from stdin\n";
        std::cerr << "Options:\n";
        std::cerr << "    -cc    extracts the cc_ images used by latex_to_html"
            " converter" << std::endl;
        exit(1);
    }
    yylex();
    //extract_images();
    return 0;
}
