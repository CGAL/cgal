/**************************************************************************
 
  cc_extract_images.C
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Extracts inline filenames from a HTML stream given on stdin,
              matches '["/][^"/<>]*.(gif|jpg|png)"'. The option '-cc' selects
              between names starting with 'cc_' prefix, and those without.
              Implemented as a finite state automaton.
  System    : C++ (g++)
  Author    : (c) 1998,2004 Lutz Kettner
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#include <stdlib.h>
#include <iostream>

#define State0   0    // start state
#define State1   1    // '"'
#define State2   2    // '"c'
#define State3   3    // '"cc'
#define State4   4    // '"cc_[^"/<>]*'
#define State5   5    // '"[^"/<>]*.'
#define State6   6    // '"[^"/<>]*.[gjp]'
#define State7   7    // '"[^"/<>]*.[gjp][ipn]'
#define State8   8    // '"[^"/<>]*.[gjp][ipn][fg]'  accepted, " starts output

const int max_buf = 32000;
char      buffer[ max_buf];
int       idx = 0;
bool      cc_option = false; // if true, select only images with 'cc_' prefix

void extract_images() {
    int  state = State0;
    int  c = std::cin.get();
    bool cc_name = false;
    while( std::cin) {
        if ( state > 0) {
            if ( idx >= max_buf) {
                std::cerr << "cc_extract_images: error: internal buffer "
                             "overflow. Reset buffer." << std::endl;
                idx = 0;
                state = State0;
            }
            buffer[idx++] = c;
        }
        switch ( state) {
        case State0:
            idx = 0;
            cc_name = false;
            if ( c == '"')
                state = State1;
            break;
        case State1:
            idx = 0;
            buffer[idx++] = c;
            if ( c == 'c')
                state = State2;
            else if ( c == '"' )
                state = State1;
            else if ( c == '<' || c == '>')
                state = State0;
            else
                state = State4;
            break;
        case State2:
            if ( c == 'c')
                state = State3;
            else if ( c == '"' )
                state = State1;
            else if ( c == '<' || c == '>')
                state = State0;
            else
                state = State4;
            break;
        case State3:
            if ( c == '_') {
                cc_name = true;
                state = State4;
            } else if ( c == '"' )
                state = State1;
            else if ( c == '<' || c == '>')
                state = State0;
            else
                state = State4;
            break;
        case State4:
            if ( c == '.')
                state = State5;
            else if ( c == '"' )
                state = State1;
            else if ( c == '<' || c == '>')
                state = State0;
            break;
        case State5:
            if ( c == 'g' || c == 'j' || c == 'p')
                state = State6;
            else if ( c == '"' )
                state = State1;
            else if ( c == '<' || c == '>')
                state = State0;
            else
                state = State4;
            break;
        case State6:
            if ( c == 'i' || c == 'p' || c == 'n')
                state = State7;
            else if ( c == '"' )
                state = State1;
            else if ( c == '<' || c == '>')
                state = State0;
            else 
                state = State4;
            break;
        case State7:
            if ( c == 'f' | c == 'g')
                state = State8;
            else if ( c == '"' )
                state = State1;
            else if ( c == '<' || c == '>')
                state = State0;
            else 
                state = State4;
            break;
        case State8:
            if ( c == '"') {
                state = State0;
                if ( cc_option == cc_name) {
                    buffer[idx-1] = '\0';  // Cancel last '"'
                    // check for proper suffix
                    if (    ( 0 == strcmp( buffer + idx-4, "gif"))
                         || ( 0 == strcmp( buffer + idx-4, "png"))
                         || ( 0 == strcmp( buffer + idx-4, "jpg"))) {
                        if ( buffer[0] == '.' && buffer[1] == '/')
                            std::cout << (buffer+2) << std::endl;
                        else
                            std::cout << buffer << std::endl;
                    } else {
                        state = State1;
                    }
                }
            } else if ( c == '/')
                state = State1;
            else if ( c == '<' || c == '>')
                state = State0;
            else 
                state = State4;
            break;
        default:
            std::cerr << "cc_extract_images: internal error: unknown state. "
                         "Reste state." << std::endl;
            state = State0;
        }
        c = std::cin.get();
    }
}


int main( int argc, char* argv[]) {
    if ( argc == 2 && 0 == strcmp( argv[1], "-cc"))
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
    extract_images();
    return 0;
}
