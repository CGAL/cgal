/**************************************************************************
 
  cc_extract_images.C
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Extracts inline filenames from a HTML stream given on stdin,
              following cc_ naming conventions: "cc_[^"]*.gif"
	      Implementierung als endlicher Automat
  System    : C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#include <stdlib.h>
#include <iostream.h>

#define State0   0    // Startzustand
#define State1   1    // '"' or '/'       erkannt
#define State2   2    // '"c'             erkannt
#define State3   3    // '"cc'            erkannt
#define State4   4    // '"cc_[^"]*'      erkannt
#define State5   5    // '"cc_[^"]*.'     erkannt
#define State6   6    // '"cc_[^"]*.g'    erkannt
#define State7   7    // '"cc_[^"]*.gi'   erkannt
#define State8   8    // '"cc_[^"]*.gif   erkannt, " starts output

const int max_buf = 32000;
char      buffer[ max_buf];
int       idx = 0;

int main( int argc, char**) {
    int c;
    int state = State0;

    if ( argc > 1) {
        cerr << "Usage: cc_extract_images" << endl;
	cerr << "Filter to extract inline image names from a HTML stream "
	        "from stdin" << endl;
        exit(1);
    }
    c = cin.get();
    while( cin) {
	if ( state > 0) {
	    if ( idx >= max_buf) {
		cerr << "cc_extract_images: error: internal buffer overflow."
		     << endl;
		exit(1);
	    }
	    buffer[idx++] = c;
	}
        switch ( state) {
	case State0:
	    if ( c == '"' || c == '/') {
		idx = 0;
		state = State1;
	    } else
		idx = 0;
	    break;
	case State1:
	    if ( c == 'c')
		state = State2;
	    else
		state = State0;
	    break;
	case State2:
	    if ( c == 'c')
		state = State3;
	    else
		state = State0;
	    break;
	case State3:
	    if ( c == '_')
		state = State4;
	    else
		state = State0;
	    break;
	case State4:
	    if ( c == '.')
		state = State5;
	    else if ( c == '"')
		state = State0;
	    break;
	case State5:
	    if ( c == 'g')
		state = State6;
	    else if ( c == '"')
		state = State0;
	    else
		state = State4;
	    break;
	case State6:
	    if ( c == 'i')
		state = State7;
	    else if ( c == '"')
		state = State0;
	    else 
		state = State4;
	    break;
	case State7:
	    if ( c == 'f')
		state = State8;
	    else if ( c == '"')
		state = State0;
	    else 
		state = State4;
	    break;
	case State8:
	    if ( c == '"') {
		buffer[idx-1] = '\0';  // Cancel last '"'
		cout << buffer << endl;
		state = State0;
	    } else 
		state = State4;
	    break;
	default:
	    cerr << "cc_extract_images: internal error: unknown state."
		 << endl;
	    exit(1);
        }
	c = cin.get();
    }
    return(0);
}
