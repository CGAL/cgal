/////////////////////////////////////////////////////////////////////////
//
// this implements sorely missing Unix pwd command on
// Windows 9* and NT 4. Used for Windows-specific installation
//
// apart from the usual pwd functionality, it prints, if given,
// at most 2 parameters from the optional parameter list
//
/////////////////////////////////////////////////////////////////////////
//
// author        : Dima Pasechnik <dima@cs.uu.nl> December 1999
//
// NB. link it with kernel32.lib

#include <iostream>
#include <windows.h>
using namespace std;

int main(int argc, char** argv) {
  char buf[256];
  int len = GetCurrentDirectory(255, buf);
  if (len) {
    if (argc > 1) cout << argv[1] << " "; // print the parameter as it is
    if (argc > 2) cout << argv[2]; // print the parameter as it is
    cout << buf;
    if (argc > 3) cout << argv[3]; // print the parameter as it is
    cout << endl;
    return 0;
  }
  else {
    int err = GetLastError();
    cerr << "Error: " << err << endl;
    return err;
  }
}

