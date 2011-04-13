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
// file          : src/CGALWin/_file.C
// package       : cgal_window (1.0)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 0.9.7
// revision_date : 23 May 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================



#include <CGAL/LEDA/basic.h>
#include <CGAL/LEDA/file.h>
#include <CGAL/LEDA/string_manip.h>


#include <cstdio>
#include <string>
#include <ctime>


#if defined(__unix__)
#include <unistd.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

#if defined(__HAS_FNMATCH_H__)
#include <fnmatch.h>
#else
inline int fnmatch(string,string,int) { return 0; }
#endif

#endif


#if defined(__win32__)

#define WIN32_EXTRA_LEAN
#define WIN32_LEAN_AND_MEAN
#include <windows.h>


#endif



#if defined(__unix__)

namespace CGAL {

using std::string;

bool is_directory(string fname)
{ struct stat stat_buf;
  if (stat(fname.c_str(),&stat_buf) != 0) return false;
  return (stat_buf.st_mode & S_IFMT) == S_IFDIR;
}

bool is_file(string fname)
{ struct stat stat_buf;
  if (stat(fname.c_str(),&stat_buf) != 0) return false;
  return (stat_buf.st_mode & S_IFMT) == S_IFREG;
}

bool is_link(string fname)
{ struct stat stat_buf;
  if (lstat(fname.c_str(),&stat_buf) != 0) return false;
  return (stat_buf.st_mode & S_IFMT) == S_IFLNK;
}


string get_directory()
{ char dir_buf[256];
  getcwd(dir_buf,256);
  return dir_buf;
}


string set_directory(string new_dir)
{ string old_dir = get_directory();
  if (is_directory(new_dir))
    chdir(new_dir.c_str());
  else
    std::cerr << "set_directory:" << new_dir << " is not a directory.";
  return old_dir;
}


static void read_directory(string dir_name, int what, std::list<string>& L, 
                                                           const char* pat=0)
{ // what == 0: all files
  // what == 1: regular files 
  // what == 2: sub-directories

 L.clear();

 if (!is_directory(dir_name))
 { 
   std::cerr << "read_directory:" << dir_name << " is not a directory."; 
   return;
  }

 DIR* dir_p = opendir(dir_name.c_str());

 dirent* dir_e;
 while ( (dir_e = readdir(dir_p)) != NULL )
 { string fname = dir_e->d_name;
   if (pat != 0 && fnmatch(pat,fname.c_str(),0)) continue;
   if (what != 0)
   { string full_name = dir_name + "/" + fname;
     if (what == 1 && !is_file(full_name)) continue;
     if (what == 2 && !is_directory(full_name)) continue;
    }
   L.push_back(fname);
  }
 closedir(dir_p);
}


string tmp_file_name() { return tmpnam(NULL); }

}

#endif



// to do ...
#if defined(__win32__)

namespace CGAL {

using namespace std;

bool is_directory(string name)
{ DWORD att = GetFileAttributes(name.c_str());
  if (att == 0xFFFFFFFF) return false;
  return (att & FILE_ATTRIBUTE_DIRECTORY) != 0; 
 }


bool is_file(string name)
{ DWORD att = GetFileAttributes(name.c_str());
  if (att == 0xFFFFFFFF) return false;
  return !(att & FILE_ATTRIBUTE_DIRECTORY); 
  //WIN32_FIND_DATA fd;
  //HANDLE ha = FindFirstFile(name.cstring(),&fd);
  //return (ha && (ha != (HANDLE)0xffffffff);
 }


bool is_link(string) { return false; }


string get_directory()
{ char dir_buf[256];
  int len = GetCurrentDirectory(256,dir_buf);
  if (dir_buf[len-1] == '\\') dir_buf[len-1] = '\0';
  return string(dir_buf);
}


string set_directory(string new_dir)
{ string old_dir = get_directory();
  if (is_directory(new_dir.c_str()))
    SetCurrentDirectory(new_dir.c_str());
  else
    std::cerr << "set_directory:" << new_dir << " is not a directory.";    
    
  return old_dir;
}



static void read_directory(string dir_name, int what, std::list<string>& L, 
                                                            const char* pat = 0)
{ 
  L.clear();

 if (!is_directory(dir_name.c_str()))
 { 
   std::cerr << "read_directory:" << dir_name << " is not a directory.";  
   return;
  }

  string cwd = set_directory(dir_name.c_str());

  if (pat == 0) pat = "*";

  WIN32_FIND_DATA fd;
  HANDLE ha = FindFirstFile(pat,&fd);

  if (ha && ha != (HANDLE)0xffffffff)
  { do { string fname = fd.cFileName;
         bool isdir = (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;
         if (what == 0 || (what == 1 && !isdir) || (what == 2 && isdir))
            L.push_back(fname);
     } while (FindNextFile(ha,&fd));
    FindClose(ha);
   }

  set_directory(cwd);
}



string tmp_file_name()
{ char name[MAX_PATH];
  char path[128];
  GetTempPath(128,path);
  time_t clock; 
  time(&clock);
  int r = ((int)clock + rand()) & 0xFFFF;
  GetTempFileName(path,"LEDA",r,name);
  return string(name);
 }

}

#endif




namespace CGAL {


using namespace std;

std::list<string> get_entries(string dir) 
{ std::list<string> L;
  read_directory(dir,0,L);
  return L;
 }

std::list<string> get_files(string dir)
{ std::list<string> L;
  read_directory(dir,1,L);
  return L;
 }


std::list<string> get_files(string dir, string pattern)
{ std::list<string> L;
  read_directory(dir,1,L,pattern.c_str());
  return L;
 }

std::list<string> get_directories(string dir)
{ std::list<string> L;
  read_directory(dir,2,L);
  return L;
 }


}

#if defined(__win32__)


namespace CGAL {


void message_box(const char* msg, const char* label)
{ std::cerr << label << ": " << msg << endl; 
//MessageBox(NULL,msg,label,MB_OK);  
}


}

#else

namespace CGAL {


void message_box(const char* msg, const char* label) 
{ std::cerr << label << ": " << msg << endl; }

}


#endif


