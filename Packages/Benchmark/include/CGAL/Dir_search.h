// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $$
// release_date  : $$
//
// file          : include/CGAL/Dir_search.h
// package       : Planar_map (5.87)
// maintainer    : Efi Fogel <efif@post.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Efi Fogel <efif@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef DIR_SEARCH_H
#define DIR_SEARCH_H

#include <sys/types.h>
#include <sys/stat.h>
#if !(defined _MSC_VER)
#include <unistd.h>
#endif
#include <list>

/*!
 */
class Dir_search {
public:
  Dir_search() {}
  Dir_search(const char * dir) { add(dir); }
  void add(const char * dir) { m_dirs.push_back(dir); }
  void add(const std::string & dir) { m_dirs.push_back(dir); }
  bool find(const char * filename, std::string & fullname)
  {
    if (filename[0] == '/') {
      fullname = filename;
      return true;
    }
    Dirs::const_iterator di;
    struct stat buf;
    for (di = m_dirs.begin(); di != m_dirs.end(); di++) {
      fullname = (*di);
      fullname.append("/");
      fullname.append(std::string(filename));
      int rc = stat(fullname.c_str(), &buf);
      if (rc < 0 || ((buf.st_mode & S_IFDIR) == S_IFDIR)) continue;
      if (rc == 0) return true;
    }
    return false;
  }

private:
  typedef std::list<std::string> Dirs;

  Dirs m_dirs;
};

#endif
