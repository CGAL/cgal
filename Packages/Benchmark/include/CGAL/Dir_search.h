// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>
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
