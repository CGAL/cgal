// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_BENCHMARK_DIR_SEARCH_HPP
#define CGAL_BENCHMARK_DIR_SEARCH_HPP

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
