#ifndef DIR_SEARCH_H
#define DIR_SEARCH_H

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
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
