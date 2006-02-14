// Copyright (c) 1999  Martin-Luther-University Halle-Wittenberg (Germany).
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
// Author(s)     : Matthias Baesken, Algorithmic Solutions


#ifndef __CGAL_WINDOW_FILE__
#define __CGAL_WINDOW_FILE__

#if defined(CGAL_USE_CGAL_HEADERS)
#include <CGAL/basic.h>
#endif

#include <CGAL/LEDA/basic.h>
#include <list>


namespace CGAL {



extern __exportF void message_box(const char* msg, const char* label);



extern __exportF std::string set_directory(std::string new_dir);
/*{\Xfuncl     sets the current working directory to |new_dir| and 
               returns the name of the old cwd. }*/

extern __exportF std::string get_directory();
/*{\Xfuncl     returns the name of the current working directory. }*/


extern __exportF std::list<std::string> get_directories(std::string dir);
/*{\Xfuncl     returns the list of names of all sub-directories in 
               directory |dir|. }*/

extern __exportF std::list<std::string> get_files(std::string dir);
/*{\Xfuncl     returns the list of names of all regular files in 
               directory |dir|. }*/

extern __exportF std::list<std::string> get_files(std::string dir, std::string pattern);
/*{\Xfuncl     returns the list of names of all regular files in 
               directory |dir| matching pattern. }*/

extern __exportF std::list<std::string> get_entries(std::string dir);
/*{\Xfuncl     returns the list of all entries (directory and files) of 
               directory |dir|. }*/


extern __exportF bool is_directory(std::string fname);
/*{\Xfuncl     returns true if |fname| is the path name of a directory
               and false otherwise. }*/

extern __exportF bool is_file(std::string fname);
/*{\Xfuncl     returns true if |fname| is the path name of a regular file
               and false otherwise. }*/

extern __exportF bool is_link(std::string fname);
/*{\Xfuncl     returns true if |fname| is the path name of a symbolic link
               and false otherwise. }*/


extern __exportF std::string tmp_file_name();
/*{\Xfuncl     returns a unique name for a temporary file. }*/

}

#endif
