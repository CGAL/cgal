// Copyright (c) 2019  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Kaimo Hu

#ifndef _CONSOLE_COLOR_H_
#define _CONSOLE_COLOR_H_

#include <iostream>

#if defined(WIN32)
#include <windows.h>
#endif

inline std::ostream& blue(std::ostream &s)
{
#if defined(WIN32)
    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hStdout, FOREGROUND_BLUE|FOREGROUND_GREEN|FOREGROUND_INTENSITY);
#else
    s << "\e[0;34m";
#endif
    return s;
}

inline std::ostream& red(std::ostream &s)
{
#if defined(WIN32)
    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hStdout, FOREGROUND_RED|FOREGROUND_INTENSITY);
#else
    s << "\e[0;31m";
#endif
    return s;
}

inline std::ostream& green(std::ostream &s)
{
#if defined(WIN32)
    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hStdout, FOREGROUND_GREEN|FOREGROUND_INTENSITY);
#else
    s << "\e[0;32m";
#endif
    return s;
}

inline std::ostream& yellow(std::ostream &s)
{
#if defined(WIN32)
    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hStdout, FOREGROUND_GREEN|FOREGROUND_RED|FOREGROUND_INTENSITY);
#else
    s << "\e[0;33m";
#endif
    return s;
}

inline std::ostream& white(std::ostream &s)
{
#if defined(WIN32)
    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hStdout, FOREGROUND_RED|FOREGROUND_GREEN|FOREGROUND_BLUE);
#else
    s << "\e[0;37m";
#endif
    return s;
}

#endif
