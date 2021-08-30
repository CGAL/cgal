// Copyright (c) 2019  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Kaimo Hu

#ifndef CGAL_POLYGON_MESH_PROCESSING_CONSOLE_COLOR_H
#define CGAL_POLYGON_MESH_PROCESSING_CONSOLE_COLOR_H

// C/C++
#include <iostream>
// CGAL
#include <CGAL/license/Polygon_mesh_processing/minimal_angle_remeshing.h>


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

#endif // CGAL_POLYGON_MESH_PROCESSING_CONSOLE_COLOR_H
