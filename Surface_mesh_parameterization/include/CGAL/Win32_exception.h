// Copyright (c) 2007-08  INRIA Sophia-Antipolis (France).
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
//
// Author(s) : Laurent Saboret

#ifndef CGAL_WIN32_EXCEPTION_H
#define CGAL_WIN32_EXCEPTION_H

#include <exception>

#include <windows.h>
#include <eh.h>


namespace CGAL {

// C++ class wrapping a Win32 structured exception.
class Win32_exception
{
// Data
private:
    unsigned int m_seNumber;
    
// Public operations
public:
    Win32_exception(unsigned int structuredExceptionNumber) 
    : m_seNumber(structuredExceptionNumber) 
    {}
    
    unsigned int getStructuredExceptionNumber() const { return m_seNumber; }
};


/// \internal
/// Class Win32_exception_handler:
/// - Translate Win32 structured exceptions to C++ exceptions.
/// - Protect application against stack overflow using _resetstkoflw()
///   (see http://msdn.microsoft.com/en-us/library/89f73td2(VS.80).aspx).
///
/// Caution: requires /EHa compilation option.
class Win32_exception_handler
{
// Data
private:
    _se_translator_function m_previous_translator;
    
// Public operations
public:

  Win32_exception_handler()
  {
    // Protect application against next stack overflow
    if (needs_stack_reset(false))
    {
      std::cerr << "Win32_exception_handler: reset stack using _resetstkoflw()\n";
      if(!_resetstkoflw()) 
      {
        std::cerr << "Win32_exception_handler: _resetstkoflw() failed, exiting!\n";
        abort(); 
      }
    }

    // Set up a function to handle win32 exceptions,
    // including stack overflow exceptions.
    m_previous_translator = _set_se_translator(translate);
  }

  ~Win32_exception_handler()
  {
    // Do not call _resetstkoflw() in this function!

    // Restore previous win32 exceptions handler
    _set_se_translator(m_previous_translator);
  }
    
// Private operations
private:

  // Translate Win32 structured exceptions to C++ exceptions
  static
  void __cdecl translate(unsigned int code, _EXCEPTION_POINTERS*)
  {
     // For stack overflow exceptions, set m_needs_stack_reset.
     // Use minimal stack space in this function.
     // Do not call _resetstkoflw() in this function!
     if (code == EXCEPTION_STACK_OVERFLOW)
     {
        std::cerr << "Win32_exception_handler: stack overflow!\n";
        needs_stack_reset(true);
     }

     // Throw our own C++  exception object
     throw Win32_exception(code);
  }

  // This method is a trick to allocate a global variable in a header file:
  // - get the previous value of m_needs_stack_reset
  // - set the new value of m_needs_stack_reset
  static
  bool needs_stack_reset(bool new_value)
  {
    // default value
    static bool m_needs_stack_reset = false;
    
    bool previous_value = m_needs_stack_reset;
    m_needs_stack_reset = new_value;
    return previous_value;
  }
};


} //namespace CGAL

#endif // CGAL_WIN32_EXCEPTION_H
