// Copyright (C) 2000 Stephen Cleary
//
// This file can be redistributed and/or modified under the terms found
//  in "copyright.html"
// This software and its documentation is provided "as is" without express or
//  implied warranty, and with no claim as to its suitability for any purpose.
//
// See http://www.boost.org for updates, documentation, and revision history.

// This file was AUTOMATICALLY GENERATED from "stdin"
//  Do NOT include directly!
//  Do NOT edit!

template <typename T0>
element_type * construct(const T0 & a0)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1>
element_type * construct(const T0 & a0, const T1 & a1)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const T0 & a0, const T1 & a1, const T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}

