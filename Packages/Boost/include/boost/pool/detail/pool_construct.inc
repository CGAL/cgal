// Copyright (C) 2000 Stephen Cleary
//
// Permission to copy, use, and distribute this software is granted, provided
// that this copyright notice appears in all copies.
// Permission to modify the code and to distribute modified code is granted,
// provided that this copyright notice appears in all copies, and a notice
// that the code was modified is included with the copyright notice.
//
// This software is provided "as is" without express or implied warranty, and
// with no claim as to its suitability for any purpose.
//
// See http://www.boost.org for updates, documentation, and revision history.

// This file was AUTOMATICALLY GENERATED from "pool_c~1.m4"
//  Do NOT include directly!
//  Do NOT edit!

template <typename T0>
element_type * construct(T0 & a0)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0); }
  catch (...) { free(ret); throw; }
  return ret;
}
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
template <typename T0>
element_type * construct(volatile T0 & a0)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0>
element_type * construct(const volatile T0 & a0)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1>
element_type * construct(T0 & a0, T1 & a1)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1>
element_type * construct(const T0 & a0, T1 & a1)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1>
element_type * construct(volatile T0 & a0, T1 & a1)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1>
element_type * construct(const volatile T0 & a0, T1 & a1)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1>
element_type * construct(T0 & a0, const T1 & a1)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1); }
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
template <typename T0, typename T1>
element_type * construct(volatile T0 & a0, const T1 & a1)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1>
element_type * construct(const volatile T0 & a0, const T1 & a1)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1>
element_type * construct(T0 & a0, volatile T1 & a1)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1>
element_type * construct(const T0 & a0, volatile T1 & a1)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1>
element_type * construct(volatile T0 & a0, volatile T1 & a1)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1>
element_type * construct(const volatile T0 & a0, volatile T1 & a1)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1>
element_type * construct(T0 & a0, const volatile T1 & a1)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1>
element_type * construct(const T0 & a0, const volatile T1 & a1)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1>
element_type * construct(volatile T0 & a0, const volatile T1 & a1)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1>
element_type * construct(const volatile T0 & a0, const volatile T1 & a1)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(T0 & a0, T1 & a1, T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const T0 & a0, T1 & a1, T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(volatile T0 & a0, T1 & a1, T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const volatile T0 & a0, T1 & a1, T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(T0 & a0, const T1 & a1, T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const T0 & a0, const T1 & a1, T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(volatile T0 & a0, const T1 & a1, T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const volatile T0 & a0, const T1 & a1, T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(T0 & a0, volatile T1 & a1, T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const T0 & a0, volatile T1 & a1, T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(volatile T0 & a0, volatile T1 & a1, T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const volatile T0 & a0, volatile T1 & a1, T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(T0 & a0, const volatile T1 & a1, T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const T0 & a0, const volatile T1 & a1, T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(volatile T0 & a0, const volatile T1 & a1, T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const volatile T0 & a0, const volatile T1 & a1, T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(T0 & a0, T1 & a1, const T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const T0 & a0, T1 & a1, const T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(volatile T0 & a0, T1 & a1, const T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const volatile T0 & a0, T1 & a1, const T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(T0 & a0, const T1 & a1, const T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
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
template <typename T0, typename T1, typename T2>
element_type * construct(volatile T0 & a0, const T1 & a1, const T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const volatile T0 & a0, const T1 & a1, const T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(T0 & a0, volatile T1 & a1, const T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const T0 & a0, volatile T1 & a1, const T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(volatile T0 & a0, volatile T1 & a1, const T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const volatile T0 & a0, volatile T1 & a1, const T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(T0 & a0, const volatile T1 & a1, const T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const T0 & a0, const volatile T1 & a1, const T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(volatile T0 & a0, const volatile T1 & a1, const T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const volatile T0 & a0, const volatile T1 & a1, const T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(T0 & a0, T1 & a1, volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const T0 & a0, T1 & a1, volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(volatile T0 & a0, T1 & a1, volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const volatile T0 & a0, T1 & a1, volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(T0 & a0, const T1 & a1, volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const T0 & a0, const T1 & a1, volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(volatile T0 & a0, const T1 & a1, volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const volatile T0 & a0, const T1 & a1, volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(T0 & a0, volatile T1 & a1, volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const T0 & a0, volatile T1 & a1, volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(volatile T0 & a0, volatile T1 & a1, volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const volatile T0 & a0, volatile T1 & a1, volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(T0 & a0, const volatile T1 & a1, volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const T0 & a0, const volatile T1 & a1, volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(volatile T0 & a0, const volatile T1 & a1, volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const volatile T0 & a0, const volatile T1 & a1, volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(T0 & a0, T1 & a1, const volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const T0 & a0, T1 & a1, const volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(volatile T0 & a0, T1 & a1, const volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const volatile T0 & a0, T1 & a1, const volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(T0 & a0, const T1 & a1, const volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const T0 & a0, const T1 & a1, const volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(volatile T0 & a0, const T1 & a1, const volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const volatile T0 & a0, const T1 & a1, const volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(T0 & a0, volatile T1 & a1, const volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const T0 & a0, volatile T1 & a1, const volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(volatile T0 & a0, volatile T1 & a1, const volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const volatile T0 & a0, volatile T1 & a1, const volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(T0 & a0, const volatile T1 & a1, const volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const T0 & a0, const volatile T1 & a1, const volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(volatile T0 & a0, const volatile T1 & a1, const volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}
template <typename T0, typename T1, typename T2>
element_type * construct(const volatile T0 & a0, const volatile T1 & a1, const volatile T2 & a2)
{
  element_type * const ret = malloc();
  if (ret == 0)
    return ret;
  try { new (ret) element_type(a0, a1, a2); }
  catch (...) { free(ret); throw; }
  return ret;
}

