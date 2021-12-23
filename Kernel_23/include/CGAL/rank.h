// Copyright (c) 2021
// GeometryFactory (France),
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_RANK_H
#define CGAL_RANK_H

namespace CGAL {

template <class RT>
int rank_11(const RT& a0)
{
  return a0!=0 ? 1 : 0;
}

template <class RT>
int rank_21(const RT& a0, const RT& a1)
{
  return (a0!=0 || a1 !=0) ? 1 : 0;
}

template <class RT>
int rank_12(const RT& a0, const RT& a1)
{
  return (a0!=0 || a1 !=0) ? 1 : 0;
}

template <class RT>
int rank_31(const RT& a0, const RT& a1, const RT& a2)
{
  return (a0!=0 || a1 !=0 || a2 !=0) ? 1 : 0;
}

template <class RT>
int rank_32(const RT& a0, const RT& b0,
            const RT& a1, const RT& b1,
            const RT& a2, const RT& b2)
{
  if (a0==0)
  {
    if (a1==0)
    {
      if (a2==0)
      {
        return rank_31<RT>(b0,b1,b2);
      }
      else
      {
        return 1 + rank_21<RT>(b0, b1);
      }
    }
    else
    {
      return 1 + rank_21<RT>(b0, a1*b2-a2*b1);
    }
  }
  else
  {
    return 1 + rank_21<RT>(a0*b1-a1*b0, a0*b2-a2*b0);
  }
}

template <class RT>
int rank_22(const RT& a0, const RT& b0,
            const RT& a1, const RT& b1)
{
  if (a0==0)
  {
    if (a1==0)
    {
      return rank_21<RT>(b0,b1);
    }
    else
      return 1 + rank_11<RT>(b0);
  }
  return 1 + rank_11<RT>(a0*b1-a1*b0);
}

template <class RT>
int rank_33(const RT& a0, const RT& b0, const RT& c0,
            const RT& a1, const RT& b1, const RT& c1,
            const RT& a2, const RT& b2, const RT& c2)
{
  if (a0==0)
  {
    if (a1==0)
    {
      if (a2==0)
      {
        return rank_32<RT>(b0, c0, b1, c1, b2, c2);
      }
      else
      {
        return 1 + rank_22<RT>(b0, c0, b1, c1);
      }
    }
    else
      return 1 + rank_22<RT>(b0, c0, a1*b2-a2*b1, a1*c2-a2*c1);
  }
  else
  {
    return 1 + rank_22<RT>(a0*b1-a1*b0, a0*c1-a1*c0, a0*b2-a2*b0, a0*c2-a2*c0);
  }
}


template <class RT>
int rank_23(const RT& a0, const RT& b0, const RT& c0,
            const RT& a1, const RT& b1, const RT& c1)
{
  if (a0==0)
  {
    if (a1==0)
    {
      return rank_22<RT>(b0, c0, b1, c1);
    }
    else
      return 1 + rank_12<RT>(b0,c0);
  }
  else
  {
    return 1 + rank_12<RT>(a0*b1-a1*b0,a0*c1-a1*c0);
  }
}


template <class RT>
int rank_34(const RT& a0, const RT& b0, const RT& c0, const RT& d0,
            const RT& a1, const RT& b1, const RT& c1, const RT& d1,
            const RT& a2, const RT& b2, const RT& c2, const RT& d2)
{
  if (a0==0)
  {
    if (a1==0)
    {
      if (a2==0)
      {
        return rank_33<RT>(b0, c0, d0, b1, c1, d1, b2, c2, d2);
      }
      else
      {
        return 1 + rank_23<RT>(b0, c0, d0, b1, c1, d1);
      }
    }
    else
      return 1 + rank_23<RT>(b0, c0, d0,a1*b2-a2*b1, a1*c2-a2*c1, a1*d2 - a2*d1);
  }
  else
  {
    return 1 + rank_23<RT>(a0*b1-a1*b0, a0*c1-a1*c0, a0*d1-a1*d0, a0*b2-a2*b0, a0*c2-a2*c0, a0*d2-a2*d0);
  }
}

} // namespace CGAL

#endif // CGAL_RANK_H
