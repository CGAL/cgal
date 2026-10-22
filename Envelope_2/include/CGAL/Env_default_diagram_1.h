// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Ron Wein   <wein@post.tau.ac.il>
//             Efi Fogel  <efifogel@gmail.com>

#ifndef CGAL_ENV_DEFAULT_DIAGRAM_1_H
#define CGAL_ENV_DEFAULT_DIAGRAM_1_H

#include <CGAL/license/Envelope_2.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/Env_default_diagram_1.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/Envelope_2/Envelope_diagram_1.h>"
#include <CGAL/Installation/internal/deprecation_warning.h>

#include <CGAL/Envelope_2/Envelope_diagram_1.h>

namespace CGAL {

template <typename Traits_, typename Allocator = CGAL_ALLOCATOR(int)>
class Env_default_diagram_1 : public Envelope_2::Envelope_diagram_1<Traits_, Allocator> {
  using Base = Envelope_2::Envelope_diagram_1<Traits_, Allocator>;

public:
  using Base::Base;
};

} //namespace CGAL

#endif
