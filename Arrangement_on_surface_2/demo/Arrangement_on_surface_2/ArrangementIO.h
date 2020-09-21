// Copyright (c) 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ahmed Essam <theartful.ae@gmail.com>

#ifndef ARRANGEMENT_DEMO_IO_H
#define ARRANGEMENT_DEMO_IO_H

#include <fstream>
#include <vector>

class ArrangementDemoTabBase;
namespace CGAL
{
    class Object;
}

struct ArrangementIO
{
  ArrangementDemoTabBase* read(std::ifstream&);
  bool write(ArrangementDemoTabBase*, std::ofstream&);
};

// removed from ArrangementDemoWindow to speed up its compilation speed
// ArrangementDemoWindow doesn't have to specialize any arrangement related
// functions and now only deals with GUI logic
// TODO: move this somewhere else?
ArrangementDemoTabBase* makeOverlayTab(const std::vector<CGAL::Object>&);

#endif // ARRANGEMENT_DEMO_IO_H
