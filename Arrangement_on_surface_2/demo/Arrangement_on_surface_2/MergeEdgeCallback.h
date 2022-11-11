// Copyright (c) 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

#ifndef MERGE_EDGE_CALLBACK_H
#define MERGE_EDGE_CALLBACK_H

#include "Callback.h"
#include <CGAL/Object.h>

namespace demo_types
{
enum class TraitsType : int;
}

/**
   Handles merging of arrangement curves selected from the scene.
*/
class MergeEdgeCallbackBase : public CGAL::Qt::Callback
{
public:
  static MergeEdgeCallbackBase*
  create(demo_types::TraitsType, CGAL::Object arr_obj, QObject* parent);

protected:
  using CGAL::Qt::Callback::Callback;
};

#endif // MERGE_EDGE_CALLBACK_H
