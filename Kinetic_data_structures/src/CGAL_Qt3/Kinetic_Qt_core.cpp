// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#include <CGAL/basic.h>

#include <CGAL/Kinetic/IO/internal/Qt_core.h>

#include "Kinetic_Qt_core.moc"

namespace CGAL { namespace Kinetic { namespace internal {

Qt_core::Qt_core() {
  //playable_=NULL;
  mode_=STOP;
}

void Qt_core::play_button() {
  CGAL_LOG(Log::SOME, "Play button pushed.\n");
  mode_=RUN;
  CGAL_KINETIC_NOTIFY(LAST_BUTTON_PRESSED);
}

void Qt_core::pause_button() {
  CGAL_LOG(Log::SOME, "Pause button pushed.\n");
  mode_=PAUSE;
  CGAL_KINETIC_NOTIFY(LAST_BUTTON_PRESSED);
}

void Qt_core::stop_button() {
  CGAL_LOG(Log::SOME, "Stop button pushed.\n");
  mode_=STOP;
  CGAL_KINETIC_NOTIFY(LAST_BUTTON_PRESSED);
}

void Qt_core::play_to_button() {
  CGAL_LOG(Log::SOME, "Play_to button pushed.\n");
  mode_=RUN_TO;
  CGAL_KINETIC_NOTIFY(LAST_BUTTON_PRESSED);
}

void Qt_core::play_through_button() {
  CGAL_LOG(Log::SOME, "Play through button pushed.\n");
  mode_= RUN_THROUGH;
  CGAL_KINETIC_NOTIFY(LAST_BUTTON_PRESSED);
}
void Qt_core::reverse_button() {
  CGAL_LOG(Log::SOME, "Reverse button pushed.\n");
  mode_=REVERSE;
  CGAL_KINETIC_NOTIFY(LAST_BUTTON_PRESSED);
}
void Qt_core::faster_button() {
  CGAL_LOG(Log::SOME, "Faster button pushed.\n");
  mode_=FASTER;
  CGAL_KINETIC_NOTIFY(LAST_BUTTON_PRESSED);
}
void Qt_core::slower_button() {
  CGAL_LOG(Log::SOME, "Slower button pushed.\n");
  mode_=SLOWER;
  CGAL_KINETIC_NOTIFY(LAST_BUTTON_PRESSED);
}
} } } //namespace CGAL::Kinetic::internal
