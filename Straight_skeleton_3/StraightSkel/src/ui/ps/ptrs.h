// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * @file   ui/ps/ptrs.h
 * @author Gernot Walzl
 * @date   2012-11-06
 */

#ifndef UI_PS_PTRS_H
#define UI_PS_PTRS_H

#include "smarter_ptr.h"

namespace ui { namespace ps {

class PSPrinter;
class SpacePSPrinter;
class PlanePSPrinter;

typedef SHARED_PTR<PSPrinter> PSPrinterSPtr;
typedef WEAK_PTR<PSPrinter> PSPrinterWPtr;
typedef SHARED_PTR<SpacePSPrinter> SpacePSPrinterSPtr;
typedef WEAK_PTR<SpacePSPrinter> SpacePSPrinterWPtr;
typedef SHARED_PTR<PlanePSPrinter> PlanePSPrinterSPtr;
typedef WEAK_PTR<PlanePSPrinter> PlanePSPrinterWPtr;

} }

#endif /* UI_PS_PTRS_H */
