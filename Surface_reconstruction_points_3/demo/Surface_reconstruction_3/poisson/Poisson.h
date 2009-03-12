// Poisson.h : main header file for the Poisson application
//
#pragma once

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"       // main symbols


// CPoissonApp:
// See Poisson.cpp for the implementation of this class
//

class CPoissonApp : public CWinApp
{
public:
	CPoissonApp();


// Overrides
public:
	virtual BOOL InitInstance();

// Implementation
	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()
};

extern CPoissonApp theApp;