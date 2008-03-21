// Gyroviz.h : main header file for the Gyroviz application
//
#pragma once

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"       // main symbols


// CGyrovizApp:
// See Gyroviz.cpp for the implementation of this class
//

class CGyrovizApp : public CWinApp
{
public:
	CGyrovizApp();


// Overrides
public:
	virtual BOOL InitInstance();

// Implementation
	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()
};

extern CGyrovizApp theApp;