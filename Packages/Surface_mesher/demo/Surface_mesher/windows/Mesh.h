// Mesh.h : main header file for the Mesh application
//
#pragma once

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "Lib/arcball/arcball.h"
#include "resource.h"

class CMeshApp : public CWinApp
{
public:
	CMeshApp();

	// arcball for copy/paste viewpoint
	CArcball m_Arcball;

// Overrides
public:
	virtual BOOL InitInstance();

// Implementation
	DECLARE_MESSAGE_MAP()
	afx_msg void OnCameraCopyviewpoint();
	afx_msg void OnCameraPasteviewpoint();
	afx_msg void OnUpdateCameraCopyviewpoint(CCmdUI *pCmdUI);
	afx_msg void OnUpdateCameraPasteviewpoint(CCmdUI *pCmdUI);
};

extern CMeshApp theApp;