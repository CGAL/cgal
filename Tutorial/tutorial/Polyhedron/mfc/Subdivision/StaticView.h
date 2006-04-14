#if !defined(AFX_STATICVIEW_H__CFA8BDC2_788D_45CC_8DC5_8EA5C01C3587__INCLUDED_)
#define AFX_STATICVIEW_H__CFA8BDC2_788D_45CC_8DC5_8EA5C01C3587__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// StaticView.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CStaticView window

class CStaticView : public CStatic
{
// Construction
public:
	CStaticView();
	BOOL CreateViewGLContext(HDC hDC);
	BOOL SetWindowPixelFormat(HDC hDC);
	void BuildListCube();
	int Init(); 

	// OpenGL specific
	HGLRC m_hGLContext;
	int m_GLPixelIndex;
	float m_xRotation;
	float m_yRotation;
	float m_zRotation;
	BOOL m_xyRotation;

	float m_xTranslation;
	float m_yTranslation;
	float m_zTranslation;

	float m_xScaling;
	float m_yScaling;
	float m_zScaling;

	float m_SpeedTranslation;
	float m_SpeedRotation;

	// Colors
	float m_ClearColorRed;
	float m_ClearColorGreen;
	float m_ClearColorBlue;

	// Animation
	float m_StepRotationX;
	float m_StepRotationY;
	float m_StepRotationZ;


// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CStaticView)
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CStaticView();

	// Generated message map functions
protected:
	//{{AFX_MSG(CStaticView)
	afx_msg void OnPaint();
	afx_msg void OnDestroy();
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnTimer(UINT nIDEvent);
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	//}}AFX_MSG

	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_STATICVIEW_H__CFA8BDC2_788D_45CC_8DC5_8EA5C01C3587__INCLUDED_)
