// MeshView.h : interface of the CMeshView class

#pragma once

#include "Lib/arcball/Arcball.h"
#include "Lib/arcball/Camera.h"
#include "Lib/arcball/Viewport.h"

class CMeshView : public CView
{
protected: // create from serialization only
	CMeshView();
	DECLARE_DYNCREATE(CMeshView)

// Attributes
public:
	CMeshDoc* GetDocument() const;

private:

	// mouse
	bool m_LeftButtonDown;
	bool m_RightButtonDown;
	CPoint m_LeftDownPos;
	CPoint m_RightDownPos;
	bool m_Moving;

	// arcball
	CCamera   m_Camera;
	CViewport m_Viewport;
	CArcball  m_Arcball;

	// OpenGL
	HGLRC m_hGLContext;
	int m_GLPixelIndex;

// Operations
public:

	// arcball
	CArcball *arcball() { return &m_Arcball; }

	// OpenGL specific
	BOOL SetWindowPixelFormat(HDC hDC);
	BOOL CreateViewGLContext(HDC hDC);

	// camera
	void InitCamera();

	// mouse 
	void HandleMouseButton(int x, int y);

  // Overrides
	public:
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);

// Implementation
public:
	virtual ~CMeshView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	DECLARE_MESSAGE_MAP()
public:
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnPaint();
	afx_msg void OnSize(UINT nType, int cx, int cy);
	virtual void OnActivateView(BOOL bActivate, CView* pActivateView, CView* pDeactiveView);
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	afx_msg void OnDestroy();
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	virtual void OnInitialUpdate();
};

#ifndef _DEBUG  // debug version in MeshView.cpp
inline CMeshDoc* CMeshView::GetDocument() const
   { return reinterpret_cast<CMeshDoc*>(m_pDocument); }
#endif

