// PoissonView.h : interface of the CPoissonView class
//

#pragma once

#include "arcball.h"

class CPoissonView : public CView
{
protected: // create from serialization only
	CPoissonView();
	DECLARE_DYNCREATE(CPoissonView)

// Attributes
public:
	CPoissonDoc* GetDocument() const;

// Operations
public:

	// mouse
	bool m_LeftButtonDown;
	bool m_MiddleButtonDown;
	bool m_RightButtonDown;
	CPoint m_prev_mouse_pos;
	bool m_moving;
	
	// keyboard
	UINT m_KeyboardFlags;

	// arcball
	Arcball *m_pArcball;
	bool first_paint;

	// OpenGL
	HGLRC m_hGLContext;
	int m_GLPixelIndex;

	// Drawing options
	bool m_view_normals;
	bool m_view_original_normals;
	bool m_view_contour;
	bool m_view_surface;
	bool m_view_delaunay_edges;
	bool m_view_points;
	bool m_view_arcball;

	// OpenGL specific
	BOOL SetWindowPixelFormat(HDC hDC);
	BOOL CreateViewGLContext(HDC hDC);

// Overrides
	public:
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
virtual BOOL PreCreateWindow(CREATESTRUCT& cs);

// Implementation
public:
	virtual ~CPoissonView();
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
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnMButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnMButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnUpdateRenderPoints(CCmdUI *pCmdUI);
	afx_msg void OnRenderPoints();
	afx_msg void OnRenderNormals();
	afx_msg void OnUpdateRenderNormals(CCmdUI *pCmdUI);
	afx_msg void OnRenderSurface();
	afx_msg void OnUpdateRenderSurface(CCmdUI *pCmdUI);
	afx_msg void OnArcballReset();
    afx_msg void OnRenderArcball();
    afx_msg void OnUpdateRenderArcball(CCmdUI *pCmdUI);
    afx_msg void OnRenderOriginalnormals();
    afx_msg void OnUpdateRenderOriginalnormals(CCmdUI *pCmdUI);
};

#ifndef _DEBUG  // debug version in PoissonView.cpp
inline CPoissonDoc* CPoissonView::GetDocument() const
   { return reinterpret_cast<CPoissonDoc*>(m_pDocument); }
#endif

