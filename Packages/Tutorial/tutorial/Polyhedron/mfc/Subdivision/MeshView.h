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
	bool m_Lighting;
	bool m_Culling;
	int m_PolygonMode;
	bool m_SmoothShading;
	bool m_UseNormals;
	bool m_FirstView;
	bool m_SuperimposeEdges;
	bool m_SuperimposeVertices;
	bool m_SuperimposeOnlyControlMesh;
	bool m_ThickerControlEdges;
	bool m_Antialiasing;
	float m_ThicknessControlEdges;
	float m_PointSize;
	bool m_DrawBoundingBox;
	bool m_DrawBoundingBoxWhenMoving;
	bool m_DrawVoronoiEdges;

	// colors
	float m_BackColor[3];
	float m_MeshColor[3];
	float m_EdgeColor[3];
	float m_VertexColor[3];
	float m_ControlEdgeColor[3];

// Operations
public:

	// arcball
	CArcball *arcball() { return &m_Arcball; }

	// OpenGL specific
	BOOL SetWindowPixelFormat(HDC hDC);
	BOOL CreateViewGLContext(HDC hDC);

	// edit color
	bool EditColor(float pColor[3]);

	// change material
	void ChangeMaterial(CString &string,bool update = true);

	// camera
	void InitCamera();
	void ViewAll(bool check_first = true);

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
	afx_msg void OnRenderLight();
	afx_msg void OnUpdateRenderLight(CCmdUI *pCmdUI);
	afx_msg void OnRenderCulling();
	afx_msg void OnUpdateRenderCulling(CCmdUI *pCmdUI);
	afx_msg void OnRenderBackcolor();
	afx_msg void OnEditCopy();
	afx_msg void OnModeFill();
	afx_msg void OnUpdateModeFill(CCmdUI *pCmdUI);
	afx_msg void OnModeWireframe();
	afx_msg void OnUpdateModeWireframe(CCmdUI *pCmdUI);
	afx_msg void OnModePoint();
	afx_msg void OnUpdateModePoint(CCmdUI *pCmdUI);
	virtual void OnInitialUpdate();
	afx_msg void OnRenderSmooth();
	afx_msg void OnUpdateRenderSmooth(CCmdUI *pCmdUI);
	afx_msg void OnViewAll();
	afx_msg void OnRenderSuperimposeEdges();
	afx_msg void OnUpdateRenderSuperimposeEdges(CCmdUI *pCmdUI);
	afx_msg void OnColorsMesh();
	afx_msg void OnColorsSuperimposededges();
	afx_msg void OnRenderAntialiasing();
	afx_msg void OnUpdateRenderAntialiasing(CCmdUI *pCmdUI);
	afx_msg void OnPredefinedmodesMesh();
	afx_msg void OnUpdatePredefinedmodesMesh(CCmdUI *pCmdUI);
	afx_msg void OnPredefinedControlmesh();
	afx_msg void OnUpdatePredefinedControlmesh(CCmdUI *pCmdUI);
	afx_msg void OnPredefinedMesh();
	afx_msg void OnUpdatePredefinedMesh(CCmdUI *pCmdUI);
	afx_msg void OnColorsThickedges();
	afx_msg void OnRenderReflectionlines();
	afx_msg void OnUpdateRenderReflectionlines(CCmdUI *pCmdUI);
	afx_msg void OnRenderVisualchooser();
	afx_msg void OnSuperimposeVertices();
	afx_msg void OnUpdateSuperimposeVertices(CCmdUI *pCmdUI);
	afx_msg void OnColorsVertices();
	afx_msg void OnMaterialPearl();
	afx_msg void OnMaterialBrass();
	afx_msg void OnMaterialBlackplastic();
	afx_msg void OnMaterialGold();
	afx_msg void OnMaterialSilver();
	afx_msg void OnMaterialJade();
	afx_msg void OnMaterialRuby();
	afx_msg void OnFileDumptops();
	afx_msg void OnUpdateFileDumptops(CCmdUI *pCmdUI);
	afx_msg void OnBoundixboxShow();
	afx_msg void OnUpdateBoundixboxShow(CCmdUI *pCmdUI);
	afx_msg void OnBoundixboxShowwhenmoving();
	afx_msg void OnUpdateBoundixboxShowwhenmoving(CCmdUI *pCmdUI);
	afx_msg void OnSpecialDrawvoronoiedges();
	afx_msg void OnUpdateSpecialDrawvoronoiedges(CCmdUI *pCmdUI);
};

#ifndef _DEBUG  // debug version in MeshView.cpp
inline CMeshDoc* CMeshView::GetDocument() const
   { return reinterpret_cast<CMeshDoc*>(m_pDocument); }
#endif

