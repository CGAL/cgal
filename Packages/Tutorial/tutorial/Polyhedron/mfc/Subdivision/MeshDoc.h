// MeshDoc.h : interface of the CMeshDoc class
//

#pragma once

// CGAL
#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include "lib/enriched_polyhedron.h"

typedef float number_type;
typedef CGAL::Simple_cartesian<number_type> Enriched_kernel;

class CMeshDoc : public CDocument
{
protected: // create from serialization only
	CMeshDoc();
	DECLARE_DYNCREATE(CMeshDoc)

// Attributes
public:

	Enriched_polyhedron<Enriched_kernel,Enriched_items> *m_pMesh;

// Operations
public:

	// status message
	void StatusMessage( char* fmt, ... );
	void ResetMeshProperties();
	void UpdateMeshProperties(bool update_component = false,
		                        bool update_boundary = false);

// Overrides
	public:
	virtual void Serialize(CArchive& ar);

// Implementation
public:
	virtual ~CMeshDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	DECLARE_MESSAGE_MAP()
public:
	virtual BOOL OnOpenDocument(LPCTSTR lpszPathName);
	virtual BOOL OnSaveDocument(LPCTSTR lpszPathName);
	afx_msg void OnSubdivisionSqrt3();
	afx_msg void OnUpdateSubdivisionSqrt3(CCmdUI *pCmdUI);
	afx_msg void OnSubdivisionSqrt3Twice();
	afx_msg void OnUpdateSubdivisionSqrt3Twice(CCmdUI *pCmdUI);
	afx_msg void OnSubdivisionQuad();
	afx_msg void OnUpdateSubdivisionQuad(CCmdUI *pCmdUI);
	afx_msg void OnSubdivisionDoo();
	afx_msg void OnUpdateSubdivisionDoo(CCmdUI *pCmdUI);
	afx_msg void OnSubdivisionCatmull();
	afx_msg void OnSubdivisionLoop();
	afx_msg void OnUpdateSubdivisionLoop(CCmdUI *pCmdUI);
	afx_msg void OnSubdivisionCatmull32868();
	afx_msg void OnUpdateSubdivisionCatmull32868(CCmdUI *pCmdUI);
};


