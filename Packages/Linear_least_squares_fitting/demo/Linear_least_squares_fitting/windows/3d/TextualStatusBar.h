// TextualStatusBar.h : header file
//

class CTextualStatusBar : public CStatusBar
{
// Data
private:
	CToolTipCtrl m_wndToolTip;

// Construction
public:
	CTextualStatusBar();

// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CTextualStatusBar)
	public:
	virtual BOOL PreTranslateMessage(MSG* pMsg);
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CTextualStatusBar();

	// Generated message map functions
protected:
	//{{AFX_MSG(CTextualStatusBar)
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnDestroy();
	afx_msg UINT OnNcHitTest(CPoint point);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////
