#pragma once


// CDialogRenderingMode dialog

class CDialogRenderingMode : public CDialog
{
	DECLARE_DYNAMIC(CDialogRenderingMode)

public:
	CDialogRenderingMode(CWnd* pParent = NULL);   // standard constructor
	virtual ~CDialogRenderingMode();
	unsigned int m_Mode;

// Dialog Data
	enum { IDD = IDD_DIALOG_RENDERING };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnClose();
	afx_msg void OnDestroy();
};
