#ifndef ARRANGEMENT_DEMO_RATIONAL_CURVE_INPUT_DIALOG_H
#define ARRANGEMENT_DEMO_RATIONAL_CURVE_INPUT_DIALOG_H

#include <QDialog>

namespace Ui
{
class RationalCurveInputDialog;
}

class RationalCurveInputDialog : public QDialog
{
    Q_OBJECT

public:
    explicit RationalCurveInputDialog(QWidget *parent = 0);
    ~RationalCurveInputDialog();
    std::string getNumeratorText();
    std::string getDenominatorText();
    Ui::RationalCurveInputDialog* getUi(){return this->ui;}

private:
    Ui::RationalCurveInputDialog *ui;
};

#endif // ALGEBRAICCURVEINPUTDIALOG_H

