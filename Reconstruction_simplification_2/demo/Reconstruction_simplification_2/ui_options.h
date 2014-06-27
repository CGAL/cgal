/********************************************************************************
** Form generated from reading UI file 'options.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_OPTIONS_H
#define UI_OPTIONS_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QDialog>
#include <QtGui/QDialogButtonBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QSpinBox>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Dialog_options
{
public:
    QDialogButtonBox *buttonBox;
    QCheckBox *use_flip_checkbox;
    QWidget *layoutWidget;
    QVBoxLayout *vboxLayout;
    QHBoxLayout *hboxLayout;
    QLabel *verbose_label;
    QSpinBox *verbose_spinbox;
    QHBoxLayout *hboxLayout1;
    QLabel *mchoice_label;
    QSpinBox *mchoice_spinbox;
    QHBoxLayout *hboxLayout2;
    QLabel *percent_label;
    QSpinBox *percent_spinbox;
    QHBoxLayout *hboxLayout3;
    QLabel *norm_tol_label;
    QDoubleSpinBox *norm_tol_spinbox;
    QHBoxLayout *hboxLayout4;
    QLabel *tang_tol_label;
    QDoubleSpinBox *tang_tol_spinbox;
    QHBoxLayout *hboxLayout5;
    QLabel *alpha_label;
    QDoubleSpinBox *alpha_spinbox;
    QHBoxLayout *hboxLayout6;
    QLabel *relocation_label;
    QSpinBox *relocation_spinbox;
    QHBoxLayout *hboxLayout7;
    QLabel *ghost_vs_solid_label;
    QDoubleSpinBox *ghost_spinbox;
    QWidget *widget;
    QVBoxLayout *vboxLayout1;
    QHBoxLayout *hboxLayout8;
    QLabel *thickness;
    QDoubleSpinBox *thickness_spinbox;
    QHBoxLayout *hboxLayout9;
    QLabel *point_size;
    QDoubleSpinBox *point_size_spinbox;
    QHBoxLayout *hboxLayout10;
    QLabel *vertex_size;
    QDoubleSpinBox *vertex_size_spinbox;

    void setupUi(QDialog *Dialog_options)
    {
        if (Dialog_options->objectName().isEmpty())
            Dialog_options->setObjectName(QString::fromUtf8("Dialog_options"));
        Dialog_options->resize(409, 379);
        buttonBox = new QDialogButtonBox(Dialog_options);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setGeometry(QRect(40, 330, 341, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::NoButton|QDialogButtonBox::Ok);
        use_flip_checkbox = new QCheckBox(Dialog_options);
        use_flip_checkbox->setObjectName(QString::fromUtf8("use_flip_checkbox"));
        use_flip_checkbox->setGeometry(QRect(220, 20, 159, 21));
        layoutWidget = new QWidget(Dialog_options);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(20, 20, 188, 304));
        vboxLayout = new QVBoxLayout(layoutWidget);
        vboxLayout->setObjectName(QString::fromUtf8("vboxLayout"));
        vboxLayout->setContentsMargins(0, 0, 0, 0);
        hboxLayout = new QHBoxLayout();
        hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
        verbose_label = new QLabel(layoutWidget);
        verbose_label->setObjectName(QString::fromUtf8("verbose_label"));

        hboxLayout->addWidget(verbose_label);

        verbose_spinbox = new QSpinBox(layoutWidget);
        verbose_spinbox->setObjectName(QString::fromUtf8("verbose_spinbox"));
        verbose_spinbox->setMaximum(2);

        hboxLayout->addWidget(verbose_spinbox);


        vboxLayout->addLayout(hboxLayout);

        hboxLayout1 = new QHBoxLayout();
        hboxLayout1->setObjectName(QString::fromUtf8("hboxLayout1"));
        mchoice_label = new QLabel(layoutWidget);
        mchoice_label->setObjectName(QString::fromUtf8("mchoice_label"));

        hboxLayout1->addWidget(mchoice_label);

        mchoice_spinbox = new QSpinBox(layoutWidget);
        mchoice_spinbox->setObjectName(QString::fromUtf8("mchoice_spinbox"));
        mchoice_spinbox->setMaximum(100);
        mchoice_spinbox->setValue(100);

        hboxLayout1->addWidget(mchoice_spinbox);


        vboxLayout->addLayout(hboxLayout1);

        hboxLayout2 = new QHBoxLayout();
        hboxLayout2->setObjectName(QString::fromUtf8("hboxLayout2"));
        percent_label = new QLabel(layoutWidget);
        percent_label->setObjectName(QString::fromUtf8("percent_label"));

        hboxLayout2->addWidget(percent_label);

        percent_spinbox = new QSpinBox(layoutWidget);
        percent_spinbox->setObjectName(QString::fromUtf8("percent_spinbox"));
        percent_spinbox->setMaximum(100);
        percent_spinbox->setValue(100);

        hboxLayout2->addWidget(percent_spinbox);


        vboxLayout->addLayout(hboxLayout2);

        hboxLayout3 = new QHBoxLayout();
        hboxLayout3->setObjectName(QString::fromUtf8("hboxLayout3"));
        norm_tol_label = new QLabel(layoutWidget);
        norm_tol_label->setObjectName(QString::fromUtf8("norm_tol_label"));

        hboxLayout3->addWidget(norm_tol_label);

        norm_tol_spinbox = new QDoubleSpinBox(layoutWidget);
        norm_tol_spinbox->setObjectName(QString::fromUtf8("norm_tol_spinbox"));
        norm_tol_spinbox->setDecimals(3);
        norm_tol_spinbox->setMaximum(100);
        norm_tol_spinbox->setSingleStep(0.05);
        norm_tol_spinbox->setValue(100);

        hboxLayout3->addWidget(norm_tol_spinbox);


        vboxLayout->addLayout(hboxLayout3);

        hboxLayout4 = new QHBoxLayout();
        hboxLayout4->setObjectName(QString::fromUtf8("hboxLayout4"));
        tang_tol_label = new QLabel(layoutWidget);
        tang_tol_label->setObjectName(QString::fromUtf8("tang_tol_label"));

        hboxLayout4->addWidget(tang_tol_label);

        tang_tol_spinbox = new QDoubleSpinBox(layoutWidget);
        tang_tol_spinbox->setObjectName(QString::fromUtf8("tang_tol_spinbox"));
        tang_tol_spinbox->setDecimals(3);
        tang_tol_spinbox->setMaximum(100);
        tang_tol_spinbox->setSingleStep(0.05);
        tang_tol_spinbox->setValue(100);

        hboxLayout4->addWidget(tang_tol_spinbox);


        vboxLayout->addLayout(hboxLayout4);

        hboxLayout5 = new QHBoxLayout();
        hboxLayout5->setObjectName(QString::fromUtf8("hboxLayout5"));
        alpha_label = new QLabel(layoutWidget);
        alpha_label->setObjectName(QString::fromUtf8("alpha_label"));

        hboxLayout5->addWidget(alpha_label);

        alpha_spinbox = new QDoubleSpinBox(layoutWidget);
        alpha_spinbox->setObjectName(QString::fromUtf8("alpha_spinbox"));
        alpha_spinbox->setDecimals(3);
        alpha_spinbox->setMaximum(100);
        alpha_spinbox->setSingleStep(0.05);
        alpha_spinbox->setValue(100);

        hboxLayout5->addWidget(alpha_spinbox);


        vboxLayout->addLayout(hboxLayout5);

        hboxLayout6 = new QHBoxLayout();
        hboxLayout6->setObjectName(QString::fromUtf8("hboxLayout6"));
        relocation_label = new QLabel(layoutWidget);
        relocation_label->setObjectName(QString::fromUtf8("relocation_label"));

        hboxLayout6->addWidget(relocation_label);

        relocation_spinbox = new QSpinBox(layoutWidget);
        relocation_spinbox->setObjectName(QString::fromUtf8("relocation_spinbox"));

        hboxLayout6->addWidget(relocation_spinbox);


        vboxLayout->addLayout(hboxLayout6);

        hboxLayout7 = new QHBoxLayout();
        hboxLayout7->setObjectName(QString::fromUtf8("hboxLayout7"));
        ghost_vs_solid_label = new QLabel(layoutWidget);
        ghost_vs_solid_label->setObjectName(QString::fromUtf8("ghost_vs_solid_label"));

        hboxLayout7->addWidget(ghost_vs_solid_label);

        ghost_spinbox = new QDoubleSpinBox(layoutWidget);
        ghost_spinbox->setObjectName(QString::fromUtf8("ghost_spinbox"));
        ghost_spinbox->setDecimals(3);
        ghost_spinbox->setMaximum(100);
        ghost_spinbox->setSingleStep(0.05);
        ghost_spinbox->setValue(100);

        hboxLayout7->addWidget(ghost_spinbox);


        vboxLayout->addLayout(hboxLayout7);

        widget = new QWidget(Dialog_options);
        widget->setObjectName(QString::fromUtf8("widget"));
        widget->setGeometry(QRect(220, 210, 174, 109));
        vboxLayout1 = new QVBoxLayout(widget);
        vboxLayout1->setObjectName(QString::fromUtf8("vboxLayout1"));
        vboxLayout1->setContentsMargins(0, 0, 0, 0);
        hboxLayout8 = new QHBoxLayout();
        hboxLayout8->setObjectName(QString::fromUtf8("hboxLayout8"));
        thickness = new QLabel(widget);
        thickness->setObjectName(QString::fromUtf8("thickness"));

        hboxLayout8->addWidget(thickness);

        thickness_spinbox = new QDoubleSpinBox(widget);
        thickness_spinbox->setObjectName(QString::fromUtf8("thickness_spinbox"));

        hboxLayout8->addWidget(thickness_spinbox);


        vboxLayout1->addLayout(hboxLayout8);

        hboxLayout9 = new QHBoxLayout();
        hboxLayout9->setObjectName(QString::fromUtf8("hboxLayout9"));
        point_size = new QLabel(widget);
        point_size->setObjectName(QString::fromUtf8("point_size"));

        hboxLayout9->addWidget(point_size);

        point_size_spinbox = new QDoubleSpinBox(widget);
        point_size_spinbox->setObjectName(QString::fromUtf8("point_size_spinbox"));

        hboxLayout9->addWidget(point_size_spinbox);


        vboxLayout1->addLayout(hboxLayout9);

        hboxLayout10 = new QHBoxLayout();
        hboxLayout10->setObjectName(QString::fromUtf8("hboxLayout10"));
        vertex_size = new QLabel(widget);
        vertex_size->setObjectName(QString::fromUtf8("vertex_size"));

        hboxLayout10->addWidget(vertex_size);

        vertex_size_spinbox = new QDoubleSpinBox(widget);
        vertex_size_spinbox->setObjectName(QString::fromUtf8("vertex_size_spinbox"));

        hboxLayout10->addWidget(vertex_size_spinbox);


        vboxLayout1->addLayout(hboxLayout10);


        retranslateUi(Dialog_options);
        QObject::connect(buttonBox, SIGNAL(accepted()), Dialog_options, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), Dialog_options, SLOT(reject()));

        QMetaObject::connectSlotsByName(Dialog_options);
    } // setupUi

    void retranslateUi(QDialog *Dialog_options)
    {
        Dialog_options->setWindowTitle(QApplication::translate("Dialog_options", "Options", 0, QApplication::UnicodeUTF8));
        use_flip_checkbox->setText(QApplication::translate("Dialog_options", "use flip", 0, QApplication::UnicodeUTF8));
        verbose_label->setText(QApplication::translate("Dialog_options", "Verbose", 0, QApplication::UnicodeUTF8));
        mchoice_label->setText(QApplication::translate("Dialog_options", "Multiply Choice #", 0, QApplication::UnicodeUTF8));
        percent_label->setText(QApplication::translate("Dialog_options", "Init with %", 0, QApplication::UnicodeUTF8));
        norm_tol_label->setText(QApplication::translate("Dialog_options", "Normal Tol", 0, QApplication::UnicodeUTF8));
        tang_tol_label->setText(QApplication::translate("Dialog_options", "Tangential Tol", 0, QApplication::UnicodeUTF8));
        alpha_label->setText(QApplication::translate("Dialog_options", "Alpha", 0, QApplication::UnicodeUTF8));
        relocation_label->setText(QApplication::translate("Dialog_options", "Relocation", 0, QApplication::UnicodeUTF8));
        ghost_vs_solid_label->setText(QApplication::translate("Dialog_options", "Ghost vs solid", 0, QApplication::UnicodeUTF8));
        thickness->setText(QApplication::translate("Dialog_options", "Line thickness", 0, QApplication::UnicodeUTF8));
        point_size->setText(QApplication::translate("Dialog_options", "Point size", 0, QApplication::UnicodeUTF8));
        vertex_size->setText(QApplication::translate("Dialog_options", "Vertex size", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class Dialog_options: public Ui_Dialog_options {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_OPTIONS_H
