#ifndef DIALOG_H
#define DIALOG_H

#include "qimagemat.h"
#include <QDialog>
#include <QImage>
#include <QToolTip>
#include <QPainter>
#include <QPen>
#include <QWheelEvent>
#include <QMessageBox>
#include <QScreen>
#include "opencv.hpp"
#include <vector>

using namespace std;
using namespace cv;

namespace Ui {
class Dialog;
}

class Dialog : public QDialog
{
    Q_OBJECT

public:
    explicit Dialog(QWidget *parent = 0);
    ~Dialog();


public:
    void setImage(QImage image);
    void windowReSize(QSize size);
    void getPoint(Point &startPoint, Point &endPoint);

private:
    Ui::Dialog *ui;
    void paintEvent(QPaintEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);

    QImage image;

    int selectStartPosx;
    int selectStartPosy;
    int imagePosx;
    int imagePosy;

    int imageCurx;
    int imageCury;
    int curAtWindowx;
    int curAtWindowy;

    Point startPoint;
    Point endPoint;
    int pointCount;

    bool mousePress;
};

#endif // DIALOG_H
