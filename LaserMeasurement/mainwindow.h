#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "QDebug"
#include "QFileDialog"
#include "QMessageBox"
#include "QValidator"
#include <QImageReader>
#include <QStringList>
#include <QStandardPaths>

#include "string"
#include <vector>
#include "fstream"
#include "io.h"
#include "math.h"

#include "opencv.hpp"

#include "qimagemat.h"
#include "dialog.h"

using namespace std;
using namespace cv;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void open1Slot();
    void open2Slot();
    void changeOpen1Slot();
    void changeAcoeSlot();
    void changeBcoeSlot();
    void changeCcoeSlot();
    void changeDcoeSlot();
    void changeCustomSet();
    void startCalibSlot();
    void Calculate();
    void CalculatePoints();
    void AutoCreateSlots();
    void markSlot();
    void measureSlot();
    void unitChangeSlot();

private:
    Ui::MainWindow *ui;

    void init();
    void setConnect();
    void findPoints();
    void cvFitPlane(const CvMat* points, float* plane);
    void measureDist();

    QImage image;
    QLabel *imageLabel;
    QString imageFileName;

    Dialog *markDialog;

    vector<string> vecFileName1Seq;
    vector<string> vecFileName2Seq;
    int imageCount;
    Size imageSize;
    vector<vector<Point3f>> chessWPointsSeq;//保存棋盘世界坐标序列（实物）
    vector<vector<Point3f>> chessCPointsSeq;//保存棋盘世界坐标对应的相机坐标（实物）
    vector<Point3f> DCPointSeq;//实物交点的摄像机坐标（这里以D点为交点）
    Size cornerSize,SqureSize, chessSize;
    vector<vector<Point2f>> cornerPointsSeq;//角点的像素坐标
    vector<Vec4f> laserLineParaSeq;
    vector<vector<Point2f>> crossPointSeq; //保存图片激光与棋盘横向直线交点坐标

    vector<vector<Vec4f>> imHLineParaSeq;//拟合好的图片横向直线方程

    Mat cameraMatrix;
    Mat distCoeffs;
    vector<Mat> rvecsMatSeq;
    vector<Mat> rotationMatrixSeq;
    vector<Mat> tvecsMatSeq;
    double Acoe,Bcoe,Ccoe,Dcoe;//直线方程Ax+By+Cz+D= 0;

    Point startPoint, endPoint;

    double distance;

    Mat tImage;
};

#endif // MAINWINDOW_H
