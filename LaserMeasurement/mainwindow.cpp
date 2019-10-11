#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    markDialog(new Dialog(this)),
    imageLabel(new QLabel),
    Acoe(0),
    Bcoe(0),
    Ccoe(0),
    Dcoe(0),
    imageCount(0)
{
    ui->setupUi(this);

    //添加imagelable控件
    imageLabel->setBackgroundRole(QPalette::Base);
    imageLabel->setScaledContents(true);

    ui->scrollArea->setWidget(imageLabel);
    ui->scrollArea->setWidgetResizable(true);
    //ui->scrollArea->setVisible(false);

    //一些设置输入框范围的初始化
    init();

    setConnect();//建立一些控件与方法的连接
}

MainWindow::~MainWindow()
{
    delete ui;
    delete imageLabel;
    delete markDialog;
}

void MainWindow::setConnect()
{
    //当棋盘格长宽选值框改变时判断值是否合理，然后改变载入图片按钮是否可选
    void(QSpinBox::*fpSpinBox)(int)=&QSpinBox::valueChanged;//获取SpinBox::valueChanged的指针
    connect(ui->rowNumSpB, fpSpinBox, this, &MainWindow::changeOpen1Slot);
    connect(ui->colNumSpB,fpSpinBox, this, &MainWindow::changeOpen1Slot);
    connect(ui->chessLenLE, &QLineEdit::textChanged, this, &MainWindow::changeOpen1Slot);

    //连接点击棋盘格载入图片时调用的方法
    connect(ui->openBtn1, &QPushButton::clicked, this, &MainWindow::open1Slot);

    //连接点击标定相机时调用的方法
    connect(ui->startCalibBtn, &QPushButton::clicked, this, &MainWindow::startCalibSlot);

    //自定义光平面方程的复选框连接的方法，将输入光平面方程的输入框置为可用或不可用
    connect(ui->customSetCB, &QCheckBox::stateChanged, this, &MainWindow::changeCustomSet);

    //打开标点测距图片
    connect(ui->openBtn2, &QPushButton::clicked, this, &MainWindow::open2Slot);

    //点击标点调用的方法
    connect(ui->markBtn, &QPushButton::clicked, this, &MainWindow::markSlot);

    //点击测距调用的方法
    connect(ui->measureBtn, &QPushButton::clicked, this, &MainWindow::measureSlot);

    //当单位复选框改变时调用改变单位方法
    void(QComboBox::*fpComboBox)(int)=&QComboBox::currentIndexChanged;
    connect(ui->unitCB, fpComboBox, this, &MainWindow::unitChangeSlot);

    //当光平面方程改变时调用的改变光平面方程的方法
    connect(ui->xCoeLE, &QLineEdit::textChanged, this, &MainWindow::changeAcoeSlot);
    connect(ui->yCoeLE, &QLineEdit::textChanged, this, &MainWindow::changeBcoeSlot);
    connect(ui->zCoeLE, &QLineEdit::textChanged, this, &MainWindow::changeCcoeSlot);
    connect(ui->constCoeLE, &QLineEdit::textChanged, this, &MainWindow::changeDcoeSlot);

    //当通过载入图片的方式计算光平面方程可点击此按钮自动生成光平面方程
    connect(ui->autoCreateBtn, &QPushButton::clicked, this, &MainWindow::AutoCreateSlots);
}

void MainWindow::open1Slot()
{
    //一个一个加载的方案
    QString fileName;
    fileName = QFileDialog::getOpenFileName(this,"载入需标定的图片", QDir::currentPath(),  "iamge(*.jpg *.bmp *.png)");
    if(fileName.isEmpty())
    {
        QMessageBox::information(this,"ERROR","载入图片失败");
        return;
    }
    //将棋盘格图片保存到序列
    vecFileName1Seq.push_back(string((const char*)fileName.toLocal8Bit()));

    fileName.clear();
    fileName = QFileDialog::getOpenFileName(this,"载入对应的带激光的图片", QDir::currentPath(),  "iamge(*.jpg *.bmp *.png)");
    if(fileName.isEmpty())
    {
        QMessageBox::information(this,"ERROR","载入图片失败");
        vecFileName1Seq.pop_back();//载入带激光的（结构光）图片失败时移除棋盘格图片，因为棋盘格图片和激光图片时要成对存在的
        return;//载入失败时跳出
    }
    //将带激光的（结构光）图片保存到序列
    vecFileName2Seq.push_back(string((const char*)fileName.toLocal8Bit()));

    Mat matImage;
    vector<Point2f> cornerPoints;
    matImage = imread(vecFileName1Seq[vecFileName1Seq.size() - 1],IMREAD_COLOR);

    cornerSize.width = ui->colNumSpB->value();//角点数
    cornerSize.height = ui->rowNumSpB->value();
    chessSize.width = ui->chessLenLE->text().toInt();//棋盘格边长
    chessSize.height = chessSize.width;

    //查找角点的像素坐标
    bool ok = findChessboardCorners(matImage, cornerSize, cornerPoints, CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_NORMALIZE_IMAGE);//提取角点

    if(!ok)
    {
        QMessageBox::information(this, "ERROR", "提取角点失败，请确认输入信息是否有误，或另换图片");
        vecFileName1Seq.pop_back();
        vecFileName2Seq.pop_back();
        return;
    }

    imageCount++;//图片计数
    //显示图片数量
    ui->imageNumLE->setText(QString::number(imageCount));

    if(imageCount == 1)//载入第一张图时获取图片信息,因为后面的图片长宽信息都一样，所以只要获取一次
    {
        imageSize.width = matImage.cols;
        imageSize.height = matImage.rows;
    }

    
    Mat matGrayImage;
    cvtColor(matImage, matGrayImage, CV_RGB2GRAY);
    //亚像素精确化
    cornerSubPix(matGrayImage, cornerPoints, Size(11, 11), Size(-1, -1), TermCriteria(CV_TERMCRIT_ITER + CV_TERMCRIT_EPS, 20, 0.01));
    cornerPointsSeq.push_back(cornerPoints);
    drawChessboardCorners(matGrayImage, cornerSize, cornerPoints, true);

    if(imageCount > 1)//当载入超过1张图片使标定按钮可用
        ui->startCalibBtn->setEnabled(true);

    vector<Point3f> chessPoints;//将棋盘世界坐标写入棋盘世界坐标序列；假设z = 0;xy在棋盘上，原点为提取棋盘格角点的第一个点
    for(int i = 0; i < cornerSize.height; i++)
        for(int j = 0; j < cornerSize.width; j++)
        {
            Point3f Point;
            Point.x = i * chessSize.width;
            Point.y = j * chessSize.height;
            Point.z = 0;
            chessPoints.push_back(Point);
        }
    chessWPointsSeq.push_back(chessPoints);



    //阈值分割,以计算激光直线方程
    Mat matImage1 = imread(vecFileName2Seq[vecFileName2Seq.size() - 1],IMREAD_GRAYSCALE);//读取带激光图片；
    Mat matThreImage;
    threshold(matImage1, matThreImage, 200, 255, THRESH_BINARY);//参数：源图像，目标图像，阈值，阈值最大值，阈值类型（这里是二进制阈值化）

    tImage = matImage1.clone();

    //提取形成直线的点
    vector<Point> laserPoints;
    for(int i = 0; i < matThreImage.rows; i++)
        for(int j = 0; j < matThreImage.cols; j++)
        {
            if(matThreImage.at<uchar>(i,j) == 255)
            {
                laserPoints.push_back(Point(j, i));
            }
        }
    Vec4f laserLinePara;
    fitLine(laserPoints, laserLinePara, DIST_L2, 0, 1e-1, 1e-2);//最小二乘法拟合直线方程
    cout << "LaserLinePara" << laserLinePara << endl;

    cv::Point point0;
    point0.x = laserLinePara[2];
    point0.y = laserLinePara[3];

    double k = laserLinePara[1] / laserLinePara[0];//计算斜率

    //计算直线的端点(y = k(x - x0) + y0)，及直线与上边缘的交点和下边缘的交点
    cv::Point point1, point2;
    point1.x = 0;
    point1.y = k * (0 - point0.x) + point0.y;
    point2.x = matImage1.rows;
    point2.y = k * (matImage.rows - point0.x) + point0.y;

    cv::line(matGrayImage, point1, point2, cv::Scalar(0, 255, 0), 2, 8, 0);//将直线画在灰度图上

    laserLineParaSeq.push_back(laserLinePara);//将方程放入方程序列，保存所有图片的直线方程

    //显示灰度图片在标签上
    QImage img;
    img = cvMat_to_QImage(matGrayImage);
    imageLabel->setPixmap(QPixmap::fromImage(img));
}

//设置输入框的输入有效性
void MainWindow::init()
{
    QDoubleValidator *validator=new QDoubleValidator(-99999999, 99999999, 10, this);

    ui->xCoeLE->setValidator(validator);
    ui->yCoeLE->setValidator(validator);
    ui->zCoeLE->setValidator(validator);
    ui->constCoeLE->setValidator(validator);

    validator->setRange(0, 999, 10);
    ui->chessLenLE->setValidator(validator);
}

//改变载入标定图片是否可用
void MainWindow::changeOpen1Slot()
{
    if(ui->rowNumSpB->value() > 2 && ui->colNumSpB->value() > 2 && ui->chessLenLE->text().toInt() > 0)
        ui->openBtn1->setEnabled(true);
    else
        ui->openBtn1->setEnabled(false);
}

//相机标定
void MainWindow::startCalibSlot()
{
    //参数：棋盘格世界坐标，角点像素坐标，图片尺寸，输出内参数矩阵，输出畸变矩阵，输出旋转向量，输出位移向量
    calibrateCamera(chessWPointsSeq, cornerPointsSeq, imageSize, cameraMatrix, distCoeffs, rvecsMatSeq, tvecsMatSeq, CV_CALIB_FIX_K3);
    QString info;
    info = "标定完成！！ 数据文件已输出.";
    ui->statusBar->showMessage(info);

    //保存存文件
    ofstream fout("caliberation_result.txt");
    fout << "相机内参数矩阵："  << endl;
    fout << cameraMatrix << endl << endl;
    fout << "畸变参数： \n";
    fout << distCoeffs << endl << endl << endl;
    Mat rotationMatrix = Mat(3, 3, CV_64F, Scalar::all(0));
    for(int i = 0; i < imageCount; i++)
    {

        fout << "第" << i + 1 << "幅图像的旋转向量" << endl;
        fout << rvecsMatSeq[i] << endl;

        //将旋转向量转换为相应的旋转矩阵
        Rodrigues(rvecsMatSeq[i], rotationMatrix);

        rotationMatrixSeq.push_back(rotationMatrix);
        fout << "第" << i + 1 << "幅图像的旋转矩阵："  << endl;
        fout << rotationMatrixSeq[i] << endl;
        cout << "RT" << rotationMatrixSeq[i] << endl;
        fout << "第" <<i + 1 << "幅图像的平移向量：" << endl;
        fout << tvecsMatSeq[i] << endl << endl;
        rotationMatrix = Mat(3, 3, CV_64F, Scalar::all(0));//这里必须重新设置为0矩阵，否则下一张图传给rotationMatrixSeq出现的错误。无法对应图片
    }

    //计算光平面方程
    Calculate();
}

//计算光平面方程
void MainWindow::Calculate()
{
    //拟合行角点形成的直线方程
    Vec4f imHLineParaT;
    vector<Vec4f> imHLineParaT2;
    vector<Point2f> points;
    for(int k = 0; k < imageCount; k++)//第k张图片
    {
        for(int i = 0; i < cornerSize.height; i++)//角点高
        {
            for(int j = 0; j < cornerSize.width; j++)//角点宽（通过高和宽（横坐标与纵坐标）确定角点位置）
                points.push_back(cornerPointsSeq[k][i * cornerSize.width + j]);//将角点放入points序列

            fitLine(points, imHLineParaT, DIST_L2, 0, 1e-2, 1e-2);//拟合横向直线方程

            cout << imHLineParaT << endl;//输出方程检查

            points.clear();//清楚数据，以免影响下一条直线方程的计算。
            imHLineParaT2.push_back(imHLineParaT);//每次计算完后都要先赋值,得到第k张图片的全部方程
        }
//        imHLineParaSeq.push_back(imHLineParaT2);//将所有的直线方程放入序列
        imHLineParaT2.clear();//赋值完后清除，否则下次会再次传入上次的值
    }

    //输出方程检查
    for(int k = 0; k < imageCount; k++)
        for(int i = 0; i < cornerSize.height; i++)
        {
            cout << imHLineParaSeq[k][i] << endl;
        }

    findPoints();//求结构光与棋盘格横向直线的交点坐标
    CalculatePoints();
    ui->autoCreateBtn->setEnabled(true);//使自动生成按钮可用


    cv::Point point0;
    point0.x = imHLineParaSeq[imageCount -1][3][2];
    point0.y = imHLineParaSeq[imageCount -1][3][3];

    double k = imHLineParaSeq[imageCount -1][3][1] / imHLineParaSeq[imageCount -1][3][0];

    //计算直线的端点(y = k(x - x0) + y0)
    cv::Point point1, point2;
    point1.x = 0;
    point1.y = k * (0 - point0.x) + point0.y;
    point2.x = tImage.cols;
    point2.y = k * (tImage.cols - point0.x) + point0.y;

    cv::line(tImage, point1, point2, cv::Scalar(0, 255, 0), 2, 8, 0);
    imshow("show", tImage);
}

void MainWindow::findPoints()//求结构光与棋盘格横向直线的交点坐标
{
    Point2f crossPoint,laserPoint0, imHPoint0;
    vector<Point2f> crossPointT;
    double laserk, imHk;

    for(int k = 0; k < imageCount; k++)//第k张图片
    {
        for(int i = 0; i < cornerSize.height; i++) //第k张图片的第i条横向直线
        {
            laserPoint0.x = laserLineParaSeq[k][2];
            laserPoint0.y = laserLineParaSeq[k][3];
            imHPoint0.x = imHLineParaSeq[k][i][2];
            imHPoint0.y = imHLineParaSeq[k][i][3];
            laserk = laserLineParaSeq[k][1] / laserLineParaSeq[k][0];
            imHk = imHLineParaSeq[k][i][1] / imHLineParaSeq[k][i][0];

            crossPoint.x = (-imHk * imHPoint0.x + imHPoint0.y + laserk * laserPoint0.x - laserPoint0.y) / (laserk - imHk);
            crossPoint.y = imHk * (crossPoint.x - imHPoint0.x) + imHPoint0.y;
            circle(tImage, crossPoint, 3, Scalar(0, 255, 0));
            crossPointT.push_back(crossPoint);//将第k张图片的的第i个交点放入crossPointT;
        }
        crossPointSeq.push_back(crossPointT)//将第k张图片的所有交点放入crossPointSeq;
    }
}

void MainWindow::CalculatePoints()//求实物相机坐标交点
{
    double leftx, lefty, ux, uy;//leftz,uz;//交比不变原理等式的左值和右上值
    Point2f A1, B1, C1, D1;
    Point3f A, B, C, D;
    double r1, r2, r4, r5, r7, r8;//r3,r6,r9;  //旋转矩阵值
    r1 = rotationMatrixSeq[0].at<double>(0,0);
    r2 = rotationMatrixSeq[0].at<double>(0,1);
    // r3 = rotationMatrixSeq[0].at<double>(0,2);//解方程时没用到
    r4 = rotationMatrixSeq[0].at<double>(1,0);
    r5 = rotationMatrixSeq[0].at<double>(1,1);
    //r6 = rotationMatrixSeq[0].at<double>(1,2);
    r7 = rotationMatrixSeq[0].at<double>(2,0);
    r8 = rotationMatrixSeq[0].at<double>(2,1);
    //r9 = rotationMatrixSeq[0].at<double>(2,2);
    double t1, t2,t3;//平移向量
    t1 = tvecsMatSeq[0].at<double>(0, 0);
    t2 = tvecsMatSeq[0].at<double>(1, 0);
    t3 = tvecsMatSeq[0].at<double>(2, 0);


    double Xw,Yw;//存储交点的世界坐标，以求得交点相机的Z轴坐标

    //生成世界坐标转换成摄像机坐标的外参矩阵Lw(即这里的T，所有图片的Lw都保存到TSeq）
    vector<Mat> TSeq;//
    Mat T;
    for(int n = 0; n < imageCount; n++)
    {
        T = Mat(4, 4, CV_64F,Scalar::all(0));
        for(int x = 0; x < 3; x++)
            for(int y = 0; y < 3; y++)
                T.at<double>(x, y) = rotationMatrixSeq[n].at<double>(x, y);
        T.at<double>(0, 3) = tvecsMatSeq[n].at<double>(0, 0);
        T.at<double>(1, 3) = tvecsMatSeq[n].at<double>(1, 0);
        T.at<double>(2, 3) = tvecsMatSeq[n].at<double>(2, 0);
        T.at<double>(3, 0) = 0;
        T.at<double>(3, 1) = 0;
        T.at<double>(3, 2) = 0;
        T.at<double>(3, 3) = 1;
        TSeq.push_back(T);
    }

    Mat world = Mat(4, 1, CV_64F, Scalar::all(0));//归一化的世界坐标向量，临时存储
    Mat world2Cam = Mat(4, 1, CV_64F, Scalar::all(0));//归一化地相机坐标向量,临时存储
    Point3f CamPoints;//保存计算好的相机坐标（实物）
    vector<Point3f> chessCPoints;//保存一张图的相机坐标（实物）

    Point3f point3fT;//临时存储某张棋盘格的世界坐标
    for(int n = 0; n < imageCount; n++)
    {
        for(int i = 0; i < cornerSize.width * cornerSize.height; i++)
        {
            point3fT = chessWPointsSeq[n][i];//提取第n张棋盘格图片的第i个角点的世界坐标
            world.at<double>(0, 0) = point3fT.x;
            world.at<double>(1, 0) = point3fT.y;
            world.at<double>(2, 0) = point3fT.z;
            world.at<double>(3, 0) = 1;

            world2Cam = TSeq[n] * world;//计算出该角点的摄像机坐标
            CamPoints.x = world2Cam.at<double>(0,0);
            CamPoints.y = world2Cam.at<double>(1,0);
            CamPoints.z = world2Cam.at<double>(2,0);
            chessCPoints.push_back(CamPoints);
        }
        chessCPointsSeq.push_back(chessCPoints);
        //cout << endl << "chessCPoints" << chessCPoints;
        chessCPoints.clear();//清除上一张图的数据
    }

    for(int n = 0; n < imageCount; n++)    //输出检查数据是否有误；发现程序关闭时才会输出第二张图片的坐标，是不是qt的bug
    {
        cout << endl << "CamePoint3f" << endl;
        cout << chessCPointsSeq[n];
    }
    cout << endl <<"end";    //再此输出个endl后发现上面的坐标全部输出，但end没有显示出来，可能是触发bug了。

    //根据交比不变原理求出结构光与棋盘格横向直线交点的摄像机坐标,即这里的D点
    for(int n = 0; n < imageCount; n++)
    {
        for(int i = 0; i < cornerSize.height; i++)
        {
            A1 = cornerPointsSeq[n][i * cornerSize.width];//A1,B1等为像素坐标点
            B1 = cornerPointsSeq[n][1 + i * cornerSize.width];
            C1 = cornerPointsSeq[n][i * cornerSize.width + cornerSize.width - 1];
            D1 = crossPointSeq[n][i];
            A = chessCPointsSeq[n][i * cornerSize.width];//ABC为图像世摄像机坐标点
            B = chessCPointsSeq[n][1 + i * cornerSize.width];
            C = chessCPointsSeq[n][i * cornerSize.width + cornerSize.width - 1];

            leftx = ((A1.x - C1.x) * (B1.x - D1.x)) / ((B1.x - C1.x) * (A1.x - D1.x));
            lefty = ((A1.y - C1.y) * (B1.y - D1.y)) / ((B1.y - C1.y) * (A1.y - D1.y));
            //leftz = leftx;
            ux = (A.x - C.x) / (B.x - C.x);
            uy = (A.y - C.y) / (B.y - C.y);
            //uz = (A.z - C.z) / (B.z - C.z);
            D.x = (leftx * A.x - ux * B.x) / (leftx - ux);
            D.y = (lefty * A.y - uy * B.y) / (lefty - uy);
            //D.z = (leftz * A.z - uz * B.y) / (leftz  -uz);

            Yw = (r4 * D.x - r1 * D.y  -r4 * t1 + r1 * t2) / (r2 * r4 - r1 * r5 );
            Xw = (D.x - r2 * Yw - t1) / r1;
            D.z = r7 * Xw + r8 * Yw + t3;//zz轴可以通过Xw和Yw求出，利用摄像机坐标变换相机坐标公式
            DCPointSeq.push_back(D);
            cout << endl << "D =";
            cout << D;
        }
    }

    //拟合光平面方程
    CvMat*points_mat = cvCreateMat(imageCount * cornerSize.height, 3 , CV_32FC1);//定义用来存储需要拟合点的矩阵
        for (int i = 0; i < imageCount * cornerSize.height; ++i)
        {
            points_mat->data.fl[i * 3 + 0] = DCPointSeq[i].x;//矩阵的值进行初始化   X的坐标值
            points_mat->data.fl[i * 3 + 1] = DCPointSeq[i].y;//  Y的坐标值
            points_mat->data.fl[i * 3 + 2] =DCPointSeq[i].z; //  Z的坐标值

        }
        float plane12[4] = { 0 };//
        cvFitPlane(points_mat, plane12);//调用方程

        cout << "A = "<< plane12[0] << endl;//直线方程Ax+By+Cz+D= 0;
        cout << "B =" << plane12[1] << endl;
        cout << "C =" << plane12[2] << endl;
        cout << "D =" << plane12[3] << endl;
        Acoe = plane12[0];
        Bcoe = plane12[1];
        Ccoe = plane12[2];
        Dcoe = plane12[3];
}

//网上查找的拟合光平面方程的函数，具体过程不理解。
void MainWindow::cvFitPlane(const CvMat* points, float* plane)
{
    // Estimate geometric centroid.
    int nrows = points->rows;
    int ncols = points->cols;
    int type = points->type;
    CvMat* centroid = cvCreateMat(1, ncols, type);
    cvSet(centroid, cvScalar(0));
    for (int c = 0; c<ncols; c++) {
        for (int r = 0; r < nrows; r++)
        {
            centroid->data.fl[c] += points->data.fl[ncols*r + c];
        }
        centroid->data.fl[c] /= nrows;
    }
    // Subtract geometric centroid from each point.
    CvMat* points2 = cvCreateMat(nrows, ncols, type);
    for (int r = 0; r<nrows; r++)
        for (int c = 0; c<ncols; c++)
            points2->data.fl[ncols*r + c] = points->data.fl[ncols*r + c] - centroid->data.fl[c];
    // Evaluate SVD of covariance matrix.
    CvMat* A = cvCreateMat(ncols, ncols, type);
    CvMat* W = cvCreateMat(ncols, ncols, type);
    CvMat* V = cvCreateMat(ncols, ncols, type);
    cvGEMM(points2, points, 1, NULL, 0, A, CV_GEMM_A_T);
    cvSVD(A, W, NULL, V, CV_SVD_V_T);
    // Assign plane coefficients by singular vector corresponding to smallest singular value.
    plane[ncols] = 0;
    for (int c = 0; c<ncols; c++) {
        plane[c] = V->data.fl[ncols*(ncols - 1) + c];
        plane[ncols] += plane[c] * centroid->data.fl[c];
    }
    // Release allocated resources.
    cvReleaseMat(&centroid);
    cvReleaseMat(&points2);
    cvReleaseMat(&A);
    cvReleaseMat(&W);
    cvReleaseMat(&V);
}

void MainWindow::AutoCreateSlots()
{
    ui->xCoeLE->setText(QString::number(Acoe));
    ui->yCoeLE->setText(QString::number(Bcoe));
    ui->zCoeLE->setText(QString::number(Ccoe));
    ui->constCoeLE->setText(QString::number(Dcoe));
}

void MainWindow::open2Slot()
{
    imageFileName = QFileDialog::getOpenFileName(this,"打开需测量的图片", QDir::currentPath(),  "iamge(*.jpg *.bmp *.png)");
    if(imageFileName.isEmpty())
    {
        QMessageBox::information(this,"ERROR","打开图片失败");
        return;
    }

    QImageReader reader(imageFileName);
    reader.setAutoTransform(true);//未理解
    image = reader.read();
    imageLabel->setPixmap(QPixmap::fromImage(image));
}

void MainWindow::changeAcoeSlot()
{
    if(!ui->xCoeLE->text().isEmpty() && !ui->yCoeLE->text().isEmpty() && !ui->zCoeLE->text().isEmpty() && !ui->constCoeLE->text().isEmpty())
    {
        ui->openBtn2->setEnabled(true);
    }
    if(ui->customSetCB->isChecked())
    {
       if(!ui->xCoeLE->text().isEmpty())
       {
           Acoe = ui->xCoeLE->text().toDouble();
       }
    }
}

void MainWindow::changeBcoeSlot()
{
    if(!ui->xCoeLE->text().isEmpty() && !ui->yCoeLE->text().isEmpty() && !ui->zCoeLE->text().isEmpty() && !ui->constCoeLE->text().isEmpty())
    {
        ui->openBtn2->setEnabled(true);
    }
    if(ui->customSetCB->isChecked())
    {
       if(!ui->yCoeLE->text().isEmpty())
       {
           Bcoe = ui->yCoeLE->text().toDouble();
       }
    }
}

void MainWindow::changeCcoeSlot()
{
    if(!ui->xCoeLE->text().isEmpty() && !ui->yCoeLE->text().isEmpty() && !ui->zCoeLE->text().isEmpty() && !ui->constCoeLE->text().isEmpty())
    {
        ui->openBtn2->setEnabled(true);
    }
    if(ui->customSetCB->isChecked())
    {
       if(!ui->zCoeLE->text().isEmpty())
       {
           Ccoe = ui->zCoeLE->text().toDouble();
       }
    }
}

void MainWindow::changeDcoeSlot()
{
    if(!ui->xCoeLE->text().isEmpty() && !ui->yCoeLE->text().isEmpty() && !ui->zCoeLE->text().isEmpty() && !ui->constCoeLE->text().isEmpty())
    {z
    }
    if(ui->customSetCB->isChecked())
    {
       if(!ui->constCoeLE->text().isEmpty())
       {
           Dcoe = ui->constCoeLE->text().toDouble();
       }
    }
}

void MainWindow::changeCustomSet()
{
    if(ui->customSetCB->isChecked())
    {
        ui->xCoeLE->setReadOnly(false);
        ui->yCoeLE->setReadOnly(false);
        ui->zCoeLE->setReadOnly(false);
        ui->constCoeLE->setReadOnly(false);
        ui->autoCreateBtn->setEnabled(false);
    }
    else
    {
        ui->xCoeLE->setReadOnly(true);
        ui->yCoeLE->setReadOnly(true);
        ui->zCoeLE->setReadOnly(true);
        ui->constCoeLE->setReadOnly(true);
        if(DCPointSeq.size() > 0)
            ui->autoCreateBtn->setEnabled(true);
    }
}

void MainWindow::markSlot()
{
    if(!image.isNull())
    {
        markDialog->setImage(image);
        markDialog->windowReSize(image.size());
        markDialog->show();
    }
    else
    {
        QMessageBox::information(this, "Error", "请先打开一张图片");
    }
}

void MainWindow::measureSlot()
{
    double u0,v0,fx,fy;
    fx = cameraMatrix.at<double>(0, 0);
    fy = cameraMatrix.at<double>(1, 1);
    u0 = cameraMatrix.at<double>(0, 2);
    v0 = cameraMatrix.at<double>(1, 2);
    markDialog->getPoint(startPoint, endPoint);//获取像素坐标
    cout << endl << "startPoint = " << startPoint;
    cout << endl << "endPoint = " << endPoint;
    double x1c, y1c, z1c, x2c, y2c, z2c;//转为摄像机坐标
    if(Acoe == 0 && Bcoe == 0 && Ccoe == 0 && Dcoe == 0)
        return;
    x1c = (-Dcoe * (startPoint.x - u0)) / (Acoe * (startPoint.x - u0) + Bcoe * (startPoint.y - v0) * (fx / fy) + Ccoe * fx);
    y1c = (-Dcoe * (startPoint.y - v0)) / (Acoe * (startPoint.x - u0) * (fy / fx) + Bcoe * (startPoint.y - v0) + Ccoe * fy);
    z1c = -Dcoe / (Acoe * (startPoint.x - u0) * (1 / fx) + Bcoe * (startPoint.y - v0) * (1 / fy) + Ccoe);
    x2c = (-Dcoe * (endPoint.x - u0)) / (Acoe * (endPoint.x - u0) + Bcoe * (endPoint.y - v0) * (fx / fy) + Ccoe * fx);
    y2c = (-Dcoe * (endPoint.y - v0)) / (Acoe * (endPoint.x - u0) * (fy / fx) + Bcoe * (endPoint.y - v0) + Ccoe * fy);
    z2c = -Dcoe / (Acoe * (endPoint.x - u0) * (1 / fx) + Bcoe * (endPoint.y - v0) * (1 / fy) + Ccoe);

    cout << endl << endl;
    cout << x1c << " " << y1c << " " << z1c;
    cout << endl << endl;
    cout << x2c << " " << y2c << " " << z2c;
    distance = sqrt(pow(x1c - x2c, 2) + pow(y1c - y2c, 2) + pow(z1c - z2c, 2));

    ui->distanceLE->setText(QString::number(distance));
}

void MainWindow::unitChangeSlot()
{
    if(ui->unitCB->currentIndex() == 0)
    {
        ui->distanceLE->setText(QString::number(distance));
    }
    if(ui->unitCB->currentIndex() == 1)
    {
        ui->distanceLE->setText(QString::number(distance / 10));
    }
    if(ui->unitCB->currentIndex() == 2)
    {
        ui->distanceLE->setText(QString::number(distance / 100));
    }
}
