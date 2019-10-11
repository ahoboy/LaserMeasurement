#include "dialog.h"
#include "ui_dialog.h"

Dialog::Dialog(QWidget *parent) :
    QDialog(parent),
    imagePosx(0),
    imagePosy(0),
    imageCurx(0),
    imageCury(0),
    curAtWindowx(0),
    curAtWindowy(0),
    pointCount(0),
    ui(new Ui::Dialog)
{
    ui->setupUi(this);

    //resize(QGuiApplication::primaryScreen()->availableSize());//初始显示大小
}

Dialog::~Dialog()
{
    delete ui;
}

void Dialog::setImage(QImage image1)
{
    image = image1.copy();
    update();
}

void Dialog::windowReSize(QSize size)
{
    this->resize(size);
}

void Dialog::mousePressEvent(QMouseEvent *event)
{
    curAtWindowx = event->x();
    curAtWindowy = event->y();

    imageCurx = (curAtWindowx-imagePosx);
    imageCury = (curAtWindowy-imagePosy);

    if(event->button()==Qt::LeftButton && pointCount == 0)
    {
        startPoint.x = imageCurx;
        startPoint.y = imageCury;
    }
    if(event->button()==Qt::LeftButton && pointCount == 1)
    {
        endPoint.x = imageCurx;
        endPoint.y = imageCury;
    }
    mousePress = true;
}

void Dialog::mouseMoveEvent(QMouseEvent *event)
{
    curAtWindowx = event->x();
    curAtWindowy = event->y();

    imageCurx = (curAtWindowx-imagePosx);
    imageCury = (curAtWindowy-imagePosy);
    QToolTip::showText(QPoint(curAtWindowx + this->x(),curAtWindowy + this->y()),"X:"+QString::number(imageCurx)+" Y:"+QString::number(imageCury),this,this->rect(),2000);

}

void Dialog::paintEvent(QPaintEvent *event)//所有的绘制事件必须要在这里执行
{
   QPainter painter(this);

   if(image.isNull())
   {
      return;
   }

  imagePosx=this->width()/2-image.width()/2;
  imagePosy=this->height()/2-image.height()/2;
  painter.drawImage(QPoint(imagePosx,imagePosy),image);

  if(mousePress == true && pointCount == 0)
  {

      painter.drawPoint(startPoint.x, startPoint.y);
      painter.drawText(QRect(this->rect().x() + startPoint.x, this->rect().y() + startPoint.y, 10, 10), "1");
  }

  if(mousePress == true && pointCount == 1)
  {
      painter.drawPoint(endPoint.x, endPoint.y);
      painter.drawText(QRect(this->rect().x() + startPoint.x, this->rect().y() + startPoint.y, 10, 10), "2");
  }
}

void Dialog::mouseReleaseEvent(QMouseEvent *event)
{
    mousePress = false;
    pointCount++;
    if(pointCount > 1)
    {
        pointCount = 0;
        this->close();
    }
}

void Dialog::getPoint(Point &startPointT, Point &endPointT)
{
    startPointT = startPoint;
    endPointT = endPoint;
}
