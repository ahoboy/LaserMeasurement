#ifndef QIMAGEMAT_H
#define QIMAGEMAT_H
#include "QImage"
#include "opencv.hpp"

QImage cvMat_to_QImage(const cv::Mat &mat );
cv::Mat QImage_to_cvMat( const QImage &image, bool inCloneImageData = true );

#endif // QIMAGEMAT_H
