#-------------------------------------------------
#
# Project created by QtCreator 2018-07-12T17:24:08
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = LaserMeasurement
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += \
        main.cpp \
        mainwindow.cpp \
    mprocess.cpp \
    qimagemat.cpp \
    dialog.cpp


HEADERS += \
        mainwindow.h \
    qimagemat.h \
    mprocess.h \
    dialog.h


FORMS += \
        mainwindow.ui \
    dialog.ui

INCLUDEPATH += D:/Programs/OpenCV/ForMingw/include\
                D:/Programs/OpenCV/ForMingw/include/opencv\
                D:/Programs/OpenCV/ForMingw/include/opencv2

LIBS +=         D:/Programs/OpenCV/ForMingw/lib/libopencv_calib3d340.dll.a\
                    D:/Programs/OpenCV/ForMingw/lib/libopencv_core340.dll.a\
                    D:/Programs/OpenCV/ForMingw/lib/libopencv_dnn340.dll.a\
                    D:/Programs/OpenCV/ForMingw/lib/libopencv_features2d340.dll.a\
                    D:/Programs/OpenCV/ForMingw/lib/libopencv_flann340.dll.a\
                    D:/Programs/OpenCV/ForMingw/lib/libopencv_highgui340.dll.a\
                    D:/Programs/OpenCV/ForMingw/lib/libopencv_imgcodecs340.dll.a\
                    D:/Programs/OpenCV/ForMingw/lib/libopencv_imgproc340.dll.a\
                    D:/Programs/OpenCV/ForMingw/lib/libopencv_ml340.dll.a\
                    D:/Programs/OpenCV/ForMingw/lib/libopencv_objdetect340.dll.a\
                    D:/Programs/OpenCV/ForMingw/lib/libopencv_photo340.dll.a\
                    D:/Programs/OpenCV/ForMingw/lib/libopencv_shape340.dll.a\
                    D:/Programs/OpenCV/ForMingw/lib/libopencv_stitching340.dll.a\
                    D:/Programs/OpenCV/ForMingw/lib/libopencv_superres340.dll.a\
                    D:/Programs/OpenCV/ForMingw/lib/libopencv_video340.dll.a\
                    D:/Programs/OpenCV/ForMingw/lib/libopencv_videoio340.dll.a\
                    D:/Programs/OpenCV/ForMingw/lib/libopencv_videostab340.dll.a

