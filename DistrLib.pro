QT -= gui

TEMPLATE = lib
CONFIG += staticlib

CONFIG += c++11

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    abstractdistribution.cpp \
    distributions/abstractdistribution.cpp \
    distributions/doublenormaldistr.cpp \
    distributions/empiricaldistr.cpp \
    distributions/extremevaluedistr.cpp \
    distributions/gammadistr.cpp \
    distributions/logisticdistr.cpp \
    distributions/lognormaldistr.cpp \
    distributions/mathematicadistr.cpp \
    distributions/moyaldistr.cpp \
    distributions/normaldistr.cpp \
    distributions/weibulldistr.cpp \
    distrlib.cpp \
    functions.cpp \
    random/mtwist.c \
    random/randistrs.c \
    rng.cpp

HEADERS += \
    distributions.h \
    distributions/abstractdistribution.h \
    distributions/doublenormaldistr.h \
    distributions/empiricaldistr.h \
    distributions/extremevaluedistr.h \
    distributions/gammadistr.h \
    distributions/logisticdistr.h \
    distributions/lognormaldistr.h \
    distributions/mathematicadistr.h \
    distributions/moyaldistr.h \
    distributions/normaldistr.h \
    distributions/weibulldistr.h \
    distrlib.h \
    functions.h \
    random/mtwist.h \
    random/randistrs.h \
    rng.h

# Default rules for deployment.
unix {
    target.path = $$[QT_INSTALL_PLUGINS]/generic
}
!isEmpty(target.path): INSTALLS += target


win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../build-alglib/release/ -lalglib
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../build-alglib/debug/ -lalglib
else:unix: LIBS += -L$$PWD/../build-alglib/ -lalglib

INCLUDEPATH += $$PWD/../alglib
DEPENDPATH += $$PWD/../alglib
