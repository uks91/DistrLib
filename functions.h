#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <QtCore>

//typedef long double qreal;

typedef QVector<qreal> array1;
typedef QVector<array1 > array2;
typedef qreal (*function)(qreal, array1);


void arrayToFile (const QString &fileName, const array1 &in);
void arrayToFile (const QString &fileName, const array2 &in);
array2 CDF(array1 series);
//array1 fileToArray (const QSttring &fileName);
array1 getMirroredData (const array1 &in, qreal margin=0, bool isMirrored=true, bool is4=false);
qreal getRoot(function func, array1 pars, qreal x0, qreal res=0);
array1 loadData (const QString &filename, int maxSize=-1);
QString number(qreal n, int d);
qreal mround (qreal x, qreal signs=3);

#endif // FUNCTIONS_H
