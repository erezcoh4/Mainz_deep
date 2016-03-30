//								      -*-c++-*-
// COLA:
//
// Copyright (c) 2001
//
// Institut für Kernphysik, Universität Mainz	tel. +49 6131 39-25802
// 55099 Mainz, Germany				fax  +49 6131 39-22964
//
// $Id: qspinSpectrometer.h 2216 2008-06-13 21:13:47Z distler $
//

#ifndef __QSPIN_SPECTROMETER_H__
#define __QSPIN_SPECTROMETER_H__

#include "Matrix3D.h"
#include "Magnets.h"
#include <cstring>

class qspinSpectrometer {
protected:
  char *base;
public:
  qspinSpectrometer(const char *basename) {
    base = new char[strlen(basename)+1];
    strcpy(base, basename);
  }
  virtual ~qspinSpectrometer() {
    delete base;
  }
  const char *getBasename() { return base; }
  /// Calculates the magnetic field (vector3D B) in the absolute Coordinate
  /// system (ACS, x downwards) at a given position (vector3D r):
  virtual int getMagField(Vector3D &B, Vector3D &r) = 0;
  virtual int showBoundary() = 0;
  virtual int checkBoundary() = 0;
  virtual int scaleField(double p_ref) = 0;
};

class qspinSpecA : public qspinSpectrometer {
protected:
  Multipol *Quad, *Sext;
  Dipol *D1, *D2;
public:
  qspinSpecA(const char *basename="qspin");
  virtual ~qspinSpecA() {
    delete Quad;
    delete Sext;
    delete D1;
    delete D2;
  }
  virtual int getMagField(Vector3D &B, Vector3D &r);
  virtual int showBoundary();
  virtual int checkBoundary();
  virtual int scaleField(double p_ref) {
    // Set spectrometer magnetic fields:
    double errect = p_ref/630.0;  // Scaling value for the max. field
    Quad->ScaleField(errect);
    Sext->ScaleField(errect);
    D1->ScaleField(errect);
    D2->ScaleField(errect);
    return 0;
  }
};

class qspinOOPS : public qspinSpectrometer {
protected:
  Dipol *Dip;
  Multipol *Quad;
public:
  qspinOOPS(const char *basename="qspinOOPS");
  virtual ~qspinOOPS() {
    delete Dip;
    delete Quad;
  }
  virtual int getMagField(Vector3D &B, Vector3D &r);
  virtual int showBoundary();
  virtual int checkBoundary();
  virtual int scaleField(double p_ref) {
    // Set spectrometer magnetic fields:
    double errect = p_ref/625.0;  // Scaling value for the max. field
    Quad->ScaleField(errect);
    Dip->ScaleField(errect);
    return 0;
  }
};
  
class qspinOHIPS : public qspinSpectrometer {
protected:
  Multipol *Q1, *Q2;
  Dipol *D1;
public:
  qspinOHIPS(const char *basename="qspinOHIPS");
  virtual ~qspinOHIPS() {
    delete Q1;
    delete Q2;
    delete D1;
  }
  virtual int getMagField(Vector3D &B, Vector3D &r);
  virtual int showBoundary();
  virtual int checkBoundary();
  virtual int scaleField(double p_ref) {
    // Set spectrometer magnetic fields:
    double errect = p_ref/269.1;  // Scaling value for the max. field
    Q1->ScaleField(errect);
    Q2->ScaleField(errect);
    D1->ScaleField(errect);
    return 0;
  }
};
  
#endif /* __QSPIN_SPECTROMETER_H__ */
