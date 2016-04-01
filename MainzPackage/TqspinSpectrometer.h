/**
 * \file TqspinSpectrometer.h
 *
 * \ingroup MainzPackage
 *
 * \brief Class def header for a class TqspinSpectrometer
 *
 * @author erezcohen
 */

/** \addtogroup MainzPackage
 
 @{*/
#ifndef TQSPINSPECTROMETER_H
#define TQSPINSPECTROMETER_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "Matrix3D.h"
#include "Magnets.h"
#include <cstring>


class TqspinSpectrometer {
protected:
    char *base;
public:
    
    TqspinSpectrometer(const char *basename) {
        base = new char[strlen(basename)+1];
        strcpy(base, basename);
    }
    virtual ~TqspinSpectrometer() {
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



class qspinSpecA : public TqspinSpectrometer {
protected:
    Multipol *Quad, *Sext;
    Dipol *D1, *D2;
    
public:
    qspinSpecA(const char *basename="qspin");
    ~qspinSpecA() { // virtual
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



#endif
/** @} */ // end of doxygen group

