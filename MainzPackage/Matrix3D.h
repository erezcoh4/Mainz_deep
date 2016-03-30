//								      -*-c++-*-
// Implementationen der Klassen Vector3D und Matix3D - Experimente mit C++
//
// V0.8, T.Pospischil, 20.05.98
//
// Copyright (c) 1998-2001
//
// Institut für Kernphysik, Universität Mainz	tel. +49 6131 39-25802
// 55099 Mainz, Germany				fax  +49 6131 39-22964
//
// $Id: Matrix3D.h 2216 2008-06-13 21:13:47Z distler $
//

#ifndef __MATRIX3D_H__
#define __MATRIX3D_H__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#ifdef HAVE_OSTREAM
#include <ostream>
#endif
#include <iosfwd>
#include <cmath>

const double deg = 3.141592653589793238 / 180.0;

/* ----------------------------- Klasse Vector3D --------------------------- */

class Vector3D{

public:
  
  double elem[3];               /* Internes Feld laeuft von 0 bis n-1 ! */
  
  /* Konstruktoren ------------ */
  inline Vector3D(){ elem[0] = 0.0; elem[1] = 0.0; elem[2] = 0.0; } 
  /* Copy-Konstruktor eines 3er Vektors */
  inline Vector3D(const Vector3D &X){
    elem[0] = X.elem[0]; elem[1] = X.elem[1]; elem[2] = X.elem[2];
  } 
  /* Konstruktor eines 3er Vektors mit 3 Werten */
  inline Vector3D(double x, double y, double z){
    elem[0] = x; elem[1] = y; elem[2] = z;
  } 

  /* Zugriff auf die Vektorelemente mit dem eckigen Klammerpaar [] */
  inline double& operator[](int i){ return elem[i-1]; }
  // Index laeuft hier von 1 bis 3 !!!
  
  /* Zuweisungs Operator */
  inline void operator=(const Vector3D &V){
    elem[0] = V.elem[0]; elem[1] = V.elem[1]; elem[2] = V.elem[2];
  }

  /* Vergleich zweier Vektoren: */
  inline int operator==(const Vector3D &a){
    if (elem[0] != a.elem[0]) return 0;
    if (elem[1] != a.elem[1]) return 0;
    if (elem[2] != a.elem[2]) return 0;
    return 1;
  }

  /* Summieren und Subtrahieren zweier gleichartiger Vektoren */
  inline Vector3D operator+(const Vector3D &a) const{ 
    return Vector3D( elem[0]+a.elem[0], elem[1]+a.elem[1], elem[2]+a.elem[2]);
  }
  inline Vector3D operator-(const Vector3D &a) const{ 
    return Vector3D( elem[0]-a.elem[0], elem[1]-a.elem[1], elem[2]-a.elem[2]);
  }
  /* Skalarprodukt zweier Vektoren */
  inline double operator*(const Vector3D &b) const{ 
    return ( elem[0]*b.elem[0] + elem[1]*b.elem[1] + elem[2]*b.elem[2] );
  }
  /* Kreuzprodukt zweier 3er Vektoren */
  inline Vector3D operator&&(const Vector3D &a) const{ 
    return Vector3D( (elem[1] * a.elem[2] - elem[2] * a.elem[1]),
		     (elem[2] * a.elem[0] - elem[0] * a.elem[2]),
		     (elem[0] * a.elem[1] - elem[1] * a.elem[0])  );
  }
  /* Multiplizieren und Dividieren eines Vektors mit einer Zahl */
  inline Vector3D operator*(const double f) const{ 
    return Vector3D( f*elem[0], f*elem[1], f*elem[2] );
  }
  inline Vector3D operator/(const double f) const{ 
    return Vector3D( elem[0]/f, elem[1]/f, elem[2]/f );
  }
  
  /* Betrag eines Vektors */
  inline double abs() const{
    return sqrt( elem[0]*elem[0] + elem[1]*elem[1] + elem[2]*elem[2] );
  }

};

// Produkt mit einer Zahl:
inline Vector3D operator*(const double f, const Vector3D &b){ 
  return Vector3D( f*b.elem[0], f*b.elem[1], f*b.elem[2] );
}
/* Winkel zwischen zwei Vektoren in rad */
inline double ang(const Vector3D &a, const Vector3D &b){ 
  double hilf=0;
  if ( (b.abs()*a.abs()) == 0 ) return 0;
  hilf = (a*b)/(a.abs()*b.abs());
  if (hilf>1) hilf = 1.0;
  return acos(hilf);     // Drehwinkel 
}
/* Ausgeben eines Vektors */
inline std::ostream& operator<<(std::ostream &s, Vector3D O){
  s << "( "; for(int i=1;i<=3;i++) s << O[i] << ' '; s << ')';
  return s;
}



/* ----------------------------- Klasse Matix3D ----------------------------- */

class Matrix3D{

 public:

  Vector3D m[3];

  /* Konstruktor einer Matrix3D */
  inline Matrix3D(){} 
  /* Copy-Konstruktor einer Matrix3D */
  inline Matrix3D(const Matrix3D & W){
    m[0] = W.m[0]; m[1] = W.m[1]; m[2] = W.m[2];
  } 
  /* Konstruktor einer Matrix3D als Drehmatrix um Drehachse D und Winkel phi */
  /* ---- !!!! ---- Drehachse muss auf 1 normiert sein ----- !!!! ----- */
  inline Matrix3D(const Vector3D &D, const double phi){
    m[0].elem[0] = D.elem[0]*D.elem[0]*(1-cos(phi))+cos(phi);
    m[0].elem[1] = D.elem[0]*D.elem[1]*(1-cos(phi))-D.elem[2]*sin(phi);
    m[0].elem[2] = D.elem[0]*D.elem[2]*(1-cos(phi))+D.elem[1]*sin(phi);
    m[1].elem[0] = D.elem[1]*D.elem[0]*(1-cos(phi))+D.elem[2]*sin(phi);
    m[1].elem[1] = D.elem[1]*D.elem[1]*(1-cos(phi))+cos(phi);
    m[1].elem[2] = D.elem[1]*D.elem[2]*(1-cos(phi))-D.elem[0]*sin(phi);
    m[2].elem[0] = D.elem[2]*D.elem[0]*(1-cos(phi))-D.elem[1]*sin(phi);
    m[2].elem[1] = D.elem[2]*D.elem[1]*(1-cos(phi))+D.elem[0]*sin(phi);
    m[2].elem[2] = D.elem[2]*D.elem[2]*(1-cos(phi))+cos(phi);
  } 
  /* Konstruktor einer Matrix3D aus 3 (Zeilen-)Vektoren */
  inline Matrix3D(const Vector3D &x, const Vector3D &y, const Vector3D &z){
    m[0] = x; m[1] = y; m[2] = z;
  } 

  /* Zugriff auf die Matrixelemente mit zwei eckigen Klammerpaaren [][] */
  // Index laeuft von 1 bis 3 !!!
  inline Vector3D& operator[](int i){ return m[i-1]; }

  /* Zuweisungs Operator */
  inline void operator=(const Matrix3D & M){
    m[0] = M.m[0]; m[1] = M.m[1]; m[2] = M.m[2]; 
  }

  /* Summieren und Subtrahieren zweier gleichartiger Matrizen */
  inline Matrix3D operator+(const Matrix3D &a) const{ 
    return Matrix3D( m[0]+a.m[0], m[1]+a.m[1], m[2]+a.m[2] );
  }
  inline Matrix3D operator-(const Matrix3D &a) const{ 
    return Matrix3D( m[0]-a.m[0], m[1]-a.m[1], m[2]-a.m[2] );
  }

  /* Multiplizieren und Dividieren einer Matrix3D mit einer Zahl */
  inline Matrix3D operator*(const double f) const {
    return Matrix3D( m[0]*f, m[1]*f, m[2]*f );
  }
  inline Matrix3D operator/(const double f) const {
    return Matrix3D( m[0]/f, m[1]/f, m[2]/f );
  }

  /* Transponieren einer Matrix3D */
  Matrix3D operator!();

};

/* Multiplizieren einer Zahl mit einer Matrix3D */
inline Matrix3D operator*(const double f, const Matrix3D &b){ 
  return Matrix3D( b.m[0]*f, b.m[1]*f, b.m[2]*f );
}

/* Multiplizieren einer Matrix3D mit einem Vector3D */
// Der Vector wird hier als Spaltenvector behandelt
inline Vector3D operator*(const Matrix3D &a, const Vector3D &u){
  return Vector3D( a.m[0]*u, a.m[1]*u, a.m[2]*u );
}
  
std::ostream& operator<<(std::ostream& s, Matrix3D O);      // Ausgeben einer Matrix3D


/* Drehen eines Vektors V0 um die Drehachse, die durch das
   Kreuzprodukt zwischen V1 und V2 definiert wird, um den Winkel zwischen
   V1 und V2.
   Rueckgabewert der Funktion ist der gedrehte Vektor, ausserdem werden
   der Drehwinkel phi und die normierte Drehachse D zurueckgegeben
   */
Vector3D
rotVec( Vector3D V0, // Vector der gedreht wird
	Vector3D V1, // Startrichtung
	Vector3D V2, // Zielrichtung
	Vector3D &D, // Drehachse
	double &phi);// Drehwinkel

#endif
