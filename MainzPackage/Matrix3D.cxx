//								      -*-c++-*-
// Implementationen der Klassen Vector3D und Matix3D - Experimente mit C++
//
// V0.5, T.Pospischil, 23.01.98
//
// Copyright (c) 1998-2001
//
// Institut für Kernphysik, Universität Mainz	tel. +49 6131 39-25802
// 55099 Mainz, Germany				fax  +49 6131 39-22964
//
// $Id: Matrix3D.cc 2216 2008-06-13 21:13:47Z distler $
//

#include "Matrix3D.h"
#include <cstdlib>

/* ----------------------------- Klasse Matrix3D -----------------------------*/

/* Transponieren einer Matrix3D */
Matrix3D
Matrix3D::operator!(){
  Matrix3D c;
  for (int i=1; i<=3; i++) for (int j=1; j<=3; j++) c[i][j] = m[j-1][i]; 
  return c;
}

/* Ausgeben einer Matrix3D */
std::ostream& operator<<(std::ostream& s, Matrix3D O){
  for(int i=1;i<=3;i++) s << O[i] << std::endl; 
  return s;
}

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
	double &phi) // Drehwinkel
{
  phi = ang(V1,V2);     // Drehwinkel

  // Start- und Zielrichung identisch => es muss nichts gedreht werden
  if (phi == 0){
    D.elem[0] = 0; D.elem[1]=0; D.elem[2]=0; // Drehachse undefiniert
    return V0; 
  }

  // Start- und Zielrichung entgegengesetzt:
  if ((phi <= -180.0) || (phi >= 180.0)){
    phi = 180.0;
    D.elem[0] = 0; D.elem[1]=0; D.elem[2]=0;  // Drehachse undefiniert
    return (-1.0*V0); 
  }

  // Berechnung der Drehachse:
  D = V1 && V2;
  if (D.abs() > 0) {
    D = D/D.abs();                               // Normierung der Drehachse
    
    // Drehmatrix:
    Matrix3D M(D,phi);
    
    // Rueckgabe des gedrehten Vektors:
    return (M * V0);
  }
  else {
    std::cerr << " ERROR in rotVec: length of rot axis eqals zero ???"
	      << std::endl;
    return V0;
  }
}
