#ifndef __VECTOR_3D_H
#define __VECTOR_3D_H
/**
 * 3D Vector helper class.
 * Author : Juan Pablo Vanegas. Git: jpvanegasc
 */
#include <iostream>
#include <cmath>


class Vector3D; // First declaration to avoid problems in namespace

Vector3D operator* (double a, Vector3D v1);

namespace vector{
    Vector3D unit_vector(Vector3D v);
    double norm2(Vector3D v1);
    double norm(Vector3D v1);
}

class Vector3D{
    double v[3];
    public:
        Vector3D(double =0.0, double =0.0, double =0.0);
        void load(double x0, double y0, double z0);
        void show(void);
        /* @return x component*/
        double x(void){return v[0];};
        /* @return y component*/
        double y(void){return v[1];};
        /* @return z component*/
        double z(void){return v[2];};
        /* @return nth component*/
        double & operator[](int n){return v[n];};
        Vector3D operator= (Vector3D v2);
        Vector3D operator+ (Vector3D v2);
        Vector3D operator+=(Vector3D v2);
        Vector3D operator- (Vector3D v2);
        Vector3D operator-=(Vector3D v2);
        Vector3D operator* (double a);
        Vector3D operator*=(double a);
        Vector3D operator/ (double a);
        double operator* (Vector3D v2);
        Vector3D operator^ (Vector3D v2);
        friend Vector3D operator* (double a, Vector3D v1);
        friend Vector3D vector::unit_vector(Vector3D v);
        friend double vector::norm2(Vector3D v1);
        friend double vector::norm(Vector3D v1);
};

#endif
