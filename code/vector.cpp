#include"vector.h"


/* Initialize vector. Defaults to zero */
Vector3D::Vector3D(double x, double y, double z){
    load(x,y,z);
}
/* Load vector values */
void Vector3D::load(double x0, double y0, double z0){
    v[0] = x0; v[1] = y0; v[2] = z0;
}
void Vector3D::show(void){
    std::cout << "(" << v[0]<< "," << v[1]<< "," << v[2]<< ")" << std::endl;
}
Vector3D Vector3D::operator=(Vector3D v2){
    v[0] = v2.v[0];
    v[1] = v2.v[1];
    v[2] = v2.v[2];
    
    return *this;
}
Vector3D Vector3D::operator+(Vector3D v2){
    Vector3D result;
    result.v[0] = v[0] + v2.v[0];
    result.v[1] = v[1] + v2.v[1];
    result.v[2] = v[2] + v2.v[2];
    
    return result;
}
Vector3D Vector3D::operator+=(Vector3D v2){
    *this = *this + v2;
    return *this;
}
Vector3D Vector3D::operator-(Vector3D v2){
    Vector3D result;
    result.v[0] = v[0] - v2.v[0];
    result.v[1] = v[1] - v2.v[1];
    result.v[2] = v[2] - v2.v[2];
    
    return result;
}
Vector3D Vector3D::operator-=(Vector3D v2){
    *this = *this - v2;
    return *this;
}
/* Vector times scalar */
Vector3D Vector3D::operator*(double a){
    Vector3D result;
    result.v[0] = v[0]*a;
    result.v[1] = v[1]*a;
    result.v[2] = v[2]*a;
    
    return result;
}
/* Vector times scalar */
Vector3D Vector3D::operator*=(double a){
    *this = (*this)*a;
    return *this;
}
/* Vector divided by scalar */
Vector3D Vector3D::operator/(double a){
    Vector3D result; double inverse = 1.0/a;
    result.v[0] = v[0]*inverse;
    result.v[1] = v[1]*inverse;
    result.v[2] = v[2]*inverse;
    
    return result;
}
/* Dot product */
double Vector3D::operator*(Vector3D v2){
    return v[0]*v2.v[0] + v[1]*v2.v[1] +v [2]*v2.v[2];;
}
/* Cross product */
Vector3D Vector3D::operator^(Vector3D v2){
    Vector3D result;
    result.v[0] = v[1]*v2.v[2]-v[2]*v2.v[1];
    result.v[1] = v[2]*v2.v[0]-v[0]*v2.v[2];
    result.v[2] = v[0]*v2.v[1]-v[1]*v2.v[0];
    
    return result;
}


Vector3D operator*(double a, Vector3D v1){
    Vector3D result;
    result = v1*a;
    return result;
}

namespace vector{
    /* Returns a unit vector with the same direction as v */
    Vector3D unit_vector(Vector3D v){
        return v/vector::norm(v);
    }
    /* @return vector norm squared */
    double norm2(Vector3D v1){
        return v1.v[0]*v1.v[0] + v1.v[1]*v1.v[1] + v1.v[2]*v1.v[2];
    }
    /* @return vector norm */
    double norm(Vector3D v1){
        return std::sqrt(vector::norm2(v1));
    }
}