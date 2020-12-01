#ifndef __N_BODY_PROBLEM_MOLECULAR_DYNAMICS_H
#define __N_BODY_PROBLEM_MOLECULAR_DYNAMICS_H

#include<iostream>
#include<cmath>
#include<fstream>

#include"vector.h"

// Geometry constants
const double square = 50.0;
const double Lx = square, Ly = square;

// Diffusion
const double GAMMA = 0.0;

// PEFRL
const double Zi = 0.1786178958448091e0;
const double Lambda = 0.2123418310626054*(-1);
const double Xi = 0.06626458266981849*(-1);

const double coef1 = (1 - 2*Lambda)/2;
const double coef2 = (1 - 2*(Xi+Zi));

// Implementation
const int Nx = 10, Ny = 10;
const int N = Nx*Ny;
const double R0 = 1.0, m0 = 1.0;
const double dt = 1.0e-3;
const double dx = Lx/(Nx+1), dy = Ly/(Ny+1);

class Body{
    private:
        Vector3D r, V, F; double m, R;
    public:
        void initialize(double x0, double y0, double Vx0, double Vy0, double m, double R);
        void add_force(Vector3D dF);
        void move_r(double dt, double coef){r += V*(dt*coef);}
        void move_v(double dt, double coef){V += F*(dt*coef/m);}
        void print(void);
        void delete_f(void){F.load(0,0,0);}
        double get_x(void){return r.x();} 
        double get_y(void){return r.y();} 
        double get_z(void){return r.z();} 
        double get_vx(void){return V.x();} 
        double get_vy(void){return V.y();} 
        double get_fx(void){return F.x();} 
        double get_fy(void){return F.y();} 

        friend class Collider;
};

class Collider{
    public:
        void calculate_force_pair(Body &molecule1, Body &molecule2);
        void calculate_all_forces(Body *molecule);
        void move_with_pefrl(Body *molecule, double dt);
};

#endif
