#include "molecular_dynamics.h"

void Body::initialize(double x0, double y0, double Vx0, double Vy0, double m0, double R0){
    r.load(x0,y0,0); V.load(Vx0,Vy0,0); m = m0; R = R0;
}
/**
 * Adds the contact force between two molecules, and if molecule is in contact with wall, 
 * calculates and adds it too. Includes diffusion. 
 * Force with wall is calculated using Hertz contact force.
 * 
 * @param dF: Force with other molecule. See Collider.calculate_force_pair for details.
 */
void Body::add_force(Vector3D dF){
    double x=r.x(), y=r.y();
    bool far_from_wall = ((R<x) && (x<(Lx-R))) && ((R<y) && (y<(Ly-R)));
    
    if(far_from_wall)
        F += dF;
    else{
        double Fx=0, Fy=0;
        double hx1 = R - x;
        double hy1 = R - y;
        
        if(hx1 > 0) Fx = 1e4*std::pow(hx1, 1.5) - GAMMA*m*V.x()*std::pow(hx1, 0.5);
        else if(hx1 < 0) {
            double hx2 = R - Lx + x;
            if(hx2 > 0) Fx = -1e4*std::pow(hx2, 1.5) - GAMMA*m*V.x()*std::pow(hx2, 0.5);
        }

        if(hy1 > 0) Fy = 1e4*std::pow(hy1, 1.5) - GAMMA*m*V.x()*std::pow(hy1, 0.5);
        else if(hy1 < 0) {
            double hy2 = R - Ly + y;
            if(hy2 > 0) Fy = -1e4*std::pow(hy2, 1.5) - GAMMA*m*V.x()*std::pow(hy2, 0.5);
        }
        
        Vector3D aux(Fx, Fy, 0);
        F += dF + aux;
    }
}
void Body::print(void){
    std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}

/**
 * Calculates contact force between a pair of molecules. 
 * Force between molecules is calculated used Hertz contact force. Includes diffusion. 
 * 
 * @params molecule1, molecule2 : Pair of molecules which force is to be calculated.
 */
void Collider::calculate_force_pair(Body &molecule1, Body &molecule2){
    Vector3D dr, dF;
    dr = molecule2.r - molecule1.r;
    
    double s = (molecule1.R + molecule2.R) - norm(dr);
    if(s>0){
        double mean_speed = (norm(molecule1.V) + norm(molecule2.V))/2;
        double F_aux = 1.0e4*std::pow(s, 1.5) - GAMMA*molecule1.m*mean_speed*std::pow(s, 0.5);
        dF = dr*F_aux;
    }

    molecule2.add_force(dF); molecule1.add_force((-1.0)*dF);
}

void Collider::calculate_all_forces(Body *molecule){
    for(int i=0; i<N; i++) molecule[i].delete_f();
    
    for(int i=0; i<N; i++)
        for(int j=i+1; j<N; j++)
            calculate_force_pair(molecule[i], molecule[j]);
}

void Collider::move_with_pefrl(Body *molecule, double dt){
    int i;
    for(i=0;i<N;i++) molecule[i].move_r(dt,Zi);
    calculate_all_forces(molecule);
    for(i=0;i<N;i++) molecule[i].move_v(dt,coef1);
    for(i=0;i<N;i++) molecule[i].move_r(dt,Xi);
    calculate_all_forces(molecule);
    for(i=0;i<N;i++) molecule[i].move_v(dt,Lambda);
    for(i=0;i<N;i++) molecule[i].move_r(dt,coef2);
    calculate_all_forces(molecule);
    for(i=0;i<N;i++) molecule[i].move_v(dt,Lambda);
    for(i=0;i<N;i++) molecule[i].move_r(dt,Xi);
    calculate_all_forces(molecule);
    for(i=0;i<N;i++) molecule[i].move_v(dt,coef1);
    for(i=0;i<N;i++) molecule[i].move_r(dt,Zi);
}