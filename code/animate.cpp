/**
 * Animation module.
 * For calculating and saving data please use main.cpp
 */
#include"molecular_dynamics.h"
#include"random64.h"
#include"animation.hpp"

int main(int argc, char *argv[]){
    Body Molecule[N];
    Collider Newton;
    CRandom ran64(1);
    double t,tdibujo;

    start_animation(argc);

    double k_T=10000;
    double x0, y0, theta, Vx0, Vy0;
    double V0 = sqrt(2.0*k_T/m0);

    for(int i=0; i<Nx; i++)
        for(int j=0; j<Ny; j++){
            x0=(i+1)*dx; y0=(j+1)*dy; 
            theta = 2.0*M_PI*ran64.r(); 
            Vx0 = V0*std::cos(theta); Vy0 = V0*std::sin(theta);
            if (i==0 && j==0){
                Molecule[i*Ny+j].initialize(x0, y0, Vx0, Vy0, m0*10, R0*2.5);
                continue;
            }
            Molecule[i*Ny+j].initialize(x0, y0, Vx0, Vy0, m0, R0);
        }

    double T_sim=8.5,T_eq=1.5;
        
    for(t=tdibujo=0.0; t<T_eq+T_sim; t+=dt,tdibujo+=dt){
        if(tdibujo>0.05 && t<T_eq){
            begin_frame(argc);
            for(int i=0; i<N; i++) Molecule[i].print();
            end_frame(argc);
            tdibujo=0.0;
        }

        Newton.move_with_pefrl(Molecule, dt);
    }
    
    return 0;
}