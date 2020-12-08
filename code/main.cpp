/**
 * Main module. 
 * This program is for calculating and saving data. For animations please use animate.cpp
 */
#include<fstream>

#include"molecular_dynamics.h"



int main(int argc, char *argv[]){
    // Initial conditions
    int N, c; double **initial = NULL;
    load_file("earth_sun.csv", initial, N, c);

    double t_min = 0, t_max = 1e3; // days
    int steps = 1000;

    // Implementation
    Body *BodySystem = new Body[N];
    Collider Newton(N);

    for(int i=0; i<N; i++) BodySystem[i] = Body(N);

    std::ofstream orbit("earth.txt"), energy("energy.txt"), angular("angular.txt");

    for(int i=0; i<N; i++)
        BodySystem[i].initialize(
            initial[i][0], initial[i][1], initial[i][3], initial[i][4], initial[i][6], initial[i][7]
        );

    double dt = (t_max - t_min)/(double) steps;
    double E0 = Newton.energy(BodySystem), L0 = Newton.angular_momentum(BodySystem[1]);

    for(int t=0; t<steps; t++){
        orbit << BodySystem[1].get_x() << '\t' << BodySystem[1].get_y() << '\n';
        energy << t << '\t' << Newton.energy(BodySystem)/E0 << '\n';
        angular << t << '\t' << Newton.angular_momentum(BodySystem[1])/L0 << '\n';
        Newton.move_with_pefrl(BodySystem, dt);
    }

    orbit.close(); energy.close(); angular.close();

    for(int i=0; i<N; i++) delete[] initial[i];
    delete[] initial;

    return 0;
}