#include<fstream>

#include"molecular_dynamics.h"


int main(int argc, char *argv[]){
    // Initial conditions
    int N, c; double **ini = NULL;
    load_file("earth_sun.csv", ini, N, c);

    double t_min = 0, t_max = 1e3; // days
    int steps = 1000;

    // Implementation
    Body *Bodies = new Body[N];
    Collider Newton(N);

    for(int i=0; i<N; i++) Bodies[i] = Body(N);

    std::ofstream orbit("earth.txt"), energy("energy.txt"), angular("angular.txt");

    for(int i=0; i<N; i++)
        Bodies[i].initialize(
            ini[i][0], ini[i][1], ini[i][2], ini[i][3], ini[i][4], ini[i][5], ini[i][6], ini[i][7]
        );

    double dt = (t_max - t_min)/(double) steps;
    double E0 = Newton.energy(Bodies), L0 = Newton.angular_momentum(Bodies[1]);

    for(int t=0; t<steps; t++){
        orbit << Bodies[1].get_x() << ',' << Bodies[1].get_y() << ',' << Bodies[i].get_z() << '\n';
        energy << t << ',' << Newton.energy(Bodies)/E0 << '\n';
        angular << t << ',' << Newton.angular_momentum(Bodies[1])/L0 << '\n';

        Newton.move_with_pefrl(Bodies, dt);
    }

    orbit.close(); energy.close(); angular.close();

    for(int i=0; i<N; i++) delete[] ini[i];
    delete[] ini;

    return 0;
}