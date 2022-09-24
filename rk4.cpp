#include <iostream>
#include <Eigen/Dense>
#include <math.h>

const int stateSize = 2;

typedef Eigen::Matrix<double,stateSize,1> stateVec;

void rk4Step(stateVec (*dz)(double, stateVec), double* z, double ti, double dt, int i);

stateVec stateDerivative(double t, stateVec z){

    stateVec dz;
    dz[0] = z[1];
    dz[1] = -z[0];

    return dz;
}


double* rk4(stateVec (*dz)(double, stateVec), double* z0, double t0, double dt, int n){

    double* z = new double[stateSize * (n + 1)]; //z0 plus n states
    double* tspan = new double(n+1);

    printf("Pointer: %p\n", z);

    tspan[0] = t0;

    for(int i = 0; i < stateSize; i++){
        z[i] = z0[i];
    }

    for(int i = 0; i < n; i++){
        tspan[i+1] = dt + tspan[i];
        rk4Step(dz, z, tspan[i], dt, i);
        //std::cout << z[(i + 1) * stateSize] << std::endl;
    }

    for(int i = 0; i < stateSize * n + 1; i = i + 2){
        std::cout << z[i] << std::endl;
    }

    return z;
    

}

void rk4Step(stateVec (*dz)(double, stateVec), double* z, double ti, double dt, int i){

    stateVec z_vec;
    for(int j = 0; j < stateSize; j++){
        z_vec[j] = z[i*stateSize+j];
    }


    stateVec k1, k2, k3, k4;

    k1 = stateDerivative(ti, z_vec);
    k2 = stateDerivative(ti + 0.5 * dt, z_vec + 0.5 * dt * k1);
    k3 = stateDerivative(ti + 0.5 * dt, z_vec + 0.5 * dt * k2);
    k4 = stateDerivative(ti + dt, z_vec + dt * k3);

    

    stateVec z1_vec = z_vec + 1.0/6.0 * (k1 + 2*k2 + 2*k3 + k4) * dt;

    //write z1_vec to z now

    for(int j = 0; j < stateSize; j++){
        z[(i+1)*stateSize + j] = z1_vec[j];
        //std::cout<<"Write to: " << (i+1)*stateSize + j <<std::endl;
    }

}

int main(){

    double z0[stateSize] = {1, .01};
    double t0 = 0;
    double dt = pow(10,-3);
    int n = 5000;

    double* z = rk4(stateDerivative, z0, t0, dt, n);

    std::cout << "Returned rk4" <<std::endl;

    // for(int i = 0; i < stateSize * n + 1; i = i + 2){
    //     std::cout << z[i] << std::endl;
    // }

    delete[] z;

}
