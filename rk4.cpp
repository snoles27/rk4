#include <iostream>
#include <strings.h>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>

const int stateSize = 2;
const int numStatesStore = 50; 

//"set startup-with-shell off" when debuggging

typedef Eigen::Matrix<double,stateSize,1> stateVec;

stateVec rk4Step(stateVec (*dz)(double, stateVec), stateVec z_vec, double ti, double dt);

stateVec stateDerivative(double t, stateVec z){

    stateVec dz;
    dz[0] = z[1];
    dz[1] = -z[0];

    return dz;
}


stateVec rk4(stateVec (*dz)(double, stateVec), stateVec z0, double t0, double dt, int n, std::string fileName){

    double* zlog = new double[stateSize * numStatesStore]; //z0
    double* tspanlog = new double[numStatesStore];
    
    stateVec current = z0;
    double currentTime = t0;
    std::ofstream file;
    file.open(fileName, std::ofstream::out | std::ofstream::trunc);

    //put first state in the storeArray
    tspanlog[0] = currentTime;
    for(int i = 0; i < stateSize; i++){
        zlog[i] = z0[i];
    }

    int j = 1; //counter for states in array

    for(int i = 0; i < n; i++){

        current = rk4Step(dz, current, currentTime, dt);
        currentTime += dt;

        //fill zlog with current state
        tspanlog[j] = currentTime;
        for(int k = 0; k < stateSize; k++){
            zlog[j*stateSize + k] = current[k];
        }

        //incriment log array
        j++;

        if (j == numStatesStore || i == (n-1))
        {
            j = 0;
            //write array to file
            for(int k = 0; k < numStatesStore; k++){
                file << tspanlog[k] << ", ";
                tspanlog[k] = 0;
                for(int l = 0; l < stateSize; l++)
                {
                    file << zlog[(k*stateSize) + l] << ", ";
                    zlog[(k*stateSize) + l] = 0.0;
                }
                file << std::endl;
            }
        }
    }
    file.close();
    return current;
}

stateVec rk4Step(stateVec (*dz)(double, stateVec), stateVec z_vec, double ti, double dt){

    stateVec k1, k2, k3, k4;
    k1 = dz(ti, z_vec);
    k2 = dz(ti + 0.5 * dt, z_vec + 0.5 * dt * k1);
    k3 = dz(ti + 0.5 * dt, z_vec + 0.5 * dt * k2);
    k4 = dz(ti + dt, z_vec + dt * k3);

    return z_vec + 1.0/6.0 * (k1 + 2*k2 + 2*k3 + k4) * dt;
}

int main(){

    stateVec z0 (1,0);
    double t0 = 0;
    double dt = pow(10,-3);
    double n = 4999;
    std::string saveFile = "data.txt";
    

    stateVec final = rk4(stateDerivative, z0, t0, dt, n, saveFile);

    std::cout << "Returned rk4" <<std::endl;


}
