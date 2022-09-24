#include <iostream>


const int n = 4;

void printAll(double* array, int size);

int main(){

    double balance[n] = {100, 200, 300, 400};
    printAll(&balance[0], n);

}

void printAll(double* array, int size){

    for(int i = 0; i < size; i++){
        std::cout << array[i] << std::endl;
    }

    std::cout<<std::endl;

    for(int i = 0; i < size; i++){
        std::cout << *(array + i) << std::endl;
    }


}