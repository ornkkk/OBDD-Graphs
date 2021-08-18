#include "obddGraph.hh"
#include <cudd.h>
#include <iostream>
#include "cuddObj.hh"
#include "obddGraph.cpp"
#include <chrono>

#include <vector>
#include <string>
#include <cmath>
#include <cstring>
#include "dddmp.h"

using namespace std::chrono;

int main(){
    Cudd mgr(0,0);
    obddGraphW G(mgr);
   
    G.AddV(6);
    G.AddE(1, 2, 11);
    G.AddE(1, 3, 14);
    G.AddE(1, 4, 15);
    G.AddE(1, 5, 10);
    G.AddE(1, 6, 10);
    G.AddE(2, 3, 9);
    G.AddE(2, 4, 16);
    G.AddE(2, 5, 15);
    G.AddE(2, 6, 9);
    G.AddE(3, 4, 11);
    G.AddE(3, 5, 13);
    G.AddE(3, 6, 6);
    G.AddE(4, 5, 8);
    G.AddE(4, 6, 10);
    G.AddE(5, 6, 9);
    //G.printAdjMat();
    //std::cout << "Done" << std::endl;
    G.makeOBDD();
    //std::cout << "Done" << std::endl;
    auto start = high_resolution_clock::now();
    //std::cout << "Done" << std::endl;
    G.makeP();
    //std::cout << "Done" << std::endl;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    
    G.makeP1();
    //std::cout << "Done" << std::endl;
    //G.printAdjOBDD();
    //G.optimizeOBDD();
    //G.printInfo();
    //G.displayOBDD();
    //int order[18]={0,1,2,3,4,5,11, 17, 10, 16, 7, 13, 9, 15, 6, 12, 8, 14};
    //G.changeOrder(order);
    //G.displayP();
    
    //BDD M = G.findMST();
    //BDD T = G.findEulerTour(M);
    BDD T = G.findTSPTour();
    //std::cout << M << std::endl;
    const std::vector<ADD> fns = {T.Add()};
    FILE* tempFile = fopen("./tempM.dot", "w");
    mgr.DumpDot(fns, NULL, NULL, tempFile);
    fclose(tempFile);
    system("python3 -m xdot ./tempM.dot");
    system("rm -rf ./tempM.dot");
    
    return 0;
}