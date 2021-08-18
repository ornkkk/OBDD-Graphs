#include "obddGraph.hh"
#include "obddGraph.cpp"
#include <iostream>
#include <fstream>

int main(){
    Cudd mgr(0,0);
    obddGraphW G(mgr);

    std::ifstream f;
    f.open("./Data/xqf131.tsp");
    if (f.fail()){
		std::cerr << "Can't open input file!\n" << std::endl;
		exit(1);		
	}
	
    std::vector<int> X;
    std::vector<int> Y;
	int i=1, N, j=0, t;
    std::string temp, t1, t2;
	while (!f.eof()){
        if(i==6){
            f >> temp >> temp >> N;
            i++;
            Y.resize(N);
            X.resize(N);
            continue;
        }
        if(i<=9){
            f.ignore(500, '\n') ;
            i++;
            continue;
        }
        else if(i == N+10){
            break;
        }
        else{
            f >> t1 >> X[j] >> Y[j];
            j++;
            i++;
        }
	}
	f.close();
    for(i=0; i<X.size(); ++i){
        std::cout << X[i] << " " << Y[i] << std::endl;
    }

    G.AddV(N);
    double dist;
    for(int j=0; j<N; ++j){
        for(int k=0; k<N; ++k){
            dist = (double)std::sqrt(std::pow((X[j]-X[k]),2) + std::pow((Y[j]-Y[k]),2));
            G.AddE(j+1, k+1, dist);
        }
    }

    G.printAdjMat();
    G.makeOBDD();
    std::cout << "Done" << std::endl;
    G.makeP();
    std::cout << "Done" << std::endl;
    G.makeP1();
    std::cout << "Done" << std::endl;
    
    
    BDD T = G.findTSPTour();
    
    const std::vector<ADD> fns = {T.Add()};
    FILE* tempFile = fopen("./tempM.dot", "w");
    mgr.DumpDot(fns, NULL, NULL, tempFile);
    fclose(tempFile);
    system("python3 -m xdot ./tempM.dot");
    system("rm -rf ./tempM.dot");
    
}