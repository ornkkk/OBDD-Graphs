#include "obddGraph.hh"
#include <iostream>
#include "cuddObj.hh"
#include <vector>
#include <string>
#include <cmath>
#include <cstring>


template<class T>
obddBase<T>::obddBase(Cudd mgr, bool directed){
    this->V = 0;
    this->E = 0;
    this->mgr = mgr;
    this->L = 0;
    this->directed = directed;
    this->D = 0;
    this->maxW = 0;
}

template<class T>
obddBase<T>::~obddBase(){
    this->adjMat.clear();
}

template<class T>
BDD obddBase<T>::CalcF(std::vector<BDD> ip){
    int i=0;
    std::vector<BDD> temp1 = ip; 
    std::vector<BDD> temp2 = ip;
    while(i!=ip.size()){
        if(ip[i] == mgr.bddOne() || ip[i] == mgr.bddZero()){
            i++;
        }
        else{
            break;
        }
    }
    if(i == ip.size()){
        bool z = vec2val(ip);
        BDD c = (z)? mgr.bddOne() : mgr.bddZero();
        return c;

    }
    temp1[i] = mgr.bddOne();
    temp2[i] = mgr.bddZero();
    return ip[i].Ite(CalcF(temp1), CalcF(temp2));
}

template<class T>
void obddBase<T>::AddV(int n){
    if(n>=1){
        for(int i=0; i<n; ++i){
            V++;
            L = ceil(log2(V));
            std::vector<T> row(V, 0);
            this->adjMat.push_back(row);
            for(auto i=0; i<V-1; ++i){
                adjMat[i].resize(V, 0);
            }
        }
    }
}

template<class T>
void obddBase<T>::makeOBDD(void){
    X.resize(L);
    Y.resize(L);
    W.resize(D);
    for(int i=0; i<L; ++i){
        X[i] = mgr.bddVar();
        mgr.pushVariableName("x"+std::to_string(i));
    }
    for(int i=0; i<L; ++i){
        Y[i] = mgr.bddVar();
        mgr.pushVariableName("y"+std::to_string(i));
    }
     for(int i=0; i<D; ++i){
        W[i] = mgr.bddVar();
        mgr.pushVariableName("d"+std::to_string(i));
    }

    std::vector<BDD> ip;
    ip.insert(ip.end(), X.begin(), X.end());
    ip.insert(ip.end(), Y.begin(), Y.end());
    ip.insert(ip.end(), W.begin(), W.end());

    adjOBDD = CalcF(ip);
}

template<class T>
void obddBase<T>::optimizeOBDD(void){
    mgr.ReduceHeap();
}

template<class T>
void obddBase<T>::changeOrder(int* order){
    mgr.ShuffleHeap(order);
}

template<class T>
std::pair<int, int> obddBase<T>::Size(void){
    return std::make_pair(V,E);
}

template<class T>
void obddBase<T>::printAdjMat(void){
    for(auto i: adjMat){
        for(auto j: i){
            std::cout << j << " ";
        }
        std::cout << std::endl;
    }
}

template<class T>
void obddBase<T>::printAdjOBDD(void){
    std::cout << adjOBDD << std::endl;
}

template<class T>
void obddBase<T>::saveOBDD(void){
    std::string loc;
    std::cout << "Enter name of the .dot file to save: ";
    std::cin >> loc;
    loc = "./"+loc+".dot";

    std::string temp[(2*L) + D];
    const char* inames[(2*L + D)];
    const char* onames[1] = {"f"};

    const std::vector<ADD> fns = {adjOBDD.Add()};
    for(int i=0; i<L; ++i){
        temp[i] = "x"+std::to_string(i);
        inames[i] = temp[i].c_str();
        temp[L+i] = "y"+std::to_string(i);
        inames[L+i] = temp[L+i].c_str();
    }
    for(int i=2*L; i<(2*L)+D; ++i){
        temp[i] = "w"+std::to_string(i-(2*L));
        inames[i] = temp[i].c_str();
    }


    FILE* t = fopen(loc.c_str(), "w");
    mgr.DumpDot(fns, inames, onames, t);
    fclose(t);

    std::string cmd = "dot -Tpdf "+loc+" -O";
    system(cmd.c_str());
}

template<class T>
void obddBase<T>::displayOBDD(void){
    //std::string temp[(2*L) + D];
    //const char* inames[(2*L + D)];
    const char* onames[1] = {"f"};

    const std::vector<ADD> fns = {adjOBDD.Add()};
    
    /*
    for(int i=0; i<L; ++i){
        temp[i] = "x"+std::to_string(i);
        inames[i] = temp[i].c_str();
        temp[L+i] = "y"+std::to_string(i);
        inames[L+i] = temp[L+i].c_str();
    }
    for(int i=2*L; i<(2*L)+D; ++i){
        temp[i] = "w"+std::to_string(i-(2*L));
        inames[i] = temp[i].c_str();
    }
    */

    FILE* t = fopen("./temp.dot", "w");
    mgr.DumpDot(fns, NULL, onames, t);
    fclose(t);
    system("python3 -m xdot ./temp.dot");
    system("rm -rf ./temp.dot");
}

template<class T>
void obddBase<T>::printInfo(void){
    std::cout << "OBDD Info: " << std::endl;
    ADD f = adjOBDD.Add();
    f.print(2*L+D,5);
    return;
}

////////////////////////////////////////////////////////////////////////////

void obddGraphW::AddE(int a, int b, int w){
    E++;
    maxW = (w>maxW)? w:maxW;
    D = ceil(log2(maxW));
    this->adjMat[a-1][b-1] = w;
    if(!directed){
        this->adjMat[b-1][a-1] = w;
    }
}

bool obddGraphW::vec2val(std::vector<BDD> xyw){
    int i = 0;
    std::string x="";
    std::string y="";
    std::string d="";
    for(i = 0; i<L; ++i){
        x += (xyw[i] == mgr.bddOne())? "1" : "0";
        y += (xyw[(L)+i] == mgr.bddOne())? "1" : "0";
    }
    for(i=2*L; i<xyw.size(); ++i){
        d += (xyw[i] == mgr.bddOne())? "1" : "0";
    }
    int u = std::stoi(x, nullptr, 2);
    int v = std::stoi(y, nullptr, 2);
    int w = std::stoi(d, nullptr, 2);

    return (u<V && v<V && w<maxW && adjMat[u][v] != 0 && adjMat[u][v] == w+1)? true:false;
}

BDD obddGraphW::findTC(BDD Xe){
    BDD R, R1, R2, R3;
    tempX2.resize(L); //z
    for(int k=0; k<L; ++k){
        tempX2[k] = mgr.bddVar();
        mgr.pushVariableName("x''"+std::to_string(k));
    }
    BDD eq = mgr.bddOne();
    for(int i=0; i<L; ++i){
        eq *= !(X[i]^Y[i]);
    }

    BDD cubeTemp = mgr.computeCube(tempX2);
    R = eq + Xe;
    do{
        R1 = R;
        R2 = R1.SwapVariables(X, tempX2);
        R3 = R1.SwapVariables(Y, tempX2);
        R = (R3*R2).ExistAbstract(cubeTemp);
    }while(R != R1);

    return R;
}

void obddGraphW::makeP(void){
    pX1.resize(L);
    pY1.resize(L);
    //pW1.resize(D);
    for(int i=0; i<L; ++i){
        pX1[i] = mgr.bddVar();
        mgr.pushVariableName("px"+std::to_string(i));
    }
    for(int i=0; i<L; ++i){
        pY1[i] = mgr.bddVar();
        mgr.pushVariableName("py"+std::to_string(i));
    }
    /*
    for(int i=0; i<D; ++i){
        pW1[i] = mgr.bddVar();
        mgr.pushVariableName("pd"+std::to_string(i));
    }
    */
    
    pX2.resize(L);
    pY2.resize(L);
    //pW2.resize(D);
    for(int i=0; i<L; ++i){
        pX2[i] = mgr.bddVar();
        mgr.pushVariableName("px'"+std::to_string(i));
    }
    for(int i=0; i<L; ++i){
        pY2[i] = mgr.bddVar();
        mgr.pushVariableName("py'"+std::to_string(i));
    }
    /*
    for(int i=0; i<D; ++i){
        pW2[i] = mgr.bddVar();
        mgr.pushVariableName("pd'"+std::to_string(i));
    }
    */
    
    
    
    std::vector<BDD> pVar;
    pVar.insert(pVar.end(), pX1.begin(), pX1.end());
    pVar.insert(pVar.end(), pY1.begin(), pY1.end());
    //pVar.insert(pVar.end(), pW1.begin(), pW1.end());
    pVar.insert(pVar.end(), pX2.begin(), pX2.end());
    pVar.insert(pVar.end(), pY2.begin(), pY2.end());
    //pVar.insert(pVar.end(), pW2.begin(), pW2.end());

    this->P = CalcP(pVar);
}

BDD obddGraphW::CalcP(std::vector<BDD> p1){
    int i=0;
    std::vector<BDD> temp1 = p1; 
    std::vector<BDD> temp2 = p1;
    while(i!=p1.size()){
        if(p1[i] == mgr.bddOne() || p1[i] == mgr.bddZero()){
            i++;
        }
        else{
            break;
        }
    }
    if(i == p1.size()){
        bool z = vec2valP(p1);
        BDD c = (z)? mgr.bddOne() : mgr.bddZero();
        return c;

    }
    temp1[i] = mgr.bddOne();
    temp2[i] = mgr.bddZero();
    return p1[i].Ite(CalcP(temp1), CalcP(temp2));
}

bool obddGraphW::vec2valP(std::vector<BDD> p1){
    int i = 0, n=2*L;
    std::string x1="";
    std::string y1="";
    //std::string d1="";
    std::string x2="";
    std::string y2="";
    //std::string d2="";
    for(i = 0; i<L; ++i){
        x1 += (p1[i] == mgr.bddOne())? "1" : "0";
        y1 += (p1[(L)+i] == mgr.bddOne())? "1" : "0";
        x2 += (p1[n+i] == mgr.bddOne())? "1" : "0";
        y2 += (p1[n+L+i] == mgr.bddOne())? "1" : "0";
    }
    /*
    for(i=2*L; i<n; ++i){
        d1 += (p1[i] == mgr.bddOne())? "1" : "0";
        d2 += (p1[n+i] == mgr.bddOne())? "1" : "0";
    }
    */

    int u1 = std::stoi(x1, nullptr, 2);
    int v1 = std::stoi(y1, nullptr, 2);
    //int w1 = std::stoi(d1, nullptr, 2);

    int u2 = std::stoi(x2, nullptr, 2);
    int v2 = std::stoi(y2, nullptr, 2);
    //int w2 = std::stoi(d2, nullptr, 2);
    //std::cout << "Done..." << std::endl;

    bool e1 = (u1<V) && (v1<V) && (adjMat[u1][v1] != 0);
    bool e2 = (u2<V) && (v2<V) && (adjMat[u2][v2] != 0);

    if((e1==true) && (e2==true)){
        int w1 = adjMat[u1][v1];
        int w2 = adjMat[u2][v2];
        if(w1 < w2){
            return true;
        }
        else if(w1 == w2){
            if(std::min(u1, v1) < std::min(u2, v2)){
                return true;
            }
            else if(std::min(u1, v1) == std::min(u2, v2)){
                if(std::max(u1, v1) < std::max(u2, v2)){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }
        else{
            return false;
        }
    }
    else{
        return false;
    }
}



BDD obddGraphW::findMST(void){
    BDD MST{mgr.bddZero()}, MST1, C, R, temp, tempCube;
    makeP();

    tempX1.resize(L);   //z
    tempY1.resize(L);   //y'
    //tempW1.resize(D);   //d'

    for(int k=0; k<L; ++k){
        tempX1[k] = mgr.bddVar();  //Z
        mgr.pushVariableName("x'"+std::to_string(k));
    }
    
    for(int k=0; k<L; ++k){
        tempY1[k] = mgr.bddVar();  //Y'
        mgr.pushVariableName("y'"+std::to_string(k));
    }

    /*
    for(int k=0; k<D; ++k){
        tempW1[k] = mgr.bddVar();  //d'
        mgr.pushVariableName("d'"+std::to_string(k));
    }
    */

    std::vector<BDD> tmp;
    tmp.insert(tmp.end(), tempX1.begin(), tempX1.end());
    tmp.insert(tmp.end(), tempY1.begin(), tempY1.end());
    //tmp.insert(tmp.end(), tempW1.begin(), tempW1.end());
    tmp.insert(tmp.end(), X.begin(), X.end());
    tmp.insert(tmp.end(), Y.begin(), Y.end());
    //tmp.insert(tmp.end(), W.begin(), W.end());

    std::vector<BDD> pVar;
    pVar.insert(pVar.end(), pX1.begin(), pX1.end());
    pVar.insert(pVar.end(), pY1.begin(), pY1.end());
    //pVar.insert(pVar.end(), pW1.begin(), pW1.end());
    pVar.insert(pVar.end(), pX2.begin(), pX2.end());
    pVar.insert(pVar.end(), pY2.begin(), pY2.end());
    //pVar.insert(pVar.end(), pW2.begin(), pW2.end());

    
    tempCube = mgr.computeCube(tempY1)*mgr.computeCube(tempX1);
    BDD Xe = adjOBDD.ExistAbstract(mgr.computeCube(W));


    do{
        R = findTC(MST);
        temp = R.SwapVariables(Y, tempX1);
        temp *= ((Xe.SwapVariables(X, tempX1)).SwapVariables(Y, tempY1));
        temp *= !((R.SwapVariables(X, tempX1)).SwapVariables(Y, tempY1));
        temp *= P.SwapVariables(pVar, tmp);
        temp = (temp).ExistAbstract(tempCube);
        C = Xe * (!R) * !temp;
        C = C + C.SwapVariables(X, Y);
        MST1 = MST;
        MST = MST1 + C;
    } while(MST != MST1);
    return MST;
}

void obddGraphW::makeP1(int functionNumber=1){
    if(functionNumber != 3){
        p1X.resize(L);
        p1Y.resize(L);
        p1Z.resize(L);
        for(int i=0; i<L; ++i){
            p1X[i] = mgr.bddVar();
            mgr.pushVariableName("p1x"+std::to_string(i));
        }
        for(int i=0; i<L; ++i){
            p1Y[i] = mgr.bddVar();
            mgr.pushVariableName("p1y"+std::to_string(i));
        }
        for(int i=0; i<L; ++i){
            p1Z[i] = mgr.bddVar();
            mgr.pushVariableName("p1z"+std::to_string(i));
        }
    
        std::vector<BDD> p1Var;
        p1Var.insert(p1Var.end(), p1X.begin(), p1X.end());
        p1Var.insert(p1Var.end(), p1Y.begin(), p1Y.end());
        p1Var.insert(p1Var.end(), p1Z.begin(), p1Z.end());

        this->P1 = CalcP1(p1Var, functionNumber);
    }
    else{
        P1 = P;
    }
}

BDD obddGraphW::CalcP1(std::vector<BDD> p1, int functionNumber){
    int i=0;
    std::vector<BDD> temp1 = p1; 
    std::vector<BDD> temp2 = p1;
    while(i!=p1.size()){
        if(p1[i] == mgr.bddOne() || p1[i] == mgr.bddZero()){
            i++;
        }
        else{
            break;
        }
    }
    if(i == p1.size()){
        bool z = vec2valP1(p1, functionNumber);
        BDD c = (z)? mgr.bddOne() : mgr.bddZero();
        return c;

    }
    temp1[i] = mgr.bddOne();
    temp2[i] = mgr.bddZero();
    return p1[i].Ite(CalcP1(temp1, functionNumber), CalcP1(temp2, functionNumber));
}

bool obddGraphW::vec2valP1(std::vector<BDD> p1, int functionNumber=1){
    int i = 0;
    //std::string x="";
    std::string y="";
    std::string z="";

    if(functionNumber==1){
        for(i = 0; i<L; ++i){
            y += (p1[(L)+i] == mgr.bddOne())? "1" : "0";
            z += (p1[(L+L)+i] == mgr.bddOne())? "1" : "0";
        }

        int v = std::stoi(y, nullptr, 2);
        int w = std::stoi(z, nullptr, 2);
        return (v<w)? true:false;
    }
    else if(functionNumber==2){
        for(i = 0; i<L; ++i){
            y += (p1[L+i] == p1[i])? "0" : "1";
            z += (p1[L+L+i] == p1[i])? "0" : "1";
        }

        int v = std::stoi(y, nullptr, 2);
        int w = std::stoi(z, nullptr, 2);
        return (v<w)? true:false;
    }
    return false;
}

BDD obddGraphW::findEulerTour(BDD &Xed){
    //BDD Xeu = Xed.ExistAbstract(mgr.computeCube(W));
    //std::cout << Xeu << std::endl;
    BDD O, tempP, temp, temp1;

    tempY2.resize(L); //y"
    tempY3.resize(L); //y"'
    Z.resize(L); //z
    for(int i=0; i<L; ++i){
        tempY2[i] = mgr.bddVar();
        mgr.pushVariableName("y''"+std::to_string(i));
        tempY3[i] = mgr.bddVar();
        mgr.pushVariableName("y'''"+std::to_string(i));
        Z[i] = mgr.bddVar();
        mgr.pushVariableName("z"+std::to_string(i));
    }

    tempP = P1.SwapVariables(p1X, X);
    tempP = tempP.SwapVariables(p1Y, Y);
    tempP = tempP.SwapVariables(p1Z, tempY1);
    
    O = Xed * Xed.SwapVariables(Y, tempY1) * tempP;

    temp = Xed.SwapVariables(Y, tempY2);
    temp *= tempP.SwapVariables(tempY1, tempY2);
    temp *= tempP.SwapVariables(Y, tempY2);
    temp = temp.ExistAbstract(mgr.computeCube(tempY2));

    O *= !temp;

    temp = Xed.SwapVariables(Y, tempY2);
    temp *= tempP.SwapVariables(Y, tempY2);
    temp = temp.ExistAbstract(mgr.computeCube(tempY2));

    temp1 = Xed.SwapVariables(Y, tempY3);
    temp1 *= tempP.SwapVariables(tempY1, tempY3);
    temp1 = temp1.ExistAbstract(mgr.computeCube(tempY3));

    O = O + (Xed * Xed.SwapVariables(Y, tempY1) * !temp * !temp1);

    BDD ET = (O.SwapVariables(X, Y)).SwapVariables(tempY1, Z);

    return ET;
}

BDD obddGraphW::findTCp(BDD Pnew){
    BDD C, C1, C2, C3;
    Pnew = Pnew.SwapVariables(Z, tempY1);

    BDD eqX{mgr.bddOne()}, eqY{mgr.bddOne()}, eqX1Y{mgr.bddOne()}, tempP;
    for(int i=0; i<L; ++i){
        eqX *= !(X[i]^tempX1[i]);
        eqY *= !(Y[i]^tempY1[i]);
        eqX1Y *= !(Y[i]^tempX1[i]);
    }


    C = (eqX * eqY) + (Pnew * eqX1Y);

    do{
        C1 = C;
        C2 = (C.SwapVariables(tempX1, tempX2)).SwapVariables(tempY1, tempY2);
        C3 = (C.SwapVariables(X, tempX2)).SwapVariables(Y, tempY2);
        C = (C2*C3).ExistAbstract(mgr.computeCube(tempX2)*mgr.computeCube(tempY2));
    } while(C != C1);

    return C;
}

BDD obddGraphW::findTSPTour(void){
    //makeP1();
    BDD tempP, temp, temp1, Pnew, SC, TSP;

    //BDD M = (findMST()).ExistAbstract(mgr.computeCube(W));
    BDD M = findMST();
    BDD ET = findEulerTour(M);

    tempP = ((P1.SwapVariables(p1X, Y)).SwapVariables(p1Y, tempX1)).SwapVariables(p1Z, X);
    temp = (M.SwapVariables(X, tempX1)*tempP).ExistAbstract(mgr.computeCube(tempX1));

    BDD RED = M * !temp;

    Pnew = ET * !RED;

    SC = findTCp(Pnew);

    temp = RED.SwapVariables(X, tempX1);
    temp *= SC.SwapVariables(Y, tempY1);

    temp1 = ((Pnew.SwapVariables(X, tempX2)).SwapVariables(Y, X)).SwapVariables(Z, tempY1);
    temp1 = temp1.ExistAbstract(mgr.computeCube(tempX2));

    temp *= !temp1;

    TSP = temp.ExistAbstract(mgr.computeCube(tempX1)*mgr.computeCube(tempY1));

    return TSP;  
}

obddGraphW::~obddGraphW(){
    adjMat.clear();
}

////////////////////////////////////////////////////////////////////////////

void obddGraphUW::AddE(int a, int b){
    E++;
    this->adjMat[a-1][b-1] = true;
    if(!directed){
        this->adjMat[b-1][a-1] = true;
    }
}

bool obddGraphUW::vec2val(std::vector<BDD> xy){
    int i = 0;
    std::string x="";
    std::string y="";
    for(i = 0; i<L; ++i){
        x += (xy[i] == mgr.bddOne())? "1" : "0";
        y += (xy[(L)+i] == mgr.bddOne())? "1" : "0";
    }
    int u = std::stoi(x, nullptr, 2);
    int v = std::stoi(y, nullptr, 2);

    return (u<V && v<V && adjMat[u][v])? true:false; 
}

obddGraphUW::~obddGraphUW(){
    adjMat.clear();
}

////////////////////////////////////////////////////////////////////////////
