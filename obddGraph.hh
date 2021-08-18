#ifndef OBDDGRAPH_HH
#define OBDDGRAPH_HH

#include <iostream>
#include "cuddObj.hh"
#include <vector>
#include <string>

class obddGraphW;
class obddGraphUW;

template<class T>
class obddBase{
protected:
    int V, E, L, D, maxW;
    bool directed;
    Cudd mgr;
    BDD adjOBDD;
    std::vector<std::vector<T>> adjMat;
    std::vector<BDD> X;
    std::vector<BDD> Y;
    std::vector<BDD> W;

    virtual bool vec2val(std::vector<BDD>)=0;
    BDD CalcF(std::vector<BDD>);
public:
    obddBase(Cudd, bool directed=false);
    void printAdjOBDD(void);
    void AddV(int n=1);
    void makeOBDD(void);
    void optimizeOBDD(void);
    void changeOrder(int* order);
    std::pair<int, int> Size(void);
    void printAdjMat(void);
    void saveOBDD(void);
    void displayOBDD(void);
    void printInfo(void);
    ~obddBase();
};

class obddGraphW: public obddBase<int>{
private:
    std::vector<BDD> pX1;
    std::vector<BDD> pY1;
    //std::vector<BDD> pW1;
    std::vector<BDD> pX2;
    std::vector<BDD> pY2;
    //std::vector<BDD> pW2;
    std::vector<BDD> p1X;
    std::vector<BDD> p1Y;
    std::vector<BDD> p1Z;
    std::vector<BDD> Z;
    std::vector<BDD> tempX1; //x'
    std::vector<BDD> tempY1; //y'
    //std::vector<BDD> tempW1; //d'
    std::vector<BDD> tempX2; //x"
    std::vector<BDD> tempY2; //y"
    std::vector<BDD> tempW2; //d"
    std::vector<BDD> tempY3; //y"'
    BDD P;
    BDD P1;
protected:
    bool vec2val(std::vector<BDD>);
    bool vec2valP1(std::vector<BDD>, int);
    BDD findTC(BDD);
    BDD CalcP(std::vector<BDD>);
    BDD CalcP1(std::vector<BDD>, int);
    bool vec2valP(std::vector<BDD>);
    BDD findTCp(BDD);
public:
    obddGraphW(Cudd mgr, bool directed=false): obddBase(mgr, directed){
        //mgr.AutodynEnable();
        //mgr.SetNextReordering(10);
    };
    void makeP(void);
    void makeP1(int);
    void AddE(int, int, int);
    BDD findMST(void);
    BDD findEulerTour(BDD&);
    BDD findTSPTour(void);
    ~obddGraphW();
};

class obddGraphUW: public obddBase<bool>{
protected:
    bool vec2val(std::vector<BDD>);
public:
    obddGraphUW(Cudd mgr, bool directed=false): obddBase(mgr, directed){};
    void AddE(int, int);
    ~obddGraphUW();
};



//#include "obddGraph.cpp"

#endif

