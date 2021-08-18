#include<iostream>
#include "cuddObj.hh"

/*
BDD CalcR(std::vector<BDD> ip){
    int i=0;
    std::vector<BDD> temp1 = ip; 
    std::vector<BDD> temp2 = ip;
    while(i!=2*L){
        if(ip[i] == mgr.bddOne() || ip[i] == mgr.bddZero()){
            i++;
        }
        else{
            break;
        }
    }
    if(i == 2*L){
        bool z = vec2valR(ip);
        BDD c = (z)? mgr.bddOne() : mgr.bddZero();
        return c;
    }

    temp1[i] = mgr.bddOne();
    temp2[i] = mgr.bddZero();
    return ip[i].Ite(CalcR(temp1), ClacR(temp2));
}

bool obddGraphW::traverseOBDD(std::vector<BDD> xyw, BDD c){
    BDD t;
    BDD e;
    ADD curr = c.Add();
    for(int i=0; i<2*L; ++i){
        if(xyw[i] == mgr.bddOne()){
            curr = BDD(mgr,Cudd_T(curr.getNode()));
        }
        else if(xyw[i] == mgr.bddOne()){
            curr = BDD(mgr,Cudd_E(curr.getNode()));
        }


    }
}

bool obddGraphW::vec2valR(std::vector<BDD> xyw){
    int i = 0;
    std::string x="";
    std::string y="";
    for(i = 0; i<L; ++i){
        x += (xyw[i] == mgr.bddOne())? "1" : "0";
        y += (xyw[(L)+i] == mgr.bddOne())? "1" : "0";
    }
    int u = std::stoi(x, nullptr, 2);
    int v = std::stoi(y, nullptr, 2);

    int n[D];
    for(int i=0;i<D;++i){
        n[i] = (2*L)+i;
    }

    BDD d = mgr.IndicesToCube(n, D);
    BDD c = adjOBDD.ExistAbstract(d);

    return (u==v || adjOBDD.ExistAbstract(d))
}

BDD obddGraphW::findTC(){
    BDD R;

}

*/