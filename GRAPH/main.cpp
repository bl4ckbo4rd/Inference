#include "problems.h"

int main (int argc, char *argv[]){

    //N        : # of Nodes
    //L        : # of samples
    //Ltot     : # of total samples in the dataset
    //T        : # of time steps for which we iterate the optization routine
    //beta     : inverse temperature
    //rho_acc  : fraction of J's that we update at each step. set it such that we uptade one link per iteration
    //eta      : learning rate used to update the couplings. set it to 1
    //f        : set to 0 to study FB (Free Boundary) 2D lattices and to 1 to study PB (Periodic Boundary) 2D lattices
    //           set to 2 to study ER graphs and to 3 to study RR graphs, 4 to study diamond lattices
    //rho      : density of positive couplings in the original graph
    //w        : seed used to extract the L random samples
    //a_size   : factor that set the initial size of the vector that contains the "good couplings". set it to ~ O(number of couplings that we expect to find).
    //a_th     : flag that set the threshold used to eliminate more couplings from the vector of "good couplings".
    //           all the couplings whose PSLgain is smaller than a_th times the largest value of the PSLgain are eliminated from the vector of good couplings
    //           and thus are not considered in the sorting
    //           this provides an extra optimization, decreasing the size of the vector of good couplings over which we optimize. set it to 0 not to use it.
    //numberMB : number of minibatch used in the gradient ascent routine
    //CLUSTER  : set to 1 to run on cluster, 0 otherwise
    
    int N, L, Ltot, T, f, w;
    double beta, rho, rho_acc, eta, a_size, a_th;
    int numberMB;
    int CLUSTER;
    
    if(argc == 15){
        int i = 1;
        N        = atoi(argv[i++]);
        L        = atoi(argv[i++]);
        Ltot     = atoi(argv[i++]);
        T        = atoi(argv[i++]);
        beta     = atof(argv[i++]);
        rho_acc  = atof(argv[i++]);
        eta      = atof(argv[i++]);
        f        = atoi(argv[i++]);
        rho      = atof(argv[i++]);
        w        = atoi(argv[i++]);
        a_size   = atof(argv[i++]);
        a_th     = atof(argv[i++]);
        numberMB = atoi(argv[i++]);
        CLUSTER  = atoi(argv[i++]);
    }
    else{
        cout << "argument: N, L, Ltot, T, beta, rho_acc, eta, f, rho, w, a_size, a_th, numberMB, CLUSTER" << endl;
        return 0;
    }
    
    
    
    f_Inference(N, L, Ltot, T, beta, rho_acc, eta, f, rho, w, a_size, a_th, numberMB, CLUSTER);
        
    return 1;
}
