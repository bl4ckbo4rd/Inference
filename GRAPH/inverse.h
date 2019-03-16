/*
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iterator>
#include <algorithm>
*/
#include "../../../../dlib/optimization.h"
#include <random>

using namespace dlib;
using namespace std;

typedef matrix<double,0,1> column_vector;                                                               //vectors used by the library dlib to perform the optimization

//this function optimizes the PSL by updating iteratively the links that give the best PSL gain
//input variables:
//N         : # of Nodes
//L         : # of samples
//T         : # of time steps for which we iterate the optization routine
//beta      : inverse temperature
//s         : sampled data
//rho_acc   : fraction of J's that we update at each step
//eta       : learning rate used to update the couplings. set it to ~ 0.5
//f         : set to 0 to study FB (Free Boundary) 2D lattices and to 1 to study PB (Periodic Boundary) 2D lattices
//rho       : density of positive couplings in the original graph
//a_size    : factor that set the initial size of the vector that contains the "good couplings"
//a_th      : flag that set the threshold used to eliminate couplings from the vector of "good couplings"
//numberMB  : number of minibatch used in the gradient ascent routine
//CLUSTER   : set to 1 to run on cluster, 0 otherwise
//test_mode : flag set to 0 by default that can be set to 1 to have a verbose version of the program

void network_setup(int, int, int, double, std::vector < std::vector < int > >&, double, double, int, double, double, double, int, int, int);

extern std::vector < std::vector <int> > a;                                                             //a is an L * N matrix
                                                                                                        //it contains the L sampled configurations, one per row.
                                                                                                        //it is defined as a global because the constructor of the class PSL,
                                                                                                        //that was the natural place where to store this data, has to be
                                                                                                        //called several times in the operator().
                                                                                                        //Making it a global variable, the operator() can access its elements.
                                                                                                        //The extern is needed in order to avoid duplicates.
                                                                                                        //In fact, this header is included both in inverse.cpp and in problems.h
                                                                                                        //and without the "extern", the matrix a would be declared in both of them.
                                                                                                        //a is filled in inverse.cpp.

//this class defined the object we would like to optimize
//we do so by using the LBSFG method in the function network setup
//we used it to perform an exact optimzation over the activated variables

class PSL
{
private:
    int L, N;                                                                                           //L: # of samples; N: # of spins
    double beta;                                                                                        //inverse temperature
    std::vector < std::vector <int> > neigh;
    std::vector <int> v_ind;
    std::vector <int> v_q;
    
    int Npar, K, K2;
    
public:
    
    PSL (int p_N, int p_L, double p_beta, int p_K, std::vector < std::vector <int> >& p_neigh, std::vector <int>& p_v_ind, std::vector <int>& p_v_q)
    {
        
        N        = p_N;
        L        = p_L;
        beta     = p_beta;
        K        = p_K;
        
        neigh.resize(N);
        
        K2 = 2*K;
        
        Npar = 0.5 * (N-1) * N;
        
        for (int r = 0; r < N; r++){
            int c = p_neigh[r].size();
            neigh[r].resize(c);
            for (int j = 0; j < c; j++){
                neigh[r][j] = p_neigh[r][j];
            }
        }
        
        v_ind.resize(K2);
        for (int q2 = 0; q2 < K2; q2++)
            v_ind[q2] = p_v_ind[q2];
        
        v_q.resize(Npar);
        for (int i = 0; i < Npar; i++)
            v_q[i]   = p_v_q[i];
        
    }
    
    
    double operator() ( const column_vector& J) const
    {
        
        int r, k, j;
        int jj, rr, ind, c;

        double PSL_r;
        double PSL = 0.;
        double h, JJ;
        
        int q2 = 0;
        
        int saved_q2 = q2;
        
        for (r = 0; r < N; r ++){
            PSL_r = 0;
            
            c = neigh[r].size();
            
            saved_q2 = q2;
            
            for (k = 0; k < L; k++){
                h = 0.;
                
                q2 = saved_q2;
                
                for (j = 0; j < c; j ++){

                    jj  =  neigh[r][j];
                    JJ  =  J(v_q[v_ind[q2]]);
                    h   += JJ * a[k][jj];
                    
                    q2++;

                }

                PSL_r += log ( 1./( 1. + exp ( -2 * beta * a[k][r] * h ) ) );

            }
            
            PSL += PSL_r;
            q2  = saved_q2 + c;
            
        }
        
        return - PSL;
     
        
    }
    
    
    
};


//this function optimizes the PSL with a gradient ascent on the activated couplings.
//it can be called in the place of the exact optimization performed with the LBSFG Newton method perfomed with the class PSL
//input variables:
//N         : # of Nodes
//L         : # of samples
//beta      : inverse temperature
//K         : # of activated J
//alpha     : learning rate in the gradiend ascent
//max_it    : maximum number of iteration in the gradient ascent "exact" optimization of the PSL
//precision : sets the presion with with we compute the maximum of the PSL in the gradient ascent
//neigh     : vector of neighbours of each site
//v_ind     : vector that associate the same ind the couple (ij) and (ji)
//v_q       : vector that takes ind and gives an index 0<q<K.


void GradientAscent(int, int, double, int, double, int, double, std::vector < std::vector <int> >&, std::vector <int>&, std::vector <int>&, column_vector &, int numberMB);

