#include "inverse.h"

std::vector < std::vector <int> > a;                                                                    //this vector was already declared in the heaader as a global vector
                                                                                                        //and we re-declare it here because it is where we need to use it.
                                                                                                        //it contains L configurations sampled from L MonteCarlo's
                                                                                                        //after a *long* waiting time, specified as a parameters of the dynamics (Nsweeps)

void network_setup(int N, int L, int T, double beta, std::vector < std::vector < int > >& s, double rho_acc, double eta, int f, double rho, double a_size, double a_th, int numberMB, int CLUSTER, int test_mode){
 
    int i, k, r, j, q, q2, rr, jj, i_tu, r_up1, r_up2;
    int ind, ind_j, ind_max;
    double denr, denj, denrr, denjj;
    double JJ;
    
    //# of total possible parameters we may need to optimize over
    int Npar = 0.5 * (N-1) * N;
    
    //# of parameters we update at each step
    int Nup  = (int) (rho_acc * Npar);
    
    //# of activated couplings
    int K;
    
    //connectivity
    int c;
    
    cout << "# of couplings we update at each step:" << endl;
    cout << Nup << endl;
    cout << "# of possible couplings" << endl;
    cout << Npar << endl;
    
    double PSL_c, BIC;
    double s2_r, s_PSL, s_BIC;
    
    double Num, Den;
    
    //first and second derivatives of the PseudoLikelihood
    std::vector <double> DP(Npar,0.);
    std::vector <double> DDP(Npar,0.);
    std::vector <double> PSLgain(Npar,0.);
    
    
    //vector of local fields
    std::vector < std::vector<double> > h(L, std::vector<double>(N,0.));
    
    //vector of couplings
    std::vector <double> J(Npar,0.);
    
    //vector of indices of the couplings
    std::vector <int> index_J(Npar,0.);
    
    //vector of indices of the activated couplings
    std::vector <int> activated_J;
    
    //vector that contains 1 on activated links and 0 otherwise
    std::vector <int> flag_J(Npar,0.);

    //vector of true couplings
    std::vector <double> J_true(Npar,0.);
    
    //vector of neighbours of each node
    std::vector < std::vector <int> > neigh;
    neigh.resize(N);
    
    //the coupling J(r,j) is contained in J(ind), where ind is
    //ind = 0.5 * ( N * (N-1) - (N-r)*(N-r-1) ) + j - (r+1)
    //map lets recover which r and j correspond to a given ind
    std::vector < std::vector <int> > map;
    map.resize(Npar);
    for (r = 0; r < N; r ++){
        for (j = r + 1; j < N; j ++){
            //index running over couplings
            ind = 0.5 * ( N * (N-1) - (N-r)*(N-r-1) ) + j - (r+1);
            map[ind].resize(2);
            map[ind][0] = r;
            map[ind][1] = j;
        }
    }

    //a contains the data and it is a global object
    a.resize(L);
    for (int l = 0; l < L; l ++){
        a[l].resize(N);
        a[l] = s[l];
    }
    
    r = 0;
    //upload the true couplings
    if (test_mode == 1){
        
        ostringstream convert1;
        convert1 << beta;
        string s_beta = convert1.str();
        
        ostringstream convert;
        convert << N;
        string size = convert.str();
        
        ostringstream convert2;                                                              // stream used for the conversion
        convert2 << rho;                                                                     // insert the textual representation of 'Number' in the characters in the stream
        string s_rho = convert2.str();
        
        
        string datafile;
        
        string path;

        if (CLUSTER)
            path = "/home/jrocchi/DATA/INFERENCE/Dynamics/";
        else
            path = "/Users/jrocchi/DATA/INFERENCE/Dynamics/";
        
        if(f == 0)
            datafile = path + "beta=" + s_beta + "/graph_" + size +"_2D-FB_" + s_rho + ".txt" ;
        
        if(f == 1)
            datafile = path + "beta=" + s_beta + "/graph_" + size +"_2D-PB_" + s_rho + ".txt" ;
        
        if(f == 2)
            datafile = path + "beta=" + s_beta + "/graph_" + size +"_ER_" + s_rho + ".txt" ;
        
        if(f == 3)
            datafile = path + "beta=" + s_beta + "/graph_" + size +"_RR_" + s_rho + ".txt" ;
        
        if(f == 4)
        datafile = path + "beta=" + s_beta + "/graph_" + size +"_Diamond_" + s_rho + ".txt" ;
        
        ifstream infile(datafile);
        
        std::vector < std::vector <int> > neighs;
        neighs.resize(N);
        
        // If we couldn't open the output file stream for reading
        if (!infile)
        {
            cerr << "Uh oh, graph.txt could not be opened for reading!" << endl;
            exit(1);
        }
        
        string temp;
        
        while (std::getline(infile, temp)) {
            
            istringstream buffer(temp);
            std::vector<int> line( (istream_iterator<int>(buffer)), istream_iterator<int>() );
            neighs[r] = line;
            
            //c is the connectivity of the node
            //the 0.5 factor comes from the structure of the file graph.txt where line i contains
            //a number of elements which is twice the number of its neighbours. see below
            //line "r" of the file graph.txt contains the list of neighbours of the node "r". near each site, there is the value of the coupling
            //e.g. :
            //1 +1 3 -1 5 +1
            //2 -1 0 +1
            //means that node 0 is connected to node 1, 3 and 5 with couplings given by +1, -1 and +1
            //and that node 1 is connected to node 2 and 0 with couplings given by -1 and +1.
            
            c = (int)(0.5 * neighs[r].size());
            
            for (int k = 0; k < c; k++){
                
                int j = neighs[r][2*k];
                
                if (j < r){
                    jj  = r;
                    rr  = j;
                }
                
                if (j > r){
                    jj  = j;
                    rr  = r;
                }
                
                ind = 0.5 * ( N * (N-1) - (N-rr)*(N-rr-1) ) + jj - (rr+1);
                J_true[ind] = neighs[r][2*k+1];
            }
            
            r++;
            
        }
        
    }
    
    //this is an auxiliary vector where I store the nodes involved in the couplings that I add in the graph
    //if at time t the coupling (i-j) is added, I save in this vector i and j and at time t+1 I will update the PSL gain on the couplings (i-k) and (j-k), for k = 1,...,N
    std::vector<int> nodes_to_update(2*Nup);

    //aN is the size of the vector good_couplings that I sort during the iteration, made by a * N couplings.
    //The size of this vector evolves during the iteration since new couplings are included at time t if they are found to be greater than the smallest element of good_couplings found at time t = 1.
    int aN = (int)(a_size * (double)N);
    
    //vect contains the first aN couplings, ordered in terms of the PSL gain
    std::vector< pair <int,double> > vect(aN);
    std::vector<int> map_ind(Npar,-1);
    
    int size_vect, cnt_ind;
    double PSL_u, mean_vect;

    for (int time = 1; time <= T; time ++){
        
        //lista dei couplings di cui devo ricalcolare il PSL gain nei tempi > 1
        std::vector<int> ind_to_update;

        //run over the whole dataset
        for (k = 0; k < L; k++){
            
            if (time == 1) {
                
                //compute the local field for each spin
                for (r = 0; r < N; r ++){
                    
                    //compute hr
                    for (j = 0; j < N; j ++){
                        if( j < r){
                            jj    = r;
                            rr    = j;
                        }
                        if( j > r){
                            jj    = j;
                            rr    = r;
                        }
                        if( j == r)
                            continue;
                        
                        ind_j = 0.5 * ( N * (N-1) - (N-rr)*(N-rr-1) ) + jj - (rr+1);
                        
                        JJ    = J[ind_j];
                        h[k][r]  += JJ * s[k][j];
                    }
                    
                }

                //compute the first and second derivative of the Pseudolikelihood for each coupling
                for (r = 0; r < N; r ++){
                    for (j = r + 1; j < N; j ++){
                        
                        //index running over couplings
                        ind = 0.5 * ( N * (N-1) - (N-r)*(N-r-1) ) + j - (r+1);
                        
                        denr  = 1 + exp( 2 * beta * s[k][r] * h[k][r] );
                        denj  = 1 + exp( 2 * beta * s[k][j] * h[k][j] );
                        
                        denrr = 2 + 2 * cosh( 2 * beta * s[k][r] * h[k][r] );
                        denjj = 2 + 2 * cosh( 2 * beta * s[k][j] * h[k][j] );
                        
                        DP[ind]  += 2 * beta * s[k][r] * s[k][j] * ( 1./denr + 1./denj );
                        DDP[ind] += - 4 * beta * beta * ( 1./denrr + 1./denjj);
                        
                    }
                }
                
            }
            if(time > 1){
                
                for (int it_u = 0; it_u < Nup; it_u+=2 ){
                
                    r_up1 = nodes_to_update[it_u];
                    r_up2 = nodes_to_update[it_u+1];
                    
                    //aggiorno i campi sui nodi r_up1 e r_up2 che sono coinvolti nell'aggiornamento del coupling fatto nell'iterazione precedente
                    
                    h[k][r_up1] = 0.;
                    h[k][r_up2] = 0.;
                    
                    //compute hr1
                    for (j = 0; j < N; j ++){
                        if( j < r_up1){
                            jj    = r_up1;
                            rr    = j;
                        }
                        if( j > r_up1){
                            jj    = j;
                            rr    = r_up1;
                        }
                        if( j == r_up1)
                            continue;
                        
                        ind_j = 0.5 * ( N * (N-1) - (N-rr)*(N-rr-1) ) + jj - (rr+1);

                        if (k == 0){
                            DP[ind_j] = 0.;
                            DDP[ind_j] = 0.;
                        }
                        
                        JJ    = J[ind_j];
                        h[k][r_up1]  += JJ * s[k][j];
                    }
                    
                    //compute hr2
                    for (j = 0; j < N; j ++){
                        if( j < r_up2){
                            jj    = r_up2;
                            rr    = j;
                        }
                        if( j > r_up2){
                            jj    = j;
                            rr    = r_up2;
                        }
                        if( j == r_up2)
                            continue;
                        
                        ind_j = 0.5 * ( N * (N-1) - (N-rr)*(N-rr-1) ) + jj - (rr+1);
                        
                        if (k == 0){
                            DP[ind_j] = 0.;
                            DDP[ind_j] = 0.;
                        }
                        
                        JJ    = J[ind_j];
                        h[k][r_up2]  += JJ * s[k][j];
                    }
                    
                    //re - compute the first and second derivative of the Pseudolikelihood for the couplings involving r_up1 and r_up2

                    //for r_up1, run over all the other nodes except r_up1, including r_up2
                    //the link r_up1 - r_up2 will be therefore excluded in the next for because it is already considered here
                    for (j = 0; j < N; j ++){
                        if ( j < r_up1 ){
                            jj    = r_up1;
                            rr    = j;
                        }
                        if ( j > r_up1 ){
                            jj    = j;
                            rr    = r_up1;
                        }
                        if(j == r_up1)
                            continue;
                        
                        ind = 0.5 * ( N * (N-1) - (N-rr)*(N-rr-1) ) + jj - (rr+1);
                        
                        //lista dei couplings di cui devo ricalcolare il PSL gain nei tempi > 1
                        if (k == 0)
                            ind_to_update.push_back(ind);
                        
                        denr  = 1 + exp( 2 * beta * s[k][rr] * h[k][rr] );
                        denj  = 1 + exp( 2 * beta * s[k][jj] * h[k][jj] );
                        
                        denrr = 2 + 2 * cosh( 2 * beta * s[k][rr] * h[k][rr] );
                        denjj = 2 + 2 * cosh( 2 * beta * s[k][jj] * h[k][jj] );
                        
                        DP[ind]  += 2 * beta * s[k][rr] * s[k][jj] * ( 1./denr + 1./denj );
                        DDP[ind] += - 4 * beta * beta * ( 1./denrr + 1./denjj);
                    }
                    
                    
                    
                    //for r_up2, run over all the other nodes except r_up2 and r_up1
                    for (j = 0; j < N; j ++){
                        if ( j < r_up2 ){
                            jj    = r_up2;
                            rr    = j;
                        }
                        if ( j > r_up2 ){
                            jj    = j;
                            rr    = r_up2;
                        }
                        if (j == r_up2 || j == r_up1)
                            continue;
                        
                        ind = 0.5 * ( N * (N-1) - (N-rr)*(N-rr-1) ) + jj - (rr+1);
                        
                        //lista dei couplings di cui devo ricalcolare il PSL gain nei tempi > 1
                        if (k == 0)
                            ind_to_update.push_back(ind);
                        
                        denr  = 1 + exp( 2 * beta * s[k][rr] * h[k][rr] );
                        denj  = 1 + exp( 2 * beta * s[k][jj] * h[k][jj] );
                        
                        denrr = 2 + 2 * cosh( 2 * beta * s[k][rr] * h[k][rr] );
                        denjj = 2 + 2 * cosh( 2 * beta * s[k][jj] * h[k][jj] );
                        
                        DP[ind]  += 2 * beta * s[k][rr] * s[k][jj] * ( 1./denr + 1./denj );
                        DDP[ind] += - 4 * beta * beta * ( 1./denrr + 1./denjj);
                    }
                
                }
                
            }
            
        }
        
        //al tempo = 1 calcolo il PSL gain per tutti i link
        if (time == 1){
            for (i = 0; i < Npar; i++){
                PSLgain[i] = - 0.5 * DP[i] * DP[i] / (L * DDP[i]);
            }
            
            
            
            for (i = 0; i < Npar; i ++)
                index_J[i] = i;
            
            //Lambda function to order the elements of J in a descendent order according to the values of PSLgain
            sort(index_J.begin(), index_J.end(), [&](const int& a, const int& b) { return (PSLgain[a] > PSLgain[b]); });
            
            //al tempo 1 ordino tutti gli PSL e poi mi prendo i primi aN couplings e li metto dentro il vettore vect
            //vect is made by the pair (ind, PSLgain)
            for (i = 0; i < aN; i ++){
                ind = index_J[i];
                map_ind[ind] = i;
                vect[i].first  = ind;
                vect[i].second = PSLgain[ind];
            }
            
            cnt_ind = aN;
            
        }
        //se al tempo t ho aggiornato il coupling (i-j),
        //al tempo t+1 ricalcolo solo il PSL gain per i link (i-k) e (j-k), con k=1,..,N, che coinvolgono i siti i e j
        if (time > 1){
            
            size_vect = vect.size();
            
            //compute the mean value of the good couplings
            mean_vect = 0.;
            for (int i=0; i<size_vect; i++)
                mean_vect += vect[i].second;
            mean_vect/=size_vect;
            
            
            //get the value of the PSL gain for the couplings that we updated in the previous time step and use it to define a threshold
            //all the links whose PSLgain is smaller than this threshold after the following updating will be excluded from the vector of good couplings
            //in order to kill the proliferation of couplings in this vector. this is done only if flag_opt is set to 1, otherwise th is useless.
            double th = a_th * PSLgain[activated_J[activated_J.size()-1]];
            
            //ri-calcolo i PSL gain per i link i-k e j-k
            int size_tu = ind_to_update.size();
            
            for (i = 0; i < size_tu; i++){
                i_tu = ind_to_update[i];
                PSLgain[i_tu] = - 0.5 * DP[i_tu] * DP[i_tu] / (L * DDP[i_tu]);
                PSL_u = PSLgain[i_tu];
                
                //check if il link di cui devo fare updating della PSL appartiene al vettore size_good_couplings.
                //if it does, we update it value
                if (map_ind[i_tu] >= 0){
                    vect[map_ind[i_tu]].second = PSL_u;
                }
                
                //if it did not but its updated value is found to be large, we insert it in the vector good_couplings
                if (map_ind[i_tu] < 0 && PSL_u > mean_vect){
                    vect.push_back( make_pair(i_tu,PSL_u) );
                    map_ind[i_tu] = cnt_ind;
                    cnt_ind ++;
                }
                
            }
            
            cout << "delete link" << endl;
            
            // erase the link that we have just updated from the good coupling
            for (int i = 0; i<Nup; i++){
                vect.erase(vect.begin());
                cout << "ind: " << activated_J[activated_J.size()-1 - i] << " PSL gain ---> " << PSLgain[activated_J[activated_J.size()-1 - i]] << endl;
                map_ind[activated_J[activated_J.size()-1 - i]] = -1;
            }
            
            if(a_th != 0){
                //erase other links that turns out to be smaller than the threshold
                //we first set  map_ind = -1 on these links
                for (i = 0; i < size_tu; i++){
                    i_tu = ind_to_update[i];
                    PSL_u = PSLgain[i_tu];
                    if (map_ind[i_tu] >= 0 && PSL_u < th){
                        cout << "ind: " << i_tu << " PSL gain ---> " << PSLgain[i_tu] << endl;
                        map_ind[i_tu] = -1;
                    }
                }
                //and then we erase them
                int j = 0;
                for (int i=0; i<vect.size(); i++){
                    
                    if(map_ind[vect[j].first] > 0)
                        j++;
                    if(map_ind[vect[j].first] < 0){
                        vect.erase(vect.begin()+j);
                    }
                }
            }
            
      
            size_vect = vect.size();

            cout << "SIZE VECTOR GOOD COUPLINGS ---------> " << size_vect << endl;
            cout << endl;
            
            if (size_vect <= 0){
                cout << "time: " << time << endl;
                cout << "NO MORE COUPLINGS IN THE VECTOR OF GOOD COUPLINGS" << endl;
                break;
            }
            
            //Lambda function to order the elements of vect in a descendent order according to the values of vect.second
            sort(vect.begin(), vect.end(), [&](const pair <int,double>& a, const pair <int,double>& b) { return a.second > b.second; });
            
            //redefine a map between the indeces of links contained inside the good couplings and N, given by their order
            for(int i=0; i<size_vect; i++){
                map_ind[vect[i].first] = i;
            }
            
            //set in index_J the value of the ordered links
            for(int i=0; i<Nup; i++)
                index_J[i] = vect[i].first;
                
        }
        
        cout << "time " << time << "-----------------------------------------------------------------------------------------------------------------" << endl;
        
        //update only the couplings whose updating gives the largest contribution the pseudo-likelihood

        for (i = 0; i < Nup; i ++){
            
            ind = index_J[i];

            cout << "update link: " << map[ind][0] << " " << map[ind][1] << endl;
            nodes_to_update[i]   = map[ind][0];
            nodes_to_update[i+1] = map[ind][1];
            
            J[ind]      = J[ind] - eta * DP[ind] / DDP[ind];

            //create the vector of neighbours of each node
            if(flag_J[ind] == 0){
                activated_J.push_back(ind);
                neigh[map[ind][0]].push_back(map[ind][1]);
                neigh[map[ind][1]].push_back(map[ind][0]);
            }

            flag_J[ind] = 1;
        }
        
        K = activated_J.size();
        
        //vector of the activated J that we use to optimize the PSL in an exact way.
        //we count (ij) and (ji) only once.
        column_vector J_start;
        J_start.set_size(K);
        
        //vector of the indices of the activated J.
        //in order to associate to (ij) and (ji) the same index, v_ind takes indices 0<q2<2K and give an index 0<q<K.
        std::vector <int> v_ind;
        v_ind.resize(2*K);
        
        //flag is a vector that allows to not count links twice.
        std::vector <int> flag ;
        flag.resize(Npar, 0);
        
        //to each ind, we associate a number 0<q<K that we use as the indeces of J_start.
        //(ij) and (ji) have different q2 but are associated to the same ind, that in turn is associated to the same index q.
        std::vector <int> v_q ;
        v_q.resize(Npar, 0);
        
        q  = 0;
        q2 = 0;
        
        
        for (r = 0; r < N; r ++){
            
            c = neigh[r].size();
            
            for (j = 0; j < c; j ++){
                
                if( neigh[r][j] < r){
                    jj    = r;
                    rr    = neigh[r][j];
                }
                if( neigh[r][j] > r){
                    jj    = neigh[r][j];
                    rr    = r;
                }
                
                ind = 0.5 * ( N * (N-1) - (N-rr)*(N-rr-1) ) + jj - (rr+1);
                
                v_ind[q2] = ind;
                q2++;
                
                if (flag[ind] == 0){
                    v_q[ind]   = q;
                    flag[ind]  = 1;
                    J_start(q) = J[ind];
                    q++;
                }
   
            }
            
        }
        
        //perform an exact optimzation over the activated variables
        
        if (time > 1){
            
            /*
            find_min_using_approximate_derivatives(lbfgs_search_strategy(2*K),
                                               objective_delta_stop_strategy(1e-1),
                                               PSL(N, L, beta, K, neigh, v_ind, v_q), J_start, 0);
            */
            
            double alpha     = 0.0001;
            double precision = 0.00001;     //set the convergence threshold in the Gradient Ascent: it is the threshold in the relative change of the PSL below which we stop the iteration
            int max_it       = 100;         //cambiare max_it non da risultati perche' in genere GradientAscent converge con pochi passi: nel caso in cui Precision e' 10^-5, al piu converge in ~25 passi con alpha=0.0001.
            
            GradientAscent(N, L, beta, K, alpha, max_it, precision, neigh, v_ind, v_q, J_start, numberMB);
            
        }

        
        //after the exact optimization, we copy the optimized vector of couplings into the vector J
        
        q  = 0;
        q2 = 0;
        
        fill(flag.begin(), flag.end(), 0);
        
        cout << "print inferred graph: " << endl;

        
        for (r = 0; r < N; r ++){
            
            c = neigh[r].size();
            
            for (j = 0; j < c; j ++){
                
                if( neigh[r][j] < r){
                    jj    = r;
                    rr    = neigh[r][j];
                }
                if( neigh[r][j] > r){
                    jj    = neigh[r][j];
                    rr    = r;
                }
                
                ind = 0.5 * ( N * (N-1) - (N-rr)*(N-rr-1) ) + jj - (rr+1);
                
                if (flag[ind] == 0){
                    
                    J[ind] = J_start(q);
                    flag[ind]  = 1;
                    q++;
                }
                
                cout << neigh[r][j] << " " << J[ind] << " ";

            }
            
            cout << endl;
            
        }
        
        /*
        if (time == T){
            cout << endl;
            cout << PSL(N, L, beta, neigh, v_ind, v_q)(J_start)/L << endl;
        }
        */
       
        
        PSL_c = 0.;
        double hh = 0.;
        
        std::vector <double> PSL_;
        PSL_.resize(N, 0.);
        
        //compute the PSL
        for (r = 0; r < N; r ++){
            
            //run over all the dataset
            for (k = 0; k < L; k ++){
                hh = 0.;
                
                //compute hr
                for (j = 0; j < N; j ++){
                    if( j < r){
                        jj    = r;
                        rr    = j;
                    }
                    if( j > r){
                        jj    = j;
                        rr    = r;
                    }
                    if( j == r)
                        continue;
                    
                    ind_j = 0.5 * ( N * (N-1) - (N-rr)*(N-rr-1) ) + jj - (rr+1);
                    JJ    = J[ind_j];
                    hh    += JJ * s[k][j];
                    
                }
                
                PSL_[r] +=log ( 1./( 1. + exp ( -2. * beta * s[k][r] * hh ) ) );
            }
            
            PSL_c += PSL_[r];
            
        }
        
        PSL_c /= L;
        
        s2_r = 0;
        
        for (r = 0; r < N; r ++)
            s2_r += pow( (PSL_[r] - PSL_c), 2.);
        
        s2_r /= L;
        
        
        int cnt_J0 = 0, cnt_J1 = 0;
        int FP = 0, TP = 0, TN = 0;
        
        if(test_mode == 1){
            
            Num = 0.;
            Den = 0.;
            int N1 = N-1;
            
            for(int r = 0; r < N1; r ++){
                
                for (int j = r+1; j < N; j ++){
                    
                    jj  = j;
                    rr  = r;
                    ind = 0.5 * ( N * (N-1) - (N-rr)*(N-rr-1) ) + jj - (rr+1);
                    
                    Num += pow( (J[ind] - J_true[ind]), 2);
                    Den += pow(J[ind], 2);
                    
                    if ( J_true[ind] == 0 ){
                        cnt_J0 ++;
                        if ( abs(J[ind]) > 0 )
                            FP ++;
                        if ( abs(J[ind]) == 0)
                            TN ++;
                    }
                   
                    if ( abs(J_true[ind]) == 1 ){
                        cnt_J1 ++;
                        if ( J[ind]*J_true[ind] > 0 )
                            TP ++;
                    }
                    
                }
                
            }
            
        }
        
        BIC = L * PSL_c - 0.5 * K * log(L);
        
        s_PSL = sqrt(s2_r / L);
        s_BIC = sqrt(L * s2_r);
        
        
        cout << "time:                 " << time << endl;
        cout << "PSL:                  " << time << " " << PSL_c  << " " << s_PSL << endl;
        cout << "BIC:                  " << time << " " << BIC  << " " << s_BIC << endl;
        cout << "# of activated links: " << K    << endl;
        
        if (test_mode == 1){
            cout << "eps:                  " << sqrt(Num/Den) << endl;
            cout << "FPR-TNR-TPR:          " << (double)FP/cnt_J0 << " " << (double)TN/cnt_J0 << " " << (double)TP/cnt_J1 << endl;
        }
        //this is the error on the reconstruction matrix, that of course we can estimate
        //only if we already know the true set of J's and that, in realistic cases,
        //cannot be computed
        
        cout <<  "-----------------------------------------------------------------------------------------------------------------" << endl;

        cout << endl;
        cout << endl;
        cout << endl;
        
        

    }
    
}

void GradientAscent(int N, int L, double beta, int K, double alpha, int max_it, double precision, std::vector < std::vector <int> >& neigh, std::vector <int>& v_ind, std::vector <int>& v_q, column_vector & J, int numberMB){

    int r, k, j;
    int jj, c;

    double JJ;
    double PSL, previous_PSL;
    
    std::vector <double> h;
    h.resize(N,0.);
    
    std::vector <double> DP;
    DP.resize(K,0.);
    
    std::vector <double> flag;
    flag.resize(K,0.);
    
    int q, q2, saved_q2;
    double denr, denjj;
    
    int sizeMB = L/numberMB;
    
    double invsizeMB = 1./sizeMB;
    
    int z = 0;

    for (int it = 0; it < max_it; it++){
    
        fill(DP.begin(), DP.end(), 0.);
        PSL = 0.;
        
        for (k = z*sizeMB; k < (z+1)*sizeMB; k++){
            
            q2 = 0;
            saved_q2 = q2;
            
            fill(h.begin(), h.end(), 0.);
            fill(flag.begin(), flag.end(), 0.);
            
            for (r = 0; r < N; r ++){
                
                c = neigh[r].size();
                saved_q2 = q2;
                q2 = saved_q2;
                
                for (j = 0; j < c; j ++){
                    
                    jj      =  neigh[r][j];
                    JJ      =  J(v_q[v_ind[q2]]);
                    h[r]    += JJ * a[k][jj];
                    
                    q2++;
                    
                }
                
                q2  = saved_q2 + c;
                
                PSL += log ( 1./( 1. + exp ( -2 * beta * a[k][r] * h[r] ) ) );
                
            }
            
            
            q  = 0;
            q2 = 0;
            saved_q2 = q2;
            
            for (r = 0; r < N; r ++){
                
                
                c = neigh[r].size();
                saved_q2 = q2;
                q2 = saved_q2;
                
                for (j = 0; j < c; j ++){
                    
                    jj      =  neigh[r][j];
                    JJ      =  J(v_q[v_ind[q2]]);
                    
                    q = v_q[v_ind[q2]];
                    
                    
                    denr  = 1 + exp( 2 * beta * a[k][r]  * h[r] );
                    denjj = 1 + exp( 2 * beta * a[k][jj] * h[jj] );
                    
                    if (flag[q] == 0){
                        DP[q]    += 2 * beta * a[k][r] * a[k][jj] * ( 1./denr + 1./denjj );
                        flag[q]  = 1;
                    }
                    
                    q2++;
                    
                }
                
                q2  = saved_q2 + c;
                
            }
            
        }
        
        z = (z+1) % (L/sizeMB);
        
        for (q = 0; q < K; q ++){
            J(q) += alpha * L * invsizeMB * DP[q];
        }
        
        if ( abs ( (PSL - previous_PSL) / PSL ) < precision ){
            cout << "IT ENDS AT " << it << endl;
            break;
        }
        
        
        previous_PSL = PSL;
        
    }
    

}

