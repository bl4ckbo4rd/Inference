#include "problems.h"


void f_Inference(int N, int L, int Ltot, int T, double beta, double rho_acc, double eta, int f, double rho, int w, double a_size, double a_th, int numberMB, int CLUSTER){


    ostringstream convert1;                                                              // stream used for the conversion      
    convert1 << beta;                                                                    // insert the textual representation of 'Number' in the characters in the stream      
    string s_beta = convert1.str();      
  
    ostringstream convert;                                                               // stream used for the conversion  
    convert << N;                                                                        // insert the textual representation of 'Number' in the characters in the stream  
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
        datafile = path + "beta=" + s_beta + "/data_" + size +"_2D-FB_" + s_rho + ".txt" ;
    
    if(f == 1)
        datafile = path + "beta=" + s_beta + "/data_" + size +"_2D-PB_" + s_rho + ".txt" ;

    if(f == 2)
        datafile = path + "beta=" + s_beta + "/data_" + size +"_ER_" + s_rho + ".txt" ;

    if(f == 3)
        datafile = path + "beta=" + s_beta + "/data_" + size +"_RR_" + s_rho + ".txt" ;
  
    if(f == 4)
    datafile = path + "beta=" + s_beta + "/data_" + size +"_Diamond_" + s_rho + ".txt" ;
    
    cout << datafile << endl;
    
     ifstream infile(datafile);     

    
    // If we couldn't open the output file stream for reading
    if (!infile)
    {
        cerr << "Uh oh, data.txt could not be opened for reading!" << endl;
        exit(1);
    }
    
    std::vector < std::vector <int> > s;
    
    s.resize(L);
    
    int i = 0;
    int j = 0;
    
    std::mt19937 gen_ber(w);
    std::bernoulli_distribution dis_ber{0.5};
   
    std::mt19937 gen_u(w);
    std::uniform_int_distribution<> dis_u{0, Ltot - 1};
    
    //flag is a vector that contains 1 in L (random) entries and 0 elsewhere
    std::vector <int> flag;
    flag.resize(Ltot,0);
    
    int cnt = 0;
    
    while(cnt < L){
        i = dis_u(gen_u);
        if(flag[i])
            continue;
        else{
            flag[i] = dis_ber(gen_ber);
            cnt += flag[i];
        }
    }
    
    string temp;
    
    i=0;
    
    while (std::getline(infile, temp)) {
        if (flag[i]){
            istringstream buffer(temp);
            std::vector<int> line( (istream_iterator<int>(buffer)), istream_iterator<int>() );
            s[j].resize(N);
            s[j] = line;
            j++;
        }
        
        i++;
        
    }
    
    int test_mode = 1;
    network_setup(N, L, T, beta, s, rho_acc, eta, f, rho, a_size, a_th, numberMB, CLUSTER, test_mode);

    
}

