#include "sp.h"

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

using namespace std;


int main_chicago(int argc, char* argv[]) //Experiment 3
//int main(int argc, char* argv[])
{	
	srand(atoi(argv[3]));
    
	SP sp;
	
    sp.read(argv[1], argv[2]);
    
    cout<<"k="<<sp.get_k()<<"\n";
    cout<<"fullk="<<sp.get_fullk()<<"\n";

    int s = rand()%(sp.get_n()-1);
    int t = rand()%(sp.get_n()-1);
                                    
    sp.set_timelimit(-1);
    
    bool highvalue = false;
    do
    {
        int pathlength = sp.solve_average(s,t).second;
        
        if (pathlength > 50)
        {
            highvalue = true;
            cout<<"nomedges;"<<pathlength<<"\n";
        }
        else
        {
            s = rand()%(sp.get_n()-1);
            t = rand()%(sp.get_n()-1);
        }
    }while(highvalue == false);
    
	vector<double> times;
	
    cout<<"path;"<<s<<";"<<t<<"\n"<<flush;
    
    double start=clock();
    sp.solve_minsummin(s,t,1);
    times.push_back((clock() - start)/CLOCKS_PER_SEC);
    
    cout<<flush;
    
    //start=clock();         
    //cout<<"lbin;"<<sp.solve_lb_in(s,t)<<"\n";
    //times.push_back((clock() - start)/CLOCKS_PER_SEC);
    //sp.set_timelimit(times.back());
    
    //cout<<flush;
    
    //start=clock();         
    //cout<<"lbout;"<<sp.solve_lb_out(s,t)<<"\n";
    //times.push_back((clock() - start)/CLOCKS_PER_SEC);
    //sp.set_timelimit(times.back());
    
    //cout<<flush;
    
    start=clock();         
    sp.solve_explainable_greedy(s,t,2);
    times.push_back((clock() - start)/CLOCKS_PER_SEC);
    sp.set_timelimit(times.back());
    
    cout<<flush;
    
    start=clock();         
    sp.solve_explainable_greedy(s,t,4);
    times.push_back((clock() - start)/CLOCKS_PER_SEC);
    sp.set_timelimit(times.back());
    
    cout<<flush;
    
    cout<<"times";
    for (int i=0; i<times.size(); ++i)
        cout<<";"<<times[i];
    cout<<"\n";
    
	
	return 0;
}




//int main_grid_ip(int argc, char* argv[]) //Experiment 1 and 2
int main(int argc, char* argv[])
{
	srand(atoi(argv[3]));
	
    for (int rep=0; rep<10; ++rep)
    {
    
        SP sp;
        
        int n = atoi(argv[1]);
        int k = atoi(argv[2]);
        sp.initgrid(n,k);
        
        int s = 0;
        int t = n*n-1;
        
        sp.set_timelimit(-1);
        vector<double> times;
        
        double start=clock();
        sp.solve_minsummin(s,t,1);
        times.push_back((clock() - start)/CLOCKS_PER_SEC);
        
        start=clock();         
        cout<<"lbin;"<<sp.solve_lb_in(s,t)<<"\n";
        times.push_back((clock() - start)/CLOCKS_PER_SEC);   
            
        start=clock();         
        cout<<"lbout;"<<sp.solve_lb_out(s,t)<<"\n";
        times.push_back((clock() - start)/CLOCKS_PER_SEC);
        
        //cout<<flush;
        
        //start=clock();         
        //sp.solve_explainable_greedy(s,t,2);
        //times.push_back((clock() - start)/CLOCKS_PER_SEC);
        //sp.set_timelimit(times.back());
        
        //start=clock();
        //sp.solve_minsummin(s,t,2);
        //times.push_back((clock() - start)/CLOCKS_PER_SEC);
        
        
        //start=clock();         
        //sp.solve_explainable_greedy(s,t,4);
        //times.push_back((clock() - start)/CLOCKS_PER_SEC);
        //sp.set_timelimit(times.back());
        
        //start=clock();
        //sp.solve_minsummin(s,t,4);
        //times.push_back((clock() - start)/CLOCKS_PER_SEC);
        
        start=clock();         
        sp.solve_explainable_greedy(s,t,2);
        times.push_back((clock() - start)/CLOCKS_PER_SEC);
        
        start=clock();         
        sp.solve_explainable_greedy(s,t,4);
        times.push_back((clock() - start)/CLOCKS_PER_SEC);
        
        start=clock();
        sp.solve_explainable_minsummin_tree(s,t,1);
        times.push_back((clock() - start)/CLOCKS_PER_SEC);
        
        start=clock();
        sp.solve_explainable_minsummin_tree(s,t,2);
        times.push_back((clock() - start)/CLOCKS_PER_SEC);
        
        
        cout<<"times";
        for (int i=0; i<times.size(); ++i)
            cout<<";"<<times[i];
        cout<<"\n"<<flush;
        
    }
    
	return 0;
}
