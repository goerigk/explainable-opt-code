#ifndef _H_SP_H
#define _H_SP_H

#include <vector>
#include <string>
#include <set>

struct Edge
{
	int s,t;
	int orig_id;
	double length;
	
	double s_lon,s_lat,e_lon,e_lat;
};

struct Graph
{
	std::vector<Edge> edges;
	
	std::vector<std::vector<int> > in_edges;
	std::vector<std::vector<int> > out_edges;
};

class SP
{
	public:
    
        //generators
        void initgrid(int _n, int _k); //create random grid graph
        void read(std::string graph_file, std::string scen_file);
        
        //solvers
        void solve_minsummin(int s, int t, int H);
        void solve_explainable_greedy(int s, int t, int H); //explainable greedy
        void solve_explainable_minsummin_tree(int s, int t, int H); //explainable IP
        
        std::pair<double,int> solve_average(int s, int t);
        double solve_lb_in(int s, int t);
        double solve_lb_out(int s, int t);
        
        //aux
        void set_timelimit(double _cplextl);
        int get_k();
		int get_fullk();
        int get_n();

    private:
		//graph data
		Graph g;
		int n,m;
	
        const int mignore = 0; //set parameter to exclude edges (if only meta-data should be used)
        const int bignore = 0; //set parameter if greedy should ignore some split values to improve computation times
    
		//scenarios
		int k;
		std::vector<std::vector<double> > raw;
		
		int fullk;
		std::vector<std::vector<double> > fullraw;
		std::vector<int> mark_fullraw;
        
        double cplextl;
        
        void solve_explainable_greedy2(int s, int t);
        void solve_explainable_greedy4(int s, int t);
        
        double date_to_val(std::string day, std::string month);
        double time_to_val(std::string time);
        
        std::pair<double, double> date_time_to_val(std::string day, std::string month, std::string time);
        
};


#endif
