#include "sp.h"
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <map>
#include <set>

#include "ilcplex/ilocplex.h"

ILOSTLBEGIN

using namespace std;


int SP::get_k()
{
	return k;
}

int SP::get_fullk()
{
	return fullk;
}

void SP::set_timelimit(double _cplextl)
{
    cplextl = _cplextl;
}


void SP::solve_minsummin(int s, int t, int H)
{
    double bigM = 0;
	for (int j=0; j<k; ++j)
    {
        double scensum = 0;
		for (int i=0; i<m; ++i)
            scensum += raw[j][i];
        bigM = max(bigM, scensum);
    }
    
    
    IloEnv env;
	IloModel model(env);
			
	vector<vector<IloNumVar> > cplexx(H);
    for (int h=0; h<H; ++h)
    {
        cplexx[h].resize(m);
        for (int i=0; i<m; ++i)
            cplexx[h][i] = IloNumVar(env, 0, 1, ILOBOOL);
    }
    
    vector<vector<IloNumVar> > cplexa(H);
    for (int h=0; h<H; ++h)
    {
        cplexa[h].resize(k);
        for (int j=0; j<k; ++j)
            cplexa[h][j] = IloNumVar(env, 0, 1, ILOBOOL);
    }
    
	vector<IloNumVar> cplexz(k);
    for (int j=0; j<k; ++j)
        cplexz[j] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
	
    
    //path constraints
    for (int h=0; h<H; ++h)
        for (int i=0; i<n; ++i)
        {
            IloExpr con(env);
            for (int j=0; j<g.out_edges[i].size(); ++j)
                con += cplexx[h][g.out_edges[i][j]];
            for (int j=0; j<g.in_edges[i].size(); ++j)
                con -= cplexx[h][g.in_edges[i][j]];
            if (i == s)
                model.add(con == 1);
            else if (i == t)
                model.add(con == -1);
            else
                model.add(con == 0);
        }
        
        
    //objective choice constraints
    for (int h=0; h<H; ++h)
        for (int j=0; j<k; ++j)
        {
            IloExpr con(env);
            for (int i=0; i<m; ++i)
                con += raw[j][i] * cplexx[h][i];
            model.add(con <= cplexz[j] + bigM*(1-cplexa[h][j]));
        }
        
    //active choice constraint
    for (int j=0; j<k; ++j)
    {
        IloExpr con(env);
        for (int h=0; h<H; ++h)
            con += cplexa[h][j];
        model.add(con == 1);
    }
    
    IloExpr obj(env);
    for (int j=0; j<k; ++j)
        obj += cplexz[j];
	model.add(IloMinimize(env, obj));
			
	IloCplex cplex(model);
	
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
    if (cplextl>-0.5)
        cplex.setParam(IloCplex::TiLim, cplextl);

	double start = clock();
	bool result = cplex.solve();			
	double time = (clock() - start)/CLOCKS_PER_SEC;

    vector<vector<double> > x(H);
    for (int h=0; h<H; ++h)
    {
        x[h].resize(m);
        for (int i=0; i<m; ++i)
            x[h][i] = round(cplex.getValue(cplexx[h][i]));
    }
    
    env.end();
    
    
    
    //evaluate
    
    double avperf = 0;
    
    for (int j=0; j<k; ++j)
    {
        double minval = bigM;
        for (int h=0; h<H; ++h)
        {
            double objval = 0;
            for (int i=0; i<m; ++i)
                objval += raw[j][i] * x[h][i];

            if (objval < minval)
                minval = objval;
        }
        
        avperf += minval;
    }
    avperf /= k;
    
    double outavperf = 0;
    for (int j=k; j<fullk; ++j)
    {
        double minval = bigM;
        for (int h=0; h<H; ++h)
        {
            double objval = 0;
            for (int i=0; i<m; ++i)
                objval += fullraw[j][i] * x[h][i];

            if (objval < minval)
                minval = objval;
        }
        
        outavperf += minval;
    }
    outavperf /= (fullk - k);
    
    cout<<"msm;"<<H<<";vals;"<<avperf<<";"<<outavperf<<"\n";        
}


void SP::solve_explainable_greedy(int s, int t, int H)
{
    if (H==2)
        solve_explainable_greedy2(s,t);
    else if (H==4)
        solve_explainable_greedy4(s,t);
    else
        cout<<"Choice of H not implemented.\n";
    
}


void SP::solve_explainable_greedy2(int s, int t)
{   
    int L=2;
    
    double bestinav = 999999;
    double bestoutav = 999999;
    
    int bestd1 = 0;
    double bestb1val = 0;
    
    vector<vector<double> > bestx;
    
    for (int d1=0; d1<m; ++d1)  
    //for (int d1=m-mignore; d1<m; ++d1) //use if only meta-data should be applied
    {       
        set<double> sb1vals;
        for (int j=0; j<k; ++j)
            sb1vals.insert(raw[j][d1]);
    
        //only single value exists
        if (sb1vals.size() == 1)
            continue;
 
        vector<double> b1vals;
        
        for (set<double>::iterator it=sb1vals.begin(); next(it)!=sb1vals.end(); ++it)
            b1vals.push_back((*it + *(next(it))) / 2);
                    
        for (int b1=0; b1<b1vals.size(); ++b1)
        {
            if (d1 < m-mignore && rand()%100 < bignore)
                continue;
            
            //solve shortest path for each class
            vector<vector<double> > avcost(2);
            avcost[0].resize(m,0);
            avcost[1].resize(m,0);
            
            vector<bool> typeexists(2,false);
            for (int j=0; j<k; ++j)
            {
                int type;
                if (raw[j][d1] <= b1vals[b1])
                    type = 0;
                else
                    type = 1;
                    
                typeexists[type] = true;
                for (int i=0; i<m; ++i)
                    avcost[type][i] += raw[j][i];
            }
                    
            for (int type=0; type<L; ++type)
                if (!typeexists[type])
                {
                    //take overall average scenario
                    for (int j=0; j<k; ++j)
                        for (int i=0; i<m; ++i)
                            avcost[type][i] += raw[j][i];
                }
            
            vector<vector<double> > curx(L);
            
            for (int type=0; type<L; ++type)
            {
                curx[type].resize(m,0);
                
                //solve shortest path
                IloEnv env;
                IloModel model(env);

                vector<IloNumVar> cplexx(m);
                for (int i=0; i<m; ++i)
                    cplexx[i] = IloNumVar(env, 0, 1, ILOBOOL);
        
                //path constraints
                for (int i=0; i<n; ++i)
                {
                    IloExpr con(env);
                    for (int j=0; j<g.out_edges[i].size(); ++j)
                        con += cplexx[g.out_edges[i][j]];
                    for (int j=0; j<g.in_edges[i].size(); ++j)
                        con -= cplexx[g.in_edges[i][j]];
                    if (i == s)
                        model.add(con == 1);
                    else if (i == t)
                        model.add(con == -1);
                    else
                        model.add(con == 0);
                }           
                
                IloExpr obj(env);
                for (int i=0; i<m-mignore; ++i)
                    obj += avcost[type][i]*cplexx[i];
                model.add(IloMinimize(env, obj));
                
                IloCplex cplex(model);
                            
                cplex.setOut(env.getNullStream());
                cplex.setParam(IloCplex::Threads, 1);
                                    
                bool result = cplex.solve();			
                
                for (int i=0; i<m-mignore; ++i)
                    if (cplex.getValue(cplexx[i]) > 0.5)
                        curx[type][i] = 1;
                
                env.end();
                
            }//end type
                    
            double curav=0;
            int choice0 = 0;
            int choice1 = 0;
            
            for (int j=0; j<k; ++j)        
            {
                //calculate decision
                int lchoice;
                if (raw[j][d1] <= b1vals[b1])
                {
                    lchoice = 0;
                    ++choice0;
                }
                else
                {
                    lchoice = 1;
                    ++choice1;
                }
                    
                double objval = 0;
                for (int i=0; i<m; ++i)
                    objval += raw[j][i] * curx[lchoice][i];
                curav += objval;                                
                }
                    
            curav /= k;
            
            if (curav < bestinav)
            {
                //update best
                bestinav = curav;
                bestx = curx;
                
                bestd1 = d1;
                bestb1val = b1vals[b1];
                        
                //evaluate out sample
                bestoutav = 0;
                for (int j=k; j<fullk; ++j)
                {
                    int lchoice;
                    if (fullraw[j][d1] <= b1vals[b1])
                        lchoice = 0;
                    else
                        lchoice = 1;
                        
                    double objval = 0;
                    for (int i=0; i<m; ++i)
                        objval += fullraw[j][i] * curx[lchoice][i];
                    bestoutav += objval;
                    
                }
            
                    bestoutav /= fullk -k;   
            }
                        
        }//end bval enumerate
    }//end dimension enumeration


    cout<<"g2d1;"<<bestd1<<";"<<bestb1val<<"\n";
        
    //evaluate bestx according to minmaxmin
    double inmm = 0;
    for (int j=0; j<k; ++j)
    {
        double minval = 999999;
        for (int l=0; l<L; ++l)
        {
            double objval = 0;
            for (int i=0; i<m; ++i)
                objval += raw[j][i] * bestx[l][i];

            minval = min(objval,minval);
        }
        
        inmm += minval;
    }
    inmm /= k;

    double outmm = 0;
    for (int j=k; j<fullk; ++j)
    {
        double minval = 999999;
        for (int l=0; l<L; ++l)
        {
            double objval = 0;
            for (int i=0; i<m; ++i)
                objval += fullraw[j][i] * bestx[l][i];

            minval = min(objval,minval);
        }
        
        outmm += minval;
    }
    outmm /= fullk - k;

    cout<<"ex;"<<L<<";vals;"<<bestinav<<";"<<bestoutav<<";"<<inmm<<";"<<outmm<<"\n";
       
}

void SP::solve_explainable_greedy4(int s, int t)
{
    //first two solutions, then four
    int L=2;
    
    double bestinav = 999999;
    double bestoutav = 999999;
    
    int bestd1 = 0;
    double bestb1val = 0;
    
    //for (int d1=m-mignore; d1<m; ++d1)
    for (int d1=0; d1<m; ++d1)
    {       
        set<double> sb1vals;
        for (int j=0; j<k; ++j)
            sb1vals.insert(raw[j][d1]);
 
        //only single value exists
        if (sb1vals.size() == 1)
            continue;
 
        vector<double> b1vals;
        for (set<double>::iterator it=sb1vals.begin(); next(it)!=sb1vals.end(); ++it)
            b1vals.push_back((*it + *(next(it))) / 2);
                    
        for (int b1=0; b1<b1vals.size(); ++b1)
        {
            if (d1 < m-mignore && rand()%100 < bignore)
                continue;
            
            //solve shortest path for each class
            vector<vector<double> > avcost(2);
            avcost[0].resize(m,0);
            avcost[1].resize(m,0);
            
            vector<bool> typeexists(2,false);
            for (int j=0; j<k; ++j)
            {
                int type;
                if (raw[j][d1] <= b1vals[b1])
                    type = 0;
                else
                    type = 1;
                    
                typeexists[type] = true;
                for (int i=0; i<m; ++i)
                    avcost[type][i] += raw[j][i];
            }
                    
            for (int type=0; type<L; ++type)
                if (!typeexists[type])
                {
                    //take overall average scenario
                    for (int j=0; j<k; ++j)
                        for (int i=0; i<m; ++i)
                            avcost[type][i] += raw[j][i];
                }
            
            vector<vector<double> > curx(L);
            
            for (int type=0; type<L; ++type)
            {
                curx[type].resize(m,0);
                            
                //solve shortest path
                IloEnv env;
                IloModel model(env);
                
                vector<IloNumVar> cplexx(m);
                for (int i=0; i<m; ++i)
                    cplexx[i] = IloNumVar(env, 0, 1, ILOBOOL);
                    
                //path constraints
                for (int i=0; i<n; ++i)
                {
                    IloExpr con(env);
                    for (int j=0; j<g.out_edges[i].size(); ++j)
                        con += cplexx[g.out_edges[i][j]];
                    for (int j=0; j<g.in_edges[i].size(); ++j)
                        con -= cplexx[g.in_edges[i][j]];
                    if (i == s)
                        model.add(con == 1);
                    else if (i == t)
                        model.add(con == -1);
                    else
                        model.add(con == 0);
                }                                
                                    
                IloExpr obj(env);
                for (int i=0; i<m; ++i)
                    obj += avcost[type][i]*cplexx[i];
                model.add(IloMinimize(env, obj));
                        
                IloCplex cplex(model);
                
                cplex.setOut(env.getNullStream());
                cplex.setParam(IloCplex::Threads, 1);
                
                bool result = cplex.solve();			
                
                for (int i=0; i<m; ++i)
                    if (cplex.getValue(cplexx[i]) > 0.5)
                        curx[type][i] = 1;
                        
                env.end();
            }//end type
                    
            double curav=0;
            for (int j=0; j<k; ++j)        
            {
                //calculate decision
                int lchoice;
                if (raw[j][d1] <= b1vals[b1])
                    lchoice = 0;
                else
                    lchoice = 1;
                    
                double objval = 0;
                for (int i=0; i<m; ++i)
                    objval += raw[j][i] * curx[lchoice][i];
                curav += objval;                                
                }
                    
            curav /= k;
                        
            if (curav < bestinav)
            {
                //update best
                bestinav = curav;
                bestd1 = d1;
                bestb1val = b1vals[b1];
                        
                //evaluate out sample
                bestoutav = 0;
                for (int j=k; j<fullk; ++j)
                {
                    int lchoice;
                    if (fullraw[j][d1] <= b1vals[b1])
                        lchoice = 0;
                    else
                        lchoice = 1;
                        
                    double objval = 0;
                    for (int i=0; i<m; ++i)
                        objval += fullraw[j][i] * curx[lchoice][i];
                    bestoutav += objval;
                    
                }
            
                    bestoutav /= fullk -k;   
            }
                        
        }//end bval enumerate
    }//end dimension enumeration

        
        
        
    //second dimension
    L=4;
    
    bestinav = 999999;
    bestoutav = 999999;
    
    int bestd2 = 0;
    double bestb2val = 0;
    
    vector<vector<double> > bestx;
    
    //for (int d2=m-mignore; d2<m; ++d2)
    for (int d2=0; d2<m; ++d2)
    {
        set<double> sb2vals;
        for (int j=0; j<k; ++j)
            sb2vals.insert(raw[j][d2]);
 
        //only single value exists
        if (sb2vals.size() == 1)
            continue;
 
        vector<double> b2vals;
        for (set<double>::iterator it=sb2vals.begin(); next(it)!=sb2vals.end(); ++it)
            b2vals.push_back( (*it + *(next(it))) / 2);
        
        for (int b2=0; b2<b2vals.size(); ++b2)
        {
            if (d2 < m-mignore && rand()%100 < bignore)
                continue;
            
            //solve shortest path for each class
            vector<vector<double> > avcost(L);
            for (int l=0; l<L; ++l)
                avcost[l].resize(m,0);
            
            vector<bool> typeexists(L,false);
            for (int j=0; j<k; ++j)
            {
                int type;
                if (raw[j][bestd1] <= bestb1val)
                {
                    if (raw[j][d2] <= b2vals[b2])
                        type = 0;
                    else
                        type = 1;
                }
                else
                {
                    if (raw[j][d2] <= b2vals[b2])
                        type = 2;
                    else
                        type = 3;
                }
                    
                typeexists[type] = true;
                for (int i=0; i<m; ++i)
                    avcost[type][i] += raw[j][i];
            }
                    
            for (int type=0; type<L; ++type)
                if (!typeexists[type])
                {
                    //take overall average scenario
                    for (int j=0; j<k; ++j)
                        for (int i=0; i<m; ++i)
                            avcost[type][i] += raw[j][i];
                }
            
            vector<vector<double> > curx(L);
            
            for (int type=0; type<L; ++type)
            {
                curx[type].resize(m,0);
                            
                //solve shortest path
                IloEnv env;
                IloModel model(env);
                
                vector<IloNumVar> cplexx(m);
                for (int i=0; i<m; ++i)
                    cplexx[i] = IloNumVar(env, 0, 1, ILOBOOL);
                    
                //path constraints
                for (int i=0; i<n; ++i)
                {
                    IloExpr con(env);
                    for (int j=0; j<g.out_edges[i].size(); ++j)
                        con += cplexx[g.out_edges[i][j]];
                    for (int j=0; j<g.in_edges[i].size(); ++j)
                        con -= cplexx[g.in_edges[i][j]];
                    if (i == s)
                        model.add(con == 1);
                    else if (i == t)
                        model.add(con == -1);
                    else
                        model.add(con == 0);
                }                                
                                    
                IloExpr obj(env);
                for (int i=0; i<m; ++i)
                    obj += avcost[type][i]*cplexx[i];
                model.add(IloMinimize(env, obj));
                        
                IloCplex cplex(model);
                
                cplex.setOut(env.getNullStream());
                cplex.setParam(IloCplex::Threads, 1);
                
                bool result = cplex.solve();			
                
                for (int i=0; i<m; ++i)
                    if (cplex.getValue(cplexx[i]) > 0.5)
                        curx[type][i] = 1;
                        
                env.end();
            }//end type
                    
            double curav=0;
            for (int j=0; j<k; ++j)        
            {
                //calculate decision
                int lchoice;
                if (raw[j][bestd1] <= bestb1val)
                {
                    if (raw[j][d2] <= b2vals[b2])
                        lchoice = 0;
                    else
                        lchoice = 1;
                }
                else
                {
                    if (raw[j][d2] <= b2vals[b2])
                        lchoice = 2;
                    else
                        lchoice = 3;
                }
                    
                double objval = 0;
                for (int i=0; i<m; ++i)
                    objval += raw[j][i] * curx[lchoice][i];
                curav += objval;                                
                }
                    
            curav /= k;
                        
            if (curav < bestinav)
            {
                //update best
                bestinav = curav;
                bestx = curx;
                
                bestd2 = d2;
                bestb2val = b2vals[b2];
                        
                //evaluate out sample
                bestoutav = 0;
                for (int j=k; j<fullk; ++j)
                {
                    int lchoice;
                    if (fullraw[j][bestd1] <= bestb1val)
                    {
                        if (fullraw[j][d2] <= b2vals[b2])
                            lchoice = 0;
                        else
                            lchoice = 1;
                    }
                    else
                    {
                        if (fullraw[j][d2] <= b2vals[b2])
                            lchoice = 2;
                        else
                            lchoice = 3;
                    }
                        
                    double objval = 0;
                    for (int i=0; i<m; ++i)
                        objval += fullraw[j][i] * curx[lchoice][i];
                    bestoutav += objval;
                    
                }
            
                    bestoutav /= fullk -k;   
            }
                        
        }//end bval enumerate
    }//end dimension enumeration        
        
        
        
    cout<<"g4d1;"<<bestd1<<";"<<bestb1val<<"\n";
    cout<<"g4d2;"<<bestd2<<";"<<bestb2val<<"\n";
        
    //evaluate bestx according to minmaxmin
    double inmm = 0;
    for (int j=0; j<k; ++j)
    {
        double minval = 999999;
        for (int l=0; l<L; ++l)
        {
            double objval = 0;
            for (int i=0; i<m; ++i)
                objval += raw[j][i] * bestx[l][i];

            minval = min(objval,minval);
        }
        
        inmm += minval;
    }
    inmm /= k;

    double outmm = 0;
    for (int j=k; j<fullk; ++j)
    {
        double minval = 999999;
        for (int l=0; l<L; ++l)
        {
            double objval = 0;
            for (int i=0; i<m; ++i)
                objval += fullraw[j][i] * bestx[l][i];

            minval = min(objval,minval);
        }
        
        outmm += minval;
    }
    outmm /= fullk - k;

    cout<<"ex;"<<L<<";vals;"<<bestinav<<";"<<bestoutav<<";"<<inmm<<";"<<outmm<<"\n";
           
}


void SP::initgrid(int _n, int _k)
{
	n = _n*_n;
	
	for (int i=0; i<_n; ++i)
	{
		for (int j=0; j<_n; ++j)
		{
			if(j+1<_n)
			{
				Edge e_1;
				e_1.s = i*_n+j;
				e_1.t = i*_n+(j+1);
				g.edges.push_back(e_1);
			}
			if(i+1<_n)
			{
				Edge e_2;
				e_2.s = i*_n+j;
				e_2.t = (i+1)*_n+j;
				g.edges.push_back(e_2);
			}
		}
	}
	
	m = g.edges.size();
	
	g.in_edges.resize(n);
	g.out_edges.resize(n);
	for (int i=0; i<m; ++i)
	{
		g.out_edges[g.edges[i].s].push_back(i);
		g.in_edges[g.edges[i].t].push_back(i);
	}
	
    //always generate 1000 test data
	k = _k;
    fullk = k+1000;
	raw.resize(k);
	fullraw.resize(fullk);

    //T is number of center scenarios
    int T=5;

    vector<vector<double> > mu(T);
    vector<vector<double> > delta(T);
    
    for (int t=0; t<T; ++t)
    {
        mu[t].resize(m);
        delta[t].resize(m);
        for (int j=0; j<m; ++j)
        {
            mu[t][j] = (rand()%201)/10.0 + 10.0; //10-30
            delta[t][j] = (rand()%101)/400.0; //0-0.25
        }
        
    }

	for (int i=0; i<k; ++i)
	{
		raw[i].resize(m);
		fullraw[i].resize(m);
		for (int j=0; j<m; ++j)
		{
			double rand_cost = 2*delta[i%T][j]*mu[i%T][j]*(rand()%101)/100.0 + (1-delta[i%T][j])*mu[i%T][j];
			raw[i][j] = rand_cost;
			fullraw[i][j] = rand_cost;
		}
	}
    
    random_shuffle(raw.begin(), raw.end());

	for (int i=k; i<fullk; ++i)
	{
		fullraw[i].resize(m);
		for (int j=0; j<m; ++j)
		{
			double rand_cost = 2*delta[i%T][j]*mu[i%T][j]*(rand()%101)/100.0 + (1-delta[i%T][j])*mu[i%T][j];
			fullraw[i][j] = rand_cost;
		}
	}
	
};


void SP::read(string graph_file, string scen_file)
{
	ifstream in(graph_file);
	string line;
	int readmode = 0;

	n = 0;
	
	while (!in.eof())
	{
		getline(in,line);
		
		// Ignore Comment and Empty Lines
		if (line=="" || line[0]=='#')
			continue;
            
		Edge e;
		
		//remove id
		size_t pos = line.find(",");
		line=line.substr(pos+1);
				
		//old id
		e.orig_id = atoi(line.c_str());
		pos = line.find(",");
		line=line.substr(pos+1);
				
		//s
		e.s = atoi(line.c_str());
		pos = line.find(",");
		line=line.substr(pos+1);
		if (e.s > n)
			n = e.s;
		
		//t
		e.t = atoi(line.c_str());
		pos = line.find(",");
		line=line.substr(pos+1);
		if (e.t > n)
			n = e.t;
		
		
		//length
		e.s_lon = atof(line.c_str());
		pos = line.find(",");
		line=line.substr(pos+1);
		
		e.s_lat = atof(line.c_str());
		pos = line.find(",");
		line=line.substr(pos+1);
		
		e.e_lon = atof(line.c_str());
		pos = line.find(",");
		line=line.substr(pos+1);
		
		e.e_lat = atof(line.c_str());
		pos = line.find(",");
		line=line.substr(pos+1);
		
		//Convert Lengths to miles
		e.length = sqrt(pow(69.0*(e.e_lon - e.s_lon),2) + pow(51.0*(e.e_lat - e.s_lat),2));
		
		g.edges.push_back(e);
	}
	
	in.close();
	
    //add artificial edges for extended data
    n += 1;
    {
        Edge e;
        e.s = n;
        e.t = n;
        
        //weekday
        g.edges.push_back(e);
        
        //time
        g.edges.push_back(e);        
    }
	
	n += 1;
	m = g.edges.size();
	
    cout<<n<<" nodes and "<<m<<" edges including extension\n"<<flush;
    
	g.in_edges.resize(n);
	g.out_edges.resize(n);
	for (int i=0; i<m; ++i)
	{
		g.out_edges[g.edges[i].s].push_back(i);
		g.in_edges[g.edges[i].t].push_back(i);
	}
	
	//read scenarios
	in.open(scen_file);
	
	
	map<int,int> ids;
	
	{
		getline(in,line);
		int coms = count(line.begin(), line.end(), ',');
		for (int i=0; i<coms; ++i)
		{
			size_t pos = line.find(",");
			line=line.substr(pos+1);
			ids[atoi(line.c_str())] = i;
		}
	}
		
	
	while (!in.eof())
	{
		getline(in,line);
		
		// Ignore Comment and Empty Lines
		if (line==""||line[0]=='#')
			continue;

		int orig_m = count(line.begin(), line.end(), ',') -2;
		
		vector<double> speed(orig_m);
		for (int i=0; i<orig_m; ++i)
		{
			//if statement differs from previous version
			if (i > 0.5)
			{
				size_t pos = line.find(",");
				line=line.substr(pos+1);
			}
			speed[i] = atof(line.c_str());
		}
		
		vector<double> scen(m);
        
        //two edges less than expected for artificial data
		for (int i=0; i<m-mignore; ++i)
		{

			//no data for edge, assume speed 20
			if (ids.count(g.edges[i].orig_id) == 0)
			{
				//cout<<"original id has no data\n";
				scen[i] = g.edges[i].length / 20;
			}
			else
			{
				double adapted_speed = speed[ids[g.edges[i].orig_id]] ;
				if(adapted_speed < 3)
				{
					adapted_speed = 3;
				}
				scen[i] = g.edges[i].length / adapted_speed;
			}
		}
        
        size_t pos = line.find(",");
        line=line.substr(pos+1);
        pos = line.find(",");
        string day=line.substr(0,pos);
        line=line.substr(pos+1);
        pos = line.find(",");
        string month=line.substr(0,pos);
        pos = line.find(",");
        line=line.substr(pos+1);

        pair<double,double> daytime = date_time_to_val(day, month, line);
        scen[m-2] = daytime.first;
        scen[m-1] = daytime.second;
		
		fullraw.push_back(scen);
	}
    
	in.close();
	
	fullk = fullraw.size();
    random_shuffle(fullraw.begin(), fullraw.end());
    k = fullk/2;
    for (int i=0; i<k; ++i)
        raw.push_back(fullraw[i]);
    
}

double SP::date_to_val(std::string day, std::string month)
{
    int iday=atoi(day.c_str());
    
    //Mr, Apr, Mai 2017
    
    if (month[1] == 'r')
    {
        //7 March 2017 = Tue
        iday = (iday+1)%7+1;
    }
    else if (month[0] == 'A')
    {
        //7 April 2017 = Fri
        iday = (iday+4)%7+1;
    }
    else
    {
        //7 May 2017 = Sun
        iday = (iday+6)%7+1;
    }
    
    return iday;
}

double SP::time_to_val(std::string time)
{
    int hour = atoi(time.c_str());
    size_t pos = time.find("_");
    time=time.substr(pos+1);
    
    int minute = atoi(time.c_str());
    pos = time.find("_");
    time=time.substr(pos+1);
    
    int second = atoi(time.c_str());
    
    return second + 60*minute + 60*60*hour;
}

int SP::get_n()
{
    return n;
}

pair<double,int> SP::solve_average(int s, int t)
{
    vector<double> avcost(m,0);
    for (int j=0; j<k; ++j)
        for (int i=0; i<m; ++i)
            avcost[i] += raw[j][i]/k;
    
    IloEnv env;
    IloModel model(env);

    vector<IloNumVar> cplexx(m);
    for (int i=0; i<m; ++i)
        cplexx[i] = IloNumVar(env, 0, 1, ILOBOOL);

    //path constraints
    for (int i=0; i<n; ++i)
    {
        IloExpr con(env);
        for (int j=0; j<g.out_edges[i].size(); ++j)
            con += cplexx[g.out_edges[i][j]];
        for (int j=0; j<g.in_edges[i].size(); ++j)
            con -= cplexx[g.in_edges[i][j]];
        if (i == s)
            model.add(con == 1);
        else if (i == t)
            model.add(con == -1);
        else
            model.add(con == 0);
    }           
    
    IloExpr obj(env);
    for (int i=0; i<m-mignore; ++i)
        obj += avcost[i]*cplexx[i];
    model.add(IloMinimize(env, obj));
    
    IloCplex cplex(model);
                
    cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::Threads, 1);
                        
    bool result = cplex.solve();			
    
    if (!result)
    {
        env.end();
        return make_pair(-1,-1);
    }
    
    double objval = cplex.getObjValue();
    
    //nr of edges
    int nredges = 0;
    for (int i=0; i<m-mignore; ++i)
        if (cplex.getValue(cplexx[i])>0.5)
            ++nredges;
    
    env.end();    
    
    return make_pair(objval,nredges);
}

pair<double,double> SP::date_time_to_val(std::string day, std::string month, std::string time)
{

    int hour = atoi(time.c_str());
    size_t pos = time.find("_");
    time=time.substr(pos+1);
    
    int minute = atoi(time.c_str());
    pos = time.find("_");
    time=time.substr(pos+1);
    
    int second = atoi(time.c_str());
    
    double time_in_sec = second + 60*minute + 60*60*hour;
    int dayminus = 0;
    
    //7 hours time difference
    if (time_in_sec < 25200)
    {
        dayminus = 1;
        time_in_sec += 86400;
    }
    time_in_sec -= 25200;


    int iday=atoi(day.c_str())-dayminus;
    
    if (month[1] == 'r')
        iday = (iday+1)%7+1;
    else if (month[0] == 'A')
        iday = (iday+4)%7+1;
    else
        iday = (iday+6)%7+1;
    
    
    return make_pair(iday, time_in_sec);
    
}

double SP::solve_lb_in(int s, int t)
{
    double sumcost = 0;
    for (int j=0; j<k; ++j)
    {
    
        IloEnv env;
        IloModel model(env);

        vector<IloNumVar> cplexx(m);
        for (int i=0; i<m; ++i)
            cplexx[i] = IloNumVar(env, 0, 1, ILOBOOL);

        //path constraints
        for (int i=0; i<n; ++i)
        {
            IloExpr con(env);
            for (int j=0; j<g.out_edges[i].size(); ++j)
                con += cplexx[g.out_edges[i][j]];
            for (int j=0; j<g.in_edges[i].size(); ++j)
                con -= cplexx[g.in_edges[i][j]];
            if (i == s)
                model.add(con == 1);
            else if (i == t)
                model.add(con == -1);
            else
                model.add(con == 0);
        }           
    
        IloExpr obj(env);
        for (int i=0; i<m-mignore; ++i)
            obj += raw[j][i]*cplexx[i];
        model.add(IloMinimize(env, obj));
        
        IloCplex cplex(model);
                
        cplex.setOut(env.getNullStream());
        cplex.setParam(IloCplex::Threads, 1);
                            
        bool result = cplex.solve();			
        
        sumcost += cplex.getObjValue()/k;
    
        env.end();    
    }
    
    return sumcost;
}

double SP::solve_lb_out(int s, int t)
{
    double sumcost = 0;
    for (int j=k; j<fullk; ++j)
    {
    
        IloEnv env;
        IloModel model(env);

        vector<IloNumVar> cplexx(m);
        for (int i=0; i<m; ++i)
            cplexx[i] = IloNumVar(env, 0, 1, ILOBOOL);

        //path constraints
        for (int i=0; i<n; ++i)
        {
            IloExpr con(env);
            for (int j=0; j<g.out_edges[i].size(); ++j)
                con += cplexx[g.out_edges[i][j]];
            for (int j=0; j<g.in_edges[i].size(); ++j)
                con -= cplexx[g.in_edges[i][j]];
            if (i == s)
                model.add(con == 1);
            else if (i == t)
                model.add(con == -1);
            else
                model.add(con == 0);
        }           
    
        IloExpr obj(env);
        for (int i=0; i<m-mignore; ++i)
            obj += fullraw[j][i]*cplexx[i];
        model.add(IloMinimize(env, obj));
        
        IloCplex cplex(model);
                
        cplex.setOut(env.getNullStream());
        cplex.setParam(IloCplex::Threads, 1);
                            
        bool result = cplex.solve();			
        
        sumcost += cplex.getObjValue()/(fullk-k);
    
        env.end();    
    }
    
    return sumcost;
}





void SP::solve_explainable_minsummin_tree(int s, int t, int H)
{
    int L = pow(2,H);
    vector<vector<int> > S(H);
    
    if (H==1)
    {
        S[0] = {0};
    }
    else if (H==2)
    {
        S[0] = {2,3};
        S[1] = {1,3};
    }
    else
    {
        cout<<"H not implemented.\n";
        exit(0);
    }
    
        
    double bestcrossval = 99999;
    double fullwcperf = 99999;
    
    vector<vector<double> > bestx(L);
    for (int l=0; l<L; ++l)
        bestx[l].resize(m);
    
    //bigM
    double bigM = 0;
    double bigE = 0;
    for (int j=0; j<k; ++j)
    {
        double scensum = 0;
        for (int i=0; i<m; ++i)
        {
            scensum += raw[j][i];
            bigE = max(bigE, raw[j][i]);
        }
        bigM = max(bigM, scensum);
    }
    ++bigE;
        
        
    IloEnv env;
    IloModel model(env);
            
    vector<vector<IloNumVar> > cplexx(L);
    for (int l=0; l<L; ++l)
    {
        cplexx[l].resize(m);
        for (int i=0; i<m; ++i)
            cplexx[l][i] = IloNumVar(env, 0, 1, ILOBOOL);
    }
    
    vector<vector<IloNumVar> > cplexa(L);
    for (int l=0; l<L; ++l)
    {
        cplexa[l].resize(k);
        for (int j=0; j<k; ++j)
            cplexa[l][j] = IloNumVar(env, 0, 1, ILOBOOL);
    }
    
    vector<IloNumVar> cplexbval(H);
    for (int h=0; h<H; ++h)
        cplexbval[h] = IloNumVar(env, -bigE, bigE, ILOFLOAT);

    vector<vector<IloNumVar> > cplexp(H);
    for (int h=0; h<H; ++h)
    {
        cplexp[h].resize(m);
        for (int i=0; i<m; ++i)
            cplexp[h][i] = IloNumVar(env, 0, 1, ILOBOOL);
    }
            
    vector<IloNumVar> cplexz(k);
    for (int j=0; j<k; ++j)
        cplexz[j] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
    
    
        
    //path constraints
    for (int l=0; l<L; ++l)
        for (int i=0; i<n; ++i)
        {
            IloExpr con(env);
            for (int j=0; j<g.out_edges[i].size(); ++j)
                con += cplexx[l][g.out_edges[i][j]];
            for (int j=0; j<g.in_edges[i].size(); ++j)
                con -= cplexx[l][g.in_edges[i][j]];
            if (i == s)
                model.add(con == 1);
            else if (i == t)
                model.add(con == -1);
            else
                model.add(con == 0);
        }
        
        
    //objective choice constraints
    for (int l=0; l<L; ++l)
        for (int j=0; j<k; ++j)
        {
            IloExpr con(env);
            for (int i=0; i<m; ++i)
                con += raw[j][i] * cplexx[l][i];
            model.add(con <= cplexz[j] + bigM*(1-cplexa[l][j]));
        }
    
    //relate choice and rule constraints
    for (int h=0; h<H; ++h)
        for (int j=0; j<k; ++j)
        {
            IloExpr ysum(env);
            for (int i=0; i<m; ++i)
                ysum += raw[j][i] * cplexp[h][i];
            
            IloExpr asum(env);
            for (int l=0; l<S[h].size(); ++l)
                asum += cplexa[S[h][l]][j];
            
            model.add(ysum <= cplexbval[h] - 0.001 + (bigE+1)*asum);
            model.add(cplexbval[h] + 0.001 - (bigE+1)*(1-asum) <= ysum);
        }

        
    //active choice constraint
    for (int j=0; j<k; ++j)
    {
        IloExpr con(env);
        for (int l=0; l<L; ++l)
            con += cplexa[l][j];
        model.add(con == 1);
    }
    
        
    //number of coefficients constraint
    for (int h=0; h<H; ++h)
    {
        IloExpr con(env);
        for (int i=0; i<m; ++i)
            con += cplexp[h][i];
        model.add(con == 1);
    }
        
    IloExpr obj(env);
    for (int j=0; j<k; ++j)
        obj += cplexz[j];
    model.add(IloMinimize(env, obj));
            
    IloCplex cplex(model);
    
    cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::Threads, 1);
    if (cplextl>-0.5)
        cplex.setParam(IloCplex::TiLim, cplextl);

    double start = clock();
    bool result = cplex.solve();			
    double time = (clock() - start)/CLOCKS_PER_SEC;
    
    vector<double> bvals(H);
    for (int h=0; h<H; ++h)
        bvals[h] = cplex.getValue(cplexbval[h]);
        
    //take midpoint of interval
    for (int h=0; h<H; ++h)
        for (int i=0; i<m; ++i)
            if (cplex.getValue(cplexp[h][i]) > 0.5)
            {
                set<double> vals;
                for (int j=0; j<k; ++j)
                    vals.insert(raw[j][i]);
                vector<double> vvals;
                for (set<double>::iterator it=vals.begin(); it!=vals.end(); ++it)
                    vvals.push_back(*it);
                for (int j=0; j<vvals.size()-1; ++j)
                    if (vvals[j] < bvals[h] && bvals[h] < vvals[j+1])
                    {
                        bvals[h] = (vvals[j] + vvals[j+1])/2;
                        cout<<"exip"<<H<<"d"<<h+1<<";"<<i<<";"<<bvals[h]<<"\n";
                        break;
                    }
            }
            
    for (int l=0; l<L; ++l)
        for (int i=0; i<m; ++i)
            bestx[l][i] = round(cplex.getValue(cplexx[l][i]));
        
    
        
    double inav = 0;

    for (int j=0; j<k; ++j)
    {
        
        int lchoice ;
        vector<double> ysum(H,0);
        for (int h=0; h<H; ++h)
            for (int i=0; i<m; ++i)
                if (cplex.getValue(cplexp[h][i]) > 0.5)
                    ysum[h] += raw[j][i];
        
        if (H==1)
        {
            if (ysum[0] <= bvals[0])
                lchoice = 1;
            else
                lchoice = 0;
        }
        else if (H==2)
        {
            if (ysum[0] <= bvals[0] && ysum[1] <= bvals[1])
                lchoice = 0;
            else if (ysum[0] > bvals[0] && ysum[1] <= bvals[1])
                lchoice = 2;
            else if (ysum[0] <= bvals[0] && ysum[1] > bvals[1])
                lchoice = 1;
            else
                lchoice = 3;
        }
    
        double objval = 0;
        for (int i=0; i<m; ++i)
            objval += raw[j][i] * bestx[lchoice][i];
        inav += objval;
    }
    inav /= k;
        
    
    double outav = 0;

    for (int j=k; j<fullk; ++j)
    {
        int lchoice ;
        vector<double> ysum(H,0);
        for (int h=0; h<H; ++h)
            for (int i=0; i<m; ++i)
                if (cplex.getValue(cplexp[h][i]) > 0.5)
                    ysum[h] += fullraw[j][i];
        
        if (H==1)
        {
            if (ysum[0] <= bvals[0])
                lchoice = 1;
            else
                lchoice = 0;
        }
        else if (H==2)
        {
            if (ysum[0] <= bvals[0] && ysum[1] <= bvals[1])
                lchoice = 0;
            else if (ysum[0] > bvals[0] && ysum[1] <= bvals[1])
                lchoice = 2;
            else if (ysum[0] <= bvals[0] && ysum[1] > bvals[1])
                lchoice = 1;
            else
                lchoice = 3;
        }
    
        double objval = 0;
        for (int i=0; i<m; ++i)
            objval += fullraw[j][i] * bestx[lchoice][i];
        outav += objval;
    }
    outav /= (fullk-k);
    
    
    
    
    //evaluate bestx according to minmaxmin
    double inmmac = 0;
    for (int j=0; j<k; ++j)
    {
        double minval = 999999;
        for (int l=0; l<L; ++l)
        {
            double objval = 0;
            for (int i=0; i<m; ++i)
                objval += raw[j][i] * bestx[l][i];

            minval = min(objval,minval);
        }
        
        inmmac += minval;
    }
    inmmac /= k;
    
    
    //evaluate bestx according to minmaxmin
    double outmmac = 0;
    for (int j=k; j<fullk; ++j)
    {
        double minval = 999999;
        for (int l=0; l<L; ++l)
        {
            double objval = 0;
            for (int i=0; i<m; ++i)
                objval += fullraw[j][i] * bestx[l][i];

            minval = min(objval,minval);
        }
        outmmac += minval;
    }
    outmmac /= (fullk-k);
    
    
    
    env.end();

    cout<<"exip;"<<L<<";vals;"<<inav<<";"<<outav<<";"<<inmmac<<";"<<outmmac<<"\n";
    
    
}
