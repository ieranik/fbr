/* Copyright 2020, Gurobi Optimization, LLC */

/* This example formulates and solves the following simple QP model:

     minimize    x^2 + x*y + y^2 + y*z + z^2 + 2 x
     subject to  x + 2 y + 3 z >= 4
                 x +   y       >= 1
                 x, y, z non-negative

   It solves it once as a continuous model, and once as an integer model.
*/

#include "gurobi_c++.h"
#include <iostream>
#include <vector>
using namespace std;

FILE* fin;
FILE* fout;

class point
{
public:
    double x, y;
    point(double xx = 0.0, double yy = 0.0)
    {
        x = xx;
        y = yy;
    }
};

int cr, nump, nume;
vector<point> locs;

void op_file()
{
    fin = fopen("input.txt", "r");
    fout = fopen("output.txt", "w");
    
    fscanf(fin,"%d",&cr);
    
    fscanf(fin, "%d", &nump);
    locs = vector<point>(nump, point());

    for (int i = 0; i < nump; i++)
    {
        fscanf(fin, "%lf %lf", &locs[i].x, &locs[i].y);
        //fprintf(fout, "%lf %lf\n", locs[i].x, locs[i].y);
    }

    fscanf(fin, "%d", &nume);
    
    
    


    //fprintf(fout, "%d\n", cr);
    //fprintf(fout, "%d\n", locs.size());
    //for (int i = 0; i < locs.size(); i++)
    //{
    //    fprintf(fout, "%f %f\n", locs[i].x, locs[i].y);
    //}

    //int edgecount = 0;
    //for (int i = 0; i < locs.size(); i++)
    //{
    //    edgecount += G[i].size();
    //}
    //edgecount /= 2;
    //edgecount += add_edges.size();

    //fprintf(fout, "%d\n", edgecount);


    //for (int i = 0; i < locs.size(); i++)
    //{
    //    for (int j = 0; j < G[i].size(); j++)
    //    {

    //        if (i < G[i][j]) fprintf(fout, "%d %d\n", i, G[i][j]);
    //    }
    //}


    //for (int i = 0; i < add_edges.size(); i++)
    //{
    //    fprintf(fout, "%d %d\n", add_edges[i].first, add_edges[i].second);
    //}







    //fin=fopen("input.txt","r");
    //fscanf(fin,"%d",&i);

}


int
main(int   argc,
     char *argv[])
{
    op_file();
    try {
        GRBEnv env = GRBEnv();

        GRBModel model = GRBModel(env);

        // Create variables

        GRBVar* x = new GRBVar[nump];
        GRBVar* y = new GRBVar[nump];
        GRBVar z = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);

        for (int i = 0; i < nump; i++)
        {
            x[i] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
            y[i] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }

        for (int i = 0; i < nume; i++)
        {
            int s, d;
            fscanf(fin, "%d %d", &s, &d);
            model.addQConstr(x[s] * x[s] - 2 * x[s] * x[d] + x[d] * x[d] + y[s] * y[s] - 2 * y[s] * y[d] + y[d] * y[d] <= cr * cr);
        }

        for (int i = 0; i < nump; i++)
        {
            model.addQConstr(x[i] * x[i] - 2 * x[i] * locs[i].x + locs[i].x * locs[i].x + y[i] * y[i] - 2 * y[i] * locs[i].y + locs[i].y * locs[i].y <= z);
        }


        
            
        

        // Set objective

        GRBQuadExpr obj = z;
        model.setObjective(obj, GRB_MINIMIZE);

        // Add constraint: x + 2 y + 3 z >= 4

        

       
        // Optimize model

        model.optimize();

        fprintf(fout, "%d %d\n", nump, cr);

        for (int i = 0; i < nump; i++)
        {
            fprintf(fout, "%lf %lf\n", locs[i].x, locs[i].y);
            fprintf(fout, "%lf %lf\n", x[i].get(GRB_DoubleAttr_X), y[i].get(GRB_DoubleAttr_X));
        }

        /*cout << x.get(GRB_StringAttr_VarName) << " "
             << x.get(GRB_DoubleAttr_X) << endl;
        cout << y.get(GRB_StringAttr_VarName) << " "
             << y.get(GRB_DoubleAttr_X) << endl;
        cout << z.get(GRB_StringAttr_VarName) << " "
             << z.get(GRB_DoubleAttr_X) << endl;*/

        cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

        // Change variable types to integer

        //x.set(GRB_CharAttr_VType, GRB_INTEGER);
        //y.set(GRB_CharAttr_VType, GRB_INTEGER);
        //z.set(GRB_CharAttr_VType, GRB_INTEGER);

        //// Optimize model

        //model.optimize();

        //cout << x.get(GRB_StringAttr_VarName) << " "
        //     << x.get(GRB_DoubleAttr_X) << endl;
        //cout << y.get(GRB_StringAttr_VarName) << " "
        //     << y.get(GRB_DoubleAttr_X) << endl;
        //cout << z.get(GRB_StringAttr_VarName) << " "
        //     << z.get(GRB_DoubleAttr_X) << endl;

        //cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    fclose(fout);
    fclose(fin);

    int var;
    cin >> var;

    return 0;
}
