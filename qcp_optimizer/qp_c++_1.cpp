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
#include <ctime>
#include <time.h>
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

double cr;
int nump, nume;
vector<point> locs;

void op_file()
{
    
    fscanf(fin, "%d", &nump);
    locs = vector<point>(nump, point());

    for (int i = 0; i < nump; i++)
    {
        fscanf(fin, "%lf %lf", &locs[i].x, &locs[i].y);
        //fprintf(fout, "%lf %lf\n", locs[i].x, locs[i].y);
    }

    fscanf(fin, "%d", &nume);

}


int
main(int   argc,
     char *argv[])
{
    fin = fopen("out256.txt", "r");
    fout = fopen("outt256.txt", "w");

    fscanf(fin, "%lf", &cr);

    clock_t s;
    clock_t e;
    clock_t t;

    t = 0;
    double sum = 0.0;
    
    

    for (int i = 0; i < 100; i++)
    {
        op_file();
        try {
            GRBEnv env = GRBEnv();

            GRBModel model = GRBModel(env);

            // Create variables

            GRBVar* x = new GRBVar[nump];
            GRBVar* y = new GRBVar[nump];

            s = clock();
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

            e = clock();

            t += e - s;

           /* fprintf(fout, "%d %d\n", nump, cr);

            for (int i = 0; i < nump; i++)
            {
                fprintf(fout, "%lf %lf\n", locs[i].x, locs[i].y);
                fprintf(fout, "%lf %lf\n", x[i].get(GRB_DoubleAttr_X), y[i].get(GRB_DoubleAttr_X));
            }*/


            double res = model.get(GRB_DoubleAttr_ObjVal);
            //cout << sqrt(res) << endl;
            fprintf(fout, "%lf\n", sqrt(res));
            sum += sqrt(res);

      


        }
        catch (GRBException e) {
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        }
        catch (...) {
            cout << "Exception during optimization" << endl;
        }




    }

    cout << "Distance: " << sum / (cr * 100) << endl;

    printf("It took me %d clicks (%f seconds).\n", t, ((float)t) / CLOCKS_PER_SEC);



    fclose(fout);
    fclose(fin);

    int var;
    cin >> var;

    return 0;
}
