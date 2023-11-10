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

#define nump 8
#define cr 250.0
#define numpair (nump * (nump - 1))
#define inf 1000000

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





int
main(int   argc,
     char *argv[])
{
    fin = fopen("in.txt", "r");
    fout = fopen("inn.txt", "w");

    vector<point> locs = vector<point>(nump, point());

    
    
    
    try {
        double sum = 0.0;
        for (int iter = 0; iter < 100; iter++)
        {
            for (int i = 0; i < nump; i++) fscanf(fin, "%lf %lf", &locs[i].x, &locs[i].y);
            GRBEnv env = GRBEnv();

            GRBModel model = GRBModel(env);

            model.set(GRB_IntParam_NonConvex, 2);

            // Create variables

            GRBVar x[nump];
            GRBVar y[nump];
            GRBVar e[nump][nump];
            GRBVar z[nump][nump][nump][nump];

            GRBVar _z = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);

            for (int i = 0; i < nump; i++)
            {
                x[i] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
                y[i] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
            }

            for (int i = 0; i < nump; i++)
            {
                for (int j = 0; j < nump; j++)
                {
                    for (int k = 0; k < nump; k++)
                    {
                        for (int l = 0; l < nump; l++)
                        {
                            z[i][j][k][l] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                        }
                    }
                }
            }

            for (int i = 0; i < nump; i++)
            {
                for (int j = 0; j < nump; j++)
                {
                    e[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);;
                }
            }

            for (int i = 0; i < nump; i++)
            {
                for (int j = i + 1; j < nump; j++)
                {
                    model.addQConstr(x[i] * x[i] - 2 * x[i] * x[j] + x[j] * x[j] + y[i] * y[i] - 2 * y[i] * y[j] + y[j] * y[j] - cr * cr <= inf * (1 - e[i][j]));
                    model.addQConstr(-inf * e[i][j] <= x[i] * x[i] - 2 * x[i] * x[j] + x[j] * x[j] + y[i] * y[i] - 2 * y[i] * y[j] + y[j] * y[j] - cr * cr);
                }
            }

            for (int i = 0; i < nump; i++)
            {
                for (int j = i + 1; j < nump; j++)
                {
                    model.addConstr(e[i][j], GRB_EQUAL, e[j][i]);
                }
            }

            //for (int i = 0; i < nump; i++) model.addConstr(e[i][i], GRB_EQUAL, 1);


            for (int s = 0; s < nump; s++)
            {
                for (int d = s + 1; d < nump; d++)
                {
                    GRBLinExpr in[nump];
                    GRBLinExpr out[nump];
                    for (int i = 0; i < nump; i++)
                    {
                        in[i] = 0;
                        out[i] = 0;
                    }

                    for (int i = 0; i < nump; i++)
                    {
                        for (int j = 0; j < nump; j++)
                        {
                            in[i] += z[s][d][j][i];
                            out[i] += z[s][d][i][j];
                        }
                    }

                    for (int i = 0; i < nump; i++)
                    {
                        for (int j = 0; j < nump; j++)
                        {
                            model.addConstr(z[s][d][i][j], GRB_LESS_EQUAL, e[i][j]);
                        }
                    }

                    for (int i = 0; i < nump; i++)
                    {
                        if (i == s)
                        {
                            model.addConstr(in[i], GRB_EQUAL, 0);
                            model.addConstr(out[i], GRB_EQUAL, 2);
                        }
                        else if (i == d)
                        {
                            model.addConstr(in[i], GRB_EQUAL, 2);
                            model.addConstr(out[i], GRB_EQUAL, 0);
                        }
                        else
                        {
                            model.addConstr(in[i], GRB_EQUAL, out[i]);
                            model.addConstr(in[i], GRB_LESS_EQUAL, 1);

                        }
                    }
                }
            }

            for (int i = 0; i < nump; i++)
            {
                model.addQConstr(x[i] * x[i] - 2 * x[i] * locs[i].x + locs[i].x * locs[i].x + y[i] * y[i] - 2 * y[i] * locs[i].y + locs[i].y * locs[i].y <= _z);
            }



            // Set objective

            GRBQuadExpr obj = _z;
            model.setObjective(obj, GRB_MINIMIZE);

            // Add constraint: x + 2 y + 3 z >= 4


            // Optimize model



            model.optimize();




            double res = model.get(GRB_DoubleAttr_ObjVal);
            cout << sqrt(res) << endl;
            fprintf(fout, "%lf\n", sqrt(res));
            sum += sqrt(res);

        }

        fprintf(fout, "%lf\n", sum);
        cout << sum << endl;

      


    }
    catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Exception during optimization" << endl;
    }




    

    



    fclose(fout);
    fclose(fin);

    int var;
    cin >> var;

    return 0;
}
