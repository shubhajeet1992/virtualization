#include<stdio.h>
#include<conio.h>
#include<stdlib.h>
#include <iostream>
#include "matrix.h"
#include "engine.h"
#include"RCPtrueLib.h"
#include"RCPApproxLib.h"
#include"RCPBreakLib.h"
#include <math.h>
#include <ilcplex/ilocplex.h>

ILOSTLBEGIN

typedef IloArray<IloNumArray>    NumMatrix;
typedef IloArray<IloNumVarArray> NumVarMatrix;

int BFS(int S, int B, double kappa[3], double beta[3], double cost[10], double BSLocationX[10], double BSLocationY[10]);
int Algorithm1(int S, int B, double kappa[3], double beta[3], double cost[10], double BSLocationX[10], double BSLocationY[10]);
int DE(int S, int B, double kappa[3], double beta[3], double cost[10], double BSLocationX[10], double BSLocationY[10]);
void Algorithm2(int S, int B, double kappa[3], double beta[3], double BSLocationX[10], double BSLocationY[10]);


int main()
{

	std::cout << "Intializing Application... \n" << std::endl;
	if (!mclInitializeApplication(NULL, 0))
	{
		std::cout << "failed" << std::endl;
		_getch();
		return -1;
	}
	else
	{
		std::cout << "Done \n" << std::endl;
	}

	printf("Initializing Libraries... \n");
	if (!RCPBreakLibInitialize())
	{
		std::cout << "RCPBreakLib initialization failed" << std::endl;
		_getch();
		return -1;
	}
	else
	{
		std::cout << "RCPBreakLib initialization done .. \n" << std::endl;
	}

	if (!RCPApproxLibInitialize())
	{
		std::cout << "RCPApproxLib initialization failed" << std::endl;
		_getch();
		return -1;
	}
	else
	{
		std::cout << "RCPApproxLib initialization done .. \n" << std::endl;
	}

	if (!RCPtrueLibInitialize())
	{
		std::cout << "RCPLib initialization failed" << std::endl;
		_getch();
		return -1;
	}
	else
	{
		std::cout << "RCPLib initialization done .. \n" << std::endl;
	}
	int S = 3;
	int B = 10;
	int infeasible = 0;
	double kappa[] = { 1000000,256000, 64000 };
	double beta[] = { 0.95, 0.9, 0.85 };
	double cost[] = {200,100,200,200,100,300,200,100,200,200};
	double BSLocationX[] = { 1.6, 0.6, 1.1, 1.4, 1.8, 2, 1.1, 0.3, 0.3, 0.6 };
	double BSLocationY[] = { 1.7, 0.6, 1.7, 0.5, 1.9, 0.7, 0.4, 0.6, 1.3, 1 };
	 // // BFS
	infeasible = BFS(S, B, kappa, beta, cost, BSLocationX, BSLocationY);
	
	// // Running Algorithm 1.. 
	 infeasible = Algorithm1(S, B, kappa, beta, cost, BSLocationX, BSLocationY);

	 // DE..
	infeasible = DE(S, B, kappa, beta, cost, BSLocationX, BSLocationY);

	
	if (infeasible == 1) {
		std::cout << "Running Algorithm 2.. \n" << std::endl;
		Algorithm2(S, B, kappa, beta, BSLocationX, BSLocationY);
	}


	std::cout << "terminating..\n" << std::endl;
	RCPtrueLibTerminate();
	RCPApproxLibTerminate();
	RCPBreakLibTerminate();
	mclTerminateApplication();
	std::cout << "done. Press any key to return." << std::endl;
	_getch();
	return 0;

}




int BFS(int S, int B, double kappa[3], double beta[3], double cost[10], double BSLocationX[10], double BSLocationY[10])
{
	int infeasible = 0;
	double rcp1 = 0, rcp2 = 0, rcp3 = 0;
	double X[] = { 0,0,0,0,0,0,0,0,0,0 };
	mwArray in_X(1, 10, mxDOUBLE_CLASS);
	mwArray BSX(1, 10, mxDOUBLE_CLASS);
	BSX.SetData(BSLocationX, 10);
	mwArray BSY(1, 10, mxDOUBLE_CLASS);
	BSY.SetData(BSLocationY, 10);
	mwArray in_Delta1(1, 10, mxDOUBLE_CLASS);
	mwArray in_Delta2(1, 10, mxDOUBLE_CLASS);
	mwArray in_Delta3(1, 10, mxDOUBLE_CLASS);
	mwArray in_Xopt(1, 10, mxDOUBLE_CLASS);
	mwArray in_Delta1opt(1, 10, mxDOUBLE_CLASS);
	mwArray in_Delta2opt(1, 10, mxDOUBLE_CLASS);
	mwArray in_Delta3opt(1, 10, mxDOUBLE_CLASS);
	mwArray in_Kappa1((mxDouble)kappa[0]);
	mwArray in_Kappa2((mxDouble)kappa[1]);
	mwArray in_Kappa3((mxDouble)kappa[2]);
	double Xopt[10];
	double Xtemp[10];
	double delta[10][3];
	int BS_select = 0;
	double Cost_new, Cost = 1800;
	std::cout << "Running Brute Force.. \n" << std::endl;
	for (int i = 0; i < 10; i++) {
		X[i] = 1;
		BS_select = BS_select + 1;
		for (int j = i; j <= 10; j++) {
			std::copy(std::begin(X), std::end(X), std::begin(Xtemp));
			Xtemp[j] = 1;
			in_X.SetData(Xtemp, 10);
			//std::cout << "Xtemp is.. "<< in_X << std::endl;
			double points[10][3][5];
			for (int k = 0; k < 10; k++) {
				for (int l = 0; l < 3; l++) {
					mwArray in_Kappa((mxDouble)kappa[l]);
					points[k][l][0] = 0;
					points[k][l][1] = 0;
					points[k][l][2] = 0;
					points[k][l][3] = 0;
					points[k][l][4] = 0;
					points[k][l][5] = 0;
					if (Xtemp[k] != 0) {
						mwArray RCP_breakpoints;
						mwArray in_BSNo((mxDouble)(k + 1));
						RCP_breakpoint(1, RCP_breakpoints, in_Kappa, in_BSNo, in_X, BSX, BSY);
						//std::cout << "break points are.. " << RCP_breakpoints << std::endl;
						double *xp = new double[5];
						RCP_breakpoints.GetData(xp, 5);
						points[k][l][1] = *xp;
						points[k][l][2] = *(xp + 1);
						points[k][l][3] = *(xp + 2);
						points[k][l][4] = *(xp + 3);
						points[k][l][5] = *(xp + 4);
					}
				}
			}
			//std::cout << "out of breakpoint loop" << std::endl;

			// // solving Problem 2 in Cplex
			IloEnv env;
			try {
				IloInt m, n;
				IloModel model(env);

				IloInt S = 3;
				IloInt B = 10;
				IloInt No_break = 6;


				NumVarMatrix x(env, B);
				NumVarMatrix y(env, B);

				for (m = 0; m < B; m++) {
					x[m] = IloNumVarArray(env, S, 0.0, 1.0, ILOFLOAT);
					y[m] = IloNumVarArray(env, S, 0.0, 1.0, ILOFLOAT);
				}

				for (m = 0; m < B; m++) {
					model.add(IloSum(x[m]) <= 1);
				}


				for (m = 0; m < B; m++) {
					for (n = 0; n < S; n++) {
						IloNum break1 = beta[n];
						IloNum break2 = beta[n] - points[m][n][1];
						IloNum break3 = beta[n] - points[m][n][2];
						IloNum break4 = beta[n] - points[m][n][3];
						IloNum break5 = beta[n] - points[m][n][4];
						IloNum break6 = beta[n] - points[m][n][5];
						IloNum slope1 = break2 / 0.20;
						IloNum slope2 = (break6 - break5) / 0.20;
						model.add(y[m][n] == IloPiecewiseLinear(x[m][n],
							slope1, IloNumArray(env, 6, 0.0, 0.20, 0.40, 0.6, 0.8, 1.0),
							IloNumArray(env, 6, break1, break2, break3, break4, break5, break6), slope2));
					}
				}

				IloExpr obj(env);
				for (m = 0; m < B; m++) {
					obj += IloSum(y[m]);
				}
				model.add(IloMinimize(env, obj));
				obj.end();

				IloCplex cplex(env);
				cplex.extract(model);
				cplex.exportModel("transport.lp");
				cplex.solve();


				if (cplex.solve()) {
					env.out() << "Solution status: " << cplex.getStatus() << std::endl;
					for (m = 0; m < B; m++) {
						env.out() << "   " << m << ": ";
						for (n = 0; n < S; n++) {
							delta[m][n] = cplex.getValue(x[m][n]);
							env.out() << cplex.getValue(x[m][n]) << "\t";
						}
						env.out() << std::endl;
					}
					env.out() << "   Objective value = " << cplex.getObjValue() << std::endl;
				}
				else {
					std::cout << "No solution" << std::endl;
					//cplex.printTime();
				}


			}
			catch (IloException& e) {
				cerr << "ERROR: " << e.getMessage() << std::endl;
			}
			catch (...) {
				cerr << "Error" << std::endl;
			}
			env.end();



			// checking feasibility
			double delta1[] = { 0,0,0,0,0,0,0,0,0,0 };
			double delta2[] = { 0,0,0,0,0,0,0,0,0,0 };
			double delta3[] = { 0,0,0,0,0,0,0,0,0,0 };
			for (int i = 0; i < BS_select; i++)
			{
				delta1[i] = delta[i][0];
			}
			for (int i = 0; i < BS_select; i++)
			{
				delta2[i] = delta[i][1];
			}
			for (int i = 0; i < BS_select; i++)
			{
				delta3[i] = delta[i][2];
			}


			in_X.SetData(Xtemp, 10);
			in_Delta1.SetData(delta1, 10);
			in_Delta2.SetData(delta2, 10);
			in_Delta3.SetData(delta3, 10);
			mwArray RCP1;
			mwArray RCP2;
			mwArray RCP3;
			RCP_true(1, RCP1, in_Kappa1, in_X, in_Delta1,BSX,BSY);
			RCP_true(1, RCP2, in_Kappa2, in_X, in_Delta2, BSX, BSY);
			RCP_true(1, RCP3, in_Kappa3, in_X, in_Delta3, BSX, BSY);
			double *xp = new double[1];
			RCP1.GetData(xp, 1);
			rcp1 = *xp;
			RCP2.GetData(xp, 1);
			rcp2 = *xp;
			RCP3.GetData(xp, 1);
			rcp3 = *xp;
			//std::cout << RCP3 << std::endl;
			//std::cout << rcp3 << std::endl;
			if (rcp1 >= beta[0] && rcp2 >= beta[1] && rcp3 >= beta[2]) {
				for (int a = 0; a < 10; a++) {
					Cost_new = cost[a] * Xtemp[a];
				}
				if (Cost_new < Cost) {
					std::copy(std::begin(Xtemp), std::end(Xtemp), std::begin(Xopt));
					Cost = Cost_new;
					in_Xopt.SetData(Xtemp, 10);
					in_Delta1opt.SetData(delta1, 10);
					in_Delta2opt.SetData(delta2, 10);
					in_Delta3opt.SetData(delta3, 10);
				}
			}
		}
	}


	if (rcp1 < beta[0] || rcp2 < beta[1] || rcp3 < beta[2]) {
		infeasible = 1;
	}

	std::cout << "Brute Force results.. \n" << std::endl;
	if (infeasible == 1) {
		std::cout << "Infeasible \n" << std::endl;
	}
	else {
		std::cout << in_Xopt << std::endl;
		std::cout << in_Delta1opt << std::endl;
		std::cout << in_Delta2opt << std::endl;
		std::cout << in_Delta3opt << std::endl;
	}

	return infeasible;
}

int Algorithm1(int S, int B, double kappa[3], double beta[3], double cost[10], double BSLocationX[10], double BSLocationY[10])
{
	int infeasible = 0;
	double rcp_new = 0, rcp = 0, rcp1 = 0, rcp2 = 0, rcp3 = 0;
	double X[] = { 0,0,0,0,0,0,0,0,0,0 };
	double deltanew[] = { 0,0,0,0,0,0,0,0,0,0 };
	double delta[10][3];
	int BS = 10, BS_select = 0;
	double Cost = 1, rcp_Cost = 0.0, rcp_Cost_new = 0.0;
	double deficit1 = 0.0, deficit2 = 0.0, deficit3 = 0.0;
	mwArray in_Kappa1((mxDouble)kappa[0]);
	mwArray in_Kappa2((mxDouble)kappa[1]);
	mwArray in_Kappa3((mxDouble)kappa[2]);
	mwArray in_Kappa((mxDouble)kappa[0]);
	mwArray BSX(1, 10, mxDOUBLE_CLASS);
	BSX.SetData(BSLocationX, 10);
	mwArray BSY(1, 10, mxDOUBLE_CLASS);
	BSY.SetData(BSLocationY, 10);
	mwArray in_X(1, 10, mxDOUBLE_CLASS);
	mwArray in_Delta1(1, 10, mxDOUBLE_CLASS);
	mwArray in_Delta2(1, 10, mxDOUBLE_CLASS);
	mwArray in_Delta3(1, 10, mxDOUBLE_CLASS);
	mwArray in_Delta(1, 10, mxDOUBLE_CLASS);
	mwArray RCP;
	std::cout << "Running Algorithm 1.. \n" << std::endl;
	while (rcp1 < beta[0] && rcp2 < beta[1] && rcp3 < beta[2] && BS_select <= BS) {

		// BS selection
		int b = 0;
		for (int i = 0; i < BS; i++) {
			double Xnew[10];
			std::copy(std::begin(X), std::end(X), std::begin(Xnew));
			if (Xnew[i] == 0) {
				Xnew[i] = 1;
				deltanew[i] = 1;
				in_X.SetData(Xnew, 10);
				in_Delta.SetData(deltanew, 10);
				RCP_Approx(1, RCP, in_Kappa, in_X, in_Delta, BSX, BSY);
				double *xp = new double[1];
				RCP.GetData(xp, 1);
				rcp_new = *xp;
				for (int a = 0; a < 10; a++) {
					Cost = cost[a] * Xnew[a];
				}
				rcp_Cost_new = rcp / Cost;
				if (rcp_Cost_new >rcp_Cost) {
					rcp_Cost = rcp_Cost_new;
					b = i;
				}
				//std::cout << in_X << std::endl;

			}
			//std::cout << b << std::endl;
		}
		X[b] = 1;
		BS_select = BS_select + 1;
		in_X.SetData(X, 10);

		// computing breakpoints
		double points[10][3][5];
		for (int i = 0; i < BS; i++) {
			for (int j = 0; j < 3; j++) {
				mwArray in_Kappa((mxDouble)kappa[j]);
				points[i][j][0] = 0;
				points[i][j][1] = 0;
				points[i][j][2] = 0;
				points[i][j][3] = 0;
				points[i][j][4] = 0;
				points[i][j][5] = 0;

				if (X[i] != 0) {
					mwArray RCP_breakpoints;
					mwArray in_BSNo((mxDouble)(i + 1));
					RCP_breakpoint(1, RCP_breakpoints, in_Kappa, in_BSNo, in_X, BSX, BSY);
					//std::cout << RCP_breakpoints << std::endl;
					double *xp = new double[5];
					RCP_breakpoints.GetData(xp, 5);
					points[i][j][1] = *xp;
					points[i][j][2] = *(xp + 1);
					points[i][j][3] = *(xp + 2);
					points[i][j][4] = *(xp + 3);
					points[i][j][5] = *(xp + 4);
				}
			}
		}

		// solving Problem 2 in Cplex
		IloEnv env;
		try {
			IloInt i, j;
			IloModel model(env);

			IloInt S = 3;
			IloInt B = BS;
			IloInt No_break = 6;


			NumVarMatrix x(env, B);
			NumVarMatrix y(env, B);

			for (i = 0; i < B; i++) {
				x[i] = IloNumVarArray(env, S, 0.0, 1.0, ILOFLOAT);
				y[i] = IloNumVarArray(env, S, 0.0, 1.0, ILOFLOAT);
			}

			for (i = 0; i < B; i++) {      // supply must meet demand
				model.add(IloSum(x[i]) <= 1);
			}


			for (i = 0; i < B; i++) {
				for (j = 0; j < S; j++) {
					IloNum break1 = beta[j];
					IloNum break2 = beta[j] - points[i][j][1];
					IloNum break3 = beta[j] - points[i][j][2];
					IloNum break4 = beta[j] - points[i][j][3];
					IloNum break5 = beta[j] - points[i][j][4];
					IloNum break6 = beta[j] - points[i][j][5];
					IloNum slope1 = break2 / 0.20;
					IloNum slope2 = (break6 - break5) / 0.20;
					model.add(y[i][j] == IloPiecewiseLinear(x[i][j],
						slope1, IloNumArray(env, 6, 0.0, 0.20, 0.40, 0.6, 0.8, 1.0),
						IloNumArray(env, 6, break1, break2, break3, break4, break5, break6), slope2));
				}
			}

			IloExpr obj(env);
			for (i = 0; i < B; i++) {
				obj += IloSum(y[i]);
			}
			model.add(IloMinimize(env, obj));
			obj.end();

			IloCplex cplex(env);
			cplex.extract(model);
			cplex.exportModel("transport.lp");
			cplex.solve();


			if (cplex.solve()) {
				env.out() << "Solution status: " << cplex.getStatus() << std::endl;
				for (i = 0; i < B; i++) {
					env.out() << "   " << i << ": ";
					for (j = 0; j < S; j++) {
						delta[i][j] = cplex.getValue(x[i][j]);
						env.out() << cplex.getValue(x[i][j]) << "\t";
					}
					env.out() << std::endl;
				}
				env.out() << "   Objective value = " << cplex.getObjValue() << std::endl;
			}
			else {
				std::cout << "No solution" << std::endl;
				//cplex.printTime();
			}


		}
		catch (IloException& e) {
			cerr << "ERROR: " << e.getMessage() << std::endl;
		}
		catch (...) {
			cerr << "Error" << std::endl;
		}
		env.end();



		// computing demand deficit
		double delta1[] = { 0,0,0,0,0,0,0,0,0,0 };
		double delta2[] = { 0,0,0,0,0,0,0,0,0,0 };
		double delta3[] = { 0,0,0,0,0,0,0,0,0,0 };
		for (int i = 0; i < BS; i++)
		{
			delta1[i] = delta[i][0];
		}
		for (int i = 0; i < BS; i++)
		{
			delta2[i] = delta[i][1];
		}
		for (int i = 0; i < BS; i++)
		{
			delta3[i] = delta[i][2];
		}


		in_X.SetData(X, 10);
		in_Delta1.SetData(delta1, 10);
		in_Delta2.SetData(delta2, 10);
		in_Delta3.SetData(delta3, 10);
		mwArray RCP1;
		mwArray RCP2;
		mwArray RCP3;
		RCP_true(1, RCP1, in_Kappa1, in_X, in_Delta1, BSX, BSY);
		RCP_true(1, RCP2, in_Kappa2, in_X, in_Delta2, BSX, BSY);
		RCP_true(1, RCP3, in_Kappa3, in_X, in_Delta3, BSX, BSY);
		double *xp = new double[1];
		RCP1.GetData(xp, 1);
		rcp1 = *xp;
		RCP2.GetData(xp, 1);
		rcp2 = *xp;
		RCP3.GetData(xp, 1);
		rcp3 = *xp;

		deficit1 = beta[0] - rcp1;
		deficit2 = beta[1] - rcp2;
		deficit3 = beta[2] - rcp3;

		// selecting SP with highest demand deficit
		if (deficit1 > deficit2 && deficit1 > deficit3) {
			mwArray in_Kappa((mxDouble)kappa[0]);
			std::copy(std::begin(delta1), std::end(delta1), std::begin(deltanew));
		}
		if (deficit2 > deficit1 && deficit2 > deficit3) {
			mwArray in_Kappa((mxDouble)kappa[1]);
			std::copy(std::begin(delta1), std::end(delta1), std::begin(deltanew));
		}
		if (deficit3 > deficit1 && deficit3 > deficit1) {
			mwArray in_Kappa((mxDouble)kappa[2]);
			std::copy(std::begin(delta1), std::end(delta1), std::begin(deltanew));
		}
		/*std::cout << rcp1 << std::endl;
		std::cout << rcp2 << std::endl;
		std::cout << rcp3 << std::endl;
		std::cout << deficit1 << std::endl;
		std::cout << deficit2 << std::endl;
		std::cout << deficit3 << std::endl;*/
	}

	std::cout << "Algorithm1 results.. \n" << std::endl;
	if (rcp1 < 0.95 || rcp2 < 0.9 || rcp3 < 0.85) {
		infeasible = 1;
		std::cout << "Infeasible \n" << std::endl;
	}
	else {
		std::cout << BS_select << std::endl;
		std::cout << in_X << std::endl;
		std::cout << in_Delta1 << std::endl;
		std::cout << in_Delta2 << std::endl;
		std::cout << in_Delta3 << std::endl;
	}
	return infeasible;
}

int DE(int S, int B, double kappa[3], double beta[3], double cost[10], double BSLocationX[10], double BSLocationY[10])
{
	int infeasible = 0;
	double rcp_new = 0, rcp = 0, rcp1 = 0, rcp2 = 0, rcp3 = 0;
	double X[] = { 0,0,0,0,0,0,0,0,0,0 };
	double deltanew[] = { 0,0,0,0,0,0,0,0,0,0 };
	double delta[10][3];
	mwArray in_Kappa1((mxDouble)kappa[0]);
	mwArray in_Kappa2((mxDouble)kappa[1]);
	mwArray in_Kappa3((mxDouble)kappa[2]);
	mwArray in_Kappa((mxDouble)kappa[0]);
	mwArray in_X(1, 10, mxDOUBLE_CLASS);
	mwArray in_Delta1(1, 10, mxDOUBLE_CLASS);
	mwArray in_Delta2(1, 10, mxDOUBLE_CLASS);
	mwArray in_Delta3(1, 10, mxDOUBLE_CLASS);
	mwArray in_Xopt(1, 10, mxDOUBLE_CLASS);
	mwArray in_Delta1opt(1, 10, mxDOUBLE_CLASS);
	mwArray in_Delta2opt(1, 10, mxDOUBLE_CLASS);
	mwArray in_Delta3opt(1, 10, mxDOUBLE_CLASS);
	mwArray BSX(1, 10, mxDOUBLE_CLASS);
	BSX.SetData(BSLocationX, 10);
	mwArray BSY(1, 10, mxDOUBLE_CLASS);
	BSY.SetData(BSLocationY, 10);
	mwArray RCP;
	int i, j, k, a;
	float r1, r2;
	double Cost = 1800, Cost_new = 0;
	int gen = 0, genmax = 1000, NP = 10000;
	float F = 0.1, CR = 0.5;           // control variables of DE  
	std::cout << "Running DE.. \n" << std::endl;
	for (j = 0; j < B; j++) {
		for (k = 0; k < S; k++) {
			delta[j][k] = 0;
		}
	}
	while ((gen < genmax)) {
		gen++;
		for (i = 0; i < NP; i++) {
			for (j = 0; j < B; j++) {
				X[j] = 0;
			}
			for (j = 0; j < B; j++) {
				for (k = 0; k < S; k++) {
					r1 = rand() % 10 + 1;
					r2 = rand() % 10 + 1;
					if (r1 / 10 <= CR) {
						delta[j][k] = F*r2 / 10;
						if (r2 > 0)
							X[j] = 1;
					}
				}
			}
			double delta1[] = { 0,0,0,0,0,0,0,0,0,0 };
			double delta2[] = { 0,0,0,0,0,0,0,0,0,0 };
			double delta3[] = { 0,0,0,0,0,0,0,0,0,0 };
			for (int i = 0; i < B; i++)
			{
				delta1[i] = delta[i][0];
			}
			for (int i = 0; i < B; i++)
			{
				delta2[i] = delta[i][1];
			}
			for (int i = 0; i < B; i++)
			{
				delta3[i] = delta[i][2];
			}
			//std::cout << "values are.. \n" << std::endl;
			in_X.SetData(X, 10);
			in_Delta1.SetData(delta1, 10);
			in_Delta2.SetData(delta2, 10);
			in_Delta3.SetData(delta3, 10);
			/*std::cout << in_X << std::endl;
			std::cout << in_Delta1 << std::endl;
			std::cout << in_Delta2 << std::endl;
			std::cout << in_Delta3 << std::endl;*/
			mwArray RCP1;
			mwArray RCP2;
			mwArray RCP3;
			RCP_true(1, RCP1, in_Kappa1, in_X, in_Delta1, BSX, BSY);
			RCP_true(1, RCP2, in_Kappa2, in_X, in_Delta2, BSX, BSY);
			RCP_true(1, RCP3, in_Kappa3, in_X, in_Delta3, BSX, BSY);
			double *xp = new double[1];
			RCP1.GetData(xp, 1);
			rcp1 = *xp;
			RCP2.GetData(xp, 1);
			rcp2 = *xp;
			RCP3.GetData(xp, 1);
			rcp3 = *xp;
			/*std::cout << "rcps are.. \n" << std::endl;
			std::cout << rcp1 << std::endl;
			std::cout << rcp2 << std::endl;
			std::cout << rcp3 << std::endl;*/
			if (rcp1 >= beta[0] && rcp2 >= beta[1] && rcp3 >= beta[2]) {
				for (a = 0; a < 10; a++) {
					Cost_new = cost[a] * X[a];
				}
				if (Cost_new < Cost) {
					Cost = Cost_new;
					in_Xopt.SetData(X, 10);
					in_Delta1opt.SetData(delta1, 10);
					in_Delta2opt.SetData(delta2, 10);
					in_Delta3opt.SetData(delta3, 10);
				}

			}
		}
	}

	std::cout << "DE results.. \n" << std::endl;
	if (rcp1 < 0.95 || rcp2 < 0.9 || rcp3 < 0.85) {
		infeasible = 1;
		std::cout << "Infeasible \n" << std::endl;
	}
	else {
		std::cout << in_Xopt << std::endl;
		std::cout << in_Delta1opt << std::endl;
		std::cout << in_Delta2opt << std::endl;
		std::cout << in_Delta3opt << std::endl;
	}
	return infeasible;
}

void Algorithm2(int S, int B, double kappa[3], double beta[3], double BSLocationX[10], double BSLocationY[10])
{
	    
		mwArray in_Kappa1((mxDouble)kappa[0]);
		mwArray in_Kappa2((mxDouble)kappa[1]);
		mwArray in_Kappa3((mxDouble)kappa[2]);
		mwArray in_Kappa((mxDouble)kappa[0]);
		mwArray in_X(1, 10, mxDOUBLE_CLASS);
		mwArray in_Delta1(1, 10, mxDOUBLE_CLASS);
		mwArray in_Delta2(1, 10, mxDOUBLE_CLASS);
		mwArray in_Delta3(1, 10, mxDOUBLE_CLASS);
		mwArray in_Xopt(1, 10, mxDOUBLE_CLASS);
		mwArray in_Delta1opt(1, 10, mxDOUBLE_CLASS);
		mwArray in_Delta2opt(1, 10, mxDOUBLE_CLASS);
		mwArray in_Delta3opt(1, 10, mxDOUBLE_CLASS);
		mwArray BSX(1, 10, mxDOUBLE_CLASS);
		BSX.SetData(BSLocationX, 10);
		mwArray BSY(1, 10, mxDOUBLE_CLASS);
		BSY.SetData(BSLocationY, 10);
		double Xopt[] = { 1,1,1,1,1,1,1,1,1,1 };
		mwArray Deltas(1, 10, mxDOUBLE_CLASS);
		in_X.SetData(Xopt, 10);
		double delta_s[10];
		double points2[10][5];
		double alpha[10];
		for (int i = 0; i < 10; i++) {
			alpha[i] = 1.0;
		}

		for (int s = 0; s < S; s++) {
			// computing breakpoints
			for (int i = 0; i < 10; i++) {
				mwArray RCP_breakpoints;
				mwArray in_BSNo((mxDouble)(i + 1));
				mwArray in_Kappa((mxDouble)kappa[s]);
				RCP_breakpoint(1, RCP_breakpoints, in_Kappa, in_BSNo, in_X, BSX, BSY);
				//std::cout << RCP_breakpoints << std::endl;
				double *xp = new double[5];
				RCP_breakpoints.GetData(xp, 5);
				points2[i][1] = *xp;
				points2[i][2] = *(xp + 1);
				points2[i][3] = *(xp + 2);
				points2[i][4] = *(xp + 3);
				points2[i][5] = *(xp + 4);
			}


			IloEnv env;
			try {
				IloInt i;
				IloModel model(env);

				IloInt S = 1;
				IloInt B = 10;
				IloInt No_break = 6;



				NumVarMatrix x(env, B);
				NumVarMatrix y(env, B);

				for (i = 0; i < B; i++) {
					IloNum a = alpha[i];
					x[i] = IloNumVarArray(env, S, 0.0, a, ILOFLOAT);
					y[i] = IloNumVarArray(env, S, 0.0, a, ILOFLOAT);
				}

				IloExpr obj(env);
				for (i = 0; i < B; i++) {
					obj += IloSum(x[i]);
				}
				model.add(IloMinimize(env, obj));

				for (i = 0; i < B; i++) {
					IloNum break1 = 0;
					IloNum break2 = points2[i][1];
					IloNum break3 = points2[i][2];
					IloNum break4 = points2[i][3];
					IloNum break5 = points2[i][4];
					IloNum break6 = points2[i][5];
					IloNum slope1 = break2 / 0.20;
					IloNum slope2 = (break6 - break5) / 0.20;
					model.add(y[i][0] == IloPiecewiseLinear(x[i][0],
						slope1, IloNumArray(env, 6, 0.0, 0.20, 0.40, 0.6, 0.8, 1.0),
						IloNumArray(env, 6, break1, break2, break3, break4, break5, break6), slope2));
				}


				IloExpr constraint1(env);
				for (i = 0; i < B; i++) {
					constraint1 += IloSum(y[i]);
				}

				model.add((constraint1) >= beta[s]);
				constraint1.end();

				IloCplex cplex(env);
				cplex.extract(model);
				cplex.exportModel("transport.lp");
				cplex.solve();


				if (cplex.solve()) {
					//env.out() << "Solution status: " << cplex.getStatus() << std::endl;
					for (i = 0; i < B; i++) {
						delta_s[i] = cplex.getValue(x[i][0]);
						alpha[i] = alpha[i] - delta_s[i];
						// env.out() << cplex.getValue(x[i][0]) << "\t";
					}
					//env.out() << " - Solution: " << std::endl;

					//env.out() << "   Cost = " << cplex.getObjValue() << std::endl;
				}
				else {
					std::cout << "No solution" << std::endl;
					//cplex.printTime();
				}


			}
			catch (IloException& e) {
				cerr << "ERROR: " << e.getMessage() << std::endl;
			}
			catch (...) {
				cerr << "Error" << std::endl;
			}
			env.end();


			Deltas.SetData(delta_s, 10);
			std::cout << "slices for SP" << s + 1 << "...\n" << std::endl;
			std::cout << Deltas << std::endl;

		}

}


