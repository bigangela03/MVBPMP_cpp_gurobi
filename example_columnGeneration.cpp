//This is an example from https://groups.google.com/g/gurobi/c/pkBNfu-iX0k

#include <vector>
#include "gurobi_c++.h"
#include <sstream>

using namespace std;

int main(int argc, char *argv[])
{
	GRBEnv* env = 0;

	try{
		//Roll width!
		int B = 9;
		//Number of demands
		vector<double> w(7);
		w = { 2, 3, 4, 5, 6, 7, 8 };
		//Width of each demands
		vector<double> q(7);
		q = { 4, 2, 6, 6, 2, 2, 2 };

		int wNumber = w.size();

		//Patterns
		vector<vector<double>> pat(wNumber);

		cout << "initial patterns \n";
		for (int i = 0; i < wNumber; i++){
			for (int j = 0; j < wNumber; j++){
				if (i == j){
					pat[i].push_back(int(B / w[i]));
					cout << pat[i][j] << "\t";
				}
				else{
					pat[i].push_back(0);
					cout << pat[i][j] << "\t";
				}
			}
			cout << "\n";
		}

		int K = pat.size();
		cout << "number of patterns= " << K << "\n";

		env = new GRBEnv();
		GRBModel Master = GRBModel(*env);

		//Variables of master problem
		vector<GRBVar> x_VarMaster;
		for (int i = 0; i < K; i++){
			ostringstream vname;
			vname << "x" << i;
			x_VarMaster.push_back(
			Master.addVar(0.0, GRB_INFINITY, 1.0, GRB_INTEGER, vname.str()));
		}
		cout << "variables added \n";
		Master.update();

		//Constraints of master problem!
		vector<GRBConstr> order_ConstMaster;
		for (int i = 0; i < wNumber; i++){
			double coef;
			GRBVar tempVar;
			vector<double> temp = pat[i];
			int sizeTemp = temp.size();
			for (int j = 0; j < sizeTemp; j++){
				if (temp[j] > 0){
					coef = temp[j];
					tempVar = x_VarMaster[j];
					break;
				}
			}
			ostringstream cname;
			cname << "order" << i;
			order_ConstMaster.push_back(Master.addConstr
				(GRBLinExpr(tempVar, coef), GRB_GREATER_EQUAL, q[i], cname.str()));
		}
		cout << "constraints added! \n";
		Master.update();

		while (true){
			GRBModel relax = Master.relax();
			relax.optimize();
			cout << "Relaxed master solved! \n";
			GRBConstr* tempCons = relax.getConstrs();
			int constrNum = relax.get(GRB_IntAttr_NumConstrs);

			//Dual variables
			vector<double> pi;
			for (int j = 0; j < constrNum; j++){
				pi.push_back(relax.getConstr(j).get(GRB_DoubleAttr_Pi));
			}
			cout << "dual vriables: \n";
			for (int j = 0; j < constrNum; j++){
				cout << pi[j] << "\t";
			}
			cout << "\n";

			cout << "Pricing Problem \n";
			GRBModel PP = GRBModel(*env);
			PP.set(GRB_IntAttr_ModelSense, -1);

			vector<GRBVar> y_VarPricing;
			for (int i = 0; i < wNumber; i++){
				ostringstream vname;
				vname << "y" << i;
				y_VarPricing.push_back(
					PP.addVar(0.0, q[i], pi[i], GRB_INTEGER, vname.str()));
			}
			PP.update();
			cout << "PP variables added! \n";
			GRBLinExpr L;
			for (int i = 0; i < wNumber; i++){
				L = y_VarPricing[i] * w[i] + L;
			}
			PP.addConstr(L, GRB_LESS_EQUAL, B, "co-PP");

			PP.update();
			cout << "PP constraints added! \n";
			PP.getEnv().set(GRB_IntParam_OutputFlag, 0); //Silent Mode
			PP.optimize();

			double objValue_PP = PP.get(GRB_DoubleAttr_ObjVal);
			cout << "********************\n";
			cout << "objective of PP problem: " << objValue_PP << "\n";
			cout << "********************\n";

			if (objValue_PP < 1.000001){
				cout << "Optimum Found!";
				break;
			}

			//New pattern
			vector<double> newPatTemp;
			for (int i = 0; i < wNumber; i++){
				double temp = int(y_VarPricing[i].get(GRB_DoubleAttr_X) + 0.5);
				newPatTemp.push_back(temp);
			}
			pat.push_back(newPatTemp);

			//Create new column
			GRBColumn col;
			for (int i = 0; i < wNumber; i++){
				if (pat[K][i] > 0){
					col.addTerm(pat[K][i], order_ConstMaster[i]);
				}
			}
			//Add new column to Master
			x_VarMaster.push_back(Master.addVar
				(0.0, GRB_INFINITY, 1.0, GRB_INTEGER, col, "x"));
			Master.update();

			K++;

		}

		Master.optimize();
	}
	catch (GRBException e)
	{
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...)
	{
		cout << "Exception during optimization" << endl;
	}


	getchar();

	return 0;
}
