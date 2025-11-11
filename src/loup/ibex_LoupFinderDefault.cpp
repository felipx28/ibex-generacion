//============================================================================
//                                  I B E X
// File        : ibex_LoupFinderDefault.cpp
// Author      : Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Jul 09, 2017
// Last update : Feb 13, 2025
//============================================================================

#include "ibex_LoupFinderDefault.h"
#include "ibex_LoupFinderInHC4.h"
#include "ibex_LoupFinderFwdBwd.h"
#include "ibex_BxpLinearRelaxArgMin.h"
#include "ibex_LoupFinderProbing.h"
#include "ibex_Solver.h"
#include <iostream>
#include <fstream>

using namespace std;

namespace ibex {

LoupFinderDefault::LoupFinderDefault(const System& sys, bool inHC4) :
	finder_probing(inHC4? (LoupFinder&) *new LoupFinderInHC4(sys) : (LoupFinder&) *new LoupFinderFwdBwd(sys)),
	finder_x_taylor(sys), sys(sys) {

}

void LoupFinderDefault::add_property(const IntervalVector& init_box, BoxProperties& prop) {
	finder_probing.add_property(init_box,prop);
	finder_x_taylor.add_property(init_box,prop);

	//--------------------------------------------------------------------------
	/* Using line search from LP relaxation minimizer seems not interesting. */
//	if (!prop[BxpLinearRelaxArgMin::get_id(finder_x_taylor.sys)]) {
//		prop.add(new BxpLinearRelaxArgMin(finder_x_taylor.sys));
//	}
	//--------------------------------------------------------------------------

}

std::pair<IntervalVector, double> LoupFinderDefault::find(const IntervalVector& box, const IntervalVector& old_loup_point, double old_loup, BoxProperties& prop) {

	pair<IntervalVector,double> p=make_pair(old_loup_point, old_loup);

	bool found=false;

	try {
		p=finder_probing.find(box,p.first,p.second,prop);
		found=true;
	} catch(NotFound&) { }

	try {


		/***************************
		 ** INICIO MODIFICACIONES **
		 ***************************/

		// std::ofstream MyFile("/home/felipe/Documents/magister/datos50.txt", std::ios::app);
		// if(!MyFile.is_open()){
		// 	std::cerr << "No se pudo abrir el archivo. Verifica la ruta y permisos." << std::endl;
		// 	exit(1);
		// }

		// // Escribe los datos
		// // 1 input
		// MyFile << "variables: " << sys.nb_var << endl;
		
		// // 2 input
		// MyFile << "restricciones: " << sys.nb_ctr << endl;

		// // 3 input
		// MyFile << "lb(f_obj): " << sys.goal->eval(box).lb() << endl;
		
		// // 4 input
		// MyFile << "ub(f_obj): " << sys.goal->eval(box).ub() << endl << endl;

		// MyFile << "-------------------------------------------------------" << endl << endl;

		// MyFile.close();

		// exit(1);

		// cout << "eval fobj: " << sys.f_ctrs[0].eval(box) << endl;
		// cout << "lb(fobj): " << sys.f_ctrs[0].eval(box).lb() << endl;
		// cout << "ub(fobj): " << sys.f_ctrs[0].eval(box).ub() << endl;
		
		//exit(1);
		//		MyFile << "inter fobj: " << sys.goal->eval(box) << endl;
		/*for(int i = 0; i < box.size(); i++) {
			cout << "ub " << box[i].ub() << endl;
			cout << "lb " << box[i].lb() << endl; 
		}*/

		//IntervalMatrix jac(0,0) ;
		//sys.ctrs_jacobian(box, jac);
		//cout << "jac " << sys.jac[0][0] << endl;
		//exit(1); //restricciones activas son las que si influencian en la decisión dentro de la caja utilizada

		//MyFile << endl; //sys.args entrega variables
		

		/***************************
		 ** FIN MODIFICACIONES **
		 ***************************/

		// TODO
		// in_x_taylor.set_inactive_ctr(entailed->norm_entailed);
		p=finder_x_taylor.find(box,p.first,p.second,prop);
		found=true;
	} catch(NotFound&) { }

	if (found) {
		//--------------------------------------------------------------------------
		/* Using line search from LP relaxation minimizer seems not interesting. */
		//	BxpLinearRelaxArgMin* argmin=(BxpLinearRelaxArgMin*) prop[BxpLinearRelaxArgMin::get_id(finder_x_taylor.sys)];
		BxpLinearRelaxArgMin* argmin=NULL;
		//--------------------------------------------------------------------------

		if (argmin && argmin->argmin()) {
			Vector loup_point = p.first.lb();
			double loup = p.second;
			LoupFinderProbing(finder_x_taylor.sys).dichotomic_line_search(loup_point, loup, *argmin->argmin(), false);
			//cout << "better loup found! " << loup << endl;
			p=make_pair(loup_point,loup);
		}
		return p;
	} else
		throw NotFound();
}

LoupFinderDefault::~LoupFinderDefault() {
	delete &finder_probing;
}

void LoupFinderDefault::enable_statistics(Statistics& stats, const std::string& prefix) {
	finder_probing.enable_statistics(stats, prefix);
	finder_x_taylor.enable_statistics(stats, prefix);
}

} /* namespace ibex */
