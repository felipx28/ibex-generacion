//                                  I B E X
// File        : ibex_Optimizer.cpp
// Author      : Gilles Chabert, Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : May 14, 2012
// Last Update : Feb 13, 2025
//============================================================================

#include "ibex_Optimizer.h"
#include "ibex_Timer.h"
#include "ibex_Function.h"
#include "ibex_NoBisectableVariableException.h"
#include "ibex_BxpOptimData.h"
#include "ibex_CovOptimData.h"
#include <stdlib.h>
#include "ibex_CellBeamSearch.h"
#include "ibex_LoupFinderDefault.h"

#include "ibex_SmearFunction.h"
#include "ibex_ExtendedSystem.h"
#include "ibex_OptimLargestFirst.h"
#include "ibex_System.h"

#include <float.h>
#include <vector>
#include <queue>
#include <stdlib.h>
#include <iomanip>

using namespace std;

namespace ibex {

/*
 * TODO: redundant with ExtendedSystem.
 */
void Optimizer::write_ext_box(const IntervalVector& box, IntervalVector& ext_box) {
	int i2=0;
	for (int i=0; i<n; i++,i2++) {
		if (i2==goal_var) i2++; // skip goal variable
		ext_box[i2]=box[i];
	}
}

void Optimizer::read_ext_box(const IntervalVector& ext_box, IntervalVector& box) {
	int i2=0;
	for (int i=0; i<n; i++,i2++) {
		if (i2==goal_var) i2++; // skip goal variable
		box[i]=ext_box[i2];
	}
}

Optimizer::Optimizer(int n, Ctc& ctc, Bsc& bsc, LoupFinder& finder,
		CellBufferOptim& buffer,
		int goal_var, double eps_x, double rel_eps_f, double abs_eps_f,
		bool enable_statistics) :
                						n(n), goal_var(goal_var),
										ctc(ctc), bsc(bsc), loup_finder(finder), buffer(buffer),
										eps_x(n, eps_x), rel_eps_f(rel_eps_f), abs_eps_f(abs_eps_f),
										trace(0), timeout(-1), extended_COV(true), anticipated_upper_bounding(true),
										status(SUCCESS),
										uplo(NEG_INFINITY), uplo_of_epsboxes(POS_INFINITY), loup(POS_INFINITY),
										loup_point(IntervalVector::empty(n)), initial_loup(POS_INFINITY), loup_changed(false),
										time(0), nb_cells(0), cov(NULL) {

	if (trace) cout.precision(12);
	
	if (enable_statistics) {
		statistics = new Statistics();
		// TODO: enable statistics for missing operators (cell buffer)
		bsc.enable_statistics(*statistics, "Bsc"); 
		ctc.enable_statistics(*statistics, "Ctc"); 
		loup_finder.enable_statistics(*statistics, "LoupFinder"); 
	} else
		statistics = NULL;
}

Optimizer::Optimizer(OptimizerConfig& config) :
	Optimizer(
		config.nb_var(), 
		config.get_ctc(), 
		config.get_bsc(), 
		config.get_loup_finder(),
		config.get_cell_buffer(),
		config.goal_var(),
		OptimizerConfig::default_eps_x, // tmp, see below
		config.get_rel_eps_f(),
		config.get_abs_eps_f(),
		config.with_statistics()) {

	(Vector&) eps_x				= config.get_eps_x();
	trace						= config.get_trace();
	timeout						= config.get_timeout();
	extended_COV				= config.with_extended_cov();
	anticipated_upper_bounding	= config.with_anticipated_upper_bounding();
}

Optimizer::~Optimizer() {
	if (cov) delete cov;
	if (statistics) delete statistics;
}

// compute the value ymax (decreasing the loup with the precision)
// the heap and the current box are contracted with y <= ymax
double Optimizer::compute_ymax() {
	if (anticipated_upper_bounding) {
		//double ymax = loup - rel_eps_f*fabs(loup); ---> wrong :the relative precision must be correct for ymax (not loup)
		double ymax = loup>0 ?
				1/(1+rel_eps_f)*loup
		:
				1/(1-rel_eps_f)*loup;

		if (loup - abs_eps_f < ymax)
			ymax = loup - abs_eps_f;
		//return ymax;
		return next_float(ymax);
	} else
		return loup;
}

bool Optimizer::update_loup(const IntervalVector& box, BoxProperties& prop) {

	try {

		pair<IntervalVector,double> p=loup_finder.find(box,loup_point,loup,prop);
		loup_point = p.first;
		loup = p.second;

		if (trace) {
			cout << "                    ";
			cout << "\033[32m loup= " << loup << "\033[0m" << endl;
//			cout << " loup point=";
//			if (loup_finder.rigorous())
//				cout << loup_point << endl;
//			else
//				cout << loup_point.lb() << endl;
		}
		return true;

	} catch(LoupFinder::NotFound&) {
		return false;
	}
}

//bool Optimizer::update_entailed_ctr(const IntervalVector& box) {
//	for (int j=0; j<m; j++) {
//		if (entailed->normalized(j)) {
//			continue;
//		}
//		Interval y=sys.ctrs[j].f.eval(box);
//		if (y.lb()>0) return false;
//		else if (y.ub()<=0) {
//			entailed->set_normalized_entailed(j);
//		}
//	}
//	return true;
//}

void Optimizer::update_uplo() {
	double new_uplo=POS_INFINITY;

	if (! buffer.empty()) {
		new_uplo= buffer.minimum();
		if (new_uplo > loup && uplo_of_epsboxes > loup) {
			cout << " loup = " << loup << " new_uplo=" << new_uplo <<  " uplo_of_epsboxes=" << uplo_of_epsboxes << endl;
			ibex_error("optimizer: new_uplo>loup (please report bug)");
		}
		if (new_uplo < uplo) {
			cout << "uplo= " << uplo << " new_uplo=" << new_uplo << endl;
			ibex_error("optimizer: new_uplo<uplo (please report bug)");
		}

		// uplo <- max(uplo, min(new_uplo, uplo_of_epsboxes))
		if (new_uplo < uplo_of_epsboxes) {
			if (new_uplo > uplo) {
				uplo = new_uplo;

				if (trace)
					cout << "\033[33m uplo= " << uplo << "\033[0m" << endl;
			}
		}
		else uplo = uplo_of_epsboxes;
	}
	else if (buffer.empty() && loup != POS_INFINITY) {
		// empty buffer : new uplo is set to ymax (loup - precision) if a loup has been found
		new_uplo=compute_ymax(); // not new_uplo=loup, because constraint y <= ymax was enforced
		//    cout << " new uplo buffer empty " << new_uplo << " uplo " << uplo << endl;

		double m = (new_uplo < uplo_of_epsboxes) ? new_uplo :  uplo_of_epsboxes;
		if (uplo < m) uplo = m; // warning: hides the field "m" of the class
		// note: we always have uplo <= uplo_of_epsboxes but we may have uplo > new_uplo, because
		// ymax is strictly lower than the loup.
	}

}

void Optimizer::update_uplo_of_epsboxes(double ymin) {

	// the current box cannot be bisected.  ymin is a lower bound of the objective on this box
	// uplo of epsboxes can only go down, but not under uplo : it is an upperbound for uplo,
	// that indicates a lowerbound for the objective in all the small boxes
	// found by the precision criterion
	assert (uplo_of_epsboxes >= uplo);
	assert(ymin >= uplo);
	if (uplo_of_epsboxes > ymin) {
		uplo_of_epsboxes = ymin;
		if (trace) {
			cout << " unprocessable tiny box: now uplo<=" << setprecision(12) <<  uplo_of_epsboxes << " uplo=" << uplo << endl;
		}
	}
}

void Optimizer::handle_cell(Cell& c) {

	contract_and_bound(c);

	if (c.box.is_empty()) {
		delete &c;
	} else {
		buffer.push(&c);
	}
}

void Optimizer::contract_and_bound(Cell& c) {

	/*======================== contract y with y<=loup ========================*/
	Interval& y=c.box[goal_var];

	double ymax;
	if (loup==POS_INFINITY) ymax = POS_INFINITY;
	// ymax is slightly increased to favour subboxes of the loup
	// TODO: useful with double heap??
	else ymax = compute_ymax()+1.e-15;

	y &= Interval(NEG_INFINITY,ymax);

	if (y.is_empty()) {
		c.box.set_empty();
		return;
	} else {
		c.prop.update(BoxEvent(c.box,BoxEvent::CONTRACT,BitSet::singleton(n+1,goal_var)));
	}

	/*================ contract x with f(x)=y and g(x)<=0 ================*/
	//cout << " [contract]  x before=" << c.box << endl;
	//cout << " [contract]  y before=" << y << endl;

	ContractContext context(c.prop);
	if (c.bisected_var!=-1) {
		context.impact.clear();
		context.impact.add(c.bisected_var);
		context.impact.add(goal_var);
	}

	ctc.contract(c.box, context);
	//cout << c.prop << endl;
	if (c.box.is_empty()) return;

	//cout << " [contract]  x after=" << c.box << endl;
	//cout << " [contract]  y after=" << y << endl;
	/*====================================================================*/

	/*========================= update loup =============================*/

	IntervalVector tmp_box(n);
	read_ext_box(c.box,tmp_box);

	c.prop.update(BoxEvent(c.box,BoxEvent::CHANGE));

	bool loup_ch=update_loup(tmp_box, c.prop);

	// update of the upper bound of y in case of a new loup found
	if (loup_ch) {
		y &= Interval(NEG_INFINITY,compute_ymax());
		c.prop.update(BoxEvent(c.box,BoxEvent::CONTRACT,BitSet::singleton(n+1,goal_var)));
	}

	//TODO: should we propagate constraints again?

	loup_changed |= loup_ch;

	if (y.is_empty()) { // fix issue #44
		c.box.set_empty();
		return;
	}

	/*====================================================================*/
	// Note: there are three different cases of "epsilon" box,
	// - NoBisectableVariableException raised by the bisector (---> see optimize(...)) which
	//   is independent from the optimizer
	// - the width of the box is less than the precision given to the optimizer ("eps_x" for
	//   the original variables and "abs_eps_f" for the goal variable)
	// - the extended box has no bisectable domains (if eps_x=0 or <1 ulp)
	if (((tmp_box.diam()-eps_x).max()<=0 && y.diam() <=abs_eps_f) || !c.box.is_bisectable()) {
		update_uplo_of_epsboxes(y.lb());
		c.box.set_empty();
		return;
	}

	// ** important: ** must be done after upper-bounding
	//kkt.contract(tmp_box);

	if (tmp_box.is_empty()) {
		c.box.set_empty();
	} else {
		// the current extended box in the cell is updated
		write_ext_box(tmp_box,c.box);
	}
}

Optimizer::Status Optimizer::optimize(const IntervalVector& init_box, double obj_init_bound) {
	start(init_box, obj_init_bound);
	return optimize();
}


Optimizer::Status Optimizer::optimize(const CovOptimData& data, double obj_init_bound) {
	start(data, obj_init_bound);
	return optimize();
}

Optimizer::Status Optimizer::optimize(const char* cov_file, double obj_init_bound) {
	CovOptimData data(cov_file);
	start(data, obj_init_bound);
	return optimize();
}

void Optimizer::start(const IntervalVector& init_box, double obj_init_bound) {

	loup=obj_init_bound; // loup => lower upperbound, uplo => upper lowerbound


	/***************************
	 ** INICIO MODIFICACIONES **
	 ***************************/

	// double search_space = 1;
	// double bigger_diam = NEG_INFINITY;
	// double lower_diam = POS_INFINITY;
	
	IntervalVector aux(init_box.size());  //crea una variable auxiliar del tamaño de la cantidad de variables

	for (int i = 0; i < init_box.size(); i++) {
		aux[i] = init_box[i]; //asigna el valor de init_box a la variable auxiliar

		//+10% 
		//aux[i] = Interval((aux[i].lb() - aux[i].diam()/20), (aux[i].ub() + aux[i].diam()/20)); // esta aumenta un 10%

		//+25% 
		//aux[i] = Interval((aux[i].lb() - aux[i].diam()/8), (aux[i].ub() + aux[i].diam()/8)); // esta aumenta un 25%

		//+50% 
		//aux[i] = Interval((aux[i].lb() - aux[i].diam()/4), (aux[i].ub() + aux[i].diam()/4)); // esta aumenta un 50%
	}
	
	// Just to initialize the "loup" for the buffer
	// TODO: replace with a set_loup function
	
	buffer.contract(loup);

	uplo=NEG_INFINITY;
	uplo_of_epsboxes=POS_INFINITY;

	nb_cells=0;

	buffer.flush();

	Cell* root=new Cell(IntervalVector(n+1));

	write_ext_box(aux, root->box);

	// add data required by the bisector
	bsc.add_property(aux, root->prop);

	// add data required by the contractor
	ctc.add_property(aux, root->prop);

	// add data required by the buffer
	buffer.add_property(aux, root->prop);

	// add data required by the loup finder
	loup_finder.add_property(aux, root->prop);

	//cout << "**** Properties ****\n" << root->prop << endl;

	loup_changed=false;
	initial_loup=obj_init_bound;

	loup_point = aux; //.set_empty();
	time=0;

	if (cov) delete cov;
	cov = new CovOptimData(extended_COV? n+1 : n, extended_COV);
	cov->data->_optim_time = 0;
	cov->data->_optim_nb_cells = 0;

	handle_cell(*root);
}

void Optimizer::start(const CovOptimData& data, double obj_init_bound) {

	loup=obj_init_bound;

	// Just to initialize the "loup" for the buffer
	// TODO: replace with a set_loup function
	buffer.contract(loup);

	uplo=data.uplo();
	loup=data.loup();
	loup_point=data.loup_point();
	uplo_of_epsboxes=POS_INFINITY;

	nb_cells=0;

	buffer.flush();

	for (size_t i=loup_point.is_empty()? 0 : 1; i<data.size(); i++) {

		IntervalVector box(n+1);

		if (data.is_extended_space())
			box = data[i];
		else {
			write_ext_box(data[i], box);
			box[goal_var] = Interval(uplo,loup);
			ctc.contract(box);
			if (box.is_empty()) continue;
		}

		Cell* cell=new Cell(box);

		// add data required by the cell buffer
		buffer.add_property(box, cell->prop);

		// add data required by the bisector
		bsc.add_property(box, cell->prop);

		// add data required by the contractor
		ctc.add_property(box, cell->prop);

		// add data required by the loup finder
		loup_finder.add_property(box, cell->prop);

		buffer.push(cell);
	}

	loup_changed=false;
	initial_loup=obj_init_bound;

	time=0;

	if (cov) delete cov;
	cov = new CovOptimData(extended_COV? n+1 : n, extended_COV);
	cov->data->_optim_time = data.time();
	cov->data->_optim_nb_cells = data.nb_cells();
}

/**
 * CALCULO DE LA FEATURE GAP REL LOUP
 * EL FEATURE ES MÁS ROBUSTO QUE EL ANTERIOR, POR LO QUE DEBIESE FUNCIONAR MEJOR
 * \brief Calculate a more robust relative gap feature based on the loup.
 */

// double calculate_improved_gap_rel_loup(double current_loup, double lb_f_obj, double ub_f_obj) {
//     double epsilon = 1e-12; // Valor más pequeño para mayor precisión
//     double gap_rel_loup;
    
//     if (!std::isfinite(current_loup) || current_loup == POS_INFINITY) {
//         // Caso 1: No hay loup válido
//         // El gap representa la incertidumbre relativa de la caja actual
//         double numerator = ub_f_obj - lb_f_obj;
//         double denominator = std::max(fabs(ub_f_obj), fabs(lb_f_obj));
        
//         if (denominator < epsilon) {
//             gap_rel_loup = 1.0; // Máxima incertidumbre cuando ambos límites son ~0
//         } else {
//             gap_rel_loup = numerator / (denominator + epsilon);
//         }
        
//         // Normalizar a [0,1] - mayor valor = menos prometedora
//         gap_rel_loup = std::min(1.0, std::max(0.0, gap_rel_loup));
        
//     } else {
//         // Caso 2: Hay loup válido
//         // Calculamos qué tan cerca está el lower bound de la mejor solución conocida
//         double numerator = current_loup - lb_f_obj;
        
//         // Si lb_f_obj > current_loup, la caja es muy prometedora (gap negativo)
//         if (numerator < 0) {
//             gap_rel_loup = 0.0; // Muy prometedora
//         } else {
//             // Normalizamos por la magnitud del loup
//             double denominator = fabs(current_loup) + epsilon;
//             gap_rel_loup = numerator / denominator;
            
//             // Limitamos a [0,1] para evitar valores extremos
//             gap_rel_loup = std::min(1.0, gap_rel_loup);
//         }
//     }
    
//     return gap_rel_loup;
// }

/**
 * Feature adicional: qué tan prometedora es la caja comparada con otras en el buffer
 * Esta feature es de prueba por ahora, solo para ver como funciona
 */
// double calculate_relative_promise(double lb_f_obj, double buffer_minimum) {
//     double epsilon = 1e-12;
    
//     if (!std::isfinite(buffer_minimum)) {
//         return 0.5; // Valor neutro si no hay referencia
//     }
    
//     double diff = lb_f_obj - buffer_minimum;
    
//     if (fabs(diff) < epsilon) {
//         return 1.0; // Es la más prometedora del buffer
//     }
    
//     // Normalizamos: valores más cercanos a buffer_minimum son más prometedores
//     double scale = std::max(fabs(buffer_minimum), 1.0);
//     return std::exp(-fabs(diff) / scale); // Decay exponencial
// }

double calculate_bounds_ratios(double lb, double ub, double epsilon, int id) {
    static double initial_lb, initial_ub, initial_obj_range;
    static bool initial_bounds_set = false;
    
    // Primera llamada: guardar valores iniciales
    if (id == 1) {
        initial_lb = lb;
        initial_ub = ub;
        initial_obj_range = ub - lb;
        initial_bounds_set = true;
    }
    
    if (!initial_bounds_set) {
        return 0.5; // Valor neutral si no hay inicialización
    }
    
    double current_obj_range = ub - lb;
    
    // Validar bounds actuales
    if (!std::isfinite(current_obj_range) || current_obj_range < epsilon) {
        return 0.5; // Valor neutral para casos inválidos
    }
    
    // Si bounds iniciales son infinitos o inválidos
    if (!std::isfinite(initial_obj_range) || initial_obj_range < epsilon) {
        return 0.5; // No podemos medir progreso, valor neutral
    }
    
    // Progreso normal: qué proporción del gap inicial queda
    double progress_ratio = current_obj_range / initial_obj_range;
    
    // Clampear entre [0, 1]
    // 1.0 = no hay progreso (gap igual al inicial)
    // 0.0 = gap cerrado completamente
    return std::max(0.0, std::min(1.0, progress_ratio));
}

double calculate_diameter_shape(double bigger_diam, double lower_diam, double epsilon) {
    // Validaciones básicas
    if (lower_diam < epsilon || !std::isfinite(lower_diam)) return 0.0;
    if (bigger_diam < epsilon || !std::isfinite(bigger_diam)) return 0.0;
    
    // Si bigger < lower, hay un error en los datos
    if (bigger_diam < lower_diam) return 0.0;
    
    double diam_ratio = bigger_diam / lower_diam;
    
    // Caja "cuadrada" (ratio cercano a 1)
    if (diam_ratio <= 1.0 + epsilon) return 0.0;
    
    // Usar log para comprimir valores grandes
    double log_diam_ratio = std::log10(diam_ratio);
    
    // Normalizar a [0, 1]
    // Asumiendo que ratios típicos están entre 1 y 1000
    // log10(1) = 0, log10(1000) = 3
    return std::max(0.0, std::min(1.0, log_diam_ratio / 3.0));
}

double calculate_diameter_progress(double current_diam, double epsilon, int id, bool is_bigger_diam) {
    static double initial_bigger_diam = 0.0, initial_lower_diam = 0.0;
    static bool bigger_set = false, lower_set = false;
    
    // Primera llamada: inicializar
    if (id == 1) {
        if (is_bigger_diam) {
            initial_bigger_diam = current_diam;
            bigger_set = true;
        } else {
            initial_lower_diam = current_diam;
            lower_set = true;
        }
        return 0.0; // Sin reducción aún
    }
    
    // Verificar que está inicializado
    if ((is_bigger_diam && !bigger_set) || (!is_bigger_diam && !lower_set)) {
        return 0.5; // Valor neutral
    }
    
    double initial_value = is_bigger_diam ? initial_bigger_diam : initial_lower_diam;
    
    // Si el valor inicial es demasiado pequeño, no podemos medir
    if (initial_value < epsilon) {
        return 0.5; // Valor neutral
    }
    
    // Validar diámetro actual
    if (!std::isfinite(current_diam) || current_diam < 0) {
        return 0.5;
    }
    
    // Medir REDUCCIÓN en lugar de PROPORCIÓN
    // reduction = 1 - (current / initial)
    // 0.0 = sin reducción
    // 1.0 = reducción completa (current → 0)
    double reduction = 1.0 - (current_diam / initial_value);
    
    // Clampear entre [0, 1]
    return std::max(0.0, std::min(1.0, reduction));
}

Optimizer::Status Optimizer::optimize() {
	Timer timer;
	timer.start();

	update_uplo();

	try {

		/*****************************************
		 ** DECLARACIÓN DE BUFFER Y LOUP FINDER **
		 *****************************************/

		CellBeamSearch * thebuffer = dynamic_cast<CellBeamSearch*>(&buffer);
		LoupFinderDefault * lfd = dynamic_cast<LoupFinderDefault*>(&loup_finder);

		//vector of cells
		queue<Cell*> aux; 

		//variable auxiliar para guardar el valor de la caja "inicial"
		IntervalVector aux_box = thebuffer->top()->box;
		double prec = 1e-7;

		/***********************************
		 ** DECLARACIÓN DE LOS BISECTORES **
		 ***********************************/
		OptimLargestFirst bisector_olf(goal_var, true, prec, 0.5);
		RoundRobin bisector_rr(prec, 0.5);

		// obtención del sistema para bisectores smear
		System system = lfd->finder_x_taylor.sys;
		
		// upper lowerbound
		double aux_uplo = uplo;
		// lower upperbound
		double aux_loup = loup;
		// región interior donde se saca el uplo y loup
		IntervalVector inner=aux_box;
		
		// simulaciones son la cantidad de datos que se generarán
		int num_sim = 1000;
		double epsilon = 1e-9;

		// se ingresa al vector la caja inicial del buffer (primer nodo a tratar)
		aux.push(thebuffer->top());
		
		// se obtiene la celda actual para el análisis
		// aqui calculamos los bounds iniciales para posteriormente calcular una proporción
		double lb_f_obj_inicial = lfd->finder_x_taylor.sys.goal->eval(aux.front()->box).lb();
		double ub_f_obj_inicial = lfd->finder_x_taylor.sys.goal->eval(aux.front()->box).ub();
		double bigger_diam_inicial = aux.front()->box.max_diam();
		double lower_diam_inicial = aux.front()->box.min_diam();		

		std::cout << "lb_f_obj_inicial: " << lb_f_obj_inicial << std::endl;
		std::cout << "ub_f_obj_inicial: " << ub_f_obj_inicial << std::endl;
		std::cout << "bigger_diam_inicial: " << bigger_diam_inicial << std::endl;
		std::cout << "lower_diam_inicial: " << lower_diam_inicial << std::endl;
		
		// estas variables son las que utilizaremos para calcular los límites dentro de las simulaciones.
		double lb_f_obj;
		double ub_f_obj;
		double bigger_diam;
		double lower_diam;
		
		// Variables finales que se utilizarán como input en la red neuronal
		int variables = n;
		BitSet active;
		double ratio_bounds;
		double ratio_bigger_diam;
		double ratio_lower_diam;
		double box_shape;
		
		// Total de restricciones del sistema
		int total_constraints = system.nb_ctr;
		
		std::ofstream InputFile("/home/felipe/Documents/magister/model2/input/prueba_nuevo_dataset/input_3_heuristicas.txt", std::ios::app);
		std::ofstream OutputFile("/home/felipe/Documents/magister/model2/output/prueba_nuevo_dataset/output_3_heuristicas.txt", std::ios::app);

		if (!InputFile.is_open()) {
			cerr << "No se pudo abrir el archivo. Comprueba la ruta y permisos." << endl;
			exit(1);
		}

		if (!OutputFile.is_open()) {
			cerr << "No se pudo abrir el archivo. Comprueba la ruta y permisos." << endl;
			exit(1);
		}

		// se ingresa la cantidad de simulaciones que se realizarán
		for (int k = 0 ; k < num_sim ; k++){

			// si el tamaño del vector es 0, se sale del ciclo
			if(aux.size() == 0) break;

			// El nodo que vas a evaluar (NO lo elimines todavía)
			Cell* nodo_original = aux.front();

			// Calcular las métricas del NODO ORIGINAL
			lb_f_obj = lfd->finder_x_taylor.sys.goal->eval(nodo_original->box).lb();
			ub_f_obj = lfd->finder_x_taylor.sys.goal->eval(nodo_original->box).ub();
			bigger_diam = nodo_original->box.max_diam();
			lower_diam = nodo_original->box.min_diam();
			active = lfd->finder_x_taylor.sys.active_ctrs(nodo_original->box);
			
			// calculamos los ratios de las 4 variables numéricas
			ratio_bounds = calculate_bounds_ratios(lb_f_obj, ub_f_obj, epsilon, k+1);
			ratio_bigger_diam = calculate_diameter_progress(bigger_diam, epsilon, k+1, true);
			ratio_lower_diam = calculate_diameter_progress(lower_diam, epsilon, k+1, false);
			box_shape = calculate_diameter_shape(bigger_diam, lower_diam, epsilon);

			// Calcular proporción de restricciones activas
			double constraint_proportion = (total_constraints > 0) ? 
				(active.size() / (double)total_constraints) : 0.0;

			InputFile << "variables: " << variables << endl;
			InputFile << "restricciones: " << constraint_proportion << endl;
			InputFile << "ratio_bounds: " << ratio_bounds << endl;
			InputFile << "ratio_bigger_diam: " << ratio_bigger_diam << endl;
			InputFile << "ratio_lower_diam: " << ratio_lower_diam << endl;
			InputFile << "box_shape: " << box_shape << endl;
			InputFile << "id: " << k+1 << endl << endl;

			/***********************************************************************
			 * Se escriben en el archivo los datos para el input de la red neuronal *
			 ***********************************************************************/
			
			// para cada técnica (LSMEAR=0, LF=1, RR=2)
			for (int i = 0 ; i < 3 ; i++){

				// CRÍTICO: COPIAR el nodo original para cada técnica
				Cell* copia_nodo = new Cell(*nodo_original);

				// se asignan los nuevos valores frontera
				uplo = aux_uplo;
				loup = aux_loup;

				// aun no se han bisectado nodos
				nb_cells = 0;

				// Empujar LA COPIA al buffer
				thebuffer->push(copia_nodo);
				
				// variable para saber si la técnica se registró o no
				bool technique_registered = false;
				
				// mientras queden nodos por revisar
				while (!thebuffer->empty()) {

					// el loup no ha cambiado
					loup_changed=false;

					// for double heap , choose randomly the buffer : top  has to be called before pop
					// celda "padre" (la que se está revisando)
					Cell *c = thebuffer->top();

					if (trace >= 2) std::cout << " current box " << c->box << endl;

					try {
						// se crea un par de celdas que corresponden a los nodos "hijos" con cada técnica
						pair<Cell*,Cell*> new_cells;

						// comienza la bisección dependiendo cada técnica
						if (i == 0) //lsmear
							new_cells=bsc.bisect(*c);
						else if (i == 1) //lf
							new_cells=bisector_olf.bisect(*c);
						else if (i == 2) //rr
							new_cells=bisector_rr.bisect(*c);

						// se elimina el nodo "padre" de la lista de nodos por revisar
						thebuffer->pop();
						
						// Se elimina el nodo "padre" de la memoria (que es una copia)
						delete c;

						nb_cells+=2;  // counting the cells handled
						
						// se manejan las celdas hijas
						handle_cell(*new_cells.first);
						handle_cell(*new_cells.second);

						// se revisa si ya no hay más nodos por revisar (deadend alcanzado)
						if(thebuffer->futurebuffer.size() == 0){
							// Registrar la técnica y su costo
							if (i == 0){
								OutputFile << "LSMEAR ";
							}
							else if (i == 1){
								OutputFile <<"LF ";
							}
							else if (i == 2){
								OutputFile << "RR ";
							}
							OutputFile << nb_cells << endl;
							technique_registered = true;

							// Procesar los nodos restantes en el buffer
							int buffer_size = thebuffer->size();
							for (int tt = 0 ; tt < buffer_size ; tt++){
								Cell *temp = thebuffer->top();
								thebuffer->pop();
								
								if(i == 0){ 
									// LSMEAR: se copian los nodos al auxiliar para futuras evaluaciones
									aux.push(temp);
								}
								else{
									// LF y RR: se eliminan los nodos generados
									delete temp;
								}
							}
							
							break; // Salir del while
						}

						if (uplo_of_epsboxes == NEG_INFINITY) {
							break;
						}
						if (loup_changed) {
							// In case of a new upper bound (loup_changed == true), all the boxes
							// with a lower bound greater than (loup - goal_prec) are removed and deleted.
							double ymax=compute_ymax();
							thebuffer->contract(ymax);

							// TODO: check if happens. What is the return code in this case?
							if (ymax <= NEG_INFINITY) {
								if (trace) std::cout << " infinite value for the minimum " << endl;
								break;
							}
						}
						update_uplo();

						if (!anticipated_upper_bounding) // useless to check precision on objective if 'true'
							if (get_obj_rel_prec()<rel_eps_f || get_obj_abs_prec()<abs_eps_f)
								break;

						if (timeout>0) timer.check(timeout);
						time = timer.get_time();

					}
					catch (NoBisectableVariableException& ) {
						update_uplo_of_epsboxes((c->box)[goal_var].lb());
						thebuffer->pop();
						delete c;
						update_uplo(); // the heap has changed -> recalculate the uplo
					}

				}

				// en caso de no registrarse la técnica en el deadend
				// se registra aquí (por término de precisión, factibilidad, etc)
				if (!technique_registered){
					if (i == 0){
						OutputFile << "LSMEAR ";
					}
					else if (i == 1){
						OutputFile <<"LF ";
					}
					else if (i == 2){
						OutputFile << "RR ";
					}
					OutputFile << nb_cells << endl;
					
					// Limpiar cualquier nodo restante en el buffer
					while(thebuffer->size() > 0){
						Cell* temp = thebuffer->top();
						thebuffer->pop();
						
						if(i == 0){
							// LSMEAR: guardar para futuras evaluaciones
							aux.push(temp);
						}
						else {
							// LF y RR: eliminar
							delete temp;
						}
					}
				}
			}
			
			OutputFile << "id: " << k+1 << endl << endl;
			
			// AHORA SÍ eliminar el nodo original después de evaluar las 3 técnicas
			aux.pop();
			delete nodo_original;
		}

		InputFile << "------------------------------------------------"  << endl;
		InputFile.close();
		
		OutputFile << "------------------------------------------------" << endl;
		OutputFile.close();
		
		// LIMPIEZA FINAL: Vaciar aux completamente para evitar fugas de memoria
		while(!aux.empty()) {
			Cell *c = aux.front();
			aux.pop();
			delete c;
		}
		
	 	timer.stop();
	 	time = timer.get_time();

		// No solution found and optimization stopped with empty buffer
		// before the required precision is reached => means infeasible problem
	 	if (uplo_of_epsboxes == NEG_INFINITY)
	 		status = UNBOUNDED_OBJ;
	 	else if (uplo_of_epsboxes == POS_INFINITY && (loup==POS_INFINITY || (loup==initial_loup && abs_eps_f==0 && rel_eps_f==0)))
	 		status = INFEASIBLE;
	 	else if (loup==initial_loup)
	 		status = NO_FEASIBLE_FOUND;
	 	else if (get_obj_rel_prec()>rel_eps_f && get_obj_abs_prec()>abs_eps_f)
	 		status = UNREACHED_PREC;
	 	else
	 		status = SUCCESS;
	}

	catch (TimeOutException& ) {
		status = TIME_OUT;
	}

	/* TODO: cannot retrieve variable names here. */
	for (int i=0; i<(extended_COV ? n+1 : n); i++)
		cov->data->_optim_var_names.push_back(string(""));

	cov->data->_optim_optimizer_status = (unsigned int) status;
	cov->data->_optim_uplo = uplo;
	cov->data->_optim_uplo_of_epsboxes = uplo_of_epsboxes;
	cov->data->_optim_loup = loup;

	cov->data->_optim_time += time;
	cov->data->_optim_nb_cells += nb_cells;
	cov->data->_optim_loup_point = loup_point;

	// for conversion between original/extended boxes
	IntervalVector tmp(extended_COV ? n+1 : n);

	// by convention, the first box has to be the loup-point.
	if (extended_COV) {
		write_ext_box(loup_point, tmp);
		tmp[goal_var] = Interval(uplo,loup);
		cov->add(tmp);
	}
	else {
		cov->add(loup_point);
	}

	while (!buffer.empty()) {
		Cell* cell=buffer.top();
		if (extended_COV) {
			cov->add(cell->box);
		} else {
			read_ext_box(cell->box,tmp);
			cov->add(tmp);
		}
		delete buffer.pop();
	}

	return status;
}
namespace {
const char* green() {
#ifndef _WIN32
	return "\033[32m";
#else
	return "";
#endif
}

const char* red(){
#ifndef _WIN32
	return "\033[31m";
#else
	return "";
#endif
}

const char* white() {
#ifndef _WIN32
	return "\033[0m";
#else
	return "";
#endif
}

}

void Optimizer::report() {

	// if (!cov || !buffer.empty()) { // not started
	// 	cout << " not started." << endl;
	// 	return;
	// }

	// switch(status) {
	// case SUCCESS:
	// 	cout << green() << " optimization successful!" << endl;
	// 	break;
	// case INFEASIBLE:
	// 	cout << red() << " infeasible problem" << endl;
	// 	break;
	// case NO_FEASIBLE_FOUND:
	// 	cout << red() << " no feasible point found (the problem may be infeasible)" << endl;
	// 	break;
	// case UNBOUNDED_OBJ:
	// 	cout << red() << " possibly unbounded objective (f*=-oo)" << endl;
	// 	break;
	// case TIME_OUT:
	// 	cout << red() << " time limit " << timeout << "s. reached " << endl;
	// 	break;
	// case UNREACHED_PREC:
	// 	cout << red() << " unreached precision" << endl;
	// 	break;
	// }
	// cout << white() <<  endl;

	// // No solution found and optimization stopped with empty buffer
	// // before the required precision is reached => means infeasible problem
	// if (status==INFEASIBLE) {
	// 	cout << " infeasible problem " << endl;
	// } else {
	// 	cout << " f* in\t[" << uplo << "," << loup << "]" << endl;
	// 	cout << "\t(best bound)" << endl << endl;

	// 	if (loup==initial_loup)
	// 		cout << " x* =\t--\n\t(no feasible point found)" << endl;
	// 	else {
	// 		if (loup_finder.rigorous())
	// 			cout << " x* in\t" << loup_point << endl;
	// 		else
	// 			cout << " x* =\t" << loup_point.lb() << endl;
	// 		cout << "\t(best feasible point)" << endl;
	// 	}
	// 	cout << endl;
	// 	double rel_prec=get_obj_rel_prec();
	// 	double abs_prec=get_obj_abs_prec();

	// 	cout << " relative precision on f*:\t" << rel_prec;
	// 	if (rel_prec <= rel_eps_f)
	// 		cout << green() << " [passed] " << white();
	// 	cout << endl;

	// 	cout << " absolute precision on f*:\t" << abs_prec;
	// 	if (abs_prec <= abs_eps_f)
	// 		cout << green() << " [passed] " << white();
	// 	cout << endl;
	// }

	// cout << " cpu time used:\t\t\t" << time << "s";
	// if (cov->time()!=time)
	// 	cout << " [total=" << cov->time() << "]";
	// cout << endl;
	// cout << " number of cells:\t\t" << nb_cells;
	// if (cov->nb_cells()!=nb_cells)
	// 	cout << " [total=" << cov->nb_cells() << "]";
	// cout << endl << endl;

	// if (statistics)
	// 	cout << "  ===== Statistics ====" << endl << endl << *statistics << endl;
	cout << nb_cells << " " << time << endl;
}



} // end namespace ibex
