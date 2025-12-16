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



// ===============================================================================
//   ESTRUCTURAS Y FUNCIONES HELPER PARA GENERACIÓN DE DATOS (DATA SCIENCE)
// ===============================================================================

// Estructura para capturar la "foto" del nodo raíz/inicial y comparar progreso
struct InitialState {
    double lb_obj;
    double ub_obj;
    double obj_range;
    double max_diam;
    double min_diam;
    bool valid_bounds; 

    InitialState(double lb, double ub, double max_d, double min_d) 
        : lb_obj(lb), ub_obj(ub), max_diam(max_d), min_diam(min_d) {
        
        obj_range = ub - lb;
        // Validamos si tenemos un rango finito y útil para comparar
        valid_bounds = std::isfinite(lb) && std::isfinite(ub) && 
                       std::isfinite(obj_range) && std::abs(obj_range) > 1e-12;
    }
};

// ===============================================================================
//   FUNCIONES HELPER CORREGIDAS (LOG SCALE)
// ===============================================================================

// Calcula LOG10 de la reducción del rango objetivo
// Retorna: 0.0 (Sin cambio) hasta -30.0 (Reducción masiva)
double calculate_bounds_log_ratio(double current_lb, double current_ub, const InitialState& init) {
    if (!init.valid_bounds) return 0.0; 

    double current_range = current_ub - current_lb;
    
    // Si el rango es inválido o negativo, retornamos 0 (sin progreso)
    if (!std::isfinite(current_range) || current_range < 0) return 0.0;

    // Evitamos división por cero y log(0)
    double ratio = current_range / (init.obj_range + 1e-100); 
    
    // Si el ratio es > 1 (raro, el rango creció), lo acotamos a 1
    if (ratio > 1.0) ratio = 1.0;
    
    // Si el ratio es extremadamente pequeño (casi cero), le ponemos un piso para no dar -inf
    if (ratio < 1e-30) ratio = 1e-30;

    return std::log10(ratio); // Retornará valores entre 0.0 y -30.0
}

// Calcula LOG10 del progreso de reducción de diámetro
double calculate_diam_log_reduction(double current_diam, double initial_diam) {
    if (initial_diam < 1e-12 || !std::isfinite(initial_diam)) return 0.0; 
    
    double ratio = current_diam / initial_diam;

    if (ratio > 1.0) ratio = 1.0;
    if (ratio < 1e-30) ratio = 1e-30; // Piso de seguridad

    return std::log10(ratio); // Retornará valores entre 0.0 y -30.0
}

// La función de Shape se mantiene igual (ya usaba log10)
double calculate_box_shape(double max_diam, double min_diam, double epsilon = 1e-9) {
    if (min_diam < epsilon || !std::isfinite(min_diam)) return 5.0; 
    if (max_diam < epsilon || !std::isfinite(max_diam)) return 0.0; 
    
    double ratio = max_diam / min_diam;
    double log_ratio = std::log10(ratio);
    
    return std::max(0.0, std::min(5.0, log_ratio));
}

// Calculamos que tan lejos está esta caja de la mejor solución encontrada
// retorna -1.0 (incertidumbre total/-inf), 0.0 caja muy cerca del óptimo
double calculate_relative_gap(double box_lb, double current_loup) {
	//caso 1: no tenemos cota inferior
	if (box_lb <= -1e15 || !std::isfinite(box_lb)) {
		return 1.0;
	}

	//caso 2: no hay solución global aún (loup infinito)
	//el gap es 1.0 porque no tenemos contra qué comparar
	if (current_loup >= 1e15 || !std::isfinite(current_loup)) {
		return 1.0;
	}

	//caso 3: calculo real
	//gap = (loup - lb) / |loup + epsilon|
	//usamos valor absoluto en el denominador para evitar cambios de signo raros
	double gap = (current_loup - box_lb) / (std::abs(current_loup) + 1.0);

	//normalizacion 
	if (gap > 1.0) gap = 1.0; // no deberia pasar si lb < loup
	if (gap < 0.0) gap = 0.0; // si lb > loup (nodo podable), el gap es 0.0

	return gap;

}

Optimizer::Status Optimizer::optimize() {
    Timer timer;
    timer.start();

    update_uplo();

    try {

        /*****************************************
         ** CONFIGURACIÓN INICIAL           **
         *****************************************/

        CellBeamSearch * thebuffer = dynamic_cast<CellBeamSearch*>(&buffer);
        LoupFinderDefault * lfd = dynamic_cast<LoupFinderDefault*>(&loup_finder);

        queue<Cell*> aux; 

        // 1. Llenar la cola auxiliar con lo que haya en el buffer inicial
        if (!thebuffer->empty()) {
            while(!thebuffer->empty()) {
                aux.push(thebuffer->pop());
            }
        } else if (!buffer.empty()) {
             // Fallback defensivo
             aux.push(thebuffer->top()); 
        }

        // 2. CONFIGURACIÓN DEL ESTADO INICIAL (REEMPLAZA A STATIC)
        // Capturamos el estado del primer nodo disponible para usarlo como línea base
        // de comparación para todos los nodos subsiguientes (progreso relativo al inicio).
        InitialState* global_root_state = nullptr;

        if (!aux.empty()) {
            Cell* root = aux.front();
            double r_lb = lfd->finder_x_taylor.sys.goal->eval(root->box).lb();
            double r_ub = lfd->finder_x_taylor.sys.goal->eval(root->box).ub();
            double r_max = root->box.max_diam();
            double r_min = root->box.min_diam();
            
            global_root_state = new InitialState(r_lb, r_ub, r_max, r_min);
        } else {
            // Caso borde: cola vacía, creamos dummy para no romper el código
            global_root_state = new InitialState(0, 0, 1, 1);
        }

        double prec = 1e-7;

        /***********************************
         ** BIS ECTORES Y ARCHIVOS      **
         ***********************************/
        OptimLargestFirst bisector_olf(goal_var, true, prec, 0.5);
        RoundRobin bisector_rr(prec, 0.5);
        System system = lfd->finder_x_taylor.sys;

		// -- Heurísticas Smear --
		SmearMax bisector_sm(system,prec);
		SmearSum bisector_ss(system,prec);
		SmearSumRelative bisector_ssr(system,prec);
        
        // Variables globales backup
        double aux_uplo = uplo;
        double aux_loup = loup;
        
        int num_sim = 1000;
        double epsilon = 1e-9;

        // Archivos
        std::ofstream InputFile("/home/felipe/Documents/magister/model2/input/prueba_nuevo_dataset/input_15_dic_ratios_depth.txt", std::ios::app);
        std::ofstream OutputFile("/home/felipe/Documents/magister/model2/output/prueba_nuevo_dataset/output_15_dic_ratios_depth.txt", std::ios::app);

        if (!InputFile.is_open() || !OutputFile.is_open()) {
            cerr << "Error abriendo archivos de texto." << endl; exit(1);
        }

        // =========================================================
        //                 BUCLE DE SIMULACIONES (K)
        // =========================================================
        for (int k = 0 ; k < num_sim ; k++){

            if(aux.empty()) break;

            // --- PASO 1: OBTENER SEMILLA (SIN SACAR DE LA COLA) ---
            Cell* seed_cell = aux.front();

            // Obtener valores actuales
            double cur_lb = lfd->finder_x_taylor.sys.goal->eval(seed_cell->box).lb();
            double cur_ub = lfd->finder_x_taylor.sys.goal->eval(seed_cell->box).ub();
            double cur_max_diam = seed_cell->box.max_diam();
            double cur_min_diam = seed_cell->box.min_diam();
            BitSet active = lfd->finder_x_taylor.sys.active_ctrs(seed_cell->box);
            int variables = n;
			double depth_radio = (double)seed_cell->depth / ((double)variables + 1e-9);

            // Logs iniciales debug
            // if (k==0) {
            //     std::cout << "Initial Root - LB: " << cur_lb << " UB: " << cur_ub << std::endl;
            // }
            
            // --- PASO 2: CÁLCULO DE FEATURES EN ESCALA LOGARÍTMICA ---
            // Nota el cambio de nombre de las variables y funciones
            double log_ratio_bounds = calculate_bounds_log_ratio(cur_lb, cur_ub, *global_root_state);
            double log_ratio_bigger = calculate_diam_log_reduction(cur_max_diam, global_root_state->max_diam);
            double log_ratio_lower  = calculate_diam_log_reduction(cur_min_diam, global_root_state->min_diam);
            
            // Box shape ya está en log, así que está bien
            double box_shape = calculate_box_shape(cur_max_diam, cur_min_diam, epsilon);

            // Escribir Input (Actualiza los nombres en el txt si quieres ser explícito)
            InputFile << "variables: " << variables << endl;
            InputFile << "restricciones: " << active.size() << endl;
			InputFile << "depth_ratio: " << depth_radio << endl;
            InputFile << "box_shape: " << box_shape << endl;
            InputFile << "log_ratio_bounds: " << log_ratio_bounds << endl;      // Valor esperado: -0.5, -10.0, -25.0, etc.
            InputFile << "log_ratio_bigger: " << log_ratio_bigger << endl;
            InputFile << "log_ratio_lower: " << log_ratio_lower << endl;
            InputFile << "id: " << k+1 << endl << endl;

            // =====================================================
            //            BUCLE DE HEURÍSTICAS (I)
            // =====================================================
            for (int i = 0 ; i < 6 ; i++){

                // 1. Restaurar bounds globales (Justicia)
                uplo = aux_uplo;
                loup = aux_loup;
                nb_cells = 0;

                // 2. Limpiar buffer
                buffer.flush();

                // 3. CLONACIÓN (Fairness): Copia independiente de la caja semilla
				Cell* current_root = new Cell(*seed_cell); 
				thebuffer->push(current_root);

				// // ================== DEBUG BLOCK INICIO ==================
				// std::cout << std::setprecision(16); // Máxima precisión para ver diferencias mínimas
				// std::cout << "\n[DEBUG CHECK] Iteracion K=" << k << " Heuristica I=" << i << std::endl;

				// // 1. Verificación de Punteros (Deben ser DIFERENTES)
				// std::cout << "  Addr Semilla: " << seed_cell << std::endl;
				// std::cout << "  Addr Clon   : " << current_root << std::endl;
				// if (seed_cell != current_root) std::cout << "  -> MEMORIA: [OK] Son objetos distintos." << std::endl;
				// else std::cout << "  -> MEMORIA: [ERROR CRITICO] Es el mismo puntero!" << std::endl;

				// // 2. Verificación de Contenido (Deben ser IDÉNTICOS)
				// // IBEX permite comparar IntervalVectors directamente, pero vamos a mirar la distancia
				// double diff = distance(seed_cell->box,current_root->box); // Distancia entre cajas
				// std::cout << "  Diferencia de Cajas: " << diff << std::endl;

				// if (diff == 0.0) std::cout << "  -> CONTENIDO: [OK] Las cajas son idénticas bit a bit." << std::endl;
				// else std::cout << "  -> CONTENIDO: [ERROR] La caja clonada es diferente." << std::endl;

				// // 3. Verificación de Estado Global (LOUP debe estar reseteado)
				// std::cout << "  Loup Actual: " << loup << " | Loup Esperado: " << aux_loup << std::endl;
				// if (loup == aux_loup) std::cout << "  -> GLOBALES: [OK] Loup reseteado." << std::endl;
				// else std::cout << "  -> GLOBALES: [ERROR] El loup está sucio." << std::endl;
				// // ================== DEBUG BLOCK FIN ==================

                bool first_iteration = true;
                
                // --- Diving Loop ---
                while (!thebuffer->empty()) {

                    loup_changed=false;
                    Cell *c = thebuffer->top();

                    try {
                        pair<Cell*,Cell*> new_cells;

                        if (i == 0) new_cells=bsc.bisect(*c);           // LSMEAR
                        if (i == 1) new_cells=bisector_olf.bisect(*c);  // LF
                        if (i == 2) new_cells=bisector_rr.bisect(*c);   // RR
						if (i == 3) new_cells=bisector_sm.bisect(*c);   // SMEAR MAX
						if (i == 4) new_cells=bisector_ss.bisect(*c);   // SMEAR SUM
						if (i == 5) new_cells=bisector_ssr.bisect(*c);  // SMEAR SUM RELATIVE

                        thebuffer->pop();

                        // Gestión de memoria nodos procesados
                        if (first_iteration){
                            first_iteration = false;
                            delete c; // Borramos el clon root
                        } else {
                            delete c; // Borramos nodos intermedios
                        } 

                        nb_cells += 2;
                        handle_cell(*new_cells.first);
                        handle_cell(*new_cells.second);

                        // --- DEADEND / CLASIFICACIÓN ---
                        if(thebuffer->futurebuffer.size() == 0){ 
                            
                            if (i == 0) OutputFile << "LSMEAR: ";
                            else if (i == 1) OutputFile <<"LF: ";
                            else if (i == 2) OutputFile << "RR: ";
							else if (i == 3) OutputFile << "SM: ";
							else if (i == 4) OutputFile << "SS: ";
							else if (i == 5) OutputFile << "SSR: ";
                            // ...
                            
                            OutputFile << nb_cells << endl;
                            
                            // Gestión diferenciada de hijos
                            int current_size = thebuffer->size();
                            if(i == 0) { 
                                // LSMEAR define el futuro: guardar hijos en aux
                                for (int tt = 0 ; tt < current_size ; tt++){
                                    Cell *survivor = thebuffer->top();
                                    thebuffer->pop();
                                    aux.push(survivor); 
                                }
                            } else {
                                // Otros: solo pruebas, borrar hijos
                                for (int tt = 0 ; tt < current_size ; tt++){
                                    Cell *trash = thebuffer->top();
                                    thebuffer->pop();
                                    delete trash;
                                }
                            }
                        } 

                        // Chequeos IBEX
                        if (uplo_of_epsboxes == NEG_INFINITY) break;
                        if (loup_changed) {
                            double ymax=compute_ymax();
                            thebuffer->contract(ymax);
                            if (ymax <= NEG_INFINITY) break;
                        }
                        update_uplo();
                        if (timeout>0) timer.check(timeout);

                    } catch (NoBisectableVariableException& ) {
                        update_uplo_of_epsboxes((c->box)[goal_var].lb());
                        thebuffer->pop();
                        if(nb_cells!=0) delete c;
                        update_uplo();
                    }
                } // Fin While Diving

            } // Fin For Heurísticas

            OutputFile << "id: " << k+1 << endl << endl;

            // --- PASO 3: LIMPIEZA FINAL DE SEMILLA ---
            aux.pop();      // Sacar de la cola
            delete seed_cell; // Liberar memoria

        } // Fin For Simulaciones

        // Limpieza de memoria de Data Science
        if (global_root_state) delete global_root_state;

        InputFile << "------------------------------------------------"  << endl;
        InputFile.close();
        OutputFile << "------------------------------------------------" << endl;
        OutputFile.close();
        
        timer.stop();
        time = timer.get_time();

        // Estado final
        if (uplo_of_epsboxes == NEG_INFINITY) status = UNBOUNDED_OBJ;
        else if (uplo_of_epsboxes == POS_INFINITY && (loup==POS_INFINITY || (loup==initial_loup && abs_eps_f==0 && rel_eps_f==0))) status = INFEASIBLE;
        else if (loup==initial_loup) status = NO_FEASIBLE_FOUND;
        else if (get_obj_rel_prec()>rel_eps_f && get_obj_abs_prec()>abs_eps_f) status = UNREACHED_PREC;
        else status = SUCCESS;
    }

    catch (TimeOutException& ) {
        status = TIME_OUT;
    }

    // Reporte final para COV (sin cambios)
    for (int i=0; i<(extended_COV ? n+1 : n); i++)
        cov->data->_optim_var_names.push_back(string(""));

    cov->data->_optim_optimizer_status = (unsigned int) status;
    cov->data->_optim_uplo = uplo;
    cov->data->_optim_uplo_of_epsboxes = uplo_of_epsboxes;
    cov->data->_optim_loup = loup;
    cov->data->_optim_time += time;
    cov->data->_optim_nb_cells += nb_cells;
    cov->data->_optim_loup_point = loup_point;

    IntervalVector tmp(extended_COV ? n+1 : n);

    if (extended_COV) {
        write_ext_box(loup_point, tmp);
        tmp[goal_var] = Interval(uplo,loup);
        cov->add(tmp);
    } else {
        cov->add(loup_point);
    }

    while (!buffer.empty()) {
        Cell* cell=buffer.top();
        if (extended_COV) cov->add(cell->box);
        else {
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
