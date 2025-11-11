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
#include "ibex_LoupFinderDefault.h"

#include "ibex_OptimLargestFirst.h"
#include "ibex_RoundRobin.h"
#include "ibex_SmearFunction.h"

#include <float.h>
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
			std::cout << "                    ";
			std::cout << "\033[32m loup= " << loup << "\033[0m" << endl;
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

	loup=obj_init_bound;

	// Just to initialize the "loup" for the buffer
	// TODO: replace with a set_loup function
	buffer.contract(loup);

	uplo=NEG_INFINITY;
	uplo_of_epsboxes=POS_INFINITY;

	nb_cells=0;

	buffer.flush();

	Cell* root=new Cell(IntervalVector(n+1));

	write_ext_box(init_box, root->box);

	// add data required by the bisector
	bsc.add_property(init_box, root->prop);

	// add data required by the contractor
	ctc.add_property(init_box, root->prop);

	// add data required by the buffer
	buffer.add_property(init_box, root->prop);

	// add data required by the loup finder
	loup_finder.add_property(init_box, root->prop);

	//cout << "**** Properties ****\n" << root->prop << endl;

	loup_changed=false;
	initial_loup=obj_init_bound;

	loup_point = init_box; //.set_empty();
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

// double calculate_bounds_ratios(double lb, double ub, double epsilon, int id) {
//     static double initial_lb, initial_ub, initial_obj_range;
//     static bool initial_bounds_set = false;
    
//     // Si es la primera llamada (id == 1), guardar valores iniciales
//     if (id == 1) {
//         initial_lb = lb;
//         initial_ub = ub;
//         initial_obj_range = initial_ub - initial_lb;
//         initial_bounds_set = true;
//     }
    
//     // Verificar que se hayan establecido los valores iniciales
//     if (!initial_bounds_set) {
//         return -5.0; // Error: no se han establecido valores iniciales
//     }
    
//     // Calcular rango actual
//     double current_obj_range = ub - lb;
    
//     // Edge cases para bounds iniciales
//     bool initial_bounds_infinite = (!std::isfinite(initial_lb) || !std::isfinite(initial_ub));
//     bool initial_bounds_invalid = (!std::isfinite(initial_obj_range) || std::abs(initial_obj_range) < epsilon);
    
//     // Edge cases para bounds actuales
//     bool current_bounds_invalid = (!std::isfinite(current_obj_range) || std::abs(current_obj_range) < epsilon);
    
//     // Caso 1: Current bounds inválidos
//     if (current_bounds_invalid) {
//         return -5.0; // Valor sentinel
//     }
    
//     // Caso 2: Bounds iniciales problemáticos - usar solo rango actual
//     if (initial_bounds_infinite || initial_bounds_invalid) {
//         double log_current_range = std::log10(current_obj_range);
//         return std::max(-5.0, std::min(5.0, log_current_range));
//     }
    
//     // Caso 3: Normal - calcular progreso relativo
//     double progress_ratio = current_obj_range / initial_obj_range;
//     return std::max(0.0, std::min(1.0, progress_ratio));
// }


// double calculate_diameter_shape(double bigger_diam, double lower_diam, double epsilon) {
//     // Edge cases
//     if (lower_diam < epsilon || !std::isfinite(lower_diam)) return 0.0;
//     if (bigger_diam < epsilon || !std::isfinite(bigger_diam)) return 0.0;
    
//     double diam_ratio = bigger_diam / lower_diam;
//     if (diam_ratio <= 1.0 + epsilon) return 0.0; // Box cuadrado
    
//     double log_diam_ratio = std::log10(diam_ratio);
//     return std::max(-3.0, std::min(3.0, log_diam_ratio));
// }

// double calculate_diameter_progress(double current_diam, double epsilon, int id, bool is_bigger_diam) {
//     static double initial_bigger_diam, initial_lower_diam;
//     static bool initial_set = false;
    
//     if (id == 1) {
//         if (is_bigger_diam) {
//             initial_bigger_diam = current_diam;
//         } else {
//             initial_lower_diam = current_diam;
//         }
//         initial_set = true;
//         return 1.0; // Valor inicial
//     }
    
//     if (!initial_set) return 0.5;
    
//     double initial_value = is_bigger_diam ? initial_bigger_diam : initial_lower_diam;
//     if (initial_value < epsilon) return 0.5;
    
//     double progress = current_diam / initial_value;
//     return std::max(0.0, std::min(1.0, progress));
// }

double calculate_sum_diam(const ibex::IntervalVector& box) {
    if (box.is_empty()) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < box.size(); ++i) {
        sum += box[i].diam();
    }
    return sum;
}

Optimizer::Status Optimizer::optimize() {
	Timer timer;
	timer.start();

	update_uplo();

	try {

        LoupFinderDefault* lfd = dynamic_cast<LoupFinderDefault*>(&loup_finder);
        const System& const_sys = lfd->finder_x_taylor.sys;  //referencia const al sistema
        System& sys = const_cast<System&>(const_sys); // obtenemos una referencia NO CONST. El sistema del optimizer no es realmente const

        std::ofstream InputFile("/home/felipe/Documents/magister/model2/input/prueba_nuevo_dataset/input_prueba_fer_bis0_sampling.txt", std::ios::app);
		std::ofstream OutputFile("/home/felipe/Documents/magister/model2/output/prueba_nuevo_dataset/output_prueba_fer_bis0_sampling.txt", std::ios::app);

        if (!InputFile.is_open() || !OutputFile.is_open()) {
			std::cerr << "No se pudo abrir el archivo. Comprueba la ruta y permisos." << std::endl;
			exit(1);
		}   

        InputFile << "variables,active_ctrs,box_shape,max_diam_grad,avg_diam_grad,max_mag_grad,avg_mag_grad" << endl;
        OutputFile << "best_heuristic_index" << endl;

        double prec = 1e-7;

        OptimLargestFirst bisector_olf(goal_var, true, prec, 0.5);
		RoundRobin bisector_rr(prec, 0.5);
        SmearMax bisector_sm(sys,prec);
		SmearSum bisector_ss(sys,prec);
		SmearSumRelative bisector_ssr(sys,prec);

        Ctc *ctc_ptr = &ctc;

        std::vector<Bsc*> biselectors;
        biselectors.push_back(&bsc);            // 0: Default (LSmear)
        biselectors.push_back(&bisector_olf);    // 1: LargestFirst
        biselectors.push_back(&bisector_rr);    // 2: RoundRobin
        biselectors.push_back(&bisector_sm);    // 3: SmearMax
        biselectors.push_back(&bisector_ss);    // 4: SmearSum
        biselectors.push_back(&bisector_ssr);   // 5: SmearSumRelative

        int samples = 500;
        int samples_generated = 0;

		std::vector<std::pair<double, int>> heuristic_results;

		int variables;
		BitSet active_ctrs;
		double bigger_diam;
		double lower_diam;
		double box_shape;
		IntervalVector goal_grad;
		double max_diam_grad = 0.0, avg_diam_grad = 0.0, max_mag_grad = 0.0, avg_mag_grad = 0.0;
		double max_reward;
		int best_heuristic_index;
		double initial_pseudo_volume;

	    while (!buffer.empty() && samples_generated < samples) {

            Cell *c = buffer.top();
            const IntervalVector& box = c->box;
			
			// Este if nos permite hacer una especie de "sampling"
			// para que los datos generados no sean correlativos
			// y puedan ser mas variados, con la idea de 
			// poder representar mejor las estapas del B&B
			if (nb_cells % 20 == 0) {
				// --- 1. CÁLCULO DE FEATURES ---
				variables = n;
				active_ctrs = sys.active_ctrs(box);
				bigger_diam = box.max_diam();
				lower_diam = box.min_diam();
				// Calculamos box_shape. Es una feature robusta y fácil.
				box_shape = (bigger_diam > 1e-12) ? lower_diam / bigger_diam : 1.0;
				
				goal_grad = sys.goal->gradient(box);
				max_diam_grad = 0.0, avg_diam_grad = 0.0, max_mag_grad = 0.0, avg_mag_grad = 0.0;
				if (!goal_grad.is_empty() && n > 0) {
					double sum_diam = 0.0, sum_mag = 0.0;
					max_diam_grad = goal_grad[0].diam();
					max_mag_grad = goal_grad[0].mag();
					for (int i = 0; i < n; ++i) {
						double cd = goal_grad[i].diam(), cm = goal_grad[i].mag();
						sum_diam += cd; sum_mag += cm;
						if (cd > max_diam_grad) max_diam_grad = cd;
						if (cm > max_mag_grad) max_mag_grad = cm;
					}
					avg_diam_grad = sum_diam / n; avg_mag_grad = sum_mag / n;
				}
				
				// --- 2. EJECUCIÓN DEL ORÁCULO ---
				max_reward = -1.0;
				best_heuristic_index = -1;
				initial_pseudo_volume = calculate_sum_diam(box); // Asegúrese de que calculate_sum_diam existe
				for (int i = 0; i < biselectors.size(); ++i) {
					try {
						Cell* temp_cell = new Cell(c->box);
						//std::cout << "caja " << temp_cell->box << " " << i<< std::endl;
						pair<Cell*, Cell*> new_cells = biselectors[i]->bisect(*temp_cell);
						ctc.contract(new_cells.first->box);
						ctc.contract(new_cells.second->box);
						double final_pseudo_volume = calculate_sum_diam(new_cells.first->box) + calculate_sum_diam(new_cells.second->box);
						double reward = initial_pseudo_volume - final_pseudo_volume;
						heuristic_results.push_back({reward, i});
						if (reward > max_reward) { max_reward = reward; best_heuristic_index = i; }
						delete temp_cell; delete new_cells.first; delete new_cells.second;
					} catch (NoBisectableVariableException&) { heuristic_results.push_back({-1.0, i}); continue; }
				}

				// 3. Ordenar el vector en orden descendente de recompensa
				std::sort(heuristic_results.begin(), heuristic_results.end(), 
					// Lambda con tipos explícitos para compatibilidad
					[](const std::pair<double, int>& a, const std::pair<double, int>& b) {
						return a.first > b.first; // Compara el primer elemento del par (la recompensa)
					}
				);

				// --- 3. ESCRITURA DE DATOS ---
				if (!heuristic_results.empty() && heuristic_results[2].first >= 0) {
		
					// Escribir features (sin cambios)
					InputFile << variables << "," << active_ctrs.size() << "," << box_shape << ","
							<< max_diam_grad << "," << avg_diam_grad << "," << max_mag_grad << "," << avg_mag_grad << endl;
					
					// Escribir el ranking en el archivo de salida
					for (size_t i = 0; i < heuristic_results.size(); ++i) {
						OutputFile << heuristic_results[i].second << (i == heuristic_results.size() - 1 ? "" : ",");
					}
					OutputFile << endl;

					samples_generated++;
				}

				heuristic_results.clear();
			}

			// --- 4. AVANCE DEL B&B ---
			buffer.pop(); // POP ÚNICO: Sacamos 'c' del buffer ahora que hemos terminado con él.
			
            try {
                
				loup_changed = false;
                pair<Cell*, Cell*> new_cells = biselectors[0]->bisect(*c);
                
                nb_cells += 2;
                handle_cell(*new_cells.first);
                handle_cell(*new_cells.second);

                if (uplo_of_epsboxes == NEG_INFINITY) { break; }
                if (loup_changed) {
                    double ymax = compute_ymax();
                    buffer.contract(ymax);
                    if (ymax <= NEG_INFINITY) { break; }
                }
                update_uplo();
                if (!anticipated_upper_bounding && (get_obj_rel_prec() < rel_eps_f || get_obj_abs_prec() < abs_eps_f)) { break; }
                if (timeout > 0) { timer.check(timeout); time = timer.get_time(); }

            } catch (NoBisectableVariableException&) {
                update_uplo_of_epsboxes((c->box)[goal_var].lb());
                update_uplo();
            }
            
            delete c; // DELETE ÚNICO: Liberamos la memoria de 'c' al final de la iteración.
        }

        InputFile.close();
        OutputFile.close();

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
		if (extended_COV)
			cov->add(cell->box);
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

    std::cout << nb_cells << " " << time << endl;
}



} // end namespace ibex