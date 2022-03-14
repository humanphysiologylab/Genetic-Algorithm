#ifndef NELMIN
#define NELMIN
#include <cmath>
#include <string>
#include "penalty.h"

template <typename F, typename Vec>
void nelmin (const F & function, Vec start, Vec & xmin, std::vector<std::pair<int, double>> & error_per_gen,
  double reqmin, const Vec & step, int check_period, int max_func_eval_number,
  int & func_eval_number, int & restart_number, int & ifault, std::vector<int> & is_mutation_applicable)
{

//****************************************************************************80
//
//  Purpose:
//
//    NELMIN minimizes a function using the Nelder-Mead algorithm.
//
//  Discussion:
//
//    This routine seeks the minimum value of a user-specified function.
//
//    Simplex function minimisation procedure due to Nelder+Mead(1965),
//    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
//    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
//    25, 97) and Hill(1978, 27, 380-2)
//
//    The function to be minimized must be defined by a function of
//    the form
//
//      function function ( x, f )
//      double function
//      double x(*)
//
//    and the name of this subroutine must be declared EXTERNAL in the
//    calling routine and passed as the argument function.
//
//    This routine does not include a termination test using the
//    fitting of a quadratic surface.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by R ONeill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Nelder, Roger Mead,
//    A simplex method for function minimization,
//    Computer Journal,
//    Volume 7, 1965, pages 308-313.
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double function ( double x[] ), the name of the routine which evaluates
//    the function to be minimized.
//
//    Input, int N, the number of variables.
//
//    Input/output, double START[N].  On input, a starting point
//    for the iteration.  On output, this data may have been overwritten.
//
//    Output, double XMIN[N], the coordinates of the point which
//    is estimated to minimize the function.
//
//    Output, double highest_value, the minimum value of the function.
//
//    Input, double REQMIN, the terminating limit for the variance
//    of function values.
//
//    Input, double STEP[N], determines the size and shape of the
//    initial simplex.  The relative magnitudes of its elements should reflect
//    the units of the variables.
//
//    Input, int check_period, the convergence check is carried out
//    every check_period iterations.
//
//    Input, int max_func_eval_number, the maximum number of function
//    evaluations.
//
//    Output, int *func_eval_number, the number of function evaluations
//    used.
//
//    Output, int *restart_number, the number of restarts.
//
//    Output, int *IFAULT, error indicator.
//    0, no errors detected.
//    1, REQMIN, N, or check_period has an illegal value.
//    2, iteration terminated because max_func_eval_number was exceeded without convergence.
//
  const int n = start.size(); //dimensionality
  const double dn = n; //just double representation of n
  const int nn = n + 1; //just n + 1 really (like a number of simplex vertices
  const double dnn = nn; //just double for nn

  const double eps = 0.001; //??
  const double rcoeff = 1.0; //reflection coefficient

  //default coeffs
  /*
  const double ccoeff = 0.5; //contraction coefficient
  const double ecoeff = 2.0; //expansion coefficient
  const double shrink_coeff = 0.5; //shrink coefficient
  */
  //fancy coeffs

  const double ccoeff = 0.75 - 1.0 / (2 * n); //contraction coefficient
  const double ecoeff = 1 + 2.0 / n; //expansion coefficient
  const double shrink_coeff = 1 - 1.0 / n; //shrink coefficient


  double initial_simplex_scale = 1;
  double eps_Scale = 0.1;
//
//  Check the input parameters.
//
  if (reqmin <= 0) {
    ifault = 1;
    return;
  }

  if (n < 1) {
    ifault = 1;
    return;
  }

  if (check_period < 1) {
    ifault = 1;
    return;
  }

  // p stores (n + 1) simplex vertices in row-major order
  Vec p(n * (n + 1));


  // function value at corresponding vertices
  Vec y(n + 1);

  func_eval_number = 0;
  restart_number = 0;


//
//  Initial or restarted loop.
//
  while (1) {
    //fill p with n+1 initial simplex vertices
    //and find function values
#pragma omp parallel for
    for (int i = 0; i < n + 1; i++) {
        Vec v = start;
        if (i != n)
            v[i] += step[i] * initial_simplex_scale;
        else
            assert(i == n);
        y[i] = function(v);

        for (int j = 0; j < n; j++)
            p[j + i * n] = v[j];        //row-major
    }
    func_eval_number += (n + 1);
//
//  The simplex construction is complete.
//
//  Find highest (highest_value at highest_index) and lowest (lowest_value at lowest_index) Y values.  highest_value = Y(highest_index) indicates
//  the vertex of the simplex to be replaced.
//
    double lowest_value = y[0];
    int lowest_index = 0;

    for (int i = 1; i < nn; i++) {
      if (y[i] < lowest_value) {
        lowest_value = y[i];
        lowest_index = i;
      }
    }
//
//  Inner loop.
//
    int jcount = check_period; //jcount tracks when it is time for a check
    while (1) {
      if (max_func_eval_number <= func_eval_number)
        break;


//first, find highest_value
      double highest_value = y[0];
      int highest_index = 0;

      for (int i = 1; i < nn; i++) {
        if (highest_value < y[i]) {
          highest_value = y[i];
          highest_index = i;
        }
      }
//
//  Calculate PBAR, the centroid of the simplex vertices
//  excepting the vertex with Y value highest_value.
//
      Vec pbar(n);
      for (int i = 0; i < n; i++) {
        double z = 0;
        for (int j = 0; j < nn; j++) {
          z = z + p[i + j * n];
        }
        z = z - p[i + highest_index * n];
        pbar[i] = z / dn;
      }
//
//  Reflection through the centroid.
//
      Vec reflection_point(n);
      for (int i = 0; i < n; i++) {
        reflection_point[i] = pbar[i] + rcoeff * (pbar[i] - p[i + highest_index * n]);
      }
      const double reflection_point_value = function(reflection_point);
      func_eval_number += 1;
//
//  Successful reflection, try also extension.
//
      if (reflection_point_value < lowest_value) {
        Vec expansion_point(n);
        for (int i = 0; i < n; i++) {
          expansion_point[i] = pbar[i] + ecoeff * (reflection_point[i] - pbar[i]);
        }
        const double expansion_point_value = function (expansion_point);
        func_eval_number += 1;

        if (reflection_point_value < expansion_point_value) {
          //reflection is better, save reflection!
          for (int i = 0; i < n; i++) {
            p[i + highest_index * n] = reflection_point[i];
          }
          y[highest_index] = reflection_point_value;
        } else {
          //expansion is better, save expansion!
          for (int i = 0; i < n; i++) {
            p[i + highest_index * n] = expansion_point[i];
          }
          y[highest_index] = expansion_point_value;
        }
      }
//
//  No extension: reflection only, outside contraction, inside contraction, or shrink
//
      else {

        bool to_shrink = 0;
        int l = 0;//number of simplex vertices which values are smaller than reflection_point_value
        for (int i = 0; i < nn; i++) {
          if (reflection_point_value < y[i]) {
            l += 1;
          }
        }

        if (l > 1) { // means that f_r < f_n
          for (int i = 0; i < n; i++) {
            p[i + highest_index * n] = reflection_point[i];
          }
          y[highest_index] = reflection_point_value;
        }
//
//  Contraction on the Y(highest_index) side of the centroid.
//
        else if ( l == 0 ) // means that f_r >= f_n+1
        {
          Vec inside_contraction(n);
          for (int i = 0; i < n; i++) {
            inside_contraction[i] = pbar[i] + ccoeff * ( p[i + highest_index * n] - pbar[i] );
          }
          const double inside_contraction_value = function (inside_contraction);
          func_eval_number += 1;

          if (y[highest_index] < inside_contraction_value) {
          //
          //  Contract the whole simplex.
          //
            to_shrink = 1;
          }
//
//  Retain contraction.
//
          else {
            for (int i = 0; i < n; i++) {
              p[i + highest_index * n] = inside_contraction[i];
            }
            y[highest_index] = inside_contraction_value;
          }
        }
//
//  Contraction on the reflection side of the centroid.
//
        else if ( l == 1 ) // means that f_n <= f_r < f_n+1
        {
          Vec outside_contraction(n);
          for (int i = 0; i < n; i++) {
            outside_contraction[i] = pbar[i] + ccoeff * ( reflection_point[i] - pbar[i] );
          }
          const double outside_contraction_value = function ( outside_contraction );
          func_eval_number += 1;
//
//  Retain reflection?
//
          if (outside_contraction_value <= reflection_point_value) {
            for (int i = 0; i < n; i++) {
              p[i + highest_index * n] = outside_contraction[i];
            }
            y[highest_index] = outside_contraction_value;
          } else { //????????? maybe here we need to shrink?????????
            to_shrink = 1;
            /* original code just saves reflection point, that's all...
            for (int i = 0; i < n; i++) {
              p[i + highest_index * n] = reflection_point[i];
            }
            y[highest_index] = reflection_point_value;
            */
          }
        } else {
            assert(0);
        }

        if (to_shrink) {
            std::cout << "SHRINK!!!!" << std::endl;
            #pragma omp parallel for
            for (int j = 0; j < nn; j++) {
                Vec vec_to_eval(n);
                for (int i = 0; i < n; i++) {
                    p[i + j * n] = (1 - shrink_coeff) * p[i + lowest_index * n] +
                                shrink_coeff * (p[i + j * n]);
                    //was p[i + j * n] = (p[i + j * n] + p[i + lowest_index * n]) * shrink_coeff; wtf???

                    vec_to_eval[i] = p[i + j * n];
                }
                y[j] = function(vec_to_eval);
            }
            func_eval_number += nn;
            //KEEP TRACK OF LOW HERE!!!!!!!!!!
            lowest_value = y[0];
            lowest_index = 0;
            for (int i = 1; i < nn; i++) {
              if (y[i] < lowest_value) {
                lowest_value = y[i];
                lowest_index = i;
              }
            }
            continue;//we will find highest value !!!!
        }
      }
//
//  Check if lowest_value improved.
//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  shrink operation is not here! it goes to the start of the loop
      if (y[highest_index] < lowest_value) {
        lowest_value = y[highest_index];
        lowest_index = highest_index;
      }
      jcount -= 1;

      if ( 0 < jcount )
        continue;

//
//  Check to see if minimum reached.
//
      if (func_eval_number <= max_func_eval_number) {
        std::cout << "e: " << lowest_value << " (" << 100.0 * func_eval_number / max_func_eval_number << " % done)" <<  std::endl;
        error_per_gen.push_back({error_per_gen.size(), lowest_value});
        jcount = check_period;

        double z = 0;
        for (int i = 0; i < nn; i++) {
          z += y[i];
        }
        const double mean = z / dnn;

        z = 0.0;
        for (int i = 0; i < nn; i++) {
          z = z + std::pow ( y[i] - mean, 2 );
        }

        const double rq = reqmin * dn;
        if ( z <= rq )
          break;
      }
    }
    //end of main inner loop

//
//  Factorial tests to check that we check_value is a local minimum.
//
    for (int i = 0; i < n; i++) {
      xmin[i] = p[i + lowest_index * n];
    }
    const double check_value = y[lowest_index];

    if (max_func_eval_number < func_eval_number) {
      ifault = 2;
      break;
    }

    ifault = 0;

#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        Vec v = xmin;
        Vec w = xmin;
        v[i] += step[i] * eps;
        w[i] -= step[i] * eps;

///@todo OpenMP
        const double a = function(v);
        const double b = function(w);
        //why just not to use v or w as a new minimum then
        if (a < check_value && a < b) {
            ifault = 2;
        }
        if (b < check_value && b < a) {
            ifault = 2;
        }
        if (a < check_value || b < check_value) {
            ifault = 2;
        }
    }
    func_eval_number += 2 * n;


    //TODO OpenMP
    for (int i = 0; i < n; i++) {
      const double del = step[i] * eps;
      Vec v = xmin;
      v[i] += del;
      double z = function(v);
      func_eval_number += 1;
      if (z < check_value) {
        ifault = 2;
        break;
      }

      Vec w = xmin;
      w[i] -= del;
      z = function(w);
      func_eval_number += 1;
      if (z < check_value) {
        ifault = 2;
        break;
      }
    }

    if (ifault == 0)
      break;
//
//  Restart the procedure.
//
    for (int i = 0; i < n; i++)
      start[i] = xmin[i];

    initial_simplex_scale = eps_Scale; //why?
    restart_number += 1;
    std::cout << "****************" <<  "RESTART" << "****************" << std::endl;
  }
}


template <typename OptimizationProblem>
std::vector<std::pair<int, double>> nelder_mead(OptimizationProblem & problem, int max_evals, double stop_crit, int check_period, double simplex_scale)
{
    int param_num = problem.get_number_parameters();
    std::vector<double> init_vector(param_num);
    int init_status = problem.initial_guess_for_optimizer(init_vector.begin());
    if (init_status == -1)
        throw("NM requires initial guess from problem");
    return nelder_mead(problem, max_evals, stop_crit, check_period, simplex_scale, init_vector);
}

template <typename OptimizationProblem>
std::vector<std::pair<int, double>> nelder_mead(OptimizationProblem & problem, int max_evals, double stop_crit, int check_period, double simplex_scale, std::vector<double> init_vector)
{
    int param_num = problem.get_number_parameters();

    std::vector<double> min_v(param_num), max_v(param_num);
    std::vector<int> is_mutation_applicable(param_num);
    int boundaries_status = problem.get_boundaries(min_v, max_v, is_mutation_applicable);

    if (boundaries_status != 0)
        throw ("nelder_mead requires boundaries");


    std::vector<double> res(param_num);

    std::vector<std::pair<int, double>> error_per_gen;

    int func_eval_number, restart_number, ifault;


    //TODO still dont know the best way to handle mutable and drifting unknowns
    std::vector<double> step(param_num);
    for (int i = 0; i < param_num; i++)
        step[i] = simplex_scale * (max_v[i] - min_v[i]);

   // nelmin( [& problem](std::vector<double> & s) { return problem.genetic_algorithm_calls(s.begin()); },
    //        init_vector, res, error_per_gen, stop_crit, step, check_period, max_evals, func_eval_number, restart_number, ifault, is_mutation_applicable);

    nelmin( [& problem](std::vector<double> & s) {
            return problem.get_objective_value(s.begin());
            },
        init_vector, res, error_per_gen, stop_crit, step, check_period, max_evals, func_eval_number, restart_number, ifault, is_mutation_applicable);


    if (ifault != 0)
        std::cout << "nelmin error code: " << ifault << std::endl;
    problem.submit_result(res);
    //problem.dump_ap(res.begin(), 10);
    return error_per_gen;
}

#endif
