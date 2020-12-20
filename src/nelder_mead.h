#ifndef NELMIN
#define NELMIN
#include <cmath>
#include <string>

template <typename F, typename Vec>
void nelmin (const F & fn, Vec start, Vec & xmin, std::vector<std::pair<int, double>> & error_per_gen,
  double reqmin, const Vec & step, int konvge, int kcount, 
  int & icount, int & numres, int & ifault)
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
//      function fn ( x, f )
//      double fn
//      double x(*)
//
//    and the name of this subroutine must be declared EXTERNAL in the
//    calling routine and passed as the argument FN.
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
//    Input, double FN ( double x[] ), the name of the routine which evaluates
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
//    Output, double YNEWLO, the minimum value of the function.
//
//    Input, double REQMIN, the terminating limit for the variance
//    of function values.
//
//    Input, double STEP[N], determines the size and shape of the
//    initial simplex.  The relative magnitudes of its elements should reflect
//    the units of the variables.
//
//    Input, int KONVGE, the convergence check is carried out 
//    every KONVGE iterations.
//
//    Input, int KCOUNT, the maximum number of function 
//    evaluations.
//
//    Output, int *ICOUNT, the number of function evaluations 
//    used.
//
//    Output, int *NUMRES, the number of restarts.
//
//    Output, int *IFAULT, error indicator.
//    0, no errors detected.
//    1, REQMIN, N, or KONVGE has an illegal value.
//    2, iteration terminated because KCOUNT was exceeded without convergence.
//
  double ynewlo;
  const int n = start.size();
  const double ccoeff = 0.5;
  const double dn = n;
  const int nn = n + 1;
  const double dnn = nn;
  const double ecoeff = 2.0;
  const double eps = 0.001;
  const double rcoeff = 1.0;
  const double rq = reqmin * dn;

  int jcount = konvge; 

  double del_simplex = 1;

//
//  Check the input parameters.
//
  if ( reqmin <= 0.0 ) {
    ifault = 1;
    return;
  }

  if ( n < 1 ) {
    ifault = 1;
    return;
  }

  if ( konvge < 1 ) {
    ifault = 1;
    return;
  }

  // p stores (n + 1) vertices in row-major order
  Vec p(n * (n + 1));
  Vec pstar(n);
  Vec p2star(n);
  Vec pbar(n);
  // function value at corresponding vertices
  Vec y(n + 1);

  icount = 0;
  numres = 0;


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
            v[i] += step[i] * del_simplex;
        else
            assert(i == n);
        y[i] = fn(v);

        for (int j = 0; j < n; j++)
            p[j + i * n] = v[j];
    }
    icount += (n + 1);
//                    
//  The simplex construction is complete.
//                    
//  Find highest (ynewlo at ihi) and lowest (ylo at ilo) Y values.  YNEWLO = Y(IHI) indicates
//  the vertex of the simplex to be replaced.
//                
    double ylo = y[0];
    int ilo = 0;

    for (int i = 1; i < nn; i++) {
      if (y[i] < ylo) {
        ylo = y[i];
        ilo = i;
      }
    }
//
//  Inner loop.
//
    while (1) {
      if ( kcount <= icount )
        break;
      
      ynewlo = y[0];
      int ihi = 0;

      for (int i = 1; i < nn; i++) {
        if (ynewlo < y[i]) {
          ynewlo = y[i];
          ihi = i;
        }
      }
//
//  Calculate PBAR, the centroid of the simplex vertices
//  excepting the vertex with Y value YNEWLO.
//
      for (int i = 0; i < n; i++) {
        double z = 0.0;
        for (int j = 0; j < nn; j++) { 
          z = z + p[i + j * n];
        }
        z = z - p[i + ihi * n];  
        pbar[i] = z / dn;
      }
//
//  Reflection through the centroid.
//
      for (int i = 0; i < n; i++) {
        pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i + ihi * n] );
      }
      const double ystar = fn ( pstar );
      icount += 1;
//
//  Successful reflection, so extension.
//
      if ( ystar < ylo ) {
        for (int i = 0; i < n; i++) {
          p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
        }
        const double y2star = fn ( p2star );
        icount += 1;
//
//  Check extension.
//
        if ( ystar < y2star ) {
          for (int i = 0; i < n; i++) {
            p[i + ihi * n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Retain extension or contraction.
//
        else
        {
          for (int i = 0; i < n; i++) {
            p[i + ihi * n] = p2star[i];
          }
          y[ihi] = y2star;
        }
      }
//
//  No extension.
//
      else
      {
        int l = 0;
        for (int i = 0; i < nn; i++) {
          if (ystar < y[i]) {
            l += 1;
          }
        }

        if ( 1 < l ) {
          for (int i = 0; i < n; i++) {
            p[i + ihi * n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Contraction on the Y(IHI) side of the centroid.
//
        else if ( l == 0 )
        {
          for (int i = 0; i < n; i++) {
            p2star[i] = pbar[i] + ccoeff * ( p[i + ihi * n] - pbar[i] );
          }
          const double y2star = fn ( p2star );
          icount += 1;
//
//  Contract the whole simplex.
//
          if ( y[ihi] < y2star )
          {
            for (int j = 0; j < nn; j++) {
              for (int i = 0; i < n; i++) {
                p[i + j * n] = ( p[i + j * n] + p[i + ilo * n] ) * 0.5;
                xmin[i] = p[i + j * n];
              }
              y[j] = fn ( xmin );
              icount += 1;
            }
            ylo = y[0];
            ilo = 0;

            for (int i = 1; i < nn; i++) {
              if ( y[i] < ylo ) {
                ylo = y[i];
                ilo = i;
              }
            }
            continue;
          }
//
//  Retain contraction.
//
          else
          {
            for (int i = 0; i < n; i++) {
              p[i + ihi * n] = p2star[i];
            }
            y[ihi] = y2star;
          }
        }
//
//  Contraction on the reflection side of the centroid.
//
        else if ( l == 1 )
        {
          for (int i = 0; i < n; i++) {
            p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
          }
          const double y2star = fn ( p2star );
          icount += 1;
//
//  Retain reflection?
//
          if ( y2star <= ystar ) {
            for (int i = 0; i < n; i++) {
              p[i + ihi * n] = p2star[i];
            }
            y[ihi] = y2star;
          } else {
            for (int i = 0; i < n; i++) {
              p[i + ihi * n] = pstar[i];
            }
            y[ihi] = ystar;
          }
        }
      }
//
//  Check if YLO improved.
//
      if ( y[ihi] < ylo )
      {
        ylo = y[ihi];
        ilo = ihi;
      }
      jcount -= 1;

      if ( 0 < jcount )
        continue;

//
//  Check to see if minimum reached.
//
      if ( icount <= kcount ) {
        error_per_gen.push_back({error_per_gen.size(), ylo});
        jcount = konvge;

        double z = 0;
        for (int i = 0; i < nn; i++) {
          z += y[i];
        }
        const double mean = z / dnn;

        z = 0.0;
        for (int i = 0; i < nn; i++) {
          z = z + std::pow ( y[i] - mean, 2 );
        }

        if ( z <= rq )
          break;
      }
    }
//
//  Factorial tests to check that YNEWLO is a local minimum.
//
    for (int i = 0; i < n; i++) {
      xmin[i] = p[i + ilo * n];
    }
    ynewlo = y[ilo];

    if (kcount < icount) {
      ifault = 2;
      break;
    }

    ifault = 0;
//TODO OpenMP
#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        Vec v = xmin;
        Vec w = xmin;
        v[i] += step[i] * eps;
        w[i] -= step[i] * eps;

        const double a = fn(v);
        const double b = fn(w);
        //why just not to use v or w as a new minimum then
        if (a < ynewlo && a < b) {
            ifault = 2;
        }
        if (b < ynewlo && b < a) {
            ifault = 2;
        }
        if (a < ynewlo || b < ynewlo) {
            ifault = 2;
        }
    }
    icount += 2 * n;

    
    for (int i = 0; i < n; i++) {
      const double del = step[i] * eps;
      xmin[i] = xmin[i] + del;
      double z;
      z = fn ( xmin );
      icount += 1;
      if ( z < ynewlo ) {
        ifault = 2;
        break;
      }
      xmin[i] = xmin[i] - del - del;
      z = fn( xmin );
      icount += 1;
      if ( z < ynewlo ) {
        ifault = 2;
        break;
      }
      xmin[i] = xmin[i] + del;
    }

    if (ifault == 0)
      break;
//
//  Restart the procedure.
//
    for (int i = 0; i < n; i++)
      start[i] = xmin[i];

    del_simplex = eps;
    numres += 1;
  }
}


template <typename InitializedRandomGenerator, typename OptimizationProblem>
std::vector<std::pair<int, double>> nelder_mead(InitializedRandomGenerator rg, OptimizationProblem & problem, int max_evals, double stop_crit, int check_period, double simplex_scale)
{
    int param_num = problem.get_number_parameters();
    std::vector<double> init_vector(param_num);
    problem.initial_guess(init_vector.begin());
    std::vector<double> min_v(param_num), max_v(param_num);
    std::vector<int> is_mutation_applicable(param_num);
    int boundaries_status = problem.get_boundaries(min_v, max_v, is_mutation_applicable);
    
    if (boundaries_status != 0)
        throw ("nelmin requires boundaries");
    
    
    for (int j = 0; j < param_num; j++) {
        if (is_mutation_applicable[j])
            init_vector[j] = std::uniform_real_distribution<double>(min_v[j], max_v[j])(rg);
    }
    std::vector<double> res(param_num);

    std::vector<std::pair<int, double>> error_per_gen;

    int icount, numres, ifault;
    
    
    //TODO still dont know the best way to handle mutable and drifting unknowns
    std::vector<double> step(param_num);
    for (int i = 0; i < param_num; i++)
        step[i] = simplex_scale * (max_v[i] - min_v[i]);
    
    nelmin( [& problem](std::vector<double> & s) { return problem.genetic_algorithm_calls(s.begin()); },
            init_vector, res, error_per_gen, stop_crit, step, check_period, max_evals, icount, numres, ifault);

    if (ifault != 0)
        throw (std::string("nelmin error code ") + std::to_string(ifault)); 
    problem.genetic_algorithm_result(res);
    return error_per_gen;
}

#endif
