#ifndef MCMC_H
#define MCMC_H

#ifdef MCMC_ENABLED

#include <BAT/BCLog.h>
#include <BAT/BCModel.h>
#include <TH1D.h>
#include <TF1.h>

#include <BAT/BCGaussianPrior.h>

#include <BAT/BCMath.h>
#include <TMath.h>

#include <string>
#include <vector>


template<typename Problem>
class MCMC_optimization_problem:
    public BCModel
{
    Problem & problem;
public:
    MCMC_optimization_problem(const std::string & name, Problem & problem, std::vector<double> initial_guess)
        : BCModel(name), problem(problem)
    {
        // Define parameters here in the constructor.
        // Also define their priors, if using built-in priors.
        // For example:
        // AddParameter("mu", -2, 1, "#mu", "[GeV]");
        // GetParameters.Back().SetPrior(new BCGaussianPrior(-1, 0.25));

        // Define observables here, too. For example:
        // AddObservable("mu_squared", 1, 4, "#mu^{2}", "[GeV^{2}]");


        //add parameters
        //addParameters return false if it cannot add parameter
        //this can happen because of the non-unique name
        //also you can add prior for each parameter right here
        //AddParameter("mu", 5.27, 5.29, "#mu", "[GeV]");
        //GetParameters().Back().SetPrior(new BCGaussianPrior(5.28, 2e-3));

        //AddParameter("sigma", 25e-3, 45e-3, "#sigma", "[GeV]");
        //GetParameters().Back().SetPrior(new BCGaussianPrior(35e-3, 3e-3));

        //AddParameter("height", 0, 10, "", "[events]");
        //GetParameters().Back().SetPriorConstant();
        //what is PriorConstant? uniform?

        int param_num = problem.get_number_parameters();
        std::vector<double> min_v(param_num), max_v(param_num);
        std::vector<int> is_mutation_applicable(param_num);
        int boundaries_status = problem.get_boundaries(min_v, max_v, is_mutation_applicable);

        if (boundaries_status != 0)
            throw("MCMC requires boundaries");
        std::vector<std::string> names = problem.get_unique_parameter_names();

        for (int i = 0; i < param_num; i++) {
            AddParameter(names[i], min_v[i], max_v[i]);
            //GetParameters().Back().SetPriorConstant();
            GetParameters().Back().SetPrior(new BCGaussianPrior(initial_guess[i], 1e-1));
            std::cout << names[i] << " " << min_v[i] << " " << max_v[i] << " " << initial_guess[i] << std::endl;
        }

        //we can also explore any functions of out model's parameters
        //we need to set limits as well
        //we don't need prior for observables
        //AddObservable("SignalYield", 900, 1100, "Y_{S}", "[events]");
        //AddObservable("Resolution",
        //    100. * GetParameter("sigma").GetLowerLimit() / GetParameter("mu").GetUpperLimit(),
        //    100. * GetParameter("sigma").GetUpperLimit() / GetParameter("mu").GetLowerLimit(),
        //    "#sigma / #mu", "[%]");
    }

    ~MCMC_optimization_problem()
    {}


    double LogLikelihood(const std::vector<double>& pars)
    {
        // return the log of the conditional probability p(data|pars).
        // This is where you define your model.
        // BCMath contains many functions you will find helpful.

        // this part determines how exactly we treat uncertainty
        // it is defined by the nature of the problem

        return problem.loglikelihood(pars);
    }


// double LogAPrioriProbability(const std::vector<double>& pars)
// {
//     // return the log of the prior probability p(pars)
//     // If you use built-in priors, leave this function commented out.
// }

/* Overload CalculateObservables if using observables
    void CalculateObservables(const std::vector<double>& pars)
    {
        // Calculate and store obvserables. For example:
        //GetObservable(0) = pow(pars[0], 2);

        // store total of number events expected
        double nu = 0;
        // loop over bins of our data
        for (int i = 1; i <= fDataHistogram.GetNbinsX(); ++i)
            // calculate expected number of events in that bin
            // and add to total expectation
            nu += pars[2] * TMath::Gaus(fDataHistogram.GetBinCenter(i), pars[0], pars[1], true);

        // store in the observable
        GetObservable(0) = nu;
        // Store sigma as percentage of mu:
        GetObservable(1) = 100. * pars[1] / pars[0];
    }
*/
};




template<typename Problem>
void mcmc(Problem & problem, const std::string & name, std::vector<double> initial_guess)
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // create new MyMod object
    MCMC_optimization_problem m(name, problem, initial_guess);

    // set precision
    m.SetPrecision(BCEngineMCMC::kQuick);//kQuick);kMedium

    BCLog::OutSummary("Model created");

    //////////////////////////////
    // perform your analysis here

    // Normalize the posterior by integrating it over the full parameter space
    // m.Normalize();

    // Write Markov Chain to a ROOT file as a TTree
    // m.WriteMarkovChain(m.GetSafeName() + "_mcmc.root", "RECREATE");

    // run MCMC, marginalizing posterior
    //m.SetProposalFunctionDof(-1);
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // run mode finding; by default using Minuit
    m.FindMode(m.GetBestFitParameters());

    // draw all marginalized distributions into a PDF file
    m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf", 2, 2);//2x2 grid

    // print summary plots
    m.PrintParameterPlot(m.GetSafeName() + "_parameters.pdf");
    m.PrintCorrelationPlot(m.GetSafeName() + "_correlation.pdf");
    m.PrintCorrelationMatrix(m.GetSafeName() + "_correlationMatrix.pdf");
    m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_update.pdf", 3, 2);

    // print results of the analysis into a text file
    m.PrintSummary();

    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();
}

#endif
#endif
