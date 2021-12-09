#ifndef OBJECTIVE_H
#define OBJECTIVE_H

#include <vector>
#include <cassert>
#include <cmath>
#include <limits>

/**
 * @brief Base objective class
 *
 * Base objective class providing dist(...) interface
 *
 */
template <typename Baseline, typename VectorOfBaselines>
class BaseObjective
{
public:
    virtual double dist(const VectorOfBaselines & a, const VectorOfBaselines & b) const
        = 0;
};

/**
 * @brief P-norm of error
 *
 * Normalized to number of baselines and baselines lenghts p-norm of error
 *
 */
template <typename Baseline, typename VectorOfBaselines, int power = 2>
class MinimizePnormError:
    public BaseObjective<Baseline, VectorOfBaselines>
{
public:
    double distBaselines(const Baseline & a, const Baseline & b) const
    {
        double res = 0;
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); i++)
            res += std::pow(std::abs(a[i] - b[i]), power);
        res /= (double) a.size();
        return res;
    }

    double dist(const VectorOfBaselines & a, const VectorOfBaselines & b) const
    {
        assert(a.size() == b.size());
        double res = 0;
        for (size_t i = 0; i < a.size(); i++) {
            res += distBaselines(a[i], b[i]);
        }
        res = std::pow(res / (double) a.size(), 1.0 / power);
        if (std::isnan(res))
            res = 1e50;
        return res;
    }
};


template<typename Baseline, typename VectorOfBaselines, int power = 2>
class LSscaleNshiftMinPnormError:
    public BaseObjective<Baseline, VectorOfBaselines>
{
public:
    template<int p>
    double sum_of_elements(const Baseline & b) const
    {
        double s = 0;
        for (const auto & e: b)
            s += std::pow(e, p);
        return s;
    }
    double inner_product(const Baseline & a, const Baseline & b) const
    {
        double s = 0;
        for (unsigned i = 0; i < a.size(); i++)
            s += a[i] * b[i];
        return s;
    }

    double LSMdistBaselines(const Baseline & g, const Baseline & f) const
    {
        //watch for the correct order: g -- experimental optical mapping AP, f -- model
        //solve least square problem for f(t) = A * g(t) + B
        assert(g.size() == f.size());

        const double c = g.size();
        const double a = sum_of_elements<2>(g);
        const double b = sum_of_elements<1>(g);
        const double e = sum_of_elements<1>(f);
        const double p = inner_product(f, g);

        const double A = (c * p  - b * e) / (a * c - b * b);
        const double B = (-b * p + a * e) / (a * c - b * b);

        const double min_peak = 10;         // mV
        const double max_peak = 60;         // mV
        const double min_rest_potential = -100;   // mV
        const double max_rest_potential = -55;   // mV

        const double peak = std::min(std::max(A + B, min_peak), max_peak);
        const double rest_potential = std::min(std::max(B, min_rest_potential), max_rest_potential);

        const double amplitude = peak - rest_potential;
        /*
        std::cout << "     A+B: " << A + B << " A: " << amplitude << " B: " << rest_potential << " p: " << peak << std::endl;
        */
        if (std::isnan(amplitude)) {
            std::cout << "         a:" << a <<
            "         b:" << b <<
            "         c:" << c <<
            "         e:" << e <<
            "         p:" << p <<
            "         A:" << A <<
            "         B:" << B << std::endl;
            for (int i = 0; i < f.size(); i++) {
                if (std::isnan(f[i])) {
                    std::cout << "f is nan" << std::endl;
                    break;
                }
            }
        }
        double res = 0;
        for (size_t i = 0; i < f.size(); i++)
            res += std::pow(std::abs(f[i] - (amplitude * g[i] + rest_potential)), power);

       // for (size_t i = 40; i < 80; i++)
         //   res += std::pow(std::abs(f[i] - (amplitude * g[i] + rest_potential)), 4);

        res /= (double) f.size();
        return res;
    }

    double dist(const VectorOfBaselines & a, const VectorOfBaselines & b) const
    {
        //watch for the correct order: a -- experimental optical mapping AP, b -- model
        assert(a.size() == b.size());
        double res = 0;
        for (size_t i = 0; i < a.size(); i++) {
            res += LSMdistBaselines(a[i], b[i]);
        }
        res = std::pow(res / (double) a.size(), 1.0/power);
        if (std::isnan(res))
            res = 1e50;
        return res;
    }
};




template<typename Baseline, typename VectorOfBaselines, int power = 2>
class LSscaleMinPnormError:
    public BaseObjective<Baseline, VectorOfBaselines>
{
public:
    template<int p>
    double sum_of_elements(const Baseline & b) const
    {
        double s = 0;
        for (const auto & e: b)
            s += std::pow(e, p);
        return s;
    }
    double inner_product(const Baseline & a, const Baseline & b) const
    {
        double s = 0;
        for (unsigned i = 0; i < a.size(); i++)
            s += a[i] * b[i];
        return s;
    }

    double LSdistBaselines(const Baseline & g, const Baseline & f) const
    {
        //watch for the correct order: g -- experimental data (not from (0,1) but wider) AP, f -- model
        //solve least square problem for f(t) = A * g(t)
        assert(g.size() == f.size());

        const double a = sum_of_elements<2>(g);
        const double p = inner_product(f, g);

        const double min_A = 1, max_A = 1.2;
        const double A = std::min(max_A, std::max(min_A, p / a));

        double res = 0;
        for (size_t i = 0; i < f.size(); i++)
            res += std::pow(std::abs(f[i] - A * g[i]), power);

        res /= (double) f.size();
        return res;
    }

    double dist(const VectorOfBaselines & a, const VectorOfBaselines & b) const
    {
        //watch for the correct order: a -- experimental optical mapping AP, b -- model
        assert(a.size() == b.size());
        double res = 0;
        for (size_t i = 0; i < a.size(); i++) {
            res += LSdistBaselines(a[i], b[i]);
        }
        res = std::pow(res / (double) a.size(), 1.0/power);
        if (std::isnan(res))
            res = 1e50;
        return res;
    }
};



template<typename Baseline, typename VectorOfBaselines, int power = 2>
class MaxInnerProduct:
    public BaseObjective<Baseline, VectorOfBaselines>
{
public:
    double innerProduct(const Baseline & a, const Baseline & b) const
    {
        double res = 0, lena = 0, lenb = 0;
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); i++) {
            res += a[i] * b[i];
            lena += a[i] * a[i];
            lenb += b[i] * b[i];
        }
        res /= sqrt(lena * lenb);
        return res;
    }

    double dist(const VectorOfBaselines & a, const VectorOfBaselines & b) const
    {
        assert(a.size() == b.size());
        double res = 0;
        for (size_t i = 0; i < a.size(); i++) {
            res += innerProduct(a[i], b[i]);
        }
        res = res / (double) a.size(); //mean value
        if (res < -1 || res > 1)
            std::cout << "ATTENTION: inner product: " << res << std::endl;
        if (std::isnan(res))
            return 1e50;//std::numeric_limits<double>::max();
        return -res;
    }
};

template<typename Baseline, typename VectorOfBaselines>
std::unique_ptr<BaseObjective<Baseline, VectorOfBaselines>> new_objective(const std::string & name)
{
    if (name == "Min2normError")
        return std::make_unique<MinimizePnormError<Baseline, VectorOfBaselines, 2>>();
    if (name == "LSscaleNshiftMin2normError")
        return std::make_unique<LSscaleNshiftMinPnormError<Baseline, VectorOfBaselines, 2>>();
    if (name == "LSscaleMin2normError")
        return std::make_unique<LSscaleMinPnormError<Baseline, VectorOfBaselines, 2>>();
    if (name == "MaxInnerProduct")
        return std::make_unique<MaxInnerProduct<Baseline, VectorOfBaselines, 2>>();
    throw std::logic_error("Unknown objective type");
}

#endif
