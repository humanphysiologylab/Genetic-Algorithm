#ifndef OBJECTIVE_H
#define OBJECTIVE_H

#include <vector>
#include <cassert>
#include <cmath>
#include <limits>
#include <iostream>
#include <memory>
#include "utils.h"
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


/**
 * @brief P-norm of error and derivatives
 *
 * Normalized to number of baselines and baselines lenghts p-norm of error
 *
 */
template <typename Baseline, typename VectorOfBaselines, int power = 2>
class MinimizePnormErrorAndDeriv:
    public BaseObjective<Baseline, VectorOfBaselines>
{
public:
    double dist(const VectorOfBaselines & ax, const VectorOfBaselines & bx) const
    {
        auto mini_dist = MinimizePnormError<Baseline, VectorOfBaselines, power>();
        VectorOfBaselines a(ax.size()), b(bx.size());
        VectorOfBaselines da(a.size()), db(b.size());
        const double ema_coef = 0.8;
        for (size_t i = 0; i < a.size(); i++) {
            a[i] = ema(ax[i], ema_coef);
            da[i] = diff(a[i]);
        }
        for (size_t i = 0; i < b.size(); i++) {
            b[i] = ema(bx[i], ema_coef);
            db[i] = diff(b[i]);
        }
        const double dfdx_weight = 10;
        return mini_dist.dist(a, b) + dfdx_weight * mini_dist.dist(da, db);
    }
};

/**
 * @brief scaled P-norm of error
 *
 * Normalized to number of baselines and baselines lenghts scaled p-norm of error
 * vector @p a is used as a reference to scale
 */
template <typename Baseline, typename VectorOfBaselines, int power = 2>
class MinimizeScaledPnormError:
    public BaseObjective<Baseline, VectorOfBaselines>
{
public:
    const double eps = 1e-1;
    double distBaselines(const Baseline & a, const Baseline & b) const
    {
        double res = 0;
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); i++)
            res += std::pow(std::abs((a[i] - b[i]) / std::max(eps, std::abs(a[i]))), power);
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


/**
 * @brief Chebyshev norm of error
 *
 * Normalized to number of baselines Chebyshev norm (aka p_inf) of error
 *
 */
template <typename Baseline, typename VectorOfBaselines>
class MinimizeChebyshevNormError:
    public BaseObjective<Baseline, VectorOfBaselines>
{
public:
    double distBaselines(const Baseline & a, const Baseline & b) const
    {
        assert(a.size() == b.size());
        double max_v = 0;
        for (size_t i = 0; i < a.size(); i++) {
            const double tmp_v = std::abs(a[i] - b[i]);
            if (max_v < tmp_v)
                max_v = tmp_v;
        }
        return max_v;
    }

    double dist(const VectorOfBaselines & a, const VectorOfBaselines & b) const
    {
        assert(a.size() == b.size());
        double res = 0;
        for (size_t i = 0; i < a.size(); i++) {
            res += distBaselines(a[i], b[i]);
        }
        res = res / (double) a.size();
        if (std::isnan(res))
            res = 1e50;
        return res;
    }
};


/**
 * @brief scaled Chebyshev norm of error
 *
 * Normalized to number of baselines scaled Chebyshev norm (aka p_inf) of error
 * vector @p a is used as a reference to scale
 */
template <typename Baseline, typename VectorOfBaselines>
class MinimizeScaledChebyshevNormError:
    public BaseObjective<Baseline, VectorOfBaselines>
{
public:
    const double eps = 1e-1;
    double distBaselines(const Baseline & a, const Baseline & b) const
    {
        assert(a.size() == b.size());
        double max_v = 0;
        for (size_t i = 0; i < a.size(); i++) {
            const double tmp_v = std::abs((a[i] - b[i]) / std::max(eps, std::abs(a[i])));
            if (max_v < tmp_v)
                max_v = tmp_v;
        }
        return max_v;
    }

    double dist(const VectorOfBaselines & a, const VectorOfBaselines & b) const
    {
        assert(a.size() == b.size());
        double res = 0;
        for (size_t i = 0; i < a.size(); i++) {
            res += distBaselines(a[i], b[i]);
        }
        res = res / (double) a.size();
        if (std::isnan(res))
            res = 1e50;
        return res;
    }
};

/**
 * @brief P-log norm of error
 *
 * Normalized to number of baselines and baselines lenghts P-log norm of error
 * vector @p a is used as a reference to scale
 */
template <typename Baseline, typename VectorOfBaselines, int power = 2>
class MinimizePlogNormError:
    public BaseObjective<Baseline, VectorOfBaselines>
{
public:
    const double eps = 1e-1;
    double distBaselines(const Baseline & a, const Baseline & b) const
    {
        assert(a.size() == b.size());
        double res = 0;
        for (size_t i = 0; i < a.size(); i++) {
            res += std::pow(std::log(1 + std::abs((a[i] - b[i]) / std::max(eps, std::abs(a[i])))), power);
        }
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

/**
 * @brief Chebyshev Log norm of error
 *
 * Normalized to number of baselines Chebyshev Log norm of error
 * vector @p a is used as a reference to scale
 */
template <typename Baseline, typename VectorOfBaselines>
class MinimizeChebyshevLogNormError:
    public BaseObjective<Baseline, VectorOfBaselines>
{
public:
    const double eps = 1e-1;
    double distBaselines(const Baseline & a, const Baseline & b) const
    {
        assert(a.size() == b.size());
        double res = 0;
        for (size_t i = 0; i < a.size(); i++) {
            const double tmp = std::log(1 + std::abs((a[i] - b[i]) / std::max(eps, std::abs(a[i]))));
            if (tmp > res)
                res = tmp;
        }
        return res;
    }

    double dist(const VectorOfBaselines & a, const VectorOfBaselines & b) const
    {
        assert(a.size() == b.size());
        double res = 0;
        for (size_t i = 0; i < a.size(); i++) {
            res += distBaselines(a[i], b[i]);
        }
        res /= (double) a.size();
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
            for (size_t i = 0; i < f.size(); i++) {
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



template<typename Baseline, typename VectorOfBaselines>
class MaxCosine:
    public BaseObjective<Baseline, VectorOfBaselines>
{
public:
    double cosine(const Baseline & a, const Baseline & b) const
    {
        double dotp = 0, sna = 0, snb = 0;
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); i++) {
            dotp += a[i] * b[i];
            sna += a[i] * a[i];
            snb += b[i] * b[i];
        }
        return dotp / std::sqrt(sna * snb);
    }

    double dist(const VectorOfBaselines & a, const VectorOfBaselines & b) const
    {
        assert(a.size() == b.size());
        double res = 0;
        for (size_t i = 0; i < a.size(); i++) {
            res += cosine(a[i], b[i]);
        }
        res = res / (double) a.size(); //mean value
        if (std::isnan(res))
            return 1e50;//std::numeric_limits<double>::max();
        return 1 - res;
    }
};

template<typename Baseline, typename VectorOfBaselines>
class MinWeirdDist:
    public BaseObjective<Baseline, VectorOfBaselines>
{
public:
    double weirdDist(const Baseline & a, const Baseline & b) const
    {
        double dotprod = 0, sna = 0, snb = 0;
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); i++) {
            dotprod += a[i] * b[i];
            sna += a[i] * a[i];
            snb += b[i] * b[i];
        }
        return std::abs(std::sqrt(snb / sna) - 1) + std::abs(dotprod / std::sqrt(sna * snb) - 1);
    }

    double dist(const VectorOfBaselines & a, const VectorOfBaselines & b) const
    {
        assert(a.size() == b.size());
        double res = 0;
        for (size_t i = 0; i < a.size(); i++) {
            res += weirdDist(a[i], b[i]);
        }
        res = res / (double) a.size(); //mean value
        if (std::isnan(res))
            return 1e50;//std::numeric_limits<double>::max();
        return res;
    }
};


/**
 * @brief Weighted 2-norm of error
 *
 * Weighted 2-norm of error to fit voltage larger than threshold more precisely.
 * It is normalized to number of baselines and baselines lenghts
 *
 */
template <typename Baseline, typename VectorOfBaselines, int power = 2>
class Weighted2normError:
    public BaseObjective<Baseline, VectorOfBaselines>
{
public:
    const double w = 100;
    const double threshold = -20;
    double distBaselines(const Baseline & a, const Baseline & b) const
    {
        double res = 0;
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); i++) {
            double alpha;
            if (a[i] > threshold || b[i] > threshold)
                alpha = w;
            else
                alpha = 1;
            res += alpha * std::pow(std::abs(a[i] - b[i]), power);
        }
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



template<typename Baseline, typename VectorOfBaselines>
std::unique_ptr<BaseObjective<Baseline, VectorOfBaselines>> new_objective(const std::string & name)
{
    if (name == "Min1normError")
        return std::make_unique<MinimizePnormError<Baseline, VectorOfBaselines, 1>>();
    if (name == "Min2normError")
        return std::make_unique<MinimizePnormError<Baseline, VectorOfBaselines, 2>>();
    if (name == "Minimize2normErrorAndDeriv")
        return std::make_unique<MinimizePnormErrorAndDeriv<Baseline, VectorOfBaselines, 2>>();
    if (name == "LSscaleNshiftMin2normError")
        return std::make_unique<LSscaleNshiftMinPnormError<Baseline, VectorOfBaselines, 2>>();
    if (name == "LSscaleMin2normError")
        return std::make_unique<LSscaleMinPnormError<Baseline, VectorOfBaselines, 2>>();
    if (name == "MaxCosine")
        return std::make_unique<MaxCosine<Baseline, VectorOfBaselines>>();
    if (name == "MinWeirdDist")
        return std::make_unique<MinWeirdDist<Baseline, VectorOfBaselines>>();
    if (name == "Weighted2normError")
        return std::make_unique<Weighted2normError<Baseline, VectorOfBaselines>>();
    if (name == "MinScaled1normError")
        return std::make_unique<MinimizeScaledPnormError<Baseline, VectorOfBaselines, 1>>();
    if (name == "MinScaled2normError")
        return std::make_unique<MinimizeScaledPnormError<Baseline, VectorOfBaselines, 2>>();
    if (name == "MinChebyshevNormError")
        return std::make_unique<MinimizeChebyshevNormError<Baseline, VectorOfBaselines>>();
    if (name == "MinScaledChebyshevNormError")
        return std::make_unique<MinimizeScaledChebyshevNormError<Baseline, VectorOfBaselines>>();
    if (name == "MinPlog2normError")
        return std::make_unique<MinimizePlogNormError<Baseline, VectorOfBaselines, 2>>();
    if (name == "MinPlog1normError")
        return std::make_unique<MinimizePlogNormError<Baseline, VectorOfBaselines, 1>>();
    if (name == "MinChebyshevLogNormError")
        return std::make_unique<MinimizeChebyshevLogNormError<Baseline, VectorOfBaselines>>();

    throw std::logic_error("Unknown objective type");
}

#endif
