#ifndef LOCATIONBOB_H
#define LOCATIONBOB_H

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <iostream>
#include <algorithm>
#include "problem.h"
#include "Cluster.h"
#include "bobyqa.h"

// 计算候选缓冲区位置的成本
double calculateCost(const ProblemData::Coordinate &candidate_location,
                     const std::vector<int> &layer1_indices,
                     const std::vector<Cluster> &Clusters_Layer1)
{
    double max_diff = 0;
    std::vector<double> maxDiss(layer1_indices.size(), 0);
    std::vector<double> manhattanDistance_S(layer1_indices.size(), 0);

    for (size_t i = 0; i < layer1_indices.size(); ++i)
    {
        int layer1_idx = layer1_indices[i];
        const Cluster &cluster1 = Clusters_Layer1[layer1_idx];

        double max_Distance = 0;
        for (const auto &ff : cluster1.ffs)
        {
            double manhattanDistance = std::pow(ff.distance, 2);
            max_Distance = std::max(manhattanDistance, max_Distance);
        }
        maxDiss[i] = max_Distance;

        manhattanDistance_S[i] = std::pow(std::abs(candidate_location.x - cluster1.GetCX()) + std::abs(candidate_location.y - cluster1.GetCY()), 2);
    }
    // 计算目标成本
    for (size_t i = 0; i < layer1_indices.size() - 1; ++i)
    {
        for (size_t j = i + 1; j < layer1_indices.size(); ++j)
        {

            double diff = std::abs((maxDiss[i] + manhattanDistance_S[i]) - (maxDiss[j] + manhattanDistance_S[j]));
            max_diff = std::max(max_diff, diff);
        }
    }
    return max_diff;
}

struct BOBYQA_DATA
{
    const std::vector<int> &var1;
    const std::vector<Cluster> &var2;
};

REAL bobyqa_obj(const INTEGER n, const REAL *x, void *data)
{
    BOBYQA_DATA *params = (BOBYQA_DATA *)data;
    const std::vector<int> &var1 = params->var1;
    const std::vector<Cluster> &var2 = params->var2;
    return calculateCost(ProblemData::Coordinate{x[0], x[1]}, var1, var2);
}

// 模拟退火算法优化缓冲区位置
ProblemData::Coordinate optimizeBufferLocation_SA(
    const Cluster &cluster2,
    const std::vector<Cluster> &Clusters_Layer1)
{
    const std::vector<int> &layer1_indices = cluster2.GetChildIndex();
    static double cost_sum = 0;
    ProblemData::Coordinate best_location = {cluster2.GetCX(), cluster2.GetCY()};

    double min_cost = calculateCost(best_location, layer1_indices, Clusters_Layer1);

    // 计算 cluster2 到其对应的第一层 buffer 的最远曼哈顿距离
    double max_manhattan_distance_x = 0.0;
    double max_manhattan_distance_y = 0;
    for (int layer1_idx : layer1_indices)
    {
        const Cluster &cluster1 = Clusters_Layer1[layer1_idx];
        double manhattan_distance_x = std::abs(cluster1.GetCX() - cluster2.GetCX()) /*+ std::abs(cluster1.GetCY() - cluster2.GetCY() )*/;
        double manhattan_distance_y = std::abs(cluster1.GetCY() - cluster2.GetCY()) /* + std::abs(cluster1.GetCX() - cluster2.GetCX() )*/;
        max_manhattan_distance_x = std::max(max_manhattan_distance_x, manhattan_distance_x);
        max_manhattan_distance_y = std::max(max_manhattan_distance_y, manhattan_distance_y);
    }
    double search_radius_x = max_manhattan_distance_x; // 更新 search_radius
    double search_radius_y = max_manhattan_distance_y;
    //std::cout << "Updated search_radius: " << search_radius_x << "  " << search_radius_y << std::endl;
    double *workspace = new double[10000];
    REAL x[2] = {cluster2.GetCX(), cluster2.GetCY()};
    // std::cout << x[0] << " " << x[1] << "\n";
    REAL xu[2] = {cluster2.GetCX() + search_radius_x, cluster2.GetCY() + search_radius_y};
    REAL xl[2] = {cluster2.GetCX() - search_radius_x, cluster2.GetCY() - search_radius_y}; // 使用更新的search_radius
    BOBYQA_DATA data = {layer1_indices, Clusters_Layer1};
    bobyqa(2, 5,
           bobyqa_obj, (void *)&data,
           x, xl, xu,
           std::min(search_radius_x, search_radius_y) * 0.5, 0.5,
           0, 1000000, workspace);
    // std::cout << x[0] << " " << x[1] << "\n";
    best_location = ProblemData::Coordinate{x[0], x[1]};
    delete workspace;

    //cost_sum += calculateCost(best_location, layer1_indices, Clusters_Layer1);
    //std::cout << "cost_sum: " << cost_sum << "\n";
    return best_location;
}

#ifndef _BOBYQA_H
#define _BOBYQA_H 1

#ifndef LOGICAL
#define LOGICAL int
#endif

#ifndef INTEGER
#define INTEGER long
#endif

#ifndef REAL
#define REAL double
#endif

/* Prototype of the objective function assumed by the BOBYQA routine.  The
   returned value is the function value at X the current variables, N is the
   number of variables and DATA is anything needed by the function (unused by
   BOBYQA itself). */

typedef REAL bobyqa_objfun(const INTEGER n, const REAL *x, void *data);

/* BOBYQA seeks the least value of a function of many variables, by applying a
   trust region method that forms quadratic models by interpolation.  There is
   usually some freedom in the interpolation conditions, which is taken up by
   minimizing the Frobenius norm of the change to the second derivative of the
   model, beginning with the zero matrix.  The values of the variables are
   constrained by upper and lower bounds.  The arguments of the subroutine are
   as follows.

   N must be set to the number of variables and must be at least two.  NPT is
   the number of interpolation conditions.  Its value must be in the interval
   [N+2,(N+1)(N+2)/2].  Choices that exceed 2*N+1 are not recommended.

   OBJFUN is provided by the user to compute the objective function value at
   the values of the variables X(1),X(2),...,X(N), which are generated
   automatically by BOBYQA in a way that satisfies the bounds given in XL and
   XU.  DATA is anything needed by the function and which is passed as is to
   OBJFUN by BOBYQA.

   Initial values of the variables must be set in X(1),X(2),...,X(N).  They
   will be changed to the values that give the least calculated F.  For
   I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper bounds,
   respectively, on X(I).  The construction of quadratic models requires XL(I)
   to be strictly less than XU(I) for each I.  Further, the contribution to a
   model from changes to the I-th variable is damaged severely by rounding
   errors if XU(I)-XL(I) is too small.

   RHOBEG and RHOEND must be set to the initial and final values of a trust
   region radius, so both must be positive with RHOEND no greater than RHOBEG.
   Typically, RHOBEG should be about one tenth of the greatest expected change
   to a variable, while RHOEND should indicate the accuracy that is required in
   the final values of the variables.  An error return occurs if any of the
   differences XU(I)-XL(I), I=1,...,N, is less than 2*RHOBEG.

   The value of IPRINT should be set to 0, 1, 2 or 3, which controls the amount
   of printing.  Specifically, there is no output if IPRINT=0 and there is
   output only at the return if IPRINT=1.  Otherwise, each new value of RHO is
   printed, with the best vector of variables so far and the corresponding
   value of the objective function.  Further, each new value of F with its
   variables are output if IPRINT=3.

   MAXFUN must be set to an upper bound on the number of calls of OBJFUN.

   The array W will be used for working space.  Its length must be at least
   (NPT+5)*(NPT+N)+3*N*(N+5)/2.  Upon successful return, the first element of W
   will be set to the function value at the solution. */

extern int bobyqa(const INTEGER n, const INTEGER npt,
                  bobyqa_objfun *objfun, void *data,
                  REAL *x, const REAL *xl, const REAL *xu,
                  const REAL rhobeg, const REAL rhoend,
                  const INTEGER iprint, const INTEGER maxfun, REAL *w);

// Possible values returned by BOBYQA

#define BOBYQA_SUCCESS (0)               // algorithm converged
#define BOBYQA_BAD_NPT (-1)              // NPT is not in the required interval
#define BOBYQA_TOO_CLOSE (-2)            // insufficient space between the bounds
#define BOBYQA_ROUNDING_ERRORS (-3)      // too much cancellation in a denominator
#define BOBYQA_TOO_MANY_EVALUATIONS (-4) // maximum number of function evaluations exceeded
#define BOBYQA_STEP_FAILED (-5)          // a trust region step has failed to reduce Q

#endif // _BOBYQA_H

#endif
