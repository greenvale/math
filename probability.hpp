#pragma once

#include "math.h"
#include <vector> 

namespace mathlib
{
    class Probability
    {
    private:
        
    public:
        static double randomRealNumber();
        //static double randomRealNumber(const double& x0, const double& x1);
        static int discreteEvent(const std::vector<double>& eventProbabilities);
    };
    
    /* ============================================================== */
    
    double Probability::randomRealNumber()
    {
        return (double) rand() / RAND_MAX;
    }
    
    int Probability::discreteEvent(const std::vector<double>& dist)
    {
        // create vector of event boundaries and cumulative variable
        std::vector<double> distBoundaries(dist.size() + 1, 0.0);
        double cumulat = 0.0;
        int numEvents = dist.size();
        
        // calculate event region boundaries
        for (int i = 0; i < numEvents; ++i)
        {
            cumulat += dist[i];
            distBoundaries[i + 1] = cumulat;
        }

        // generate random number
        double randomNum = randomRealNumber();

        // identify the region index that the random number lies in
        for (int i = 0; i < numEvents; ++i)
        {
            if ((randomNum >= distBoundaries[i]) && (randomNum <= distBoundaries[i + 1]))
            {
                return i;
            }
        }
        return -1;
    }
}
