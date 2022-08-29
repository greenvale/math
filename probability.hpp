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
        static int randomDiscreteEvent(const std::vector<double>& eventProbabilities);
    };
    
    /* ============================================================== */
    
    float Probability::randomRealNumber()
    {
        return (double) rand() / RAND_MAX;
    }
    
    unsigned int Probability::randomDiscreteEvent(const std::vector<double>& eventProbabilities)
    {
        // create vector of event boundaries and cumulative variable
        std::vector<double> eventBoundaries(eventProbabilities.size() + 1, 0.0);
        double cumulat = 0.0;

        // calculate event region boundaries
        for (unsigned int i = 0; i < eventProbabilities.size(); ++i)
        {
            cumulat += prob[i];
            eventBoundaries[i + 1] = cumulat;
        }

        // generate random number
        double randomNum = randomNumber();

        // identify the region index that the random number lies in
        for (unsigned int i = 0; i < eventProbabilities.size(); ++i)
        {
            if ((randomNum >= eventBoundaries[i]) && (randomNum <= eventBoundaries[i + 1]))
            {
                return i;
            }
        }
    }
}
}
