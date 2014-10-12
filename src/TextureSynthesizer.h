#include <gsl/gsl_multimin.h>
#include "Signal.h"

namespace TextureSynthesis
{
    struct OptimizationData
    {
    };

    class TextureSynthesizer
    {
    public:
        TextureSynthesizer(const Signal& targetSignal);
        void synthesize(Signal& outSignal);

    private:
        Signal _targetSignal;
        std::vector<std::complex<double>> _targetStats;
        double _downsampleRate, _stepSize, _tolerance;
        int _numSubbands;

        static double distanceFromTarget(const gsl_vector *v, void *params);
        static void gradient(const gsl_vector *v, void *params, gsl_vector *df);
        static void gradAndDist(const gsl_vector *v, void *params, double *f, gsl_vector *df);
    };
}
