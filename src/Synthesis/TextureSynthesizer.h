#include <gsl/gsl_multimin.h>
#include "Filtering/Signal.h"
#include "Synthesis/StatsGenerator.h"

namespace TextureSynthesis
{
    class TextureSynthesizer;

    /* Container for all data necessary for each iteration of gradient descent optimization
     */
    struct OptimizationData
    {
        /* Initialize all data necessary for gradient descent optimization
         */
        void initialize(TextureSynthesizer& synthesizer, const Signal& initialSignal,
                int downsampleSize, double downsampleRate);

        /* Envelopes and signals in current iteration
         */
        std::vector<Signal> cochlearEnvelopes;
        std::vector<std::vector<Signal>> modulationSignals;

        /* Statistics in current iteration
         */
        std::vector<std::complex<double>> currentStats;

        /* Target statistics
         */
        std::vector<std::complex<double>> targetStats;

        /* Objects for computing statistics
         */
        TextureFilterer textureFilterer;
        StatsGenerator statsGenerator;
    };

    /* Synthesizes textures similar to a target signal
     */
    class TextureSynthesizer
    {
    public:
        /* Initialize synthesizer with target signal
         */
        TextureSynthesizer(const Signal& targetSignal);
        
        /* Fill outSignal with texture similar to targetSignal
         */
        void synthesize(Signal& outSignal);

        friend class OptimizationData;

    private:
        /* Target signal for synthesizing
         */
        Signal _targetSignal;

        /* Gradient descent data
         */
        double _stepSize, _tolerance;
        OptimizationData _curOptimizationData;

        /* Gradient descent functions
         */
        static double distanceFromTarget(const gsl_vector *v, void *params);
        static double distanceFromTarget(OptimizationData *data);
        static void gradient(const gsl_vector *v, void *params, gsl_vector *df);
        static void gradAndDist(const gsl_vector *v, void *params, double *f, gsl_vector *df);
    };
}
