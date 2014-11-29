#include "Synthesis/TextureSynthesizer.h"

#include "Synthesis/TextureFilterer.h"
#include "Filtering/Downsampler.h"
#include <complex>

using namespace TextureSynthesis;

void OptimizationData::initialize(TextureSynthesizer& synthesizer,
        const Signal& initialSignal, int downsampleSize, double downsampleRate)
{
    // Initialize envelopes and modulation signals
    for(int i = 0; i < TextureFilterer::numCochlearEnvelopes(); i++)
    {
        cochlearEnvelopes.push_back(Signal(downsampleSize, downsampleRate));

        std::vector<Signal> someModulationSignals;
        for(int j = 0; j < TextureFilterer::numModulationSignals(); j++)
            someModulationSignals.push_back(Signal(downsampleSize, downsampleRate));
        modulationSignals.push_back(someModulationSignals);
    }

    // Initialize target statistics from target signal
    textureFilterer.auditoryFilter(synthesizer._targetSignal, cochlearEnvelopes,
            modulationSignals, true);
    statsGenerator.computeStatistics(cochlearEnvelopes, modulationSignals,
            targetStats);

    // Initialize temporary statistics
    textureFilterer.auditoryFilter(initialSignal, cochlearEnvelopes,
            modulationSignals, false);
    statsGenerator.computeStatistics(cochlearEnvelopes, modulationSignals,
            currentStats);
}

TextureSynthesizer::TextureSynthesizer(const Signal& targetSignal) :
    _targetSignal(targetSignal),
    _curOptimizationData(),
    _stepSize(1000),
    _tolerance(10000)
{ }

void TextureSynthesizer::synthesize(Signal& outSignal)
{
    // Calculate the size and sample rate of outSignal's envelopes after downsampling
    std::pair<int, double> downsampleSizeAndRate =
        Downsampler::newSizeAndRate(outSignal.size(), outSignal.sampleRate,
                TextureFilterer::targetDownsampleRate());
    int downsampleSize = downsampleSizeAndRate.first;
    double downsampleRate = downsampleSizeAndRate.second;

    // Resize and randomly initialize outSignal
    outSignal.resize(pow(2, (int)(log(outSignal.size()) / log(2))));
    srand(time(nullptr));
    for(int i = 0; i < outSignal.size(); i++)
        outSignal[i] = 2 * (rand() / (double)RAND_MAX) - 1;

    // Initialize the optimization data for this synthesize call
    _curOptimizationData.initialize(*this, outSignal, downsampleSize, downsampleRate);

    // Initialize the parameters of our gsl minimizer function
    gsl_multimin_function_fdf statmin;
    statmin.n = downsampleSize * TextureFilterer::numCochlearEnvelopes();
    statmin.f = &distanceFromTarget;
    statmin.df = &gradient;
    statmin.fdf = &gradAndDist;
    statmin.params = &_curOptimizationData;

    // Set the initial guess for gradient descent
    gsl_vector *init = gsl_vector_alloc(statmin.n);
    for(int i = 0; i < statmin.n; i++)
        gsl_vector_set(init, i, std::real(_curOptimizationData.cochlearEnvelopes
                    [i / downsampleSize][i % downsampleSize]));

    // Initialize the gsl minimizer
    gsl_multimin_fdfminimizer *s;
    s = gsl_multimin_fdfminimizer_alloc (gsl_multimin_fdfminimizer_conjugate_fr, statmin.n);
    gsl_multimin_fdfminimizer_set(s, &statmin, init, _stepSize, _tolerance);

    // Run minimization loop
    int iter = 0;
    int status;
    do
    {
        iter++;
        std::cout << "iter: " << iter << "\n";
        status = gsl_multimin_fdfminimizer_iterate(s);
    }
    while (status == GSL_CONTINUE && iter < 5);
    std::cout << "End status: " << status << "\n";

    // Convert results into output signal
    _curOptimizationData.textureFilterer.recombine(_curOptimizationData.cochlearEnvelopes,
            outSignal);

    // Cleanup
    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(init);
}

double TextureSynthesizer::distanceFromTarget(const gsl_vector *v, void *params)
{
    // Grab cochlear envelopes from optimization data
    OptimizationData *data = (OptimizationData*)params;
    std::vector<Signal> *cochlearEnvelopes = &(data->cochlearEnvelopes);

    // Copy modified minimization arguments into optimization data
    int numEnvelopes = cochlearEnvelopes->size();
    int envelopeSize = cochlearEnvelopes->at(0).size();
    for(int i = 0; i < numEnvelopes; i++)
        for(int j = 0; j < envelopeSize; j++)
            (*cochlearEnvelopes)[i][j] = gsl_vector_get(v, envelopeSize * i + j);

    // Shell out to other overload
    return distanceFromTarget(data);
}

double TextureSynthesizer::distanceFromTarget(OptimizationData *data)
{
    // Compute new modulation signals and statistics from updated envelopes
    data->textureFilterer.modulationFilter(data->cochlearEnvelopes, data->modulationSignals);
    data->statsGenerator.computeStatistics(data->cochlearEnvelopes, data->modulationSignals,
            data->currentStats);

    // Sum up L2 distance between current stats and target stats
    double distance = 0;
    std::complex<double> difference;
    int numStats = data->targetStats.size();
    for(int i = 0; i < numStats; i++)
    {
        difference = data->currentStats.at(i) - data->targetStats.at(i);
        distance += pow(std::real(difference), 2) + pow(std::imag(difference), 2);
    }

    return distance;
}

double TextureSynthesizer::partialDerivative(OptimizationData *data, int curEnvelope,
        int curSample)
{
    // TODO: Compute partial derivatives!
    data->statsGenerator.computePartials(data->cochlearEnvelopes, data->modulationSignals,
            curEnvelope, curSample, data->partialDerivatives);

    int numPartials = data->partialDerivatives.size();
    double grad = 0;
    std::complex<double> term;
    for(int i = 0; i < numPartials; i++)
    {
        term = data->currentStats[i] * data->partialDerivatives[i];
        grad += 2 * (std::real(term) + std::imag(term));
    }

    return grad;
}

static int numGradients = 0;    // for output while debugging

void TextureSynthesizer::gradient(const gsl_vector *v, void *params, gsl_vector *df)
{
    numGradients++;

    // Get optimizationdata from void *params
    OptimizationData *data = (OptimizationData*)params;
    int numEnvelopes = data->cochlearEnvelopes.size();
    int envelopeSize = data->cochlearEnvelopes[0].size();

    // Copy modified minimization arguments into optimization data
    for(int i = 0; i < numEnvelopes; i++)
        for(int j = 0; j < envelopeSize; j++)
            (data->cochlearEnvelopes)[i][j] = gsl_vector_get(v, i * envelopeSize + j);

    // Calculate initial objective function value
    double initialDistance = distanceFromTarget(data);

    // Compute partial derivates
    for(int i = 0; i < numEnvelopes; i++)
        for(int j = 0; j < envelopeSize; j++)
            gsl_vector_set(df, i * envelopeSize + j, partialDerivative(data, i, j));
}

void TextureSynthesizer::gradAndDist(const gsl_vector *v, void *params, double *f,
        gsl_vector *df)
{
    *f = distanceFromTarget(v, params);
    gradient(v, params, df);
}
