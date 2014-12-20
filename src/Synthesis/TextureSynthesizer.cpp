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
            modulationSignals, false);
    statsGenerator.computeStatistics(cochlearEnvelopes, modulationSignals,
            targetStats);

    // Initialize temporary statistics
    textureFilterer.auditoryFilter(initialSignal, cochlearEnvelopes,
            modulationSignals, true);
    statsGenerator.computeStatistics(cochlearEnvelopes, modulationSignals,
            currentStats);
}

TextureSynthesizer::TextureSynthesizer(const Signal& targetSignal) :
    _targetSignal(targetSignal),
    _curOptimizationData(),
    _stepSize(30),
    _tolerance(.1)
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
    statmin.f = &distance;
    statmin.df = &gradient;
    statmin.fdf = &gradAndDist;
    statmin.params = &_curOptimizationData;

    // Set the initial guess for gradient descent
    gsl_vector *init = gsl_vector_alloc(statmin.n);
    for(int i = 0; i < TextureFilterer::numCochlearEnvelopes(); i++)
        for(int j = 0; j < downsampleSize; j++)
            gsl_vector_set(init, downsampleSize * i + j,
                    std::real(_curOptimizationData.cochlearEnvelopes[i][j]));

    // Initialize the gsl minimizer
    gsl_multimin_fdfminimizer *s;
    s = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, statmin.n);
    gsl_multimin_fdfminimizer_set(s, &statmin, init, _stepSize, _tolerance);

    // Run minimization loop
    int iter = 0;
    int status;
    do
    {
        iter++;
        std::cout << "iter: " << iter << "\n";
        status = gsl_multimin_fdfminimizer_iterate(s);
        printf("status: %s\n", gsl_strerror(status));
    }
    while(iter < 50);

    // TODO: figure out where to recombine more
 
    // Convert into output signal
    _curOptimizationData.textureFilterer.recombine(_curOptimizationData.cochlearEnvelopes,
           outSignal); 

    // TODO: compare spectrums of input and output
    // then try running just on that frequency (filter input)

    // Cleanup
    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(init);
}

double TextureSynthesizer::distance(const gsl_vector *v, void *params)
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

    // Compute new modulation signals and statistics from updated envelopes
    data->textureFilterer.modulationFilter(data->cochlearEnvelopes, data->modulationSignals);
    data->statsGenerator.computeStatistics(data->cochlearEnvelopes, data->modulationSignals,
            data->currentStats);

    // Compute distance gradient from jacobian of statistics function
    int numStats = data->targetStats.size();

    // Sum up L2 distance between current stats and target stats
    double distance = 0;
    for(int i = 0; i < numStats; i++)
        distance += pow(data->currentStats[i] - data->targetStats[i], 2);

    std::cout << "distance: " << distance << "\n";
    return distance;
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


    data->textureFilterer.modulationFilter(data->cochlearEnvelopes, data->modulationSignals);
    data->statsGenerator.computeStatistics(data->cochlearEnvelopes, data->modulationSignals,
            data->currentStats, data->textureFilterer.modulationBank(), data->jacobian);

    int numStats = data->targetStats.size();
    double gradient;
    for(int env = 0; env < numEnvelopes; env++)
        for(int sample = 0; sample < envelopeSize; sample++)
        {
            gradient = 0;
            for(int stat = 0; stat < numStats; stat++)
                gradient += 2 * (data->currentStats[stat] - data->targetStats[stat])
                    * data->jacobian[stat][env][sample];
            gsl_vector_set(df, env * envelopeSize + sample, gradient);
        }
}

void TextureSynthesizer::gradAndDist(const gsl_vector *v, void *params, double *f,
        gsl_vector *df)
{
    *f = distance(v, params);
    gradient(v, params, df);
}
