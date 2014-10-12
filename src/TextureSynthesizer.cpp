#include "TextureSynthesizer.h"

#include <complex>
#include "StatsGenerator.h"

using namespace TextureSynthesis;

TextureSynthesizer::TextureSynthesizer(const Signal& targetSignal)
    : _targetSignal(targetSignal),
    _targetStats(),
    _downsampleRate(400),
    _numSubbands(32),
    _stepSize(1000),
    _tolerance(10000)
{
    StatsGenerator::computeStatistics(targetSignal, _targetStats);
}

void TextureSynthesizer::synthesize(Signal& outSignal)
{
    std::vector<std::complex<double>> *outArray = &(outSignal._signal);

    gsl_multimin_function_fdf statmin;
    statmin.n = outArray->size() * _downsampleRate / _targetSignal._sampleRate * _numSubbands;
    statmin.f = &distanceFromTarget;
    statmin.df = &gradient;
    statmin.fdf = &gradAndDist;
    statmin.params = nullptr;

    gsl_vector *init = gsl_vector_alloc(outArray->size());
    srand(time(nullptr));
    for(int i = 0; i < statmin.n; i++)
        init->data[i] = 2 * rand() / RAND_MAX - 1;

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    T = gsl_multimin_fdfminimizer_conjugate_fr;
    s = gsl_multimin_fdfminimizer_alloc (T, outArray->size());
    gsl_multimin_fdfminimizer_set(s, &statmin, init, _stepSize, _tolerance);

    int iter = 0;
    int status;
    do
    {
        iter++;
        std::cout << "iter: " << iter << "\n";
        status = gsl_multimin_fdfminimizer_iterate(s);

        if (status)
            break;

        status = gsl_multimin_test_gradient (s->gradient, pow(2, 10));
    }
    while (status == GSL_CONTINUE && iter < 5);

    for(int i = 0; i < outArray->size(); i++)
        (*outArray)[i] = gsl_vector_get(s->x, i);

    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (init);
}

// TODO: Write these functions!
double TextureSynthesizer::distanceFromTarget(const gsl_vector *v, void *params)
{
    OptimizationData *data = (OptimizationData*)params;

    double distance = 0;



//    Signal signal(outLen, sampleRate);
//    for(int i = 0; i < outLen; i++)
//    {
//        signal._signal[i] = v->data[i];
//    }
//
//    std::vector<std::complex<double>> statistics;
//    StatsGenerator::computeStatistics(signal, statistics);
//
//    double distance = 0;
//    int numStats = targetStats.size();
//    for(int i = 0; i < numStats; i++)
//    {
//        distance += std::norm(statistics.at(i) - targetStats.at(i));
//    }

    return distance;
}

void TextureSynthesizer::gradient(const gsl_vector *v, void *params, gsl_vector *df)
{
    OptimizationData *data = (OptimizationData*)params;

//    // Turn v (subband envelopes) into signal
//    std::vector<Signal> upsampledEnvelopes;
//    Signal signal(outLen, sampleRate);
//
//    double d1, d2;
//    d1 = distanceFromTarget(v, params);
//
//    gsl_vector *copy = gsl_vector_alloc(sourceLen);
//    memcpy(copy->data, v->data, sourceLen);
//
//    for(int i = 0; i < sourceLen; i++)
//    {
//        double oldVal = copy->data[i];
//        double delta = copy->data[i] - oldVal;
//
//        d2 = distanceFromTarget(copy, params);
//        df->data[i] = (d2 - d1) / delta;
//
//        copy->data[i] = oldVal;
//    }
//
//    gsl_vector_free(copy);
}

void TextureSynthesizer::gradAndDist(const gsl_vector *v, void *params, double *f,
        gsl_vector *df)
{
    *f = distanceFromTarget(v, params);
    gradient(v, params, df);
}
