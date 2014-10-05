#include "aquila/aquila.h"
#include "FilterBank.h"
#include "LowpassFilter.h"
#include "CochlearFilter.h"
#include "AudioDevice.h"
#include "AuditoryModel.h"

#include <complex>
#include <iostream>
#include <cmath>
#include <gsl/gsl_multimin.h>

using namespace TextureSynthesis;

static std::vector<std::complex<double>> targetStats;
static int sourceLen = pow(2, 8);
static double sampleRate = 1000;

double distanceFromTarget(const gsl_vector *v, void *params)
{
    Signal signal(sourceLen, sampleRate);
    for(int i = 0; i < sourceLen; i++)
        signal._signal[i] = v->data[i];

    std::vector<std::complex<double>> statistics;
    computeStatistics(signal, statistics);

    double distance = 0;
    int numStats = targetStats.size();
    for(int i = 0; i < numStats; i++)
    {
        distance += std::norm(statistics.at(i) - targetStats.at(i));
    }

    return distance;
}

void gradient(const gsl_vector *v, void *params, gsl_vector *df)
{
    double change = 1.05;
    double d1, d2;
    d1 = distanceFromTarget(v, params);

    gsl_vector *copy = gsl_vector_alloc(sourceLen);
    memcpy(copy->data, v->data, sourceLen);

    for(int i = 0; i < sourceLen; i++)
    {
        double oldVal = copy->data[i];
        copy->data[i] *= change;
        double delta = copy->data[i] - oldVal;

        d2 = distanceFromTarget(copy, params);
        df->data[i] = (d2 - d1) / delta;

        copy->data[i] = oldVal;
    }

    gsl_vector_free(copy);
}

void fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df)
{
    *f = distanceFromTarget(v, params);
    gradient(v, params, df);
}

int main()
{
    Aquila::SineGenerator generator(sampleRate);
    generator.setFrequency(300).setAmplitude(1).generate(sourceLen);
    computeStatistics(Signal(generator), targetStats);

    // GSL minimization

    gsl_multimin_function_fdf statmin;
    statmin.n = sourceLen;
    statmin.f = &distanceFromTarget;
    statmin.df = &gradient;
    statmin.fdf = &fdf;
    statmin.params = nullptr;

    gsl_vector *init = gsl_vector_alloc(sourceLen);
    srand(time(nullptr));
    for(int i = 0; i < sourceLen; i++)
        init->data[i] = 2 * rand() / RAND_MAX - 1;

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    T = gsl_multimin_fdfminimizer_conjugate_fr;
    s = gsl_multimin_fdfminimizer_alloc (T, sourceLen);
    gsl_multimin_fdfminimizer_set(s, &statmin, init, 0.1, 0.05);

    int iter = 0;
    int status;
    do
    {
        iter++;
        std::cout << iter << "\n";
        status = gsl_multimin_fdfminimizer_iterate(s);

        if (status)
            break;

        status = gsl_multimin_test_gradient (s->gradient, 1e-3);
    }
    while (status == GSL_CONTINUE && iter < 5);

    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (init);

    return 0;
}
