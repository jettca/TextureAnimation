#pragma once

#include <vector>
#include "aquila/aquila.h"

#include "Filtering/Filter.h"
#include "Filtering/Signal.h"

namespace TextureSynthesis
{
    class FilterBank
    {
    public:
        FilterBank();

        void addFilter(std::shared_ptr<Filter> filter);
        void apply(const Signal& signal, std::vector<Signal>& outSignals);
        void apply(const Aquila::SpectrumType& spectrum, double sampleRate,
                std::vector<Signal>& outSignals);
        void apply(Signal& signal, int index);

    private:
        std::vector<std::shared_ptr<Filter>> _filters;
    };
}
