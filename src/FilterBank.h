#pragma once

#include <vector>
#include "aquila/aquila.h"

#include "Filter.h"
#include "Signal.h"

namespace TextureSynthesis
{
    class FilterBank
    {
    public:
        FilterBank();

        void addFilter(std::shared_ptr<Filter> filter);
        std::vector<Signal> apply(Aquila::SignalSource signal);
        std::vector<Signal> apply(Aquila::SpectrumType spectrum, double sampleRate);

    private:
        std::vector<std::shared_ptr<Filter>> _filters;
    };
}
