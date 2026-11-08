#include "utils/performance.h"
#include <cassert>
#include <cmath>
#include <sstream>
#include <stdexcept>

using namespace performance;
static const Measurement& result(Section s) { return measurements()[static_cast<std::size_t>(s)]; }

int main()
{
    configure(Mode::Off);
    {
        Scope run(Section::Run);
        count(Counter::Steps);
        measure(Section::Riemann, [] {}, true);
    }
    for (const auto& item : measurements()) assert(item.calls == 0);
    assert(counters[0] == 0);
    std::ostringstream output;
    report(output, true);
    assert(output.str().empty());

    configure(Mode::Summary);
    {
        Scope run(Section::Run);
        Scope hydro(Section::Hydro);
        bool reset_rejected = false;
        try { configure(Mode::Off); }
        catch (const std::logic_error&) { reset_rejected = true; }
        assert(reset_rejected && mode == Mode::Summary);
        measure(Section::Riemann, [] {}, true);
        hydro.stop();
        hydro.stop(); // Explicit stop followed by destruction must not double-count.
        count(Counter::Steps, 3);
    }
    assert(result(Section::Riemann).calls == 0);
    assert(result(Section::Hydro).calls == 1);
    assert(counters[0] == 3);
    assert(result(Section::Run).inclusive_seconds >= result(Section::Hydro).inclusive_seconds);
    double sum_self = 0;
    for (const auto& item : measurements()) sum_self += item.self_seconds;
    assert(std::abs(sum_self - result(Section::Run).inclusive_seconds) < 1e-9);

    configure(Mode::Detail);
    try
    {
        Scope run(Section::Run);
        measure(Section::Riemann, [] { throw std::runtime_error("trial failure"); }, true);
    }
    catch (const std::runtime_error&) {}
    assert(result(Section::Riemann).calls == 1);
    assert(result(Section::Run).calls == 1);
    output.precision(3);
    report(output, false);
    assert(output.precision() == 3);
    assert(output.str().find("failed; partial timings") != std::string::npos);
    assert(output.str().find("EVERY face") != std::string::npos);
    std::ostringstream json;
    write_json(json, false);
    assert(json.str().find("\"status\": \"failed\"") != std::string::npos);

    // Repeated/nested same-category scopes must still partition self time.
    configure(Mode::Detail);
    {
        Scope run(Section::Run);
        Scope a(Section::Hydro);
        Scope b(Section::Hydro);
    }
    assert(result(Section::Hydro).calls == 2);
    sum_self = 0;
    for (const auto& item : measurements()) sum_self += item.self_seconds;
    assert(std::abs(sum_self - result(Section::Run).inclusive_seconds) < 1e-9);
    configure(Mode::Off);
}
