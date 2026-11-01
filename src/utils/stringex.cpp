#include "utils.h"
namespace utils
{
    bool iequals(const std::string& a, const std::string& b)
    {
        return a.size() == b.size() &&
            std::equal(a.begin(), a.end(), b.begin(),
                [](unsigned char x, unsigned char y)
                {
                    return std::tolower(x) == std::tolower(y);
                });
    }
}