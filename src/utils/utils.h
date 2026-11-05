/**
 * @file utils.h
 * @brief Mathematical and string utilities.
 */
#pragma once
#include <algorithm>
#include <cctype>
#include <string>
namespace utils
{
    /**
     * @brief Compare two strings case-insensitively using character conversion.
     * @param[in] a First string.
     * @param[in] b Second string.
     * @return True if lengths match and corresponding lowercased bytes match.
     * @note Uses std::tolower with unsigned-byte inputs and the current C locale; this is not Unicode case folding.
     */
    bool iequals(const std::string& a, const std::string& b);
}