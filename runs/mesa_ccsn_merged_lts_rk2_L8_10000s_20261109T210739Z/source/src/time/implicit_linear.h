#pragma once
#include <cstddef>
#include <functional>
#include <vector>

// Private linear algebra for the implicit integrator; no external dependency.
namespace hd_implicit::detail
{
using Vector = std::vector<double>;
using LinearAction = std::function<void(const Vector&, Vector&)>;

struct SparseMatrix
{
    // Sorted scalar CSR, including a structural diagonal in every row.
    std::vector<std::size_t> row, column, diagonal;
    Vector value;
    std::size_t size() const { return diagonal.size(); }
    void multiply(const Vector& x, Vector& y) const;
};

// Right-preconditioned restarted GMRES. ILU(0) keeps exactly the CSR pattern;
// failure to factor falls back to diagonal scaling. Always verify the true
// unpreconditioned residual before accepting, including at happy breakdown.
bool solve_gmres(const SparseMatrix& matrix, const LinearAction& action,
                 const Vector& b, Vector& x, double tolerance, int& iterations,
                 int restart = 30, int max_iterations = 120);
}
