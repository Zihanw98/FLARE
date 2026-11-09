#include "implicit_linear.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace hd_implicit::detail
{
void SparseMatrix::multiply(const Vector& x, Vector& y) const
{
    if (x.size() != size()) throw std::invalid_argument("Sparse matrix vector size mismatch");
    y.assign(size(), 0.0);
    for (std::size_t i = 0; i < size(); ++i)
        for (auto p = row[i]; p < row[i+1]; ++p)
            y[i] += value[p] * x[column[p]];
}

namespace
{
double dot(const Vector& a, const Vector& b)
{
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}
double norm(const Vector& a) { return std::sqrt(dot(a, a)); }

bool factor_ilu0(const SparseMatrix& A, Vector& lu)
{
    lu = A.value;
    for (std::size_t i = 0; i < A.size(); ++i)
    {
        double scale = 0;
        for (auto p = A.row[i]; p < A.row[i+1]; ++p) scale = std::max(scale, std::abs(lu[p]));
        for (auto p = A.row[i]; p < A.diagonal[i]; ++p)
        {
            const auto j = A.column[p];
            lu[p] /= lu[A.diagonal[j]];
            // Merge sorted upper rows; drop fill outside the original pattern.
            auto target = p + 1;
            for (auto q = A.diagonal[j]+1; q < A.row[j+1]; ++q)
            {
                while (target < A.row[i+1] && A.column[target] < A.column[q]) ++target;
                if (target == A.row[i+1]) break;
                if (A.column[target] == A.column[q]) lu[target] -= lu[p] * lu[q];
            }
        }
        if (!std::isfinite(lu[A.diagonal[i]])
            || std::abs(lu[A.diagonal[i]]) <= 1e-14 * std::max(scale, 1e-100)) return false;
    }
    return std::all_of(lu.begin(), lu.end(), [](double v) { return std::isfinite(v); });
}
}

bool solve_gmres(const SparseMatrix& A, const LinearAction& action,
                 const Vector& b, Vector& x, double relative_tolerance, int& iterations, int restart, int max_iterations)
{
    const auto n = A.size();
    if (!n || b.size() != n || restart < 1 || max_iterations < 1)
        throw std::invalid_argument("Invalid GMRES dimensions or iteration limits");
    const int m = static_cast<int>(std::min<std::size_t>(restart, n));
    Vector lu;
    const bool ilu_ok = factor_ilu0(A, lu);
    auto precondition = [&](const Vector& r, Vector& z) {
        z = r;
        if (!ilu_ok)
        {
            for (std::size_t i = 0; i < n; ++i)
            {
                const double d = A.value[A.diagonal[i]];
                if (std::isfinite(d) && std::abs(d) > 1e-100) z[i] /= d;
            }
            return;
        }
        for (std::size_t i = 0; i < n; ++i)
            for (auto p = A.row[i]; p < A.diagonal[i]; ++p) z[i] -= lu[p]*z[A.column[p]];
        for (std::size_t end = n; end > 0; --end)
        {
            const auto i = end-1;
            for (auto p = A.diagonal[i]+1; p < A.row[i+1]; ++p) z[i] -= lu[p]*z[A.column[p]];
            z[i] /= lu[A.diagonal[i]];
        }
    };
    // Restart bounds the Krylov basis storage independently of iteration count.
    std::vector<Vector> V(m+1, Vector(n));
    Vector H((m+1)*m), cosine(m), sine(m), g(m+1), y(m), r(n), w(n), z(n), combination(n);
    x.assign(n, 0.0);
    const double bnorm = norm(b);
    if (!std::isfinite(bnorm)) return false;
    if (bnorm == 0) return true;
    const double tolerance = relative_tolerance * bnorm;
    iterations = 0;
    while (iterations < max_iterations)
    {
        action(x, w);
        for (std::size_t i = 0; i < n; ++i) r[i] = b[i]-w[i];
        const double beta = norm(r);
        if (!std::isfinite(beta)) return false;
        if (beta <= tolerance) return true;
        for (std::size_t i = 0; i < n; ++i) V[0][i] = r[i]/beta;
        std::fill(H.begin(), H.end(), 0.0);
        std::fill(g.begin(), g.end(), 0.0);
        g[0] = beta;
        int used = 0;
        bool breakdown = false;
        for (int j = 0; j < m && iterations < max_iterations; ++j)
        {
            precondition(V[j], z);
            action(z, w);
            const double before = norm(w);
            if (!std::isfinite(before)) return false;
            // Twice modified Gram-Schmidt limits loss of orthogonality.
            for (int pass = 0; pass < 2; ++pass)
                for (int i = 0; i <= j; ++i)
                {
                    const double h = dot(w, V[i]);
                    H[i*m+j] += h;
                    for (std::size_t k = 0; k < n; ++k) w[k] -= h*V[i][k];
                }
            const double next = norm(w);
            H[(j+1)*m+j] = next;
            breakdown = next <= 32*std::numeric_limits<double>::epsilon()*before;
            if (!breakdown)
                for (std::size_t k = 0; k < n; ++k) V[j+1][k] = w[k]/next;
            for (int i = 0; i < j; ++i)
            {
                const double a = H[i*m+j], bnext = H[(i+1)*m+j];
                H[i*m+j] = cosine[i]*a + sine[i]*bnext;
                H[(i+1)*m+j] = -sine[i]*a + cosine[i]*bnext;
            }
            const double h = std::hypot(H[j*m+j], H[(j+1)*m+j]);
            if (!std::isfinite(h) || h == 0) return false;
            cosine[j] = H[j*m+j]/h;
            sine[j] = H[(j+1)*m+j]/h;
            H[j*m+j] = h;
            H[(j+1)*m+j] = 0;
            g[j+1] = -sine[j]*g[j];
            g[j] *= cosine[j];
            used = j+1;
            ++iterations;
            if (std::abs(g[j+1]) <= tolerance || breakdown) break;
        }
        for (int i = used-1; i >= 0; --i)
        {
            double rhs = g[i];
            for (int j = i+1; j < used; ++j) rhs -= H[i*m+j]*y[j];
            y[i] = rhs/H[i*m+i];
        }
        std::fill(combination.begin(), combination.end(), 0.0);
        for (int j = 0; j < used; ++j)
            for (std::size_t k = 0; k < n; ++k) combination[k] += y[j]*V[j][k];
        precondition(combination, z);
        for (std::size_t k = 0; k < n; ++k) x[k] += z[k];
        action(x, w);
        for (std::size_t k = 0; k < n; ++k) r[k] = b[k]-w[k];
        if (norm(r) <= tolerance) return true;
        if (breakdown) return false;
    }
    return false;
}
}
