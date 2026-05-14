#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <sqlite3.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <csignal>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;
using Poly = std::vector<unsigned long>;
using FactorId = std::uint32_t;
using BranchKey = std::vector<std::uint32_t>;

namespace {

volatile std::sig_atomic_t stop_requested_flag = 0;

void request_stop(int) {
    stop_requested_flag = 1;
}

bool stop_requested() {
    return stop_requested_flag != 0;
}

struct VectorHash {
    std::size_t operator()(const BranchKey& values) const {
        std::size_t seed = values.size();
        for (auto value : values) {
            seed ^= static_cast<std::size_t>(value) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

struct Options {
    unsigned long p = 0;
    int genus = -1;
    int genus_start = -1;
    int genus_end = -1;
    int genus_step = 1;
    int max_sparsity = -1;
    std::uint64_t limit = 0;
    int progress_interval = 1000;
    std::string mode = "enumerate";
    unsigned long random_seed = 1;
    int random_max_factors = 0;
    std::uint64_t irreducible_memory_budget_mb = 1024;
    std::filesystem::path out;
    std::filesystem::path out_dir = std::filesystem::path("results");
};

struct BranchCandidate {
    unsigned long leading = 1;
    std::vector<Poly> factors;
    bool infinity_branch = false;
    std::vector<int> pattern;
    std::vector<FactorId> factor_ids;
};

struct FactorTransform {
    bool computed = false;
    bool valid = false;
    unsigned long scalar = 0;
    FactorId factor_id = 0;
};

struct BranchCanonicalInfo {
    BranchKey canonical_key;
    std::uint64_t orbit_size = 1;
};

struct Stats {
    std::uint64_t processed = 0;
    std::uint64_t sparse = 0;
    std::uint64_t duplicate = 0;
    std::uint64_t rejected_hasse_witt = 0;
    std::uint64_t rejected_exact = 0;
    std::uint64_t canonicalized = 0;
    cpp_int total_presentations = -1;
    cpp_int sparse_presentations = 0;
};

struct BranchDivisorTypeState {
    int model = 0;
    int total_degree = 0;
    bool infinity_branch = false;
    std::vector<int> pattern;
    int factor_count = 0;
    int forced_mod2_sparsity = 0;
    std::uint64_t leading_count = 1;
    cpp_int presentations = 0;
    std::uint64_t attempts = 0;
    std::uint64_t hasse_witt_passes = 0;
    std::uint64_t exact_rejections = 0;
    std::uint64_t sparse_hits = 0;
};

enum class CandidateOutcome {
    HasseWittRejected,
    Duplicate,
    ExactRejected,
    Sparse,
};

unsigned long mod_pow(unsigned long base, unsigned long exponent, unsigned long p) {
    unsigned long result = 1 % p;
    base %= p;
    while (exponent > 0) {
        if (exponent & 1UL) {
            result = static_cast<unsigned long>((static_cast<unsigned long long>(result) * base) % p);
        }
        base = static_cast<unsigned long>((static_cast<unsigned long long>(base) * base) % p);
        exponent >>= 1UL;
    }
    return result;
}

unsigned long mod_inv(unsigned long value, unsigned long p) {
    if (value % p == 0) {
        throw std::runtime_error("division by zero modulo p");
    }
    return mod_pow(value, p - 2, p);
}

void trim(Poly& poly) {
    while (poly.size() > 1 && poly.back() == 0) {
        poly.pop_back();
    }
}

std::string json_poly(const Poly& poly) {
    std::ostringstream out;
    out << "[";
    for (std::size_t i = 0; i < poly.size(); ++i) {
        if (i) out << ",";
        out << poly[i];
    }
    out << "]";
    return out.str();
}

std::string json_ints(const std::vector<int>& values) {
    std::ostringstream out;
    out << "[";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) out << ",";
        out << values[i];
    }
    out << "]";
    return out.str();
}

std::string json_factorization_partition(const std::vector<int>& degree_counts) {
    std::vector<int> partition;
    for (int degree = 1; degree < static_cast<int>(degree_counts.size()); ++degree) {
        int multiplicity = degree_counts[static_cast<std::size_t>(degree)];
        for (int i = 0; i < multiplicity; ++i) partition.push_back(degree);
    }
    return json_ints(partition);
}

std::string json_branch_divisor_type(bool infinity_branch, const std::vector<int>& degree_counts) {
    std::vector<int> branch_type;
    if (infinity_branch) branch_type.push_back(0);
    for (int degree = 1; degree < static_cast<int>(degree_counts.size()); ++degree) {
        int multiplicity = degree_counts[static_cast<std::size_t>(degree)];
        for (int i = 0; i < multiplicity; ++i) branch_type.push_back(degree);
    }
    return json_ints(branch_type);
}

std::string json_polys(std::vector<Poly> polynomials) {
    std::sort(polynomials.begin(), polynomials.end(), [](const Poly& a, const Poly& b) {
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    });
    std::ostringstream out;
    out << "[";
    for (std::size_t i = 0; i < polynomials.size(); ++i) {
        if (i) out << ",";
        out << json_poly(polynomials[i]);
    }
    out << "]";
    return out.str();
}

std::string json_cpp_ints(const std::vector<cpp_int>& values) {
    std::ostringstream out;
    out << "[";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) out << ",";
        out << values[i];
    }
    out << "]";
    return out.str();
}

std::string cpp_int_to_string(const cpp_int& value) {
    std::ostringstream out;
    out << value;
    return out.str();
}

std::string key_poly(const Poly& poly) {
    std::ostringstream out;
    out << poly.size() << ":";
    for (auto c : poly) {
        out << c << ",";
    }
    return out.str();
}

std::string key_int_vector(const std::vector<unsigned long>& values) {
    std::ostringstream out;
    for (auto value : values) {
        out << value << ",";
    }
    return out.str();
}

cpp_int pow_cpp_int(unsigned long base, int exponent) {
    cpp_int result = 1;
    for (int i = 0; i < exponent; ++i) result *= base;
    return result;
}

int mobius(int n) {
    int primes = 0;
    for (int d = 2; d * d <= n; ++d) {
        if (n % d != 0) continue;
        int count = 0;
        while (n % d == 0) {
            n /= d;
            ++count;
        }
        if (count > 1) return 0;
        ++primes;
    }
    if (n > 1) ++primes;
    return primes % 2 == 0 ? 1 : -1;
}

std::uint64_t irreducible_count_formula(unsigned long p, int degree) {
    cpp_int total = 0;
    for (int divisor = 1; divisor <= degree; ++divisor) {
        if (degree % divisor != 0) continue;
        int mu = mobius(divisor);
        if (mu == 0) continue;
        cpp_int term = pow_cpp_int(p, degree / divisor);
        total += mu > 0 ? term : -term;
    }
    total /= degree;
    if (total < 0 || total > std::numeric_limits<std::uint64_t>::max()) {
        throw std::runtime_error("irreducible count does not fit in uint64_t");
    }
    return total.convert_to<std::uint64_t>();
}

std::uint64_t saturating_add(std::uint64_t lhs, std::uint64_t rhs) {
    if (lhs > std::numeric_limits<std::uint64_t>::max() - rhs) {
        return std::numeric_limits<std::uint64_t>::max();
    }
    return lhs + rhs;
}

std::uint64_t saturating_mul(std::uint64_t lhs, std::uint64_t rhs) {
    if (rhs != 0 && lhs > std::numeric_limits<std::uint64_t>::max() / rhs) {
        return std::numeric_limits<std::uint64_t>::max();
    }
    return lhs * rhs;
}

unsigned long smallest_nonsquare(unsigned long p) {
    std::set<unsigned long> squares;
    for (unsigned long x = 1; x < p; ++x) {
        squares.insert((x * x) % p);
    }
    for (unsigned long x = 2; x < p; ++x) {
        if (!squares.count(x)) return x;
    }
    throw std::runtime_error("no nonsquare in prime field");
}

Poly flint_to_poly(const nmod_poly_t src) {
    slong len = nmod_poly_length(src);
    if (len <= 0) return Poly{0};
    Poly out(static_cast<std::size_t>(len));
    for (slong i = 0; i < len; ++i) {
        out[static_cast<std::size_t>(i)] = nmod_poly_get_coeff_ui(src, i);
    }
    trim(out);
    return out;
}

void poly_to_flint(nmod_poly_t dst, const Poly& poly, unsigned long p) {
    nmod_poly_zero(dst);
    for (std::size_t i = 0; i < poly.size(); ++i) {
        nmod_poly_set_coeff_ui(dst, static_cast<slong>(i), poly[i] % p);
    }
}

bool is_irreducible_flint(const Poly& poly, unsigned long p) {
    nmod_poly_t f;
    nmod_poly_init(f, p);
    poly_to_flint(f, poly, p);
    int ok = nmod_poly_is_irreducible(f);
    nmod_poly_clear(f);
    return ok != 0;
}

Poly multiply_flint(const Poly& a, const Poly& b, unsigned long p) {
    nmod_poly_t fa, fb, fc;
    nmod_poly_init(fa, p);
    nmod_poly_init(fb, p);
    nmod_poly_init(fc, p);
    poly_to_flint(fa, a, p);
    poly_to_flint(fb, b, p);
    nmod_poly_mul(fc, fa, fb);
    Poly out = flint_to_poly(fc);
    nmod_poly_clear(fa);
    nmod_poly_clear(fb);
    nmod_poly_clear(fc);
    return out;
}

Poly pow_flint(const Poly& a, unsigned long exponent, unsigned long p) {
    nmod_poly_t fa, fb;
    nmod_poly_init(fa, p);
    nmod_poly_init(fb, p);
    poly_to_flint(fa, a, p);
    nmod_poly_pow(fb, fa, exponent);
    Poly out = flint_to_poly(fb);
    nmod_poly_clear(fa);
    nmod_poly_clear(fb);
    return out;
}

Poly monic_poly(Poly poly, unsigned long p) {
    trim(poly);
    if (poly.empty() || (poly.size() == 1 && poly[0] == 0)) return Poly{0};
    unsigned long inverse = mod_inv(poly.back(), p);
    for (auto& coefficient : poly) coefficient = (coefficient * inverse) % p;
    trim(poly);
    return poly;
}

Poly product_tree_flint(std::vector<Poly> factors, unsigned long p) {
    if (factors.empty()) return Poly{1};
    while (factors.size() > 1) {
        std::vector<Poly> next;
        next.reserve((factors.size() + 1) / 2);
        for (std::size_t i = 0; i < factors.size(); i += 2) {
            if (i + 1 == factors.size()) {
                next.push_back(std::move(factors[i]));
            } else {
                next.push_back(multiply_flint(factors[i], factors[i + 1], p));
            }
        }
        factors = std::move(next);
    }
    return factors[0];
}

Poly expand_candidate(const BranchCandidate& candidate, unsigned long p) {
    Poly product_poly = product_tree_flint(candidate.factors, p);
    if (candidate.leading != 1) {
        for (auto& coefficient : product_poly) {
            coefficient = (coefficient * candidate.leading) % p;
        }
    }
    trim(product_poly);
    return product_poly;
}

std::vector<std::vector<unsigned long>> binomial_table(int degree, unsigned long p) {
    std::vector<std::vector<unsigned long>> table(degree + 1);
    table[0] = {1};
    for (int n = 1; n <= degree; ++n) {
        table[n].assign(n + 1, 0);
        table[n][0] = 1;
        table[n][n] = 1;
        for (int k = 1; k < n; ++k) {
            table[n][k] = (table[n - 1][k - 1] + table[n - 1][k]) % p;
        }
    }
    return table;
}

std::vector<unsigned long> linear_power(
    unsigned long alpha,
    unsigned long beta,
    int exponent,
    unsigned long p,
    const std::vector<std::vector<unsigned long>>& binom
) {
    std::vector<unsigned long> result(exponent + 1, 0);
    for (int k = 0; k <= exponent; ++k) {
        result[k] = (((binom[exponent][k] * mod_pow(alpha, static_cast<unsigned long>(k), p)) % p)
            * mod_pow(beta, static_cast<unsigned long>(exponent - k), p)) % p;
    }
    return result;
}

using Matrix4 = std::array<unsigned long, 4>;

std::vector<Matrix4> pgl2_representatives(unsigned long p) {
    std::map<std::array<unsigned long, 4>, Matrix4> reps;
    for (unsigned long a = 0; a < p; ++a)
        for (unsigned long b = 0; b < p; ++b)
            for (unsigned long c = 0; c < p; ++c)
                for (unsigned long d = 0; d < p; ++d) {
                    unsigned long det = (a * d + p - (b * c) % p) % p;
                    if (det == 0) continue;
                    Matrix4 entries{a, b, c, d};
                    unsigned long first = 0;
                    for (auto entry : entries) {
                        if (entry != 0) {
                            first = entry;
                            break;
                        }
                    }
                    unsigned long scale = mod_inv(first, p);
                    Matrix4 rep{};
                    for (int i = 0; i < 4; ++i) {
                        rep[i] = (entries[i] * scale) % p;
                    }
                    reps[rep] = rep;
                }
    std::vector<Matrix4> out;
    for (const auto& item : reps) out.push_back(item.second);
    return out;
}

std::vector<std::vector<unsigned long>> action_matrix(
    const Matrix4& matrix,
    int degree,
    unsigned long p,
    const std::vector<std::vector<unsigned long>>& binom
) {
    auto [a, b, c, d] = matrix;
    std::vector<std::vector<unsigned long>> powers_ab;
    std::vector<std::vector<unsigned long>> powers_cd;
    for (int e = 0; e <= degree; ++e) {
        powers_ab.push_back(linear_power(a, b, e, p, binom));
        powers_cd.push_back(linear_power(c, d, e, p, binom));
    }

    std::vector<std::vector<unsigned long>> columns(degree + 1, std::vector<unsigned long>(degree + 1, 0));
    for (int i = 0; i <= degree; ++i) {
        const auto& left = powers_ab[i];
        const auto& right = powers_cd[degree - i];
        for (int lx = 0; lx < static_cast<int>(left.size()); ++lx) {
            if (left[lx] == 0) continue;
            for (int rx = 0; rx < static_cast<int>(right.size()); ++rx) {
                if (right[rx] == 0) continue;
                int x_degree = lx + rx;
                columns[i][x_degree] = (columns[i][x_degree] + left[lx] * right[rx]) % p;
            }
        }
    }
    return columns;
}

std::vector<unsigned long> apply_action(
    const std::vector<unsigned long>& binary_form,
    const std::vector<std::vector<unsigned long>>& action,
    unsigned long p
) {
    std::vector<unsigned long> result(binary_form.size(), 0);
    for (std::size_t i = 0; i < binary_form.size(); ++i) {
        if (binary_form[i] == 0) continue;
        for (std::size_t j = 0; j < binary_form.size(); ++j) {
            result[j] = (result[j] + binary_form[i] * action[i][j]) % p;
        }
    }
    return result;
}

std::vector<unsigned long> normalize_square_scalar(const std::vector<unsigned long>& form, unsigned long p) {
    std::set<unsigned long> squares;
    for (unsigned long x = 1; x < p; ++x) squares.insert((x * x) % p);
    std::optional<std::vector<unsigned long>> best;
    for (auto scalar : squares) {
        std::vector<unsigned long> scaled = form;
        for (auto& value : scaled) value = (value * scalar) % p;
        if (!best || scaled < *best) best = std::move(scaled);
    }
    if (!best) throw std::runtime_error("could not normalize binary form");
    return *best;
}

class FiniteExtensionInt {
public:
    unsigned long p;
    int degree;
    std::uint64_t size;
    Poly modulus;
    std::optional<std::vector<int8_t>> character_table;

    FiniteExtensionInt(unsigned long prime, int extension_degree)
        : p(prime), degree(extension_degree), size(1), modulus(find_modulus(prime, extension_degree)) {
        for (int i = 0; i < degree; ++i) {
            if (size > std::numeric_limits<std::uint64_t>::max() / p) {
                throw std::runtime_error("extension size overflow");
            }
            size *= p;
        }
    }

    static Poly find_modulus(unsigned long p, int degree) {
        if (degree == 1) return Poly{0, 1};
        std::uint64_t count = 1;
        for (int i = 0; i < degree; ++i) count *= p;
        for (std::uint64_t index = 0; index < count; ++index) {
            std::uint64_t n = index;
            Poly poly(static_cast<std::size_t>(degree + 1), 0);
            for (int i = 0; i < degree; ++i) {
                poly[static_cast<std::size_t>(i)] = n % p;
                n /= p;
            }
            poly[static_cast<std::size_t>(degree)] = 1;
            if (is_irreducible_flint(poly, p)) return poly;
        }
        throw std::runtime_error("failed to find extension modulus");
    }

    unsigned long coefficient(std::uint64_t value, int index) const {
        for (int i = 0; i < index; ++i) value /= p;
        return value % p;
    }

    std::uint64_t add(std::uint64_t lhs, std::uint64_t rhs) const {
        std::uint64_t result = 0;
        std::uint64_t place = 1;
        for (int i = 0; i < degree; ++i) {
            result += ((lhs % p) + (rhs % p)) % p * place;
            lhs /= p;
            rhs /= p;
            place *= p;
        }
        return result;
    }

    std::uint64_t constant(unsigned long value) const {
        return value % p;
    }

    std::uint64_t multiply(std::uint64_t lhs, std::uint64_t rhs) const {
        std::vector<unsigned long> product(static_cast<std::size_t>(2 * degree - 1), 0);
        std::uint64_t left = lhs;
        for (int i = 0; i < degree; ++i) {
            unsigned long lc = left % p;
            left /= p;
            if (lc == 0) continue;
            std::uint64_t right = rhs;
            for (int j = 0; j < degree; ++j) {
                unsigned long rc = right % p;
                right /= p;
                if (rc == 0) continue;
                product[static_cast<std::size_t>(i + j)] = (product[static_cast<std::size_t>(i + j)] + lc * rc) % p;
            }
        }
        for (int d = 2 * degree - 2; d >= degree; --d) {
            unsigned long coeff = product[static_cast<std::size_t>(d)];
            if (coeff == 0) continue;
            for (int j = 0; j < degree; ++j) {
                int idx = d - degree + j;
                product[static_cast<std::size_t>(idx)] =
                    (product[static_cast<std::size_t>(idx)] + p - (coeff * modulus[static_cast<std::size_t>(j)]) % p) % p;
            }
        }
        std::uint64_t result = 0;
        std::uint64_t place = 1;
        for (int i = 0; i < degree; ++i) {
            result += product[static_cast<std::size_t>(i)] * place;
            place *= p;
        }
        return result;
    }

    std::uint64_t evaluate(const Poly& poly, std::uint64_t x) const {
        std::uint64_t result = 0;
        for (auto it = poly.rbegin(); it != poly.rend(); ++it) {
            result = add(multiply(result, x), constant(*it));
        }
        return result;
    }

    const std::vector<int8_t>& quadratic_characters() {
        if (!character_table) {
            std::vector<int8_t> chars(static_cast<std::size_t>(size), -1);
            chars[0] = 0;
            for (std::uint64_t x = 1; x < size; ++x) {
                chars[static_cast<std::size_t>(multiply(x, x))] = 1;
            }
            character_table = std::move(chars);
        }
        return *character_table;
    }

    int character(std::uint64_t value) {
        return quadratic_characters()[static_cast<std::size_t>(value)];
    }
};

class SqliteWriter {
public:
    sqlite3* db = nullptr;
    std::filesystem::path db_path;
    bool in_transaction = false;
    int pending_writes = 0;
    static constexpr int write_batch_size = 100;

    explicit SqliteWriter(const std::filesystem::path& path) : db_path(path) {
        std::filesystem::create_directories(path.parent_path());
        if (sqlite3_open(path.string().c_str(), &db) != SQLITE_OK) {
            throw std::runtime_error("failed to open sqlite output");
        }
        exec("PRAGMA journal_mode=WAL;");
        exec("PRAGMA synchronous=NORMAL;");
        exec(
            "CREATE TABLE IF NOT EXISTS sparse_curves ("
            "canonical_key TEXT PRIMARY KEY,"
            "coefficients TEXT NOT NULL,"
            "branch_factors TEXT,"
            "branch_infinity_branch INTEGER,"
            "branch_leading_coefficient INTEGER,"
            "branch_factorization_pattern TEXT,"
            "branch_divisor_type TEXT,"
            "lpoly TEXT NOT NULL,"
            "sparsity INTEGER NOT NULL,"
            "rational_branch_count INTEGER NOT NULL"
            ");"
            "CREATE TABLE IF NOT EXISTS enumeration_summary ("
            "id INTEGER PRIMARY KEY CHECK (id = 1),"
            "prime INTEGER NOT NULL,"
            "genus INTEGER NOT NULL,"
            "max_sparsity INTEGER,"
            "enumeration_mode TEXT NOT NULL,"
            "irreducible_memory_budget_mb INTEGER,"
            "limit_count INTEGER,"
            "total_coefficient_vectors TEXT,"
            "processed INTEGER NOT NULL,"
            "sparse_presentations TEXT,"
            "sparse_isomorphism_classes INTEGER NOT NULL,"
            "canonicalized_isomorphism_classes INTEGER NOT NULL,"
            "elapsed_seconds REAL NOT NULL,"
            "status_counts TEXT NOT NULL"
            ");"
            "CREATE TABLE IF NOT EXISTS enumeration_progress ("
            "position INTEGER PRIMARY KEY,"
            "processed INTEGER NOT NULL,"
            "sparse_presentations TEXT,"
            "sparse_isomorphism_classes INTEGER NOT NULL,"
            "canonicalized_isomorphism_classes INTEGER NOT NULL,"
            "elapsed_seconds REAL NOT NULL"
            ");"
        );
        ensure_column("sparse_curves", "branch_factors", "TEXT");
        ensure_column("sparse_curves", "branch_infinity_branch", "INTEGER");
        ensure_column("sparse_curves", "branch_leading_coefficient", "INTEGER");
        ensure_column("sparse_curves", "branch_factorization_pattern", "TEXT");
        ensure_column("sparse_curves", "branch_divisor_type", "TEXT");
        ensure_column("enumeration_summary", "total_coefficient_vectors", "TEXT");
        ensure_column("enumeration_summary", "sparse_presentations", "TEXT");
        ensure_column("enumeration_summary", "irreducible_memory_budget_mb", "INTEGER");
        ensure_column("enumeration_progress", "sparse_presentations", "TEXT");
        begin_transaction();
    }

    ~SqliteWriter() {
        try {
            close();
        } catch (...) {
        }
    }

    void close() {
        if (!db) return;
        if (in_transaction) commit(true);
        exec("PRAGMA wal_checkpoint(TRUNCATE);");
        exec("PRAGMA journal_mode=DELETE;");
        sqlite3_close(db);
        db = nullptr;

        std::error_code ignored;
        std::filesystem::remove(db_path.string() + "-wal", ignored);
        std::filesystem::remove(db_path.string() + "-shm", ignored);
    }

    void exec(const std::string& sql) {
        char* error = nullptr;
        if (sqlite3_exec(db, sql.c_str(), nullptr, nullptr, &error) != SQLITE_OK) {
            std::string message = error ? error : "sqlite error";
            sqlite3_free(error);
            throw std::runtime_error(message);
        }
    }

    bool column_exists(const std::string& table, const std::string& column) {
        std::string sql = "PRAGMA table_info(" + table + ")";
        sqlite3_stmt* stmt = nullptr;
        if (sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr) != SQLITE_OK) {
            throw std::runtime_error("failed to inspect sqlite schema");
        }
        bool found = false;
        while (sqlite3_step(stmt) == SQLITE_ROW) {
            const unsigned char* name = sqlite3_column_text(stmt, 1);
            if (name && column == reinterpret_cast<const char*>(name)) {
                found = true;
                break;
            }
        }
        sqlite3_finalize(stmt);
        return found;
    }

    void ensure_column(const std::string& table, const std::string& column, const std::string& type) {
        if (!column_exists(table, column)) {
            exec("ALTER TABLE " + table + " ADD COLUMN " + column + " " + type + ";");
        }
    }

    void begin_transaction() {
        if (!in_transaction) {
            exec("BEGIN TRANSACTION;");
            in_transaction = true;
            pending_writes = 0;
        }
    }

    void commit(bool force = false) {
        if (!in_transaction) return;
        if (!force && pending_writes < write_batch_size) return;
        exec("COMMIT;");
        in_transaction = false;
        pending_writes = 0;
        if (!force) begin_transaction();
    }

    void mark_write() {
        ++pending_writes;
        commit(false);
    }

    static std::string quote(const std::string& value) {
        std::string out = "'";
        for (char c : value) {
            if (c == '\'') out += "''";
            else out += c;
        }
        out += "'";
        return out;
    }

    void insert_sparse(
        const std::string& canonical_key,
        const Poly& coefficients,
        const BranchCandidate& candidate,
        const std::vector<cpp_int>& lpoly,
        int sparsity,
        int rational_branch_count
    ) {
        std::ostringstream sql;
        sql << "INSERT OR REPLACE INTO sparse_curves ("
            << "canonical_key,"
            << "coefficients,"
            << "branch_factors,"
            << "branch_infinity_branch,"
            << "branch_leading_coefficient,"
            << "branch_factorization_pattern,"
            << "branch_divisor_type,"
            << "lpoly,"
            << "sparsity,"
            << "rational_branch_count"
            << ") VALUES ("
            << quote(canonical_key) << ","
            << quote(json_poly(coefficients)) << ","
            << quote(json_polys(candidate.factors)) << ","
            << (candidate.infinity_branch ? 1 : 0) << ","
            << candidate.leading << ","
            << quote(json_factorization_partition(candidate.pattern)) << ","
            << quote(json_branch_divisor_type(candidate.infinity_branch, candidate.pattern)) << ","
            << quote(json_cpp_ints(lpoly)) << ","
            << sparsity << ","
            << rational_branch_count << ");";
        exec(sql.str());
        mark_write();
    }

    void write_progress(const Stats& stats, double elapsed) {
        std::ostringstream sql;
        sql << "INSERT OR REPLACE INTO enumeration_progress ("
            << "position,"
            << "processed,"
            << "sparse_presentations,"
            << "sparse_isomorphism_classes,"
            << "canonicalized_isomorphism_classes,"
            << "elapsed_seconds"
            << ") VALUES ("
            << stats.processed << ","
            << stats.processed << ","
            << quote(cpp_int_to_string(stats.sparse_presentations)) << ","
            << stats.sparse << ","
            << stats.canonicalized << ","
            << elapsed << ");";
        exec(sql.str());
        mark_write();
    }

    void write_summary(const Options& opts, const Stats& stats, double elapsed) {
        std::string status_label = stop_requested()
            ? "interrupted"
            : (stats.total_presentations >= 0 && cpp_int(stats.processed) >= stats.total_presentations ? "complete" : "partial");
        std::ostringstream status;
        status << "{\"duplicate\":" << stats.duplicate
               << ",\"rejected_exact\":" << stats.rejected_exact
               << ",\"rejected_hasse_witt_uncanonicalized\":" << stats.rejected_hasse_witt
               << ",\"sparse\":" << stats.sparse
               << ",\"run_status\":\"" << status_label << "\"}";
        std::ostringstream sql;
        sql << "INSERT OR REPLACE INTO enumeration_summary ("
            << "id,"
            << "prime,"
            << "genus,"
            << "max_sparsity,"
            << "enumeration_mode,"
            << "irreducible_memory_budget_mb,"
            << "limit_count,"
            << "total_coefficient_vectors,"
            << "processed,"
            << "sparse_presentations,"
            << "sparse_isomorphism_classes,"
            << "canonicalized_isomorphism_classes,"
            << "elapsed_seconds,"
            << "status_counts"
            << ") VALUES ("
            << "1,"
            << opts.p << ","
            << opts.genus << ","
            << (opts.max_sparsity < 0 ? "NULL" : std::to_string(opts.max_sparsity)) << ","
            << quote(opts.mode) << ","
            << opts.irreducible_memory_budget_mb << ","
            << (opts.limit == 0 ? "NULL" : std::to_string(opts.limit)) << ","
            << quote(cpp_int_to_string(stats.total_presentations)) << ","
            << stats.processed << ","
            << quote(cpp_int_to_string(stats.sparse_presentations)) << ","
            << stats.sparse << ","
            << stats.canonicalized << ","
            << elapsed << ","
            << quote(status.str()) << ");";
        exec(sql.str());
        mark_write();
        commit(true);
    }
};

bool for_combinations(
    int n,
    int k,
    int start,
    std::vector<int>& current,
    const std::function<bool(const std::vector<int>&)>& callback
);

class Context {
public:
    explicit Context(const Options& options)
        : opts(options),
          D(2 * options.genus + 2),
          binom(binomial_table(D, options.p)),
          pgl2(pgl2_representatives(options.p)) {
        for (const auto& matrix : pgl2) {
            actions.push_back(action_matrix(matrix, D, opts.p, binom));
        }
        load_irreducibles();
    }

    Options opts;
    int D;
    std::vector<std::vector<unsigned long>> binom;
    std::vector<Matrix4> pgl2;
    std::vector<std::vector<std::vector<unsigned long>>> actions;
    std::unordered_map<int, std::vector<std::vector<std::vector<unsigned long>>>> factor_actions;
    std::vector<std::vector<Poly>> irreducibles;
    std::vector<std::vector<FactorId>> irreducible_ids;
    std::vector<std::uint64_t> irreducible_counts;
    std::vector<bool> irreducible_materialized;
    std::vector<Poly> factors_by_id;
    std::unordered_map<std::string, FactorId> factor_id_by_key;
    std::vector<std::vector<FactorTransform>> factor_transform_cache;
    std::unordered_set<std::string> seen_canonical;
    std::unordered_set<BranchKey, VectorHash> exact_branch_result;
    std::unordered_set<BranchKey, VectorHash> seen_branch_canonical;
    std::unordered_map<BranchKey, BranchKey, VectorHash> branch_orbit_cache;
    std::map<int, FiniteExtensionInt> extensions;
    std::unordered_map<std::string, std::vector<int8_t>> factor_character_cache;

    bool is_square(unsigned long value) const {
        value %= opts.p;
        if (value == 0) return false;
        return mod_pow(value, (opts.p - 1) / 2, opts.p) == 1;
    }

    unsigned long leading_square_class(unsigned long value) const {
        return is_square(value) ? 1 : smallest_nonsquare(opts.p);
    }

    FactorId register_factor(const Poly& factor, bool cache_transforms = false) {
        std::string key = key_poly(factor);
        auto it = factor_id_by_key.find(key);
        if (it != factor_id_by_key.end()) return it->second;
        FactorId id = static_cast<FactorId>(factors_by_id.size());
        factor_id_by_key.emplace(std::move(key), id);
        factors_by_id.push_back(factor);
        factor_transform_cache.emplace_back(cache_transforms ? pgl2.size() : 0);
        return id;
    }

    FactorId factor_id(const Poly& factor) {
        return register_factor(factor);
    }

    std::vector<FactorId> candidate_factor_ids(const BranchCandidate& candidate) {
        if (candidate.factor_ids.size() == candidate.factors.size()) {
            return candidate.factor_ids;
        }
        std::vector<FactorId> ids;
        ids.reserve(candidate.factors.size());
        for (const auto& factor : candidate.factors) ids.push_back(factor_id(factor));
        return ids;
    }

    const std::vector<std::vector<std::vector<unsigned long>>>& action_matrices_for_degree(int degree) {
        auto it = factor_actions.find(degree);
        if (it != factor_actions.end()) return it->second;
        auto degree_binom = binomial_table(degree, opts.p);
        std::vector<std::vector<std::vector<unsigned long>>> matrices;
        matrices.reserve(pgl2.size());
        for (const auto& matrix : pgl2) {
            matrices.push_back(action_matrix(matrix, degree, opts.p, degree_binom));
        }
        auto inserted = factor_actions.emplace(degree, std::move(matrices));
        return inserted.first->second;
    }

    FactorTransform transform_factor(FactorId id, std::size_t matrix_index) {
        std::size_t row_index = static_cast<std::size_t>(id);
        bool cache_enabled = !factor_transform_cache[row_index].empty();
        if (cache_enabled && factor_transform_cache[row_index][matrix_index].computed) {
            return factor_transform_cache[row_index][matrix_index];
        }

        Poly factor = factors_by_id[static_cast<std::size_t>(id)];
        int degree = static_cast<int>(factor.size()) - 1;
        auto transformed = apply_action(factor, action_matrices_for_degree(degree)[matrix_index], opts.p);
        trim(transformed);
        FactorTransform result;
        result.computed = true;
        if (transformed.empty() || (transformed.size() == 1 && transformed[0] == 0)) {
            if (cache_enabled) factor_transform_cache[row_index][matrix_index] = result;
            return result;
        }
        if (transformed.size() == 1) {
            if (cache_enabled) factor_transform_cache[row_index][matrix_index] = result;
            return result;
        }
        result.scalar = transformed.back();
        Poly monic = monic_poly(std::move(transformed), opts.p);
        if (monic.size() <= 1) {
            if (cache_enabled) factor_transform_cache[row_index][matrix_index] = result;
            return result;
        }
        result.factor_id = register_factor(monic);
        result.valid = true;
        if (cache_enabled) factor_transform_cache[row_index][matrix_index] = result;
        return result;
    }

    BranchKey branch_key_from_ids(
        unsigned long leading,
        const std::vector<FactorId>& ids,
        bool infinity_branch
    ) const {
        std::vector<FactorId> sorted_ids = ids;
        std::sort(sorted_ids.begin(), sorted_ids.end());
        BranchKey key;
        key.reserve(sorted_ids.size() + 2);
        key.push_back(infinity_branch ? 1U : 0U);
        key.push_back(static_cast<std::uint32_t>(leading_square_class(leading)));
        for (FactorId id : sorted_ids) key.push_back(id);
        return key;
    }

    BranchKey branch_key(const BranchCandidate& candidate) {
        return branch_key_from_ids(candidate.leading, candidate_factor_ids(candidate), candidate.infinity_branch);
    }

    bool is_enumerated_branch_key(const BranchKey& key) const {
        if (key.size() < 2) return false;
        bool infinity_branch = key[0] != 0;
        if (infinity_branch) {
            return key[1] == 1;
        }
        for (std::size_t i = 2; i < key.size(); ++i) {
            const Poly& factor = factors_by_id[static_cast<std::size_t>(key[i])];
            if (factor.size() == 2) return false;
        }
        return true;
    }

    std::optional<BranchKey> transformed_branch_key(
        const BranchCandidate& candidate,
        const std::vector<FactorId>& ids,
        std::vector<FactorId>& transformed_ids,
        std::size_t matrix_index
    ) {
        unsigned long leading = candidate.leading % opts.p;
        bool infinity_branch = false;
        transformed_ids.clear();
        transformed_ids.reserve(ids.size() + (candidate.infinity_branch ? 1 : 0));
        int finite_degree = 0;
        auto [a, b, c, d] = pgl2[matrix_index];

        for (FactorId id : ids) {
            const auto& transform = transform_factor(id, matrix_index);
            if (!transform.valid) {
                const Poly& factor = factors_by_id[static_cast<std::size_t>(id)];
                if (factor.size() == 2) {
                    unsigned long root = (opts.p - factor[0] % opts.p) % opts.p;
                    unsigned long denominator = (c * root + d) % opts.p;
                    if (denominator == 0) {
                        infinity_branch = true;
                        continue;
                    }
                }
                return std::nullopt;
            }
            leading = (leading * transform.scalar) % opts.p;
            transformed_ids.push_back(transform.factor_id);
            finite_degree += static_cast<int>(factors_by_id[static_cast<std::size_t>(transform.factor_id)].size()) - 1;
        }

        if (candidate.infinity_branch) {
            if (c == 0) {
                infinity_branch = true;
                leading = (leading * d) % opts.p;
            } else {
                Poly finite_factor{(d * mod_inv(c, opts.p)) % opts.p, 1};
                FactorId id = register_factor(finite_factor);
                transformed_ids.push_back(id);
                finite_degree += 1;
                leading = (leading * c) % opts.p;
            }
        }

        if (infinity_branch) {
            if (finite_degree != 2 * opts.genus + 1) return std::nullopt;
            if (!is_square(leading)) return std::nullopt;
            leading = 1;
        } else {
            if (finite_degree != 2 * opts.genus + 2) return std::nullopt;
            leading = leading_square_class(leading);
        }

        return branch_key_from_ids(leading, transformed_ids, infinity_branch);
    }

    std::optional<BranchKey> transformed_branch_key(
        const BranchCandidate& candidate,
        std::size_t matrix_index
    ) {
        std::vector<FactorId> ids = candidate_factor_ids(candidate);
        std::vector<FactorId> transformed_ids;
        return transformed_branch_key(candidate, ids, transformed_ids, matrix_index);
    }

    std::uint64_t normalized_branch_orbit_size(
        const BranchCandidate& candidate,
        const BranchKey& candidate_key
    ) {
        std::vector<FactorId> ids = candidate_factor_ids(candidate);
        std::vector<FactorId> transformed_ids;
        transformed_ids.reserve(ids.size() + (candidate.infinity_branch ? 1 : 0));
        std::unordered_set<BranchKey, VectorHash> orbit;
        for (std::size_t matrix_index = 0; matrix_index < pgl2.size(); ++matrix_index) {
            auto key = transformed_branch_key(candidate, ids, transformed_ids, matrix_index);
            if (key && is_enumerated_branch_key(*key)) orbit.insert(std::move(*key));
        }
        if (orbit.empty()) orbit.insert(candidate_key);
        return static_cast<std::uint64_t>(orbit.size());
    }

    std::optional<BranchCanonicalInfo> factorized_canonical_branch_info(
        const BranchCandidate& candidate,
        const BranchKey& candidate_key,
        bool store_full_orbit
    ) {
        auto cached = branch_orbit_cache.find(candidate_key);
        if (cached != branch_orbit_cache.end()) return BranchCanonicalInfo{cached->second, 1};

        std::vector<FactorId> ids = candidate_factor_ids(candidate);
        std::vector<FactorId> transformed_ids;
        transformed_ids.reserve(ids.size() + (candidate.infinity_branch ? 1 : 0));
        std::optional<BranchKey> best;
        std::unordered_set<BranchKey, VectorHash> orbit;
        std::unordered_set<BranchKey, VectorHash> enumerated_orbit;
        orbit.reserve(pgl2.size());
        enumerated_orbit.reserve(pgl2.size());
        for (std::size_t matrix_index = 0; matrix_index < pgl2.size(); ++matrix_index) {
            auto key = transformed_branch_key(candidate, ids, transformed_ids, matrix_index);
            if (!key) return std::nullopt;
            if (!best || *key < *best) best = *key;
            if (is_enumerated_branch_key(*key)) enumerated_orbit.insert(*key);
            orbit.insert(std::move(*key));
        }
        if (!best) return std::nullopt;
        if (store_full_orbit) {
            for (const auto& key : orbit) branch_orbit_cache.emplace(key, *best);
        } else {
            branch_orbit_cache.emplace(candidate_key, *best);
            branch_orbit_cache.emplace(*best, *best);
        }
        if (enumerated_orbit.empty()) enumerated_orbit.insert(candidate_key);
        return BranchCanonicalInfo{*best, static_cast<std::uint64_t>(enumerated_orbit.size())};
    }

    std::uint64_t memory_budget_bytes() const {
        constexpr std::uint64_t mib = 1024ULL * 1024ULL;
        return saturating_mul(opts.irreducible_memory_budget_mb, mib);
    }

    std::uint64_t estimated_irreducible_degree_bytes(int degree, std::uint64_t count) const {
        std::uint64_t coefficient_bytes = saturating_mul(static_cast<std::uint64_t>(degree + 1), sizeof(unsigned long));
        std::uint64_t poly_copy_bytes = saturating_add(sizeof(Poly), coefficient_bytes);
        std::uint64_t transform_bytes = saturating_add(
            sizeof(std::vector<FactorTransform>),
            saturating_mul(static_cast<std::uint64_t>(pgl2.size()), sizeof(FactorTransform))
        );
        std::uint64_t map_bytes = 96;
        std::uint64_t per_factor = 0;
        per_factor = saturating_add(per_factor, saturating_mul(2, poly_copy_bytes));
        per_factor = saturating_add(per_factor, sizeof(FactorId));
        per_factor = saturating_add(per_factor, transform_bytes);
        per_factor = saturating_add(per_factor, map_bytes);
        return saturating_mul(count, per_factor);
    }

    std::uint64_t monic_polynomial_count(int degree) const {
        std::uint64_t total = 1;
        for (int i = 0; i < degree; ++i) {
            if (total > std::numeric_limits<std::uint64_t>::max() / opts.p) {
                throw std::runtime_error("monic polynomial count overflow");
            }
            total *= opts.p;
        }
        return total;
    }

    Poly monic_polynomial_from_index(int degree, std::uint64_t index) const {
        Poly poly(static_cast<std::size_t>(degree + 1), 0);
        for (int i = 0; i < degree; ++i) {
            poly[static_cast<std::size_t>(i)] = index % opts.p;
            index /= opts.p;
        }
        poly[static_cast<std::size_t>(degree)] = 1;
        return poly;
    }

    std::uint64_t irreducible_count(int degree) const {
        return irreducible_counts[static_cast<std::size_t>(degree)];
    }

    void load_irreducibles() {
        irreducibles.assign(static_cast<std::size_t>(D + 1), {});
        irreducible_ids.assign(static_cast<std::size_t>(D + 1), {});
        irreducible_counts.assign(static_cast<std::size_t>(D + 1), 0);
        irreducible_materialized.assign(static_cast<std::size_t>(D + 1), false);
        std::uint64_t budget = memory_budget_bytes();
        std::uint64_t used = 0;
        std::cout << "irreducible_load: 0/" << D
                  << " budget_mb=" << opts.irreducible_memory_budget_mb << std::endl;
        for (int degree = 1; degree <= D; ++degree) {
            if (stop_requested()) break;
            std::uint64_t count = irreducible_count_formula(opts.p, degree);
            irreducible_counts[static_cast<std::size_t>(degree)] = count;
            std::uint64_t estimated_bytes = estimated_irreducible_degree_bytes(degree, count);
            auto started = std::chrono::steady_clock::now();
            bool materialize = estimated_bytes <= budget - std::min(used, budget);
            if (materialize) {
                irreducibles[static_cast<std::size_t>(degree)].reserve(static_cast<std::size_t>(count));
                irreducible_ids[static_cast<std::size_t>(degree)].reserve(static_cast<std::size_t>(count));
                std::uint64_t total = monic_polynomial_count(degree);
                for (std::uint64_t index = 0; index < total; ++index) {
                    if (stop_requested()) break;
                    Poly poly = monic_polynomial_from_index(degree, index);
                    if (is_irreducible_flint(poly, opts.p)) {
                        FactorId id = register_factor(poly, true);
                        irreducibles[static_cast<std::size_t>(degree)].push_back(std::move(poly));
                        irreducible_ids[static_cast<std::size_t>(degree)].push_back(id);
                    }
                }
                if (stop_requested()) {
                    irreducibles[static_cast<std::size_t>(degree)].clear();
                    irreducible_ids[static_cast<std::size_t>(degree)].clear();
                    materialize = false;
                } else {
                    irreducible_materialized[static_cast<std::size_t>(degree)] = true;
                    used = saturating_add(used, estimated_bytes);
                }
            }
            double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
            std::cout << "irreducible_load: " << degree << "/" << D
                      << " degree=" << degree
                      << " mode=" << (materialize ? "memory" : "stream")
                      << " source=flint"
                      << " count=" << count
                      << " estimated_bytes=" << estimated_bytes
                      << " used_budget_bytes=" << used
                      << " seconds=" << seconds << std::endl;
        }
    }

    bool stream_irreducible_combinations_rec(
        int degree,
        int multiplicity,
        std::uint64_t start_index,
        std::vector<Poly>& selected,
        std::vector<FactorId>& selected_ids,
        const std::function<bool(const std::vector<Poly>&, const std::vector<FactorId>&)>& callback
    ) {
        if (static_cast<int>(selected.size()) == multiplicity) {
            return callback(selected, selected_ids);
        }
        std::uint64_t total = monic_polynomial_count(degree);
        for (std::uint64_t index = start_index; index < total; ++index) {
            Poly poly = monic_polynomial_from_index(degree, index);
            if (!is_irreducible_flint(poly, opts.p)) continue;
            FactorId id = register_factor(poly, false);
            selected.push_back(std::move(poly));
            selected_ids.push_back(id);
            if (!stream_irreducible_combinations_rec(degree, multiplicity, index + 1, selected, selected_ids, callback)) {
                return false;
            }
            selected.pop_back();
            selected_ids.pop_back();
        }
        return true;
    }

    bool for_irreducible_combinations(
        int degree,
        int multiplicity,
        const std::function<bool(const std::vector<Poly>&, const std::vector<FactorId>&)>& callback
    ) {
        if (multiplicity > static_cast<int>(irreducible_count(degree))) return true;
        if (irreducible_materialized[static_cast<std::size_t>(degree)]) {
            const auto& pool = irreducibles[static_cast<std::size_t>(degree)];
            const auto& pool_ids = irreducible_ids[static_cast<std::size_t>(degree)];
            std::vector<int> combo;
            bool keep_going = true;
            keep_going = for_combinations(static_cast<int>(pool.size()), multiplicity, 0, combo, [&](const std::vector<int>& indices) {
                if (!keep_going || stop_requested()) return false;
                std::vector<Poly> factors;
                std::vector<FactorId> ids;
                factors.reserve(indices.size());
                ids.reserve(indices.size());
                for (int idx : indices) {
                    factors.push_back(pool[static_cast<std::size_t>(idx)]);
                    ids.push_back(pool_ids[static_cast<std::size_t>(idx)]);
                }
                keep_going = callback(factors, ids);
                return keep_going && !stop_requested();
            });
            return keep_going;
        }

        std::vector<Poly> selected;
        std::vector<FactorId> selected_ids;
        selected.reserve(static_cast<std::size_t>(multiplicity));
        selected_ids.reserve(static_cast<std::size_t>(multiplicity));
        return stream_irreducible_combinations_rec(degree, multiplicity, 0, selected, selected_ids, callback);
    }

    std::pair<Poly, FactorId> random_irreducible_factor(int degree, std::mt19937_64& rng) {
        if (irreducible_materialized[static_cast<std::size_t>(degree)]) {
            const auto& pool = irreducibles[static_cast<std::size_t>(degree)];
            const auto& pool_ids = irreducible_ids[static_cast<std::size_t>(degree)];
            std::uniform_int_distribution<std::size_t> dist(0, pool.size() - 1);
            std::size_t index = dist(rng);
            return {pool[index], pool_ids[index]};
        }
        std::uniform_int_distribution<unsigned long> coeff_dist(0, opts.p - 1);
        while (true) {
            Poly poly(static_cast<std::size_t>(degree + 1), 0);
            for (int i = 0; i < degree; ++i) {
                poly[static_cast<std::size_t>(i)] = coeff_dist(rng);
            }
            poly[static_cast<std::size_t>(degree)] = 1;
            if (is_irreducible_flint(poly, opts.p)) {
                FactorId id = register_factor(poly, false);
                return {std::move(poly), id};
            }
        }
    }

    FiniteExtensionInt& extension(int degree) {
        auto it = extensions.find(degree);
        if (it == extensions.end()) {
            it = extensions.emplace(degree, FiniteExtensionInt(opts.p, degree)).first;
        }
        return it->second;
    }

    std::vector<int> lpoly_mod_p(const BranchCandidate& candidate) {
        unsigned long exponent = (opts.p - 1) / 2;
        std::vector<Poly> powers;
        for (const auto& factor : candidate.factors) {
            powers.push_back(pow_flint(factor, exponent, opts.p));
        }
        Poly h = product_tree_flint(std::move(powers), opts.p);
        unsigned long scalar = mod_pow(candidate.leading, exponent, opts.p);
        for (auto& coefficient : h) coefficient = (coefficient * scalar) % opts.p;

        nmod_mat_t H;
        nmod_poly_t charpoly;
        nmod_mat_init(H, opts.genus, opts.genus, opts.p);
        nmod_poly_init(charpoly, opts.p);
        for (int i = 1; i <= opts.genus; ++i) {
            for (int j = 1; j <= opts.genus; ++j) {
                int idx = static_cast<int>(opts.p) * i - j;
                unsigned long value = idx >= 0 && idx < static_cast<int>(h.size()) ? h[static_cast<std::size_t>(idx)] : 0;
                nmod_mat_entry(H, i - 1, j - 1) = value;
            }
        }
        nmod_mat_charpoly_danilevsky(charpoly, H);
        std::vector<int> out;
        for (int i = 1; i <= opts.genus; ++i) {
            out.push_back(static_cast<int>(nmod_poly_get_coeff_ui(charpoly, opts.genus - i)));
        }
        nmod_poly_clear(charpoly);
        nmod_mat_clear(H);
        return out;
    }

    std::string canonical_key(const Poly& affine) {
        std::vector<unsigned long> binary(static_cast<std::size_t>(D + 1), 0);
        for (std::size_t i = 0; i < affine.size() && i < binary.size(); ++i) {
            binary[i] = affine[i] % opts.p;
        }
        std::optional<std::vector<unsigned long>> best;
        for (const auto& action : actions) {
            auto transformed = normalize_square_scalar(apply_action(binary, action, opts.p), opts.p);
            if (!best || transformed < *best) best = std::move(transformed);
        }
        if (!best) throw std::runtime_error("failed to canonicalize");
        return key_int_vector(*best);
    }

    int rational_branch_count(const BranchCandidate& candidate) const {
        int count = candidate.infinity_branch ? 1 : 0;
        for (const auto& factor : candidate.factors) {
            if (factor.size() == 2) ++count;
        }
        return count;
    }

    int leading_character(const BranchCandidate& candidate, int extension_degree) {
        if (extension_degree == 1) {
            if (candidate.leading % opts.p == 0) return 0;
            return mod_pow(candidate.leading, (opts.p - 1) / 2, opts.p) == 1 ? 1 : -1;
        }
        auto& ext = extension(extension_degree);
        return ext.character(ext.constant(candidate.leading));
    }

    const std::vector<int8_t>& factor_characters(const Poly& factor, int extension_degree) {
        std::string cache_key = std::to_string(extension_degree) + "|" + key_poly(factor);
        auto it = factor_character_cache.find(cache_key);
        if (it != factor_character_cache.end()) return it->second;

        auto& ext = extension(extension_degree);
        std::vector<int8_t> chars(static_cast<std::size_t>(ext.size));
        for (std::uint64_t x = 0; x < ext.size; ++x) {
            chars[static_cast<std::size_t>(x)] = static_cast<int8_t>(ext.character(ext.evaluate(factor, x)));
        }
        auto inserted = factor_character_cache.emplace(std::move(cache_key), std::move(chars));
        return inserted.first->second;
    }

    cpp_int point_count(const BranchCandidate& candidate, int extension_degree) {
        auto& ext = extension(extension_degree);
        std::vector<int8_t> product(static_cast<std::size_t>(ext.size), 1);
        for (const auto& factor : candidate.factors) {
            const auto& chars = factor_characters(factor, extension_degree);
            for (std::size_t i = 0; i < product.size(); ++i) {
                product[i] = static_cast<int8_t>(product[i] * chars[i]);
            }
        }
        long long sum = 0;
        for (auto value : product) sum += value;
        int lead_chi = leading_character(candidate, extension_degree);
        int infinity = candidate.infinity_branch ? 1 : (lead_chi == 1 ? 2 : 0);
        cpp_int count = infinity;
        count += ext.size;
        count += lead_chi * sum;
        return count;
    }

    std::optional<std::vector<cpp_int>> exact_lpoly(const BranchCandidate& candidate) {
        std::vector<cpp_int> power_sums(static_cast<std::size_t>(opts.genus + 1));
        std::vector<cpp_int> coefficients(static_cast<std::size_t>(opts.genus + 1));
        coefficients[0] = 1;
        cpp_int q = 1;
        int sparsity = 0;
        for (int k = 1; k <= opts.genus; ++k) {
            q *= opts.p;
            power_sums[static_cast<std::size_t>(k)] = q + 1 - point_count(candidate, k);
            cpp_int total = 0;
            for (int i = 1; i <= k; ++i) {
                total += coefficients[static_cast<std::size_t>(k - i)] * power_sums[static_cast<std::size_t>(i)];
            }
            if (total % k != 0) {
                throw std::runtime_error("Newton identity produced a nonintegral coefficient");
            }
            coefficients[static_cast<std::size_t>(k)] = -total / k;
            if (k < opts.genus && coefficients[static_cast<std::size_t>(k)] != 0) {
                ++sparsity;
                if (opts.max_sparsity >= 0 && sparsity > opts.max_sparsity) {
                    return std::nullopt;
                }
            }
        }
        std::vector<cpp_int> out;
        for (int k = 1; k <= opts.genus; ++k) out.push_back(coefficients[static_cast<std::size_t>(k)]);
        return out;
    }
};

void generate_patterns_rec(
    int remaining,
    int smallest,
    int minimum,
    std::vector<int>& counts,
    std::vector<std::vector<int>>& out
) {
    if (remaining == 0) {
        out.push_back(counts);
        return;
    }
    for (int degree = std::max(smallest, minimum); degree <= remaining; ++degree) {
        counts[static_cast<std::size_t>(degree)]++;
        generate_patterns_rec(remaining - degree, degree, minimum, counts, out);
        counts[static_cast<std::size_t>(degree)]--;
    }
}

std::vector<std::vector<int>> factorization_patterns(int total_degree, bool skip_linear) {
    std::vector<std::vector<int>> out;
    std::vector<int> counts(static_cast<std::size_t>(total_degree + 1), 0);
    generate_patterns_rec(total_degree, skip_linear ? 2 : 1, skip_linear ? 2 : 1, counts, out);
    return out;
}

std::vector<std::uint8_t> divide_by_one_plus_t_mod2(const std::vector<std::uint8_t>& poly) {
    if (poly.size() < 2) {
        throw std::runtime_error("mod-2 branch divisor quotient is too small");
    }
    std::vector<std::uint8_t> quotient(poly.size() - 1, 0);
    quotient[0] = poly[0] & 1U;
    for (std::size_t i = 1; i < quotient.size(); ++i) {
        quotient[i] = static_cast<std::uint8_t>((poly[i] ^ quotient[i - 1]) & 1U);
    }
    if (((poly.back() ^ quotient.back()) & 1U) != 0) {
        throw std::runtime_error("branch divisor type did not produce a divisible mod-2 L-polynomial");
    }
    return quotient;
}

int branch_type_forced_mod2_sparsity(bool infinity_branch, const std::vector<int>& pattern, int genus) {
    int total_degree = infinity_branch ? 1 : 0;
    for (int degree = 1; degree < static_cast<int>(pattern.size()); ++degree) {
        total_degree += degree * pattern[static_cast<std::size_t>(degree)];
    }

    std::vector<std::uint8_t> product(static_cast<std::size_t>(total_degree + 1), 0);
    product[0] = 1;
    int current_degree = 0;
    auto multiply_by_orbit = [&](int degree) {
        if (degree <= 0) return;
        for (int i = current_degree; i >= 0; --i) {
            product[static_cast<std::size_t>(i + degree)] ^= product[static_cast<std::size_t>(i)];
        }
        current_degree += degree;
    };

    if (infinity_branch) multiply_by_orbit(1);
    for (int degree = 1; degree < static_cast<int>(pattern.size()); ++degree) {
        int multiplicity = pattern[static_cast<std::size_t>(degree)];
        for (int i = 0; i < multiplicity; ++i) multiply_by_orbit(degree);
    }

    std::vector<std::uint8_t> quotient = divide_by_one_plus_t_mod2(divide_by_one_plus_t_mod2(product));
    int sparsity = 0;
    int max_degree = std::min(genus - 1, static_cast<int>(quotient.size()) - 1);
    for (int degree = 1; degree <= max_degree; ++degree) {
        if (quotient[static_cast<std::size_t>(degree)] != 0) ++sparsity;
    }
    return sparsity;
}

cpp_int binomial_count(std::uint64_t n, int k) {
    if (k < 0 || static_cast<std::uint64_t>(k) > n) return 0;
    std::uint64_t kk = static_cast<std::uint64_t>(k);
    kk = std::min(kk, n - kk);
    cpp_int result = 1;
    for (std::uint64_t i = 1; i <= kk; ++i) {
        result *= n - kk + i;
        result /= i;
    }
    return result;
}

cpp_int branch_pattern_presentation_count(
    const Context& ctx,
    const std::vector<int>& pattern,
    std::uint64_t leading_count
) {
    cpp_int total = leading_count;
    for (int degree = 1; degree < static_cast<int>(pattern.size()); ++degree) {
        int multiplicity = pattern[static_cast<std::size_t>(degree)];
        if (multiplicity == 0) continue;
        total *= binomial_count(ctx.irreducible_count(degree), multiplicity);
    }
    return total;
}

int pattern_factor_count(const std::vector<int>& pattern) {
    int total = 0;
    for (int multiplicity : pattern) total += multiplicity;
    return total;
}

std::vector<BranchDivisorTypeState> build_branch_divisor_type_states(const Context& ctx, bool apply_random_factor_cap) {
    if (apply_random_factor_cap && ctx.opts.random_max_factors < 0) {
        throw std::runtime_error("--random-max-factors must be nonnegative");
    }

    std::vector<BranchDivisorTypeState> states;
    for (int model = 0; model < 2; ++model) {
        int total_degree = model == 0 ? 2 * ctx.opts.genus + 1 : 2 * ctx.opts.genus + 2;
        bool skip_linear = model == 1;
        bool infinity_branch = model == 0;
        std::uint64_t leading_count = model == 0 ? 1 : 2;
        for (const auto& pattern : factorization_patterns(total_degree, skip_linear)) {
            int factors = pattern_factor_count(pattern);
            if (factors < 1) continue;
            if (apply_random_factor_cap && ctx.opts.random_max_factors > 0 && factors > ctx.opts.random_max_factors) continue;
            int forced_mod2_sparsity = branch_type_forced_mod2_sparsity(infinity_branch, pattern, ctx.opts.genus);
            if (ctx.opts.max_sparsity >= 0 && forced_mod2_sparsity > ctx.opts.max_sparsity) continue;
            cpp_int presentations = branch_pattern_presentation_count(ctx, pattern, leading_count);
            if (presentations == 0) continue;
            BranchDivisorTypeState state;
            state.model = model;
            state.total_degree = total_degree;
            state.infinity_branch = infinity_branch;
            state.pattern = pattern;
            state.factor_count = factors;
            state.forced_mod2_sparsity = forced_mod2_sparsity;
            state.leading_count = leading_count;
            state.presentations = std::move(presentations);
            states.push_back(std::move(state));
        }
    }
    if (states.empty()) {
        throw std::runtime_error("no feasible branch divisor types after feasibility and mod-2 filters");
    }
    return states;
}

cpp_int total_branch_divisor_presentations(const std::vector<BranchDivisorTypeState>& states) {
    cpp_int total = 0;
    for (const auto& state : states) total += state.presentations;
    return total;
}

std::string random_pattern_label(const BranchDivisorTypeState& state) {
    return json_branch_divisor_type(state.infinity_branch, state.pattern);
}

std::string random_pattern_progress_summary(const std::vector<BranchDivisorTypeState>& states) {
    const BranchDivisorTypeState* best = nullptr;
    for (const auto& state : states) {
        if (state.attempts == 0) continue;
        if (!best) {
            best = &state;
            continue;
        }
        double lhs_sparse = static_cast<double>(state.sparse_hits) / static_cast<double>(state.attempts);
        double rhs_sparse = static_cast<double>(best->sparse_hits) / static_cast<double>(best->attempts);
        double lhs_hw = static_cast<double>(state.hasse_witt_passes) / static_cast<double>(state.attempts);
        double rhs_hw = static_cast<double>(best->hasse_witt_passes) / static_cast<double>(best->attempts);
        if (lhs_sparse > rhs_sparse ||
            (lhs_sparse == rhs_sparse && lhs_hw > rhs_hw) ||
            (lhs_sparse == rhs_sparse && lhs_hw == rhs_hw && state.sparse_hits > best->sparse_hits)) {
            best = &state;
        }
    }
    std::ostringstream out;
    out << states.size();
    if (best) {
        out << " best=" << random_pattern_label(*best)
            << " attempts=" << best->attempts
            << " hw_pass=" << best->hasse_witt_passes
            << " sparse=" << best->sparse_hits;
    }
    return out.str();
}

std::size_t choose_random_pattern(
    const std::vector<BranchDivisorTypeState>& states,
    std::uint64_t total_attempts,
    std::mt19937_64& rng
) {
    std::vector<double> weights;
    weights.reserve(states.size());
    double log_scale = std::log(static_cast<double>(total_attempts + states.size() + 2));
    for (const auto& state : states) {
        double attempts = static_cast<double>(state.attempts);
        double sparse_hits = static_cast<double>(state.sparse_hits);
        double hasse_witt_pass_rate = state.attempts == 0
            ? 0.0
            : static_cast<double>(state.hasse_witt_passes) / attempts;
        double factor_prior = 1.0 / static_cast<double>(
            state.factor_count * state.factor_count * state.factor_count
        );
        double mod2_prior = 1.0 / static_cast<double>(
            (state.forced_mod2_sparsity + 1) * (state.forced_mod2_sparsity + 1)
        );
        double smoothed_hit_rate = (sparse_hits + 0.25) / (attempts + 8.0);
        double exploration = std::sqrt(log_scale / (attempts + 1.0));
        double hit_boost = state.sparse_hits == 0
            ? 1.0
            : 1.0 + std::min(20.0, 4.0 * std::sqrt(sparse_hits));
        double hasse_witt_boost = state.sparse_hits == 0
            ? 1.0 + std::min(10.0, 0.25 * std::sqrt(static_cast<double>(state.hasse_witt_passes)))
            : 1.0;
        weights.push_back(
            factor_prior *
            mod2_prior *
            (0.02 + smoothed_hit_rate + 30.0 * hasse_witt_pass_rate + 0.20 * exploration) *
            hit_boost *
            hasse_witt_boost
        );
    }
    std::discrete_distribution<std::size_t> dist(weights.begin(), weights.end());
    return dist(rng);
}

int sparsity_without_last(const std::vector<int>& coeffs) {
    int sparsity = 0;
    for (std::size_t i = 0; i + 1 < coeffs.size(); ++i) {
        if (coeffs[i] != 0) ++sparsity;
    }
    return sparsity;
}

int sparsity_without_last(const std::vector<cpp_int>& coeffs) {
    int sparsity = 0;
    for (std::size_t i = 0; i + 1 < coeffs.size(); ++i) {
        if (coeffs[i] != 0) ++sparsity;
    }
    return sparsity;
}

CandidateOutcome process_candidate(Context& ctx, SqliteWriter& writer, const BranchCandidate& candidate, Stats& stats) {
    std::vector<int> modp = ctx.lpoly_mod_p(candidate);
    if (ctx.opts.max_sparsity >= 0 && sparsity_without_last(modp) > ctx.opts.max_sparsity) {
        ++stats.rejected_hasse_witt;
        return CandidateOutcome::HasseWittRejected;
    }

    BranchKey branch_key = ctx.branch_key(candidate);
    if (ctx.exact_branch_result.count(branch_key)) {
        ++stats.duplicate;
        return CandidateOutcome::Duplicate;
    }

    if (ctx.branch_orbit_cache.find(branch_key) != ctx.branch_orbit_cache.end()) {
        ++stats.duplicate;
        ctx.exact_branch_result.insert(std::move(branch_key));
        return CandidateOutcome::Duplicate;
    }

    std::optional<std::vector<cpp_int>> exact;
    bool sparse_first = ctx.opts.mode == "random";
    if (sparse_first) {
        exact = ctx.exact_lpoly(candidate);
        if (!exact) {
            ++stats.rejected_exact;
            ctx.exact_branch_result.insert(std::move(branch_key));
            return CandidateOutcome::ExactRejected;
        }
    }

    std::uint64_t orbit_size = 1;
    std::optional<BranchCanonicalInfo> canonical_branch =
        ctx.factorized_canonical_branch_info(candidate, branch_key, ctx.opts.mode == "enumerate");
    if (canonical_branch) {
        orbit_size = canonical_branch->orbit_size;
        auto inserted_branch = ctx.seen_branch_canonical.insert(canonical_branch->canonical_key);
        if (!inserted_branch.second) {
            ++stats.duplicate;
            ctx.exact_branch_result.insert(std::move(branch_key));
            return CandidateOutcome::Duplicate;
        }
    } else {
        orbit_size = ctx.normalized_branch_orbit_size(candidate, branch_key);
    }

    Poly affine = expand_candidate(candidate, ctx.opts.p);
    std::string canonical = ctx.canonical_key(affine);
    if (ctx.seen_canonical.count(canonical)) {
        ++stats.duplicate;
        ctx.exact_branch_result.insert(std::move(branch_key));
        return CandidateOutcome::Duplicate;
    }
    ctx.seen_canonical.insert(canonical);
    ++stats.canonicalized;

    if (!exact) {
        exact = ctx.exact_lpoly(candidate);
        if (!exact) {
            ++stats.rejected_exact;
            ctx.exact_branch_result.insert(std::move(branch_key));
            return CandidateOutcome::ExactRejected;
        }
    }

    int sparsity = sparsity_without_last(*exact);
    ++stats.sparse;
    stats.sparse_presentations += orbit_size;
    writer.insert_sparse(canonical, affine, candidate, *exact, sparsity, ctx.rational_branch_count(candidate));
    ctx.exact_branch_result.insert(std::move(branch_key));
    return CandidateOutcome::Sparse;
}

std::string run_status(const Stats& stats) {
    if (stop_requested()) return "interrupted";
    if (stats.total_presentations >= 0 && cpp_int(stats.processed) >= stats.total_presentations) {
        return "complete";
    }
    return "partial";
}

std::string format_duration(double seconds) {
    if (!std::isfinite(seconds) || seconds < 0) return "?";
    std::uint64_t total = static_cast<std::uint64_t>(seconds + 0.5);
    std::uint64_t days = total / 86400;
    total %= 86400;
    std::uint64_t hours = total / 3600;
    total %= 3600;
    std::uint64_t minutes = total / 60;
    std::uint64_t secs = total % 60;
    std::ostringstream out;
    out << days << "d "
        << std::setw(2) << std::setfill('0') << hours << "h "
        << std::setw(2) << std::setfill('0') << minutes << "m "
        << std::setw(2) << std::setfill('0') << secs << "s";
    return out.str();
}

bool finite_positive_total(const cpp_int& total) {
    return total > 0;
}

std::string progress_percent_label(std::uint64_t processed, const cpp_int& total) {
    if (!finite_positive_total(total)) return "?";
    long double total_value = total.convert_to<long double>();
    if (!std::isfinite(static_cast<double>(total_value)) || total_value <= 0) return "?";
    long double percent = 100.0L * static_cast<long double>(processed) / total_value;
    if (percent > 100.0L) percent = 100.0L;
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << static_cast<double>(percent) << "%";
    return out.str();
}

std::string estimated_remaining_label(std::uint64_t processed, const cpp_int& total, double elapsed) {
    if (processed == 0 || !finite_positive_total(total) || !std::isfinite(elapsed) || elapsed < 0) return "?";
    long double total_value = total.convert_to<long double>();
    if (!std::isfinite(static_cast<double>(total_value)) || total_value <= 0) return "?";
    long double remaining = (total_value / static_cast<long double>(processed) - 1.0L) * static_cast<long double>(elapsed);
    if (remaining < 0) remaining = 0;
    return format_duration(static_cast<double>(remaining));
}

void print_progress_block(
    const Context& ctx,
    const Stats& stats,
    const std::string& total_label,
    double elapsed,
    const std::string& random_pattern_summary = ""
) {
    std::cout << "prime: " << ctx.opts.p << "\n"
              << "genus: " << ctx.opts.genus << "\n"
              << "progress: " << stats.processed << "/" << total_label << "\n"
              << "progress_percent: " << progress_percent_label(stats.processed, stats.total_presentations) << "\n"
              << "elapsed: " << format_duration(elapsed) << "\n"
              << "estimated_remaining: " << estimated_remaining_label(stats.processed, stats.total_presentations, elapsed) << "\n"
              << "sparse_presentations: " << stats.sparse_presentations << "\n"
              << "sparse_isomorphism_classes: " << stats.sparse << "\n"
              << "canonicalized_isomorphism_classes: " << stats.canonicalized << "\n";
    if (!random_pattern_summary.empty()) {
        std::cout << "random_patterns: " << random_pattern_summary << "\n";
    }
    std::cout << "-\n";
}

bool for_combinations(
    int n,
    int k,
    int start,
    std::vector<int>& current,
    const std::function<bool(const std::vector<int>&)>& callback
) {
    if (stop_requested()) return false;
    if (static_cast<int>(current.size()) == k) {
        return callback(current);
    }
    for (int i = start; i <= n - (k - static_cast<int>(current.size())); ++i) {
        current.push_back(i);
        bool keep_going = for_combinations(n, k, i + 1, current, callback);
        current.pop_back();
        if (!keep_going) return false;
    }
    return true;
}

bool enumerate_factor_choices_rec(
    Context& ctx,
    SqliteWriter& writer,
    const std::vector<std::pair<int, int>>& active,
    std::size_t index,
    std::vector<Poly>& selected,
    std::vector<FactorId>& selected_ids,
    const std::vector<unsigned long>& leading_coefficients,
    bool infinity_branch,
    const std::vector<int>& pattern,
    Stats& stats,
    const std::string& total_label,
    const std::chrono::steady_clock::time_point& started
) {
    if (stop_requested()) return false;
    if (ctx.opts.limit && stats.processed >= ctx.opts.limit) return false;
    if (index == active.size()) {
        for (unsigned long leading : leading_coefficients) {
            if (stop_requested()) return false;
            if (ctx.opts.limit && stats.processed >= ctx.opts.limit) return false;
            BranchCandidate candidate{leading, selected, infinity_branch, pattern, selected_ids};
            process_candidate(ctx, writer, candidate, stats);
            ++stats.processed;
            if (ctx.opts.progress_interval > 0 && stats.processed % static_cast<std::uint64_t>(ctx.opts.progress_interval) == 0) {
                double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
                print_progress_block(ctx, stats, total_label, elapsed);
                writer.write_progress(stats, elapsed);
            }
        }
        return true;
    }

    auto [degree, multiplicity] = active[index];
    bool keep_going = true;
    keep_going = ctx.for_irreducible_combinations(degree, multiplicity, [&](const std::vector<Poly>& factors, const std::vector<FactorId>& ids) {
        std::size_t old_size = selected.size();
        std::size_t old_id_size = selected_ids.size();
        selected.insert(selected.end(), factors.begin(), factors.end());
        selected_ids.insert(selected_ids.end(), ids.begin(), ids.end());
        keep_going = enumerate_factor_choices_rec(
            ctx,
            writer,
            active,
            index + 1,
            selected,
            selected_ids,
            leading_coefficients,
            infinity_branch,
            pattern,
            stats,
            total_label,
            started
        );
        selected.resize(old_size);
        selected_ids.resize(old_id_size);
        return keep_going;
    });
    return keep_going;
}

void enumerate_mode(Context& ctx, SqliteWriter& writer, Stats& stats) {
    auto started = std::chrono::steady_clock::now();
    std::vector<BranchDivisorTypeState> types = build_branch_divisor_type_states(ctx, false);
    stats.total_presentations = total_branch_divisor_presentations(types);
    std::string total_label = cpp_int_to_string(stats.total_presentations);
    std::cout << "branch_divisor_types: " << types.size()
              << " mod2_filter=" << (ctx.opts.max_sparsity >= 0 ? "on" : "off")
              << std::endl;
    print_progress_block(ctx, stats, total_label, 0.0);

    for (const auto& type : types) {
        if (stop_requested()) return;
        std::vector<unsigned long> leadings = type.model == 0
            ? std::vector<unsigned long>{1}
            : std::vector<unsigned long>{1, smallest_nonsquare(ctx.opts.p)};
        std::vector<std::pair<int, int>> active;
        for (int degree = 1; degree < static_cast<int>(type.pattern.size()); ++degree) {
            if (type.pattern[static_cast<std::size_t>(degree)] > 0) {
                active.emplace_back(degree, type.pattern[static_cast<std::size_t>(degree)]);
            }
        }
        std::vector<Poly> selected;
        std::vector<FactorId> selected_ids;
        if (!enumerate_factor_choices_rec(
            ctx,
            writer,
            active,
            0,
            selected,
            selected_ids,
            leadings,
            type.infinity_branch,
            type.pattern,
            stats,
            total_label,
            started
        )) {
            return;
        }
    }
}

BranchCandidate random_candidate_from_pattern(
    Context& ctx,
    const BranchDivisorTypeState& state,
    std::mt19937_64& rng
) {
    std::vector<Poly> factors;
    std::vector<FactorId> factor_ids;
    factors.reserve(static_cast<std::size_t>(state.factor_count));
    factor_ids.reserve(static_cast<std::size_t>(state.factor_count));

    for (int degree = 1; degree < static_cast<int>(state.pattern.size()); ++degree) {
        int multiplicity = state.pattern[static_cast<std::size_t>(degree)];
        if (multiplicity == 0) continue;
        if (ctx.irreducible_materialized[static_cast<std::size_t>(degree)]) {
            const auto& pool = ctx.irreducibles[static_cast<std::size_t>(degree)];
            const auto& pool_ids = ctx.irreducible_ids[static_cast<std::size_t>(degree)];
            std::unordered_set<std::size_t> selected_indices;
            std::uniform_int_distribution<std::size_t> dist(0, pool.size() - 1);
            while (static_cast<int>(selected_indices.size()) < multiplicity) {
                selected_indices.insert(dist(rng));
            }
            for (std::size_t index : selected_indices) {
                factors.push_back(pool[index]);
                factor_ids.push_back(pool_ids[index]);
            }
        } else {
            std::unordered_set<std::string> seen;
            while (static_cast<int>(seen.size()) < multiplicity) {
                auto [factor, id] = ctx.random_irreducible_factor(degree, rng);
                std::string key = key_poly(factor);
                if (!seen.insert(key).second) continue;
                factors.push_back(std::move(factor));
                factor_ids.push_back(id);
            }
        }
    }

    std::vector<std::size_t> order(factors.size());
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), rng);
    std::vector<Poly> shuffled_factors;
    std::vector<FactorId> shuffled_factor_ids;
    shuffled_factors.reserve(factors.size());
    shuffled_factor_ids.reserve(factor_ids.size());
    for (std::size_t index : order) {
        shuffled_factors.push_back(std::move(factors[index]));
        shuffled_factor_ids.push_back(factor_ids[index]);
    }

    unsigned long leading = 1;
    if (!state.infinity_branch) {
        std::uniform_int_distribution<int> leading_dist(0, static_cast<int>(state.leading_count) - 1);
        leading = leading_dist(rng) == 0 ? 1 : smallest_nonsquare(ctx.opts.p);
    }
    return BranchCandidate{leading, std::move(shuffled_factors), state.infinity_branch, state.pattern, std::move(shuffled_factor_ids)};
}

void random_mode(Context& ctx, SqliteWriter& writer, Stats& stats) {
    auto started = std::chrono::steady_clock::now();
    std::mt19937_64 rng(ctx.opts.random_seed);
    std::vector<BranchDivisorTypeState> patterns = build_branch_divisor_type_states(ctx, true);
    std::cout << "branch_divisor_types: " << patterns.size()
              << " mod2_filter=" << (ctx.opts.max_sparsity >= 0 ? "on" : "off")
              << std::endl;
    std::uint64_t steps = ctx.opts.limit;
    if (!steps) steps = std::numeric_limits<std::uint64_t>::max();
    stats.total_presentations = steps == std::numeric_limits<std::uint64_t>::max()
        ? cpp_int(-1)
        : cpp_int(steps);
    std::string total_label = steps == std::numeric_limits<std::uint64_t>::max()
        ? "?"
        : std::to_string(steps);
    print_progress_block(ctx, stats, total_label, 0.0, random_pattern_progress_summary(patterns));
    for (std::uint64_t step = 0; step < steps && !stop_requested(); ++step) {
        std::size_t pattern_index = choose_random_pattern(patterns, stats.processed, rng);
        BranchDivisorTypeState& pattern_state = patterns[pattern_index];
        BranchCandidate candidate = random_candidate_from_pattern(ctx, pattern_state, rng);
        CandidateOutcome outcome = process_candidate(ctx, writer, candidate, stats);
        ++pattern_state.attempts;
        if (outcome != CandidateOutcome::HasseWittRejected) {
            ++pattern_state.hasse_witt_passes;
        }
        if (outcome == CandidateOutcome::ExactRejected) {
            ++pattern_state.exact_rejections;
        } else if (outcome == CandidateOutcome::Sparse) {
            ++pattern_state.sparse_hits;
        }
        ++stats.processed;
        if (ctx.opts.progress_interval > 0 && stats.processed % static_cast<std::uint64_t>(ctx.opts.progress_interval) == 0) {
            double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
            print_progress_block(ctx, stats, total_label, elapsed, random_pattern_progress_summary(patterns));
            writer.write_progress(stats, elapsed);
        }
    }
}

std::filesystem::path default_out_path(const Options& opts) {
    std::string sparsity = opts.max_sparsity < 0 ? "all" : "s_" + std::to_string(opts.max_sparsity);
    std::string suffix = opts.mode == "enumerate" ? "cpp" : opts.mode + "_cpp";
    return std::filesystem::path("results") /
        ("p" + std::to_string(opts.p) + "_g" + std::to_string(opts.genus) + "_" + sparsity + "_" + suffix + ".sqlite");
}

Options parse_args(int argc, char** argv) {
    Options opts;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error(std::string("missing value for ") + name);
            return argv[++i];
        };
        if (arg == "--p") opts.p = std::stoul(need("--p"));
        else if (arg == "--genus" || arg == "-g") opts.genus = std::stoi(need("--genus"));
        else if (arg == "--genus-start") opts.genus_start = std::stoi(need("--genus-start"));
        else if (arg == "--genus-end") opts.genus_end = std::stoi(need("--genus-end"));
        else if (arg == "--genus-step") opts.genus_step = std::stoi(need("--genus-step"));
        else if (arg == "--max-sparsity") opts.max_sparsity = std::stoi(need("--max-sparsity"));
        else if (arg == "--limit") opts.limit = std::stoull(need("--limit"));
        else if (arg == "--progress-interval") opts.progress_interval = std::stoi(need("--progress-interval"));
        else if (arg == "--enumeration-mode") opts.mode = need("--enumeration-mode");
        else if (arg == "--random-seed") opts.random_seed = std::stoul(need("--random-seed"));
        else if (arg == "--random-max-factors") opts.random_max_factors = std::stoi(need("--random-max-factors"));
        else if (arg == "--irreducible-memory-budget-mb") opts.irreducible_memory_budget_mb = std::stoull(need("--irreducible-memory-budget-mb"));
        else if (arg == "--out") opts.out = need("--out");
        else if (arg == "--out-dir") opts.out_dir = need("--out-dir");
        else {
            throw std::runtime_error("unknown argument: " + arg);
        }
    }
    if (opts.p < 3 || opts.p % 2 == 0) throw std::runtime_error("--p must be an odd prime");
    if (opts.mode != "enumerate" && opts.mode != "random") {
        throw std::runtime_error("--enumeration-mode must be enumerate or random");
    }
    return opts;
}

void run_one(Options opts) {
    if (opts.genus < 1) throw std::runtime_error("--genus must be positive");
    if (opts.out.empty()) opts.out = default_out_path(opts);
    auto started = std::chrono::steady_clock::now();
    Stats stats;
    SqliteWriter writer(opts.out);
    {
        Context ctx(opts);
        if (!stop_requested()) {
            if (opts.mode == "enumerate") enumerate_mode(ctx, writer, stats);
            else random_mode(ctx, writer, stats);
        }
    }
    double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    writer.write_summary(opts, stats, elapsed);
    writer.close();
    std::cout << "output: " << opts.out << "\n"
              << "stats: processed=" << stats.processed
              << " sparse_presentations=" << stats.sparse_presentations
              << " sparse=" << stats.sparse
              << " duplicate=" << stats.duplicate
              << " rejected_hasse_witt=" << stats.rejected_hasse_witt
              << " rejected_exact=" << stats.rejected_exact << "\n"
              << "run_status: " << run_status(stats) << "\n"
              << "elapsed_seconds: " << elapsed << std::endl;
}

}  // namespace

int main(int argc, char** argv) {
    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);
    try {
        Options opts = parse_args(argc, argv);
        if (opts.genus_start >= 0 || opts.genus_end >= 0) {
            if (opts.genus_start < 1 || opts.genus_end < opts.genus_start) {
                throw std::runtime_error("invalid genus batch range");
            }
            for (int g = opts.genus_start; g <= opts.genus_end; g += opts.genus_step) {
                Options one = opts;
                one.genus = g;
                std::string sparsity = one.max_sparsity < 0 ? "all" : "s_" + std::to_string(one.max_sparsity);
                std::string suffix = one.mode == "enumerate" ? "cpp" : one.mode + "_cpp";
                one.out = one.out_dir / ("p" + std::to_string(one.p) + "_g" + std::to_string(g) + "_" + sparsity + "_" + suffix + ".sqlite");
                std::cout << "batch: p=" << one.p << " genus=" << g << " out=" << one.out << std::endl;
                run_one(one);
                if (stop_requested()) return 130;
            }
        } else {
            run_one(opts);
            if (stop_requested()) return 130;
        }
        return 0;
    } catch (const std::exception& exc) {
        std::cerr << "error: " << exc.what() << std::endl;
        return 1;
    }
}
