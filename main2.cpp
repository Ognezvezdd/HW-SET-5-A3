#include <iostream>
#include <stack>
#include <deque>
#include <random>
#include <map>
#include <cmath>
#include <set>
#include <unordered_set>
#include <sstream>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <set>
#include <algorithm>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
using namespace std;

using namespace std;
#define all(a) a.begin(), a.end()
#define pb push_back
#define get(v) for (int & iq : v) cin >> iq
#define give(vv) for (int & iqq : vv) cout << iqq << " "
#define vi vector <int>
#define pii pair <int, int>
#define SOLVE int t; cin >> t; while (t--) {solve();}
typedef __int128 lll;
typedef long long ll;
typedef long double ld;
#define int ll
ll inf = 1e9 + 7, mod = 1e6 + 3;
int N = 1e6;

class RandomStreamGen {
public:
    explicit RandomStreamGen(uint64_t seed = std::random_device{}()) : rng(seed) {
        alphabet =
                "abcdefghijklmnopqrstuvwxyz"
                "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                "0123456789"
                "-";
    }

    string nextString(size_t maxLen = 30) {
        uniform_int_distribution<int> lenDist(1, (int) maxLen);
        uniform_int_distribution<int> chDist(0, (int) alphabet.size() - 1);

        int len = lenDist(rng);
        string s;
        s.reserve(len);
        for (int i = 0; i < len; i++) s.push_back(alphabet[chDist(rng)]);
        return s;
    }

    vector<string> generateStream(size_t N, size_t maxLen = 30) {
        vector<string> S;
        S.reserve(N);
        for (size_t i = 0; i < N; i++) S.push_back(nextString(maxLen));
        return S;
    }

    static vector<double> makeTimeFractions(double step = 0.05) {
        vector<double> t;
        for (double x = step; x <= 1.0 + 1e-12; x += step) t.push_back(min(1.0, x));
        return t;
    }

private:
    mt19937_64 rng;
    string alphabet;
};

class HashFuncGen {
public:
    explicit HashFuncGen(uint64_t seed = std::random_device{}()) : seed(seed) {
    }

    uint32_t hash32(const string &s) const {
        uint64_t h = 1469598103934665603ULL ^ seed;
        for (unsigned char c: s) {
            h ^= (uint64_t) c;
            h *= 1099511628211ULL;
        }

        h ^= (h >> 33);
        h *= 0xff51afd7ed558ccdULL;
        h ^= (h >> 33);
        h *= 0xc4ceb9fe1a85ec53ULL;
        h ^= (h >> 33);

        return (uint32_t) (h & 0xFFFFFFFFu);
    }

private:
    uint64_t seed;
};

class HyperLogLog {
public:
    explicit HyperLogLog(int B) : B(B), m(1u << B), reg(m, 0) {
        if (B < 4 || B > 16) throw runtime_error("B should be in [4..16]");
    }

    void reset() { fill(reg.begin(), reg.end(), 0); }

    void add(uint32_t x) {
        uint32_t idx = x >> (32 - B);
        uint32_t w = (x << B);

        int rank = rho(w, 32 - B);
        reg[idx] = max(reg[idx], (uint8_t) rank);
    }

    double estimate() const {
        double alpha = alphaM(m);

        double sum = 0.0;
        int zeros = 0;
        for (uint8_t r: reg) {
            sum += ldexp(1.0, -(int) r);
            if (r == 0) zeros++;
        }

        double raw = alpha * (double) m * (double) m / sum;

        if (raw <= 2.5 * m && zeros > 0) {
            return (double) m * log((double) m / (double) zeros);
        }

        return raw;
    }

    uint32_t registersCount() const { return m; }

private:
    int B;
    uint32_t m;
    vector<uint8_t> reg;

    static int rho(uint32_t w, int bits) {
        if (w == 0) return bits + 1;

        int lz = __builtin_clz(w);
        int rank = lz + 1;

        if (rank > bits + 1) rank = bits + 1;
        return rank;
    }

    static double alphaM(uint32_t m) {
        if (m == 16) return 0.673;
        if (m == 32) return 0.697;
        if (m == 64) return 0.709;
        return 0.7213 / (1.0 + 1.079 / (double) m);
    }
};

static vector<size_t> timeToPrefixSizes(size_t N, const vector<double> &times) {
    vector<size_t> ks;
    ks.reserve(times.size());
    for (double t: times) ks.push_back((size_t) floor(N * t));
    for (auto &k: ks) k = min(k, N);
    for (size_t i = 1; i < ks.size(); i++) ks[i] = max(ks[i], ks[i - 1]);
    return ks;
}

signed main() {
#ifdef _LOCAL
    freopen("input.txt", "r", stdin);
#endif
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    const int B = 14;
    const size_t N = 200000;
    const int STREAMS = 30;
    const double STEP = 0.05;

    auto times = RandomStreamGen::makeTimeFractions(STEP);
    auto ks = timeToPrefixSizes(N, times);

    vector<double> one_F0(times.size()), one_N(times.size());

    vector<vector<double> > all_estimates(times.size());

    for (int s = 0; s < STREAMS; s++) {
        RandomStreamGen gen(12345 + s * 999);
        HashFuncGen hgen(77777 + s * 131);
        HyperLogLog hll(B);

        auto S = gen.generateStream(N);

        unordered_set<string> exactSet;
        exactSet.reserve(N * 2);

        size_t prev = 0;
        for (size_t j = 0; j < times.size(); j++) {
            size_t k = ks[j];

            for (size_t i = prev; i < k; i++) {
                exactSet.insert(S[i]);
                hll.add(hgen.hash32(S[i]));
            }
            prev = k;

            size_t F0 = exactSet.size();
            double Nt = hll.estimate();

            all_estimates[j].push_back(Nt);

            if (s == 0) {
                one_F0[j] = (double) F0;
                one_N[j] = Nt;
            }
        }

        cerr << "stream " << (s + 1) << "/" << STREAMS << " done\n";
    }

    vector<double> meanN(times.size()), stdN(times.size());
    for (size_t j = 0; j < times.size(); j++) {
        auto &v = all_estimates[j];
        double mean = accumulate(v.begin(), v.end(), 0.0) / (double) v.size();
        meanN[j] = mean;

        double var = 0.0;
        for (double x: v) var += (x - mean) * (x - mean);
        var /= (double) v.size();
        stdN[j] = sqrt(var);
    }

    // graph1.csv
    {
        ofstream out("graph1.csv");
        out << "t,F0,N\n";
        for (size_t j = 0; j < times.size(); j++) {
            out << times[j] << "," << one_F0[j] << "," << one_N[j] << "\n";
        }
    }

    // graph2.csv
    {
        ofstream out("graph2.csv");
        out << "t,meanN,stdN,upper,lower\n";
        for (size_t j = 0; j < times.size(); j++) {
            double upper = meanN[j] + stdN[j];
            double lower = meanN[j] - stdN[j];
            out << times[j] << "," << meanN[j] << "," << stdN[j] << ","
                    << upper << "," << lower << "\n";
        }
    }

    uint32_t m = (1u << B);
    double rse104 = 1.04 / sqrt((double) m);
    double rse132 = 1.32 / sqrt((double) m);
    cout << "Saved graph1.csv and graph2.csv\n";
    cout << "B=" << B << ", m=" << m << "\n";
    cout << "Theoretical RSE ~ 1.04/sqrt(m) = " << rse104 * 100.0 << "%\n";
    cout << "Wider bound     ~ 1.32/sqrt(m) = " << rse132 * 100.0 << "%\n";
    cout << "Registers memory ~ " << (m * sizeof(uint8_t)) << " bytes (+vector overhead)\n";

    return 0;
}
