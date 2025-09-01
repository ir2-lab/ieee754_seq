#include "ieee754_seq.h"

#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

template <class T>
void test()
{
    typedef ieee754_seq<T, 2, 0, 1> seq_t;

    cout << left;
    int w = 10;
    for (auto i = seq_t::begin(); i != i.end(); ++i) {
        cout << i << ' ';
        cout << setw(w) << *i << ' ';
        cout << setw(w) << i.log2v() << endl;
    }
};

template <class T>
void test_lin_interp()
{
    typedef ieee754_seq<T, 3, 0, 2> seq_t;
    typedef typename seq_t::log_lin_interp interp_t;

    vector<T> y(seq_t::count);
    for (int i = 0; i < y.size(); ++i) {
        T v = (seq_t::at(i) - seq_t::min) / (seq_t::max - seq_t::min);
        y[i] = v * v;
    }

    interp_t interp(y);

    cout << left;
    int w = 10;
    for (int i = 0; i < y.size() - 1; ++i) {
        T x = 0.5 * (seq_t::at(i) + seq_t::at(i + 1));
        T v = (x - seq_t::min) / (seq_t::max - seq_t::min);
        v *= v;
        cout << setw(w) << x << ' ';
        cout << setw(w) << v << ' ';
        cout << setw(w) << interp(x) << endl;
    }
};

template <class T>
void test_log_interp()
{
    typedef ieee754_seq<T, 3, 0, 2> seq_t;
    typedef typename seq_t::log_log_interp interp_t;

    vector<T> y(seq_t::count);
    for (int i = 0; i < y.size(); ++i) {
        T v = (seq_t::at(i) - seq_t::min) / (seq_t::max - seq_t::min);
        y[i] = std::exp2(v * v);
    }

    interp_t interp(y);

    cout << left;
    int w = 10;
    for (int i = 0; i < y.size() - 1; ++i) {
        T x = 0.5 * (seq_t::at(i) + seq_t::at(i + 1));
        T v = (x - seq_t::min) / (seq_t::max - seq_t::min);
        v = std::exp2(v * v);
        cout << setw(w) << x << ' ';
        cout << setw(w) << v << ' ';
        cout << setw(w) << interp(x) << endl;
    }
};

int main()
{
    cout << "float sequence" << endl;
    cout << "==============" << endl;
    test<float>();
    cout << endl;

    cout << "double sequence" << endl;
    cout << "===============" << endl;
    test<double>();
    cout << endl;

    cout << "float lin interp" << endl;
    cout << "==============" << endl;
    test_lin_interp<float>();
    cout << endl;

    cout << "float log interp" << endl;
    cout << "==============" << endl;
    test_log_interp<float>();
    cout << endl;

    return 0;
}
