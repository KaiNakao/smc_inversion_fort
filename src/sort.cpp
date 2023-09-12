#include <algorithm>

extern "C" void qsort_(int *arr, int *index, int *n) {
    std::sort(index, index+(*n),[&](int a, int b) {
        return arr[a - 1] > arr[b - 1];
        });
}