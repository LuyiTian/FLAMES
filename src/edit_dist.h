#include <string>
#include <vector>
#include <map>
#include <bitset>
#include <algorithm>
#include <stdlib.h>

namespace scutil{
    int hamming_distance(const std::string &A, const std::string &B);
    double edit_distance(const std::string& A, const std::string& B);
    unsigned int edit_distance1(const int64_t *a, const unsigned int asize, const int64_t *b, const unsigned int bsize);

}

