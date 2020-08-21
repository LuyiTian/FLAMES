#include "edit_dist.h"

int scutil::hamming_distance(const std::string &A, const std::string &B)
{
    int dist = 0;
    for (int i = 0; i < B.length(); ++i)
    {
        dist += (A[i] != B[i]);
    }
    return dist;
}


double scutil::edit_distance(const std::string& A, const std::string& B)
{
    int NA = A.size();
    int NB = B.size();
    double x, y, z;

    std::vector<std::vector<double>> M(NA + 1, std::vector<double>(NB + 1));

    for (int i = 0; i <= NA; ++i)
        M[i][0] = i;

    for (int i = 0; i <= NB; ++i)
        M[0][i] = i;

    for (int a = 1; a <= NA; ++a)
    {
        for (int b = 1; b <= NB; ++b)
        {
            x = M[a - 1][b] + 1.01;
            y = M[a][b - 1] + 1.001;
            z = M[a - 1][b - 1] + (A[a - 1] == B[b - 1] ? 0 : 1);

            M[a][b] = std::min(std::min(x, y), z);
        }
    }


    return M[A.size()][B.size()];
}

// https://github.com/aflc/editdistance/blob/master/editdistance/_editdistance.cpp

template<typename T, typename TVALUE>
unsigned int edit_distance_bpv(T &cmap, int64_t const *vec, size_t const &vecsize, unsigned int const &tmax, unsigned int const &tlen) {
    int D = tmax * 64 + tlen;
    TVALUE D0, HP, HN, VP, VN;
    uint64_t top = (1LL << (tlen - 1));
    uint64_t lmb = (1LL << 63);

    for(size_t i = 0; i <= tmax; ++i) {
        VP[i] = 0;
        VN[i] = 0;
    }
    for(size_t i = 0; i < tmax; ++i) VP[i] = ~0;
    for(size_t i = 0; i < tlen; ++i) VP[tmax] |= (1LL << i);
    for(size_t i = 0; i < vecsize; ++i) {
        TVALUE &PM = cmap[vec[i]];
        for(int r = 0; r <= tmax; ++r) {
            uint64_t X = PM[r];
            if(r > 0 && (HN[r - 1] & lmb)) X |= 1LL;
            D0[r] = (((X & VP[r]) + VP[r]) ^ VP[r]) | X | VN[r];
            HP[r] = VN[r] | ~(D0[r] | VP[r]);
            HN[r] = D0[r] & VP[r];
            X = (HP[r] << 1LL);
            if(r == 0 || HP[r - 1] & lmb) X |= 1LL;
            VP[r] = (HN[r] << 1LL) | ~(D0[r] | X);
            if(r > 0 && (HN[r - 1] & lmb)) VP[r] |= 1LL;
            VN[r] = D0[r] & X;
        }
        if(HP[tmax] & top) ++D;
        else if(HN[tmax] & top) --D;
    }
    return D;
}


template <size_t N>
struct varr {
    uint64_t arr_[N];
    uint64_t & operator[](size_t const &i) {
        return arr_[i];
    }
};

template<size_t N>
unsigned int edit_distance_map_(int64_t const *a, size_t const asize, int64_t const *b, size_t const bsize) {
    typedef std::map<int64_t, varr<N> > cmap_v;
    cmap_v cmap;
    unsigned int tmax = (asize - 1) >> 6;
    unsigned int tlen = asize - tmax * 64;
    for(size_t i = 0; i < tmax; ++i) {
        for(size_t j = 0; j < 64; ++j) cmap[a[i * 64 + j]][i] |= (1LL << j);
    }
    for(size_t i = 0; i < tlen; ++i) cmap[a[tmax * 64 + i]][tmax] |= (1LL << i);
    return edit_distance_bpv<cmap_v, typename cmap_v::mapped_type>(cmap, b, bsize, tmax, tlen);
}


unsigned int scutil::edit_distance1(const int64_t *a, const unsigned int asize, const int64_t *b, const unsigned int bsize) 
{
    if(asize == 0) return bsize;
    else if(bsize == 0) return asize;
    int64_t const *ap, *bp;
    unsigned int const *asizep, *bsizep;
    if(asize < bsize) ap = b, bp = a, asizep = &bsize, bsizep = &asize;
    else ap = a, bp = b, asizep = &asize, bsizep = &bsize;

    size_t vsize = ((*asizep - 1) >> 6) + 1;  // usually the barcode won't exceed 128bp

    if(vsize == 1) return edit_distance_map_<1>(ap, *asizep, bp, *bsizep);
    else if(vsize == 2) return edit_distance_map_<2>(ap, *asizep, bp, *bsizep);
}