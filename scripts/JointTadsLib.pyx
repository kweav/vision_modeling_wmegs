# distutils: language = c++
#cython: language_level=3

import cython
cimport numpy as np
import numpy

ctypedef np.float32_t DTYPE_t
ctypedef np.float64_t DTYPE_64_t
ctypedef np.int32_t DTYPE_int_t
ctypedef np.int64_t DTYPE_int64_t
ctypedef np.uint32_t DTYPE_uint_t
ctypedef np.int8_t DTYPE_int8_t
cdef double Inf = numpy.inf

cdef extern from "math.h":
    double exp(double x) nogil
    double log(double x) nogil
    double log10(double x) nogil
    double sqrt(double x) nogil
    double pow(double x, double x) nogil
    double abs(double x) nogil
    double round(double x) nogil
    double floor(double x) nogil
    double ceil(double x) nogil


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_indicator_tads(
        np.ndarray[DTYPE_int_t, ndim=2] joint_EP not None,
        np.ndarray[DTYPE_int_t, ndim=1] coords not None,
        np.ndarray[DTYPE_t, ndim=2] RNA not None,
        np.ndarray[DTYPE_t, ndim=2] cStateBetas not None,
        np.ndarray[DTYPE_t, ndim=2] iStateBetas not None,
        np.ndarray[DTYPE_t, ndim=2] pStateBetas not None,
        np.ndarray[DTYPE_t, ndim=1] scores not None,
        np.ndarray[DTYPE_int_t, ndim=1] paths not None,
        np.ndarray[DTYPE_t, ndim=1] cStateSum not None,
        np.ndarray[DTYPE_int_t, ndim=1] TSS_list not None,
        np.ndarray[DTYPE_int_t, ndim=2] tads not None,
        int max_size):
    cdef long long int i, j, k, l, list_N, cindex, tindex, start, end, tadN
    cdef double mse, pred
    cdef long long int jointN = joint_EP.shape[0]
    cdef long long int tssN = RNA.shape[0]
    cdef long long int cellN = RNA.shape[1]
    with nogil:
        paths[0] = 0
        scores[0] = 0
        for i in range(1, jointN + 1):
            scores[i] = Inf
            paths[i] = -1
        for i in range(jointN):
            start = coords[i]
            list_N = 0
            for l in range(cellN):
                cStateSum[l] = 0
            if scores[i] == Inf:
                j = i - 1
                while j >= 0:
                    if joint_EP[j, 1] == 1:
                        break
                    if scores[j] < scores[i]:
                        scores[i] = scores[j]
                        paths[i] = j
                    j -= 1
            j = i + 1
            while j <= jointN and coords[j - 1] - coords[i] < max_size:
                if joint_EP[j - 1, 1] == 0:
                    cindex = joint_EP[j - 1, 0]
                    for l in range(cellN):
                        cStateSum[l] += cStateBetas[cindex, l]
                else:
                    TSS_list[list_N] = joint_EP[j - 1, 0]
                    list_N += 1
                if list_N > 0:
                    mse = scores[i]
                    for k in range(list_N):
                        tindex = TSS_list[k]
                        for l in range(cellN):
                            pred = (cStateSum[l] * iStateBetas[tindex, l] +
                                    pStateBetas[tindex, l])
                            mse += pow(RNA[tindex, l] - pred, 2)
                    if mse < scores[j]:
                        scores[j] = mse
                        paths[j] = i
                j += 1
        if scores[jointN] == Inf:
            j = jointN - 1
            while j >= 0:
                if joint_EP[j, 1] == 1:
                    break
                if scores[j] < scores[jointN]:
                    scores[jointN] = scores[j]
                    paths[jointN] = j
                j -= 1
        tadN = 0
        start = jointN
        while start > 0:
            end = start
            start = paths[end]
            list_N = 0
            for i in range(start, end):
                if joint_EP[i, 1] == 1:
                    list_N += 1
            if list_N > 0:
                tads[tadN, 0] = start
                tads[tadN, 1] = end
                tadN += 1
    return tadN

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_new_tads(
        np.ndarray[DTYPE_int_t, ndim=2] joint_EP not None,
        np.ndarray[DTYPE_int_t, ndim=1] coords not None,
        np.ndarray[DTYPE_t, ndim=2] RNA not None,
        np.ndarray[DTYPE_t, ndim=2] cStateBetas not None,
        np.ndarray[DTYPE_t, ndim=2] pStateBetas not None,
        np.ndarray[DTYPE_t, ndim=1] scores not None,
        np.ndarray[DTYPE_int_t, ndim=1] paths not None,
        np.ndarray[DTYPE_t, ndim=1] cStateSum not None,
        np.ndarray[DTYPE_int_t, ndim=1] TSS_list not None,
        np.ndarray[DTYPE_int_t, ndim=2] tads not None,
        int max_size):
    cdef long long int i, j, k, l, list_N, cindex, tindex, start, end, tadN
    cdef double mse, pred
    cdef long long int jointN = joint_EP.shape[0]
    cdef long long int tssN = RNA.shape[0]
    cdef long long int cellN = RNA.shape[1]
    with nogil:
        paths[0] = 0
        scores[0] = 0
        for i in range(1, jointN + 1):
            scores[i] = Inf
            paths[i] = -1
        for i in range(jointN):
            start = coords[i]
            list_N = 0
            for l in range(cellN):
                cStateSum[l] = 0
            if scores[i] == Inf:
                j = i - 1
                while j >= 0:
                    if joint_EP[j, 1] == 1:
                        break
                    if scores[j] < scores[i]:
                        scores[i] = scores[j]
                        paths[i] = j
                    j -= 1
            j = i + 1
            while j <= jointN and coords[j - 1] - coords[i] < max_size:
                if joint_EP[j - 1, 1] == 0:
                    cindex = joint_EP[j - 1, 0]
                    for l in range(cellN):
                        cStateSum[l] += cStateBetas[cindex, l]
                else:
                    TSS_list[list_N] = joint_EP[j - 1, 0]
                    list_N += 1
                if list_N > 0:
                    mse = scores[i]
                    for k in range(list_N):
                        tindex = TSS_list[k]
                        for l in range(cellN):
                            pred = (1 + cStateSum[l]) * pStateBetas[tindex, l]
                            mse += pow(RNA[tindex, l] - pred, 2)
                    if mse < scores[j]:
                        scores[j] = mse
                        paths[j] = i
                j += 1
        if scores[jointN] == Inf:
            j = jointN - 1
            while j >= 0:
                if joint_EP[j, 1] == 1:
                    break
                if scores[j] < scores[jointN]:
                    scores[jointN] = scores[j]
                    paths[jointN] = j
                j -= 1
        tadN = 0
        start = jointN
        while start > 0:
            end = start
            start = paths[end]
            list_N = 0
            for i in range(start, end):
                if joint_EP[i, 1] == 1:
                    list_N += 1
            if list_N > 0:
                tads[tadN, 0] = start
                tads[tadN, 1] = end
                tadN += 1
    return tadN


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_new_tads1(
        np.ndarray[DTYPE_int_t, ndim=2] joint_EP not None,
        np.ndarray[DTYPE_int_t, ndim=1] coords not None,
        np.ndarray[DTYPE_t, ndim=2] RNA not None,
        np.ndarray[DTYPE_t, ndim=2] cStateBetas not None,
        np.ndarray[DTYPE_t, ndim=2] pStateBetas not None,
        np.ndarray[DTYPE_t, ndim=1] scores not None,
        np.ndarray[DTYPE_int_t, ndim=1] paths not None,
        np.ndarray[DTYPE_t, ndim=1] cStateSum not None,
        np.ndarray[DTYPE_int_t, ndim=1] TSS_list not None,
        np.ndarray[DTYPE_int_t, ndim=2] tads not None,
        int max_size):
    cdef long long int i, j, k, l, list_N, cindex, tindex, start, end, tadN
    cdef double mse, pred
    cdef long long int jointN = joint_EP.shape[0]
    cdef long long int tssN = RNA.shape[0]
    cdef long long int cellN = RNA.shape[1]
    with nogil:
        paths[0] = 0
        scores[0] = 0
        for i in range(1, jointN + 1):
            scores[i] = Inf
            paths[i] = -1
        for i in range(jointN):
            start = coords[i]
            list_N = 0
            for l in range(cellN):
                cStateSum[l] = 0
            if scores[i] == Inf:
                j = i - 1
                while j >= 0:
                    if joint_EP[j, 1] == 1:
                        break
                    if scores[j] < scores[i]:
                        scores[i] = scores[j]
                        paths[i] = j
                    j -= 1
            j = i + 1
            while j <= jointN and coords[j - 1] - coords[i] < max_size:
                if joint_EP[j - 1, 1] == 0:
                    cindex = joint_EP[j - 1, 0]
                    for l in range(cellN):
                        cStateSum[l] += cStateBetas[cindex, l]
                else:
                    TSS_list[list_N] = joint_EP[j - 1, 0]
                    list_N += 1
                if list_N > 0:
                    mse = scores[i]
                    for k in range(list_N):
                        tindex = TSS_list[k]
                        for l in range(cellN):
                            pred = cStateSum[l] + pStateBetas[tindex, l]
                            mse += pow(RNA[tindex, l] - pred, 2)
                    if mse < scores[j]:
                        scores[j] = mse
                        paths[j] = i
                j += 1
        if scores[jointN] == Inf:
            j = jointN - 1
            while j >= 0:
                if joint_EP[j, 1] == 1:
                    break
                if scores[j] < scores[jointN]:
                    scores[jointN] = scores[j]
                    paths[jointN] = j
                j -= 1
        tadN = 0
        start = jointN
        while start > 0:
            end = start
            start = paths[end]
            list_N = 0
            for i in range(start, end):
                if joint_EP[i, 1] == 1:
                    list_N += 1
            if list_N > 0:
                tads[tadN, 0] = start
                tads[tadN, 1] = end
                tadN += 1
    return tadN


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_new_tads2(
        np.ndarray[DTYPE_int_t, ndim=2] joint_EP not None,
        np.ndarray[DTYPE_int_t, ndim=1] coords not None,
        np.ndarray[DTYPE_t, ndim=2] RNA not None,
        np.ndarray[DTYPE_t, ndim=2] cStateBetas not None,
        np.ndarray[DTYPE_t, ndim=2] pStateBetas not None,
        np.ndarray[DTYPE_t, ndim=1] scores not None,
        np.ndarray[DTYPE_int_t, ndim=1] paths not None,
        np.ndarray[DTYPE_t, ndim=1] cStateSum not None,
        np.ndarray[DTYPE_int_t, ndim=1] TSS_list not None,
        np.ndarray[DTYPE_int_t, ndim=2] tads not None,
        int max_size):
    cdef long long int i, j, k, l, list_N, cindex, tindex, start, end, tadN
    cdef double mse, pred
    cdef long long int jointN = joint_EP.shape[0]
    cdef long long int tssN = RNA.shape[0]
    cdef long long int cellN = RNA.shape[1]
    with nogil:
        paths[0] = 0
        scores[0] = 0
        for i in range(1, jointN + 1):
            scores[i] = Inf
            paths[i] = -1
        for i in range(jointN):
            start = coords[i]
            list_N = 0
            for l in range(cellN):
                cStateSum[l] = 0
            if scores[i] == Inf:
                j = i - 1
                while j >= 0:
                    if joint_EP[j, 1] == 1:
                        break
                    if scores[j] < scores[i]:
                        scores[i] = scores[j]
                        paths[i] = j
                    j -= 1
            j = i + 1
            while j <= jointN and coords[j - 1] - coords[i] < max_size:
                if joint_EP[j - 1, 1] == 0:
                    cindex = joint_EP[j - 1, 0]
                    for l in range(cellN):
                        cStateSum[l] += cStateBetas[cindex, l]
                else:
                    TSS_list[list_N] = joint_EP[j - 1, 0]
                    list_N += 1
                if list_N > 0:
                    mse = scores[i]
                    for k in range(list_N):
                        tindex = TSS_list[k]
                        for l in range(cellN):
                            pred = (cStateSum[l] + 1) * pStateBetas[tindex, l]
                            mse += pow(RNA[tindex, l] - pred, 2)
                    if mse < scores[j]:
                        scores[j] = mse
                        paths[j] = i
                j += 1
        if scores[jointN] == Inf:
            j = jointN - 1
            while j >= 0:
                if joint_EP[j, 1] == 1:
                    break
                if scores[j] < scores[jointN]:
                    scores[jointN] = scores[j]
                    paths[jointN] = j
                j -= 1
        tadN = 0
        start = jointN
        while start > 0:
            end = start
            start = paths[end]
            list_N = 0
            for i in range(start, end):
                if joint_EP[i, 1] == 1:
                    list_N += 1
            if list_N > 0:
                tads[tadN, 0] = start
                tads[tadN, 1] = end
                tadN += 1
    return tadN


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_cre_tads(
        np.ndarray[DTYPE_t, ndim=2] RNA not None,
        np.ndarray[DTYPE_t, ndim=2] cStateBetas not None,
        np.ndarray[DTYPE_int_t, ndim=2] TSS_ranges not None,
        np.ndarray[DTYPE_int_t, ndim=1] TSS_mapping not None,
        np.ndarray[DTYPE_t, ndim=1] scores not None,
        np.ndarray[DTYPE_int_t, ndim=1] paths not None,
        np.ndarray[DTYPE_t, ndim=1] cStateSum not None,
        np.ndarray[DTYPE_int_t, ndim=1] TSS_list not None,
        np.ndarray[DTYPE_int_t, ndim=2] tads not None):
    cdef long long int i, j, j2, k, l, m, mid, listN, cindex, tindex, start, end, tadN
    cdef double mse, pred
    cdef long long int creN = cStateBetas.shape[0]
    cdef long long int tssN = RNA.shape[0]
    cdef long long int cellN = RNA.shape[1]
    with nogil:
        paths[0] = 0
        scores[0] = 0
        for i in range(1, creN + 1):
            scores[i] = Inf
            paths[i] = -1
        for i in range(tssN):
            mid = TSS_ranges[i, 1]
            for j2 in range(mid - TSS_ranges[i, 0] + 1):
                j = mid - j2
                if scores[j] == Inf:
                    k = j - 1
                    while k >= 0:
                        if TSS_mapping[k] >= 0:
                            break
                        if scores[k] < scores[j]:
                            scores[j] = scores[k]
                            paths[j] = k
                        k -= 1
                for l in range(cellN):
                    cStateSum[l] = 0
                for k in range(j, mid):
                    for l in range(cellN):
                        cStateSum[l] += cStateBetas[k, l]
                listN = 0
                for k in range(mid, TSS_ranges[i, 2]):
                    if TSS_mapping[k] >= 0:
                        TSS_list[listN] = k
                        listN += 1
                    else:
                        for l in range(cellN):
                            cStateSum[l] += cStateBetas[k, l]
                    mse = scores[j]
                    for m in range(listN):
                        cindex = TSS_list[m]
                        tindex = TSS_mapping[cindex]
                        for l in range(cellN):
                            pred = cStateSum[l] + cStateBetas[cindex, l]
                            mse += pow(RNA[tindex, l] - pred, 2)
                    if mse < scores[k + 1]:
                        scores[k + 1] = mse
                        paths[k + 1] = j
        k = creN
        j = k - 1
        while TSS_mapping[j] == -1:
            if scores[j] < scores[k]:
                k = j
            j -= 1
        tadN = 0
        start = k
        while start > 0:
            end = start
            start = paths[end]
            listN = 0
            for i in range(start, end):
                if TSS_mapping[i] >= 0:
                    listN += 1
                    break
            if listN > 0:
                tads[tadN, 0] = start
                tads[tadN, 1] = end
                tadN += 1
    return tadN


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_new_tads_single(
        np.ndarray[DTYPE_int_t, ndim=2] joint_EP not None,
        np.ndarray[DTYPE_int_t, ndim=1] TSSs not None,
        np.ndarray[DTYPE_t, ndim=1] RNAs not None,
        np.ndarray[DTYPE_t, ndim=1] CREs not None,
        np.ndarray[DTYPE_t, ndim=1] scores not None,
        np.ndarray[DTYPE_int_t, ndim=1] paths not None,
        np.ndarray[DTYPE_t, ndim=1] predicted not None,
        np.ndarray[DTYPE_int_t, ndim=1] TSS_list not None,
        np.ndarray[DTYPE_int_t, ndim=2] tads not None):
    cdef long long int i, j, j_rev, k, m, tss, list_N, cindex
    cdef long long int tindex, prev_path, start, end, tadN
    cdef double prev_best, new_score, score
    cdef long long int jointN = joint_EP.shape[0]
    cdef long long int tssN = TSSs.shape[0]
    with nogil:
        for i in range(tssN):
            tss = TSSs[i]
            predicted[0] = 0
            for j_rev in range(tss + 1 - joint_EP[tss, 2]):
                j = tss - j_rev
                if joint_EP[j, 1] == 0:
                    cindex = joint_EP[j, 0]
                    predicted[0] += CREs[cindex]
                if scores[j] == Inf:
                    k = j - 1
                    while k > 0 and scores[k] == Inf:
                        k -= 1
                    prev_best = scores[k]
                    prev_path = k
                else:
                    prev_best = scores[j]
                    prev_path = j
                predicted[1] = predicted[0]
                list_N = 0
                for k in range(tss, joint_EP[j, 3]):
                    new_score = prev_best
                    if joint_EP[k, 1] == 0:
                        cindex = joint_EP[k, 0]
                        predicted[1] += CREs[cindex]
                    else:
                        TSS_list[list_N] = joint_EP[k, 0]
                        list_N += 1
                    for m in range(list_N):
                        tindex = TSS_list[m]
                        new_score += pow(predicted[1] - RNAs[tindex], 2)
                    if new_score < scores[k + 1]:
                        scores[k + 1] = new_score
                        paths[k + 1] = prev_path
        end = TSSs[tssN - 1] + 1
        i = end
        score = scores[end]
        while i < jointN and scores[i + 1] < Inf:
            if scores[i + 1] < score:
                score = scores[i + 1]
                end = i
            i += 1
        tadN = 0
        while end > 0:
            start = paths[end]
            for i in range(start, end):
                if joint_EP[i, 1] == 1:
                    tads[tadN, 0] = start
                    tads[tadN, 1] = end
                    tadN += 1
                    break
            end = start
    return tadN


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def beta_gradient(
        np.ndarray[DTYPE_t, ndim=2] betas not None,
        np.ndarray[DTYPE_t, ndim=1] predicted not None,
        np.ndarray[DTYPE_t, ndim=1] observed not None,
        np.ndarray[DTYPE_t, ndim=2] pfeatures not None,
        np.ndarray[DTYPE_t, ndim=2] cfeatures not None,
        np.ndarray[DTYPE_t, ndim=2] pBetaStates not None,
        np.ndarray[DTYPE_t, ndim=2] cBetaStates not None,
        np.ndarray[DTYPE_t, ndim=2] BetaGrad not None,
        double learning_rate):
    # En = sum_i(pBi * pSin) * (1 + sum(cBi * cSin))
    # Cost_n = (En - On)^2 = En^2 - 2EnOn + On^2
    # Cost_n = (sum_i(pBi * pSin) * (1 + sum_i(cBi * cSin)) - On)^2
    # dCost_n/dpBi = 2(En - On) * pSin * (1 + sum(cBi * cSin))
    # dCost_n/dcBi = 2(En - On) * cSin * sum(pBi * pSin)
    cdef long long int i, j
    cdef double mse, dCdE, p, o, dCdE_csum, dCdE_psum
    cdef long long int stateN = BetaGrad.shape[1]
    cdef long long int tssN = predicted.shape[0]
    with nogil:
        for j in range(stateN):
            BetaGrad[0, j] = 0
            BetaGrad[1, j] = 0
        for i in range(tssN):
            p = predicted[i]
            o = observed[i]
            if p != o:
                dCdE = 2 * (p - o)
                dCdE_csum = dCdE * cBetaStates[i, stateN]
                dCdE_psum = dCdE * pBetaStates[i, stateN]
                for j in range(stateN):
                    BetaGrad[0, j] += pBetaStates[i, j] * dCdE_csum
                    BetaGrad[1, j] += cBetaStates[i, j] * dCdE_psum
        for j in range(stateN):
            BetaGrad[0, j] /= tssN
            BetaGrad[1, j] /= tssN
            betas[0, j] += learning_rate * BetaGrad[0, j]
            betas[1, j] += learning_rate * BetaGrad[1, j]
        mse = 0
        for i in range(tssN):
            pBetaStates[i, stateN] = 0
            cBetaStates[i, stateN] = 1
            for j in range(stateN):
                pBetaStates[i, j] = pfeatures[i, j] * betas[0, j]
                cBetaStates[i, j] = cfeatures[i, j] * betas[1, j]
                pBetaStates[i, stateN] += pBetaStates[i, j]
                cBetaStates[i, stateN] += cBetaStates[i, j]
            predicted[i] = pBetaStates[i, stateN] * cBetaStates[i, stateN]
            mse += pow(observed[i] - predicted[i], 2)
    return mse


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def beta_gradient2(
        np.ndarray[DTYPE_t, ndim=2] betas not None,
        np.ndarray[DTYPE_t, ndim=1] predicted not None,
        np.ndarray[DTYPE_t, ndim=1] observed not None,
        np.ndarray[DTYPE_t, ndim=2] pfeatures not None,
        np.ndarray[DTYPE_t, ndim=2] cfeatures not None,
        np.ndarray[DTYPE_t, ndim=2] pBetaStates not None,
        np.ndarray[DTYPE_t, ndim=2] cBetaStates not None,
        np.ndarray[DTYPE_t, ndim=2] BetaGrad not None,
        double learning_rate,
        int update):
    # En = sum_i(pBi * pSin) * (1 + sum(cBi * cSin))
    # Cost_n = (En - On)^2 = En^2 - 2EnOn + On^2
    # dEn/dpBi = pSin * (1 + sum(cBi * cSin))
    # dEn/dcBi = cSin * sum(pBi * pSin)
    cdef long long int i, j
    cdef double mse, dCdE, p, o, dCdE_csum, dCdE_psum
    cdef long long int stateN = BetaGrad.shape[1]
    cdef long long int tssN = predicted.shape[0]
    with nogil:
        for j in range(stateN):
            BetaGrad[0, j] = 0
            BetaGrad[1, j] = 0
        for i in range(tssN):
            p = predicted[i]
            o = observed[i]
            if p != o:
                dCdE = 2 * (p - o)
                dCdE_csum = dCdE * cBetaStates[i, stateN]
                dCdE_psum = dCdE * pBetaStates[i, stateN]
                for j in range(stateN):
                    BetaGrad[0, j] += pBetaStates[i, j] * dCdE_csum
                    BetaGrad[1, j] += cBetaStates[i, j] * dCdE_psum
        for j in range(stateN):
            BetaGrad[0, j] /= tssN
            BetaGrad[1, j] /= tssN
            if update == 1:
                betas[0, j] += learning_rate * BetaGrad[0, j]
            betas[1, j] += learning_rate * BetaGrad[1, j]
        mse = 0
        for i in range(tssN):
            pBetaStates[i, stateN] = 0
            cBetaStates[i, stateN] = 1
            for j in range(stateN):
                pBetaStates[i, j] = pfeatures[i, j] * betas[0, j]
                cBetaStates[i, j] = cfeatures[i, j] * betas[1, j]
                pBetaStates[i, stateN] += pBetaStates[i, j]
                cBetaStates[i, stateN] += cBetaStates[i, j]
            predicted[i] = pBetaStates[i, stateN] * cBetaStates[i, stateN]
            mse += pow(observed[i] - predicted[i], 2)
    return mse


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def beta_indicator_gradient(
        np.ndarray[DTYPE_t, ndim=2] betas not None,
        np.ndarray[DTYPE_t, ndim=1] predicted not None,
        np.ndarray[DTYPE_t, ndim=1] observed not None,
        np.ndarray[DTYPE_t, ndim=2] pfeatures not None,
        np.ndarray[DTYPE_t, ndim=2] cfeatures not None,
        np.ndarray[DTYPE_t, ndim=2] pBetaStates not None,
        np.ndarray[DTYPE_t, ndim=2] iBetaStates not None,
        np.ndarray[DTYPE_t, ndim=2] cBetaStates not None,
        np.ndarray[DTYPE_t, ndim=2] BetaGrad not None,
        double learning_rate,
        int update):
    # En = sum_i(pBi * pSin) + sum_i(IBi * pSin) * sum_i(cBi * cSin)
    # Cost_n = (En - On)^2 = En^2 - 2EnOn + On^2
    # dCost_n/dEn = 2(En - On)
    # dEn/dpBi = pSin
    # dEn/dIBi = pSin * sum_i(cBi * cSin)
    # dEn/dcBi = cSin * sum_i(IBi * pSin)
    cdef long long int i, j
    cdef double mse, dCdE, p, o, dCdE_csum, dCdE_Isum
    cdef long long int stateN = BetaGrad.shape[1]
    cdef long long int tssN = predicted.shape[0]
    with nogil:
        for j in range(stateN):
            BetaGrad[0, j] = 0
            BetaGrad[1, j] = 0
            BetaGrad[2, j] = 0
        for i in range(tssN):
            p = predicted[i]
            o = observed[i]
            if p != o:
                dCdE = 2 * (p - o)
                dCdE_csum = dCdE * cBetaStates[i, stateN]
                dCdE_Isum = dCdE * iBetaStates[i, stateN]
                for j in range(stateN):
                    BetaGrad[0, j] += pfeatures[i, j] * dCdE
                    BetaGrad[1, j] += pfeatures[i, j] * dCdE_csum
                    BetaGrad[2, j] += cfeatures[i, j] * dCdE_Isum
        for j in range(stateN):
            BetaGrad[0, j] /= tssN
            BetaGrad[1, j] /= tssN
            BetaGrad[2, j] /= tssN
            if update > 1:
                betas[0, j] += learning_rate * BetaGrad[0, j]
            if update > 0:
                betas[1, j] += learning_rate * BetaGrad[1, j]
            betas[2, j] += learning_rate * BetaGrad[2, j]
        mse = 0
        for i in range(tssN):
            pBetaStates[i, stateN] = 0
            iBetaStates[i, stateN] = 0
            cBetaStates[i, stateN] = 0
            for j in range(stateN):
                pBetaStates[i, j] = pfeatures[i, j] * betas[0, j]
                iBetaStates[i, j] = pfeatures[i, j] * betas[1, j]
                cBetaStates[i, j] = cfeatures[i, j] * betas[2, j]
                pBetaStates[i, stateN] += pBetaStates[i, j]
                iBetaStates[i, stateN] += iBetaStates[i, j]
                cBetaStates[i, stateN] += cBetaStates[i, j]
            predicted[i] = (pBetaStates[i, stateN] +
                            iBetaStates[i, stateN] *
                            cBetaStates[i, stateN])
            mse += pow(observed[i] - predicted[i], 2)
    return mse

from libc.stdio cimport printf

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_distance_tads(
        np.ndarray[DTYPE_int_t, ndim=2] joint_EP not None,
        np.ndarray[DTYPE_t, ndim=2] RNA not None,
        np.ndarray[DTYPE_t, ndim=2] cStateBetas not None,
        np.ndarray[DTYPE_t, ndim=2] pStateBetas not None,
        np.ndarray[DTYPE_int_t, ndim=1] TSSs not None,
        np.ndarray[DTYPE_int_t, ndim=2] distances not None,
        np.ndarray[DTYPE_int_t, ndim=1] dist_starts not None,
        np.ndarray[DTYPE_t, ndim=1] scores not None,
        np.ndarray[DTYPE_int_t, ndim=1] paths not None,
        np.ndarray[DTYPE_t, ndim=2] predicted not None,
        np.ndarray[DTYPE_int_t, ndim=1] TSS_list not None,
        np.ndarray[DTYPE_int_t, ndim=1] CRE_list not None,
        np.ndarray[DTYPE_t, ndim=3] distCREs not None,
        np.ndarray[DTYPE_int_t, ndim=2] tads not None,
        int max_size,
        double alpha):
    cdef long long int i, j, k, l, m, i_inv, j_inv, cIndex, tIndex, dIndex
    cdef long long int offset, tPos, cPos, start, end, tadN
    cdef double mse, alphadist
    cdef long long int jointN = joint_EP.shape[0]
    cdef long long int tssN = RNA.shape[0]
    cdef long long int cellN = RNA.shape[1]
    cdef long long int maxTSS = predicted.shape[0]
    cdef long long int maxCRE = distances.shape[1]
    with nogil:
        # Initialize the path and score arrays
        paths[0] = 0
        scores[0] = 0
        for i in range(1, jointN + 1):
            scores[i] = Inf
            paths[i] = -1
        # Find distances raised to the alpha power
        for i in range(tssN):
            offset = dist_starts[i]
            for j in range(maxCRE):
                if distances[i, j] > 0:
                    alphadist = pow(distances[i, j], alpha)
                    for m in range(cellN):
                        distCREs[i, j, m] = alphadist * cStateBetas[j + offset, m]
                else:
                    break
        # Try each TSS as the first in all possible TADs including it
        for i_inv in range(tssN):
            i = TSSs[i_inv]
            # For each TSS, try each CRE in range upstream as the TAD start
            # until we reach the next TSS
            for j_inv in range(i - joint_EP[i, 3] + 1):
                j = i - j_inv
                # If the new TAD start point hasn't been a TAD endpoint,
                # find the best scoring upstream endpoint occurring after 
                # the next upstream TSS
                k = j - 1
                while k >= 0:
                    if joint_EP[k, 1] == 1:
                        break
                    if scores[k] < scores[j]:
                        scores[j] = scores[k]
                        paths[j] = k
                    k -= 1
                # Fill in the list of upstream CREs included in this TAD
                cPos = 0
                for k in range(j, i):
                    CRE_list[cPos] = joint_EP[k, 0]
                    cPos += 1
                # Test each downstream endpoint
                tPos = 0
                for k in range(i, joint_EP[j, 4]):
                    if joint_EP[k, 1] == 1:
                        # If the next element is a TSS, add it to the TSS list
                        tIndex = joint_EP[k, 0]
                        TSS_list[tPos] = tIndex
                        # Add an initial expr prediction based on the promoter
                        for m in range(cellN):
                            predicted[tPos, m] = pStateBetas[tIndex, m]
                        # Add each CRE in CRE list to the prediction
                        for l in range(cPos):
                            cIndex = CRE_list[l]
                            dIndex = cIndex - dist_starts[tIndex]
                            for m in range(cellN):
                                predicted[tPos, m] += distCREs[tIndex, dIndex, m]
                        tPos += 1
                    else:
                        # If the next element is a CRE, add it to the CRE list
                        cIndex = joint_EP[k, 0]
                        CRE_list[cPos] = cIndex
                        # For each TSS, add this CRE to the expr prediction
                        for l in range(tPos):
                            tIndex = TSS_list[l]
                            dIndex = cIndex - dist_starts[tIndex]
                            for m in range(cellN):
                                predicted[l, m] += distCREs[tIndex, dIndex, m]
                        cPos += 1
                    # Find the total error in expr prediction with this TAD
                    mse = scores[j]
                    for l in range(tPos):
                        tIndex = TSS_list[l]
                        for m in range(cellN):
                            mse += pow(RNA[tIndex, m] - predicted[l, m], 2)
                    # If error is less than another TAD ending here, replace
                    # the score and path with this one
                    if mse < scores[k + 1]:
                        scores[k + 1] = mse
                        paths[k + 1] = j
        if joint_EP[jointN - 1, 1] == 1:
            start = jointN
        else:
            mse = scores[jointN]
            start = jointN
            j = start - 1
            while j >= 0:
                if joint_EP[j, 1] == 1:
                    break
                if scores[j] < mse:
                    mse = scores[j]
                    start = j
                j -= 1
        tadN = 0
        while start > 0:
            end = start
            start = paths[end]
            tPos = 0
            for i in range(start, end):
                if joint_EP[i, 1] == 1:
                    tPos += 1
            if tPos > 0:
                tads[tadN, 0] = start
                tads[tadN, 1] = end
                tadN += 1
    return tadN


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def predict_distance_expression(
        np.ndarray[DTYPE_t, ndim=2]     RNA not None,
        np.ndarray[DTYPE_int_t, ndim=2] joint_EP not None,
        np.ndarray[DTYPE_int_t, ndim=2] tads not None,
        np.ndarray[DTYPE_int_t, ndim=2] pstates not None,
        np.ndarray[DTYPE_int_t, ndim=2] cstates not None,
        np.ndarray[DTYPE_t, ndim=2]     distAlphas not None,
        np.ndarray[DTYPE_int_t, ndim=1] dist_starts not None,
        np.ndarray[DTYPE_t, ndim=2]     betas not None,
        np.ndarray[DTYPE_t, ndim=2]     predicted not None,
        np.ndarray[DTYPE_int_t, ndim=1] TSS_list not None,
        np.ndarray[DTYPE_int_t, ndim=1] CRE_list not None):
    cdef long long int i, j, k, l, cIndex, tIndex
    cdef long long int tPos, cPos, tad_start, tad_end
    cdef double mse, distAlpha
    cdef long long int tadN = tads.shape[0]
    cdef long long int cellN = RNA.shape[1]
    cdef long long int stateN = betas.shape[1]
    with nogil:
        mse = 0
        for i in range(tadN):
            tad_start = tads[i, 0]
            tad_end = tads[i, 1]
            tPos = 0
            cPos = 0
            for j in range(tad_start, tad_end):
                if joint_EP[j, 1] == 1:
                    TSS_list[tPos] = joint_EP[j, 0]
                    tPos += 1
                else:
                    CRE_list[cPos] = joint_EP[j, 0]
                    cPos += 1
            for j in range(tPos):
                tIndex = TSS_list[j]
                for l in range(cellN):
                    predicted[tIndex, l] = betas[0, pstates[tIndex, l]]
                for k in range(cPos):
                    cIndex = CRE_list[k]
                    distAlpha = distAlphas[tIndex, cIndex - dist_starts[tIndex]]
                    for l in range(cellN):
                        predicted[tIndex, l] += betas[1, cstates[cIndex, l]] * distAlpha
                for l in range(cellN):
                    mse += pow(RNA[tIndex, l] - predicted[tIndex, l], 2)
    return mse


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def fit_distance_model(
        np.ndarray[DTYPE_t, ndim=2]     RNA not None,
        np.ndarray[DTYPE_int_t, ndim=2] joint_EP not None,
        np.ndarray[DTYPE_int_t, ndim=2] tads not None,
        np.ndarray[DTYPE_int_t, ndim=2] pstates not None,
        np.ndarray[DTYPE_int_t, ndim=2] cstates not None,
        np.ndarray[DTYPE_t, ndim=2]     ln_dist not None,
        np.ndarray[DTYPE_t, ndim=2]     distAlphas not None,
        np.ndarray[DTYPE_int_t, ndim=1] dist_starts not None,
        np.ndarray[DTYPE_t, ndim=2]     betas not None,
        np.ndarray[DTYPE_t, ndim=2]     predicted not None,
        np.ndarray[DTYPE_int_t, ndim=1] TSS_list not None,
        np.ndarray[DTYPE_int_t, ndim=1] CRE_list not None,
        np.ndarray[DTYPE_t, ndim=1]     alpha_gradient not None,
        np.ndarray[DTYPE_t, ndim=2]     beta_gradient not None):
    cdef long long int i, j, k, l, cIndex, tIndex, dIndex
    cdef long long int tPos, cPos, tad_start, tad_end, cstate
    cdef double mse, distAlpha, dC_dPred
    cdef long long int tadN = tads.shape[0]
    cdef long long int cellN = RNA.shape[1]
    with nogil:
        # Pred = sum_i(P_i * pbeta_i + sum_j(C_ij * cbeta_i * dist_j^alpha))
        # C = (Pred - expr)^2
        # dC/dPred = 2 * (pred - expr)
        # dPred/dpbeta_i = P_i
        # dPred/dcbeta_i = C_ij * dist^alpha
        # dPred/dalpha = C_ij * cbeta_i * dist_j^alpha * ln(dist_j)
        mse = 0
        for i in range(tadN):
            tad_start = tads[i, 0]
            tad_end = tads[i, 1]
            tPos = 0
            cPos = 0
            for j in range(tad_start, tad_end):
                if joint_EP[j, 1] == 1:
                    TSS_list[tPos] = joint_EP[j, 0]
                    tPos += 1
                else:
                    CRE_list[cPos] = joint_EP[j, 0]
                    cPos += 1
            for j in range(tPos):
                tIndex = TSS_list[j]
                for k in range(cPos):
                    cIndex = CRE_list[k]
                    dIndex = cIndex - dist_starts[tIndex]
                    distAlpha = distAlphas[tIndex, dIndex]
                    for l in range(cellN):
                        dC_dPred = 2 * (predicted[tIndex, l] - RNA[tIndex, l])
                        cstate = cstates[cIndex, l]
                        beta_gradient[1, cstate] += distAlpha * dC_dPred
                        alpha_gradient[0] += betas[1, cstate] * distAlpha * ln_dist[tIndex, dIndex]
                        mse += dC_dPred * dC_dPred / 4
    return mse
























