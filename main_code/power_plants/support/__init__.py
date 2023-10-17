import numpy as np

def TURTON_PEC_calculator(S, K_1, K_2, K_3, B_1=1., B_2=1., F_M=1., P=1., C_1=0., C_2=0., C_3=0.):

    CEPCI_turton = 397
    CEPCI_curr = 567.5
    dollar_euro = 0.92

    CP_0 = np.power(10, (K_1 + K_2 * np.log10(S) + K_3 * (np.log10(S) ** 2)))
    F_P = np.power(10, (C_1 + C_2 * np.log10(P) + C_3 * (np.log10(P) ** 2)))
    C_BM = CP_0 * (B_1 + B_2 * F_M * F_P)

    return C_BM * CEPCI_curr / CEPCI_turton * dollar_euro
