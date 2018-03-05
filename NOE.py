import pandas as pd
import numpy as np
import multiprocessing as mp


def mT(A):
    return (A * m_p)


def reduced(m1, m2):
    return (m1 * m2 / (m1 + m2))


def r(mH, A):
    return (4 * A * m_p * mH / (A * m_p + mH) ** 2)


def F(eR, A):
    if eR == 0:
        return (1)
    else:
        qF = np.sqrt(2 * mT(A) * eR)
        cF = 1.23 * A ** (1 / 3) - 0.6
        rF = np.sqrt(cF ** 2 + 7 * ((np.pi * aF) ** 2) / 3 - 5 * sF ** 2)
        qrF = qF * rF / 0.197
        return (3 * np.exp(-(qF * sF / 0.197) ** 2 / 2) * (np.sin(qrF) - np.cos(qrF) * qrF) / qrF ** 3)
"""
def F(eR, A):
    return(1)
"""

def si(mH, si0, A, eR):
    return (si0 * (A * reduced(A * m_p, mH) * F(eR, A) / reduced(m_p, mH)) ** 2)


def f(A):
    if A == Ox:
        return (0.465)
    elif A == Si:
        return (0.289)
    elif A == Al:
        return (0.089)
    else:
        return (0.048)


def lambdainv(mH, si0, A, eR):
    return (5.62e+23 * rhoE * si(mH, si0, A, eR) * f(A) / 0.891 / A / m_p)


def lambdaeff(lambdainvSi, lambdainvOx, lambdainvAl, lambdainvFe):
    return ((lambdainvSi + lambdainvOx + lambdainvAl + lambdainvFe) ** -1)


def pA(lambdainv, leff):
    return (lambdainv * leff)


def ler(mH, A, randomCos):
    return (1 - r(mH, A) * (1 - randomCos) / 2)


def lambdadis(leff, delta):
    return (-leff * (1 + delta) * np.log(np.random.random_sample()))


def rCos():
    return (2 * np.random.random_sample() - 1)


def weight(x, delta, leff):
    return ((1 + delta) * np.exp(- delta * x / (1 + delta) / leff))


def phi():
    return (2 * np.pi * np.random.random_sample())


def r01():
    return (np.random.random_sample())


def a_sel(ra, pASi, pAOxSi, pAFe):
    if ra >= 0 and ra < pASi:
        return (Si)
    elif ra >= pASi and ra < pAOxSi:
        return (Ox)
    elif ra >= pAOxSi and ra < 1 - pAFe:
        return (Al)
    else:
        return (Fe)


def diffusion(i):
    v_ini = v[i]
    CosTheta = r1 = r01()
    eRec = 0
    l_inv_Si = lambdainv(mH, sigmap, Si, eRec)
    l_inv_Ox = lambdainv(mH, sigmap, Ox, eRec)
    l_inv_Al = lambdainv(mH, sigmap, Al, eRec)
    l_inv_Fe = lambdainv(mH, sigmap, Fe, eRec)
    leff = lambdaeff(l_inv_Si, l_inv_Ox, l_inv_Al, l_inv_Fe)
    l = lambdadis(leff, delta)
    w = wi = weight(l, delta, leff)
    ztot = l * CosTheta
    E0 = 0.5 * mH * v_ini ** 2
    p = 0
    totalleft = 1
    Ef_a = E0 * totalleft
    if ztot >= d and Ef_a >= Emin:
        return [Ef_a, CosTheta, w]
    if Ef_a < Emin:
        return [1, w]
    while Ef_a >= Emin and ztot < d and ztot > 0:
        p += 1
        ra = r01()
        l_inv_Si = lambdainv(mH, sigmap, Si, eRec)
        l_inv_Ox = lambdainv(mH, sigmap, Ox, eRec)
        l_inv_Al = lambdainv(mH, sigmap, Al, eRec)
        l_inv_Fe = lambdainv(mH, sigmap, Fe, eRec)
        leff = lambdaeff(l_inv_Si, l_inv_Ox, l_inv_Al, l_inv_Fe)
        pAOx = pA(l_inv_Ox, leff)
        pASi = pA(l_inv_Si, leff)
        pAFe = pA(l_inv_Fe, leff)
        pAOxSi = pAOx + pASi
        A = a_sel(ra, pASi, pAOxSi, pAFe)
        mHmA = mH / mT(A)
        l = lambdadis(leff, delta)
        wi = weight(l, delta, leff)
        CosXiCM = rCos()
        CosXiLab = (mHmA + CosXiCM) / np.sqrt(1 + mHmA * (2 * CosXiCM + mHmA))
        CosTheta = r1 * CosXiLab - np.sqrt(1 - r1 ** 2) * np.sqrt(1 - CosXiLab ** 2) * np.cos(phi())
        Ef_b = E0 * totalleft
        left = ler(mH, A, CosXiCM)
        totalleft *= left
        Ef_a = E0 * totalleft
        eRec = Ef_b - Ef_a
        w *= wi
        z = l * CosTheta
        ztot += z
        if ztot < 0:
            return [0, w]
        if ztot >= d and Ef_a >= Emin:
            return [Ef_a, CosTheta, w]
        if Ef_a < Emin:
            return [1, w]
        r1 = CosTheta


def diffusion_Pb(i):
    save_ini = s[i // rep]
    r1_Pb = save_ini[1]
    eRec = 0
    l_Pb = (5.62e+23 * rhoPb * si(mH, sigmap, Pb, eRec) / Pb / m_p) ** -1
    l2 = lambdadis(l_Pb, delta)
    wiPb = wiPbSum = weight(l2, delta, l_Pb)
    w_Pb = save_ini[2] * wiPb
    ztot_Pb = l2 * r1_Pb
    E0_Pb = Efa_Pb = save_ini[0]
    p_Pb = 0
    totalleft = 1
    eRecDAMIC = E0_Pb * r(mH, Si) * (1 - rCos()) / 2
    if ztot_Pb >= d_Pb and Efa_Pb >= Emin:
        return save_ini + [Efa_Pb, eRecDAMIC, w_Pb]
    if Efa_Pb < Emin:
        return [1, w_Pb]
    while Efa_Pb >= Emin and ztot_Pb < d_Pb and ztot_Pb > 0:
        p_Pb += 1
        l_Pb = (5.62e+23 * rhoPb * si(mH, sigmap, Pb, eRec) / Pb / m_p) ** -1
        l2 = lambdadis(l_Pb, delta)
        mHmPb = mH / mT(Pb)
        CosXiCM_Pb = rCos()
        wiPb = weight(l2, delta, l_Pb)
        CosXiLab_Pb = (mHmPb + CosXiCM_Pb) / np.sqrt(1 + mHmPb * (2 * CosXiCM_Pb + mHmPb))
        CosTheta_Pb = r1_Pb * CosXiLab_Pb - np.sqrt(1 - r1_Pb ** 2) * np.sqrt(1 - CosXiLab_Pb ** 2) * np.cos(phi())
        wiPbSum += wiPb
        Efb_Pb = E0_Pb * totalleft
        left = ler(mH, Pb, CosXiCM_Pb)
        totalleft *= left
        Efa_Pb = E0_Pb * totalleft
        eRec = Efb_Pb - Efa_Pb
        w_Pb *= wiPb
        z = l2 * CosTheta_Pb
        ztot_Pb += z
        if ztot_Pb < 0:
            return [0, w_Pb]
        eRecDAMIC = Efa_Pb * r(mH, Si) * (1 - rCos()) / 2
        if ztot_Pb >= d_Pb and Efa_Pb >= Emin:
            return save_ini + [Efa_Pb, eRecDAMIC, w_Pb]
        if Efa_Pb < Emin:
            return [1, w_Pb]
        r1_Pb = CosTheta_Pb


if __name__ == '__main__':
    print('Loading the velocity distribution on the Earth surface and setting the parameters ...')
    v_df = pd.read_csv('vi.csv')
    v = [v_row[0] for v_row in v_df.values]
    v_len = len(v)
    delta = float(input("Set path length modification parameter, delta = "))
    n_cores = int(input("# of cores for multiprocessing = "))
    nj = float(input("# of particles at the Earth's surface to be simulated = "))
    mH = float(input("DM mass in GeV = "))
    sigmap = float(input("DM-nucleon cross_section in ubarn = ")) * 1e-30
    rep = int(input("Repetition factor from top to the bottom of the lead shield = "))
    # atomic mass numbers
    Si = 28
    Ox = 16
    Al = 27
    Fe = 56
    Pb = 207
    Cu = 63
    # mass densities in gr/cm^3
    rhoPb = 11.34
    rhoCu = 8.96
    rhoE = 2.7
    # Nuclear form factor parameters
    aF = 0.52
    sF = 0.9
    # proton mass in GeV
    m_p = 0.938
    # DM local mass density in GeV/cm^3
    rhoDM = 0.3
    # DAMIC exposure in Kg*days
    e = 0.107
    # DAMIC detector depth in cm
    d = 350 * 30.48
    # lead shield thickness
    d_Pb = 6 * 2.54
    # Nuclear recoil energy threshold of the detector
    E_th = 5.5e-7
    # output resetting
    nElost = nUp = nElost_Pb = nUp_Pb = ws = ws_Pb = 0
    # Min energy that a DM particle needs to have to potentially trigger the DAMIC detector in GeV
    Emin = E_th / (1 - ler(mH, Si, -1))
    pool = mp.Pool(n_cores)
    s = []
    for i in range(int(nj / v_len)):
        save = pool.map(diffusion, range(v_len))
        for j in range(v_len):
            if save[j][0] == 0:
                nUp += save[j][1]
            elif save[j][0] == 1:
                nElost += save[j][1]
            else:
                ws += save[j][2]
                s.append(save[j])
        if i % 10 == 0 and i != 0:
            print(str(int(i*1000000))+' particles diffusion has been simulated')
    attenuation_factor_Earth = ws / nj
    print('Number of particles that are deflected back to atmosphere = ', nUp)
    print('Number of particles that lost a large fraction of their energy and cannot trigger the detector = ', nElost)
    print('Number of particles that reached the lead shield', len(s))
    print("Earth attenuation factor = ", attenuation_factor_Earth)
    n_Pb = rep * len(s)
    pool_Pb = mp.Pool(n_cores)
    s_Pb = []
    save_Pb = pool_Pb.map(diffusion_Pb, range(n_Pb))
    for j in range(n_Pb):
        if save_Pb[j][0] == 0:
            nUp_Pb += save_Pb[j][1]
        elif save_Pb[j][0] == 1:
            nElost_Pb += save_Pb[j][1]
        else:
            ws_Pb += save_Pb[j][5]
            s_Pb.append(save_Pb[j])
    sRec = []
    for k in range(len(s_Pb)):
        if s_Pb[k][4] >= E_th:
            sRec.append(s_Pb[k])
    nUp = rep * nUp
    nElost = rep * nElost
    ws = rep * ws
    nj *= rep
    attenuation_factor = ws_Pb / nj
    print('Number of capable DM particles = ', len(s_Pb))
    print("Total attenuation factor = ", attenuation_factor)
    print('Number of successful DM particles = ', len(sRec))
    # factor 2 in calculation of total number of events is due the Earth-shielding of DM particles entering the Earth from below the horizon
    print("Expected total number of events by DAMIC = " + format(1.46e+42 * e * 0.3 * attenuation_factor * sum(
        si(mH, sigmap, Si, sRec[i][4]) * np.sqrt(2 * sRec[i][3] / mH) * sRec[i][5] for i in
        range(len(sRec))) / ws_Pb / Si / mH / 2, '.5f'))
