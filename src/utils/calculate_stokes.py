import numpy as np
from collections import OrderedDict


def get_stokes_3(I_i, sigmas, t, e, phi):
    if I_i.__len__() != 3:
        raise ValueError("Intensity vector is not a 3-vector")
    elif sigmas.__len__() != 3:
        raise ValueError("Uncertainties vector is not a 3-vector")
    elif t.__len__() != 3:
        raise ValueError("Throughput vector is not a 3-vector")
    elif e.__len__() != 3:
        raise ValueError("Efficiency vector is not a 3-vector")
    elif phi.__len__() != 3:
        raise ValueError("Angle vector is not a 3-vector")

    I_1, I_2, I_3 = I_i
    I_1 /= 0.5 * t[0]
    I_2 /= 0.5 * t[1]
    I_3 /= 0.5 * t[2]

    B_1 = (
        e[1] * e[2] * np.sin(np.radians(2.0 * phi[2] - 2.0 * phi[1])),
        e[0] * e[2] * np.sin(np.radians(2.0 * phi[0] - 2.0 * phi[2])),
        e[0] * e[1] * np.sin(np.radians(2.0 * phi[1] - 2.0 * phi[0])),
    )
    B_2 = (
        e[1] * np.sin(np.radians(2.0 * phi[1])) - e[2] *
        np.sin(np.radians(2.0 * phi[2])),
        e[2] * np.sin(np.radians(2.0 * phi[2])) - e[0] *
        np.sin(np.radians(2.0 * phi[0])),
        e[0] * np.sin(np.radians(2.0 * phi[0])) - e[1] *
        np.sin(np.radians(2.0 * phi[1])),
    )
    B_3 = (
        e[2] * np.cos(np.radians(2.0 * phi[2])) - e[1] *
        np.cos(np.radians(2.0 * phi[1])),
        e[0] * np.cos(np.radians(2.0 * phi[0])) - e[2] *
        np.cos(np.radians(2.0 * phi[2])),
        e[1] * np.cos(np.radians(2.0 * phi[1])) - e[0] *
        np.cos(np.radians(2.0 * phi[0])),
    )

    B = np.array([B_1, B_2, B_3])

    Omega_1 = e[0] * e[1] * np.sin(np.radians(2.0 * phi[1] - 2.0 * phi[0]))
    Omega_2 = e[1] * e[2] * np.sin(np.radians(2.0 * phi[2] - 2.0 * phi[1]))
    Omega_3 = e[0] * e[2] * np.sin(np.radians(2.0 * phi[0] - 2.0 * phi[2]))

    Omega = Omega_1 + Omega_2 + Omega_3

    B /= Omega

    I, Q, U = np.dot(B, np.array([I_1, I_2, I_3]))

    A = np.copy(B)

    A[0] /= 0.5 * t[0]
    A[1] /= 0.5 * t[1]
    A[2] /= 0.5 * t[2]

    sigma_I = (A[0][0] ** 2.0 * sigmas[0] ** 2.0) + (A[0][1] **
                                                     2.0 * sigmas[1] ** 2.0) + (A[0][2] ** 2.0 * sigmas[2] ** 2.0)
    sigma_Q = (A[1][0] ** 2.0 * sigmas[0] ** 2.0) + (A[1][1] **
                                                     2.0 * sigmas[1] ** 2.0) + (A[1][2] ** 2.0 * sigmas[2] ** 2.0)
    sigma_U = (A[2][0] ** 2.0 * sigmas[0] ** 2.0) + (A[2][1] **
                                                     2.0 * sigmas[1] ** 2.0) + (A[2][2] ** 2.0 * sigmas[2] ** 2.0)

    sigma_IQ = (A[0][0] * A[1][0] * sigmas[0] ** 2.0) + (A[0][1] * A[1]
                                                         [1] * sigmas[1] ** 2.0) + (A[0][2] * A[1][2] * sigmas[2] ** 2.0)
    sigma_IU = (A[0][0] * A[2][0] * sigmas[0] ** 2.0) + (A[0][1] * A[2]
                                                         [1] * sigmas[1] ** 2.0) + (A[0][2] * A[2][2] * sigmas[2] ** 2.0)
    sigma_QU = (A[1][0] * A[2][0] * sigmas[0] ** 2.0) + (A[1][1] * A[2]
                                                         [1] * sigmas[1] ** 2.0) + (A[1][2] * A[2][2] * sigmas[2] ** 2.0)

    return (I, Q, U), (sigma_I, sigma_Q, sigma_U, sigma_IQ, sigma_IU, sigma_QU)


def get_stokes(i, sigmas, t, e, phi):

    # t = [1.0, 1.0, 1.0, 1.0]
    # e = [1.0, 1.0, 1.0, 1.0]
    # phi = [0., 90., 45., 135.]

    if i.__len__() < 3:
        raise ValueError("Intensity vector is not a 3-vector")
    elif sigmas.__len__() < 3:
        raise ValueError("Uncertainties vector is not a 3-vector")
    elif t.__len__() < 3:
        raise ValueError("Throughput vector is not a 3-vector")
    elif e.__len__() < 3:
        raise ValueError("Efficiency vector is not a 3-vector")
    elif phi.__len__() < 3:
        raise ValueError("Angle vector is not a 3-vector")

    phi_tmp = []
    for ppp in phi:
        while ppp >= 180.0:
            ppp -= 180.0
        while ppp < 0:
            ppp += 180.0
        if ppp not in phi_tmp:
            phi_tmp.append(ppp)
    if phi_tmp.__len__() < 3:
        raise ValueError("Less than three different angles")

    i_1 = np.sum([ii * tt / (ss * ss) for ii, tt, ss in zip(i, t, sigmas)])
    i_2 = np.sum([ii * tt * ee * np.cos(np.radians(2.0 * pp)) / (ss * ss)
                  for ii, tt, ee, pp, ss in zip(i, t, e, phi, sigmas)])
    i_3 = np.sum([ii * tt * ee * np.sin(np.radians(2.0 * pp)) / (ss * ss)
                  for ii, tt, ee, pp, ss in zip(i, t, e, phi, sigmas)])

    I_k = (i_1, i_2, i_3)

    t_1 = np.sum([tt * tt / (ss * ss) for tt, ss in zip(t, sigmas)])
    t_2 = np.sum([tt * tt * ee * np.cos(np.radians(2.0 * pp)) / (ss * ss)
                  for tt, ee, pp, ss in zip(t, e, phi, sigmas)])
    t_3 = np.sum([tt * tt * ee * np.sin(np.radians(2.0 * pp)) / (ss * ss)
                  for tt, ee, pp, ss in zip(t, e, phi, sigmas)])

    t_k = (t_1, t_2, t_3)

    ee_1_1 = t_1
    ee_1_2 = t_2
    ee_1_3 = t_3
    e_1 = (1.0 / ee_1_1) * np.sqrt((ee_1_2) ** 2.0 + (ee_1_3) ** 2.0)

    ee_2_1 = t_2
    ee_2_2 = np.sum([tt * tt * ee * ee * np.cos(np.radians(2.0 * pp)) * np.cos(
        np.radians(2.0 * pp)) / (ss * ss) for tt, ee, pp, ss in zip(t, e, phi, sigmas)])
    ee_2_3 = np.sum([tt * tt * ee * ee * np.sin(np.radians(2.0 * pp)) * np.cos(
        np.radians(2.0 * pp)) / (ss * ss) for tt, ee, pp, ss in zip(t, e, phi, sigmas)])
    e_2 = (1.0 / ee_2_1) * np.sqrt((ee_2_2) ** 2.0 + (ee_2_3) ** 2.0)

    ee_3_1 = t_3
    ee_3_2 = np.sum([tt * tt * ee * ee * np.sin(np.radians(2.0 * pp)) * np.cos(
        np.radians(2.0 * pp)) / (ss * ss) for tt, ee, pp, ss in zip(t, e, phi, sigmas)])
    ee_3_3 = np.sum([tt * tt * ee * ee * np.sin(np.radians(2.0 * pp)) * np.sin(
        np.radians(2.0 * pp)) / (ss * ss) for tt, ee, pp, ss in zip(t, e, phi, sigmas)])
    e_3 = (1.0 / ee_3_1) * np.sqrt((ee_3_2) ** 2.0 + (ee_3_3) ** 2.0)

    e_k = (e_1, e_2, e_3)

    phi_1 = 0.5 * np.degrees(np.arctan2(ee_1_3, ee_1_2))
    phi_2 = 0.5 * np.degrees(np.arctan2(ee_3_2, ee_2_2))
    phi_3 = 0.5 * np.degrees(np.arctan2(ee_3_3, ee_2_3))
    phi_k = (phi_1, phi_2, phi_3)

    d2h_di2 = np.sum([tt * tt / (4.0 * ss * ss) for tt, ss in zip(t, sigmas)])
    d2h_dq2 = np.sum([tt * tt * ee * ee * np.cos(np.radians(2.0 * pp)) * np.cos(
        np.radians(2.0 * pp)) / (4.0 * ss * ss) for tt, ee, ss, pp in zip(t, e, sigmas, phi)])
    d2h_du2 = np.sum([tt * tt * ee * ee * np.sin(np.radians(2.0 * pp)) * np.sin(
        np.radians(2.0 * pp)) / (4.0 * ss * ss) for tt, ee, ss, pp in zip(t, e, sigmas, phi)])

    d2h_didq = np.sum([tt * tt * ee * np.cos(np.radians(2.0 * pp)) /
                       (4.0 * ss * ss) for tt, ee, ss, pp in zip(t, e, sigmas, phi)])
    d2h_dqdi = d2h_didq

    d2h_didu = np.sum([tt * tt * ee * np.sin(np.radians(2.0 * pp)) /
                       (4.0 * ss * ss) for tt, ee, ss, pp in zip(t, e, sigmas, phi)])
    d2h_dudi = d2h_didu

    d2h_dqdu = np.sum([tt * tt * ee * ee * np.sin(np.radians(2.0 * pp)) * np.cos(
        np.radians(2.0 * pp)) / (4.0 * ss * ss) for tt, ee, ss, pp in zip(t, e, sigmas, phi)])
    d2h_dudq = d2h_dqdu

    r_1 = (d2h_di2, d2h_didq, d2h_didu)
    r_2 = (d2h_dqdi, d2h_dq2, d2h_dqdu)
    r_3 = (d2h_dudi, d2h_dudq, d2h_du2)

    C = np.array((r_1, r_2, r_3)).astype(np.float32)
    # inverse
    C = np.linalg.inv(C)

    stokes, sigmas = get_stokes_3(I_k, (0, 0, 0), t_k, e_k, phi_k)

    I, Q, U = stokes
    sigma_I = C[0][0]
    sigma_Q = C[1][1]
    sigma_U = C[2][2]
    sigma_IQ = C[0][1]
    sigma_IU = C[0][2]
    sigma_QU = C[1][2]

    #PA and PD

    PA = 0.5 * np.degrees(np.arctan2(U, Q))
    PD = np.sqrt(Q * Q + U * U) / I * 100.0

    dpa_du = 1.0 / (2.0 * Q * (1.0 + (U * U) / (Q * Q)))
    dpa_dq = -U / (2.0 * (Q * Q + U * U))

    sigma_PA = dpa_du ** 2 * sigma_U + dpa_dq ** 2 * \
        sigma_Q + dpa_du * dpa_dq * sigma_QU
    sigma_PA = np.degrees(np.sqrt(abs(sigma_PA)))

    dpd_di = -np.sqrt(Q * Q + U * U) / (I * I)
    dpd_dq = Q / (I * np.sqrt(Q * Q + U * U))
    dpd_du = U / (I * np.sqrt(Q * Q + U * U))

    sigma_PD = dpd_dq ** 2 * sigma_Q + dpd_du ** 2 * sigma_U + dpd_di ** 2 * sigma_I +\
        2.0 * sigma_IQ * dpd_dq * dpd_di + 2.0 * sigma_IU * dpd_di * dpd_du +\
        2.0 * sigma_QU * dpd_dq * dpd_du
    sigma_PD = np.sqrt(abs(sigma_PD)) * 100.0

    sigma_I = np.sqrt(abs(sigma_I))
    sigma_Q = np.sqrt(abs(sigma_Q))
    sigma_U = np.sqrt(abs(sigma_U))
    sigma_IQ = np.sqrt(abs(sigma_IQ))
    sigma_IU = np.sqrt(abs(sigma_IU))
    sigma_QU = np.sqrt(abs(sigma_QU))

    stokes = OrderedDict()
    stokes['I'] = I
    stokes['Q'] = Q
    stokes['U'] = U
    stokes['sigma_I'] = sigma_I
    stokes['sigma_Q'] = sigma_Q
    stokes['sigma_U'] = sigma_U
    stokes['sigma_IQ'] = sigma_IQ
    stokes['sigma_IU'] = sigma_IU
    stokes['sigma_QU'] = sigma_QU

    stokes['PD'] = PD
    stokes['sigma_PD'] = sigma_PD

    stokes['PA'] = abs(PA)
    stokes['sigma_PA'] = sigma_PA

    u = U / I
    q = Q / I
    PD2 = np.sqrt(u ** 2 + q ** 2) * 100.0

    du_dU = 1.0 / I
    du_dI = -U / (I ** 2)

    dq_dQ = 1.0 / I
    dq_dI = -Q / (I ** 2)

    dPD2_du = u / np.sqrt(u ** 2 + q ** 2)
    dPD2_dq = q / np.sqrt(u ** 2 + q ** 2)

    j1 = (1, 0, 0)
    j2 = (q, I, 0)
    j3 = (u, 0, I)

    J = np.array((j1, j2, j3)).astype(np.float32)

    C2 = J.T.dot(C).dot(J)

    sigma_qu = np.sqrt(abs(C2[1][2]))

    sigma_u = (du_dU * sigma_U) ** 2 + (du_dI * sigma_I) ** 2 + \
        2.0 * du_dU * du_dI * sigma_IU ** 2
    sigma_u = np.sqrt(abs(sigma_u))
    sigma_q = (dq_dQ * sigma_Q) ** 2 + (dq_dI * sigma_I) ** 2 + \
        2.0 * dq_dQ * dq_dI * sigma_IQ ** 2
    sigma_q = np.sqrt(abs(sigma_q))
    # + 2.0 * dPD2_du * dPD2_dq * sigma_qu ** 2
    sigma_PD2 = (dPD2_du * sigma_u) ** 2 + (dPD2_dq * sigma_q) ** 2
    sigma_PD2 = np.sqrt(abs(sigma_PD2)) * 100.0

    # stokes['PD2'] = PD2
    # stokes['q'] = q
    # stokes['u'] = u
    # stokes['sigma_PD2'] = sigma_PD2
    # stokes['sigma_q'] = sigma_q
    # stokes['sigma_u'] = sigma_u

    return stokes
