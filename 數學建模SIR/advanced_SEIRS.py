"""about this module:
Parameters:
| beta
| sigma  潛伏期(E)至症狀出現前(但有傳染力)(I_pre)速率=1/潛伏期(E)時間
| lambda_ 症狀出現前(但有傳染力)(I_pre)至有症狀(I_sym)或無症狀(I_asym)速率=1/症狀出現前(但有傳染力)(I_pre)時間
| a  受感染者(I_pre)出現症狀的機率
| gamma_a 無症狀(I_asym)速率至康復(R)速率=1/無症狀(I_asym)時間
| h  有症狀者(I_sym)住院的機率
| eta  有症狀(I_sym)至住院(H)速率=1/將會住院者在有症狀(I_sym)的時間
| f_s  有症狀(I_sym)且不會住院者死亡機率
| mu_s  有症狀(I_sym)至不住院死亡(F)速率=1/不住院且會死亡者在有症狀(I_sym)的時間
| gamma_s  有症狀(I_sym)至不住院康復(R)速率=1/不住院且會康復者在有症狀(I_sym)的時間
| f_h  住院者(H)死亡機率
| mu_h 住院者(H)至死亡(F)速率=1/將會死亡者在住院(H)的時間
| gamma_h  住院者(H)至康復(R)速率=1/將會康復者在住院(H)的時間
| xi   康復(R)重回至易感染(S)速率=1/已康復期(R)時間
| 
| Q_sigma  已隔離者中，潛伏期(Q_E)至症狀出現前(但有傳染力)(Q_I_pre)速率=1/潛伏期(Q_E)時間
| Q_lambda_  已隔離者中，症狀出現前(但有傳染力)(Q_I_pre)至有症狀(Q_I_sym)或無症狀(Q_I_asym)速率=1/症狀出現前(但有傳染力)(Q_I_pre)時間
| Q_gamma_a  已隔離者中，無症狀(Q_I_asym)速率至康復(Q_R)速率=1/無症狀(Q_I_asym)時間
| Q_eta  已隔離者中，有症狀(Q_I_sym)至住院(H)速率=1/將會住院者在有症狀(Q_I_sym)的時間
| Q_mu_s  已隔離者中，有症狀(Q_I_sym)至不住院死亡(F)速率=1/不住院且會死亡者在有症狀(Q_I_sym)的時間
| Q_gamma_s  已隔離者中，有症狀(I_sym)至不住院康復(Q_R)速率=1/不住院且會康復者在有症狀(Q_I_sym)的時間
| rho  已隔離的康復者(Q_R)解放回康復者(R)的速率=1/康復者隔離(Q_R)時間
| 
| theta_E
| psi_E
| theta_I_pre
| psi_I_pre
| theta_I_asym
| psi_I_asym
| theta_I_sym
| psi_I_sym

Compartments:
| S
| E
| I_pre
| I_asym
| I_sym
| H
| F
| R
| Q_E
| Q_I_pre
| Q_I_asym
| Q_I_sym
| Q_R

Others:
| N
| day
| interval
"""
from pylab import *
from scipy import interpolate
import numpy as np
import time

# Parameters
beta = 0.6
sigma = 0.6  # 潛伏期(E)至症狀出現前(但有傳染力)(I_pre)速率=1/潛伏期(E)時間
lambda_ = 0.1  # 症狀出現前(但有傳染力)(I_pre)至有症狀(I_sym)或無症狀(I_asym)速率=1/症狀出現前(但有傳染力)(I_pre)時間
a = 0.5  # 受感染者(I_pre)出現症狀的機率
gamma_a = 0.2  # 無症狀(I_asym)速率至康復(R)速率=1/無症狀(I_asym)時間
h = 0.5  # 有症狀者(I_sym)住院的機率
eta = 0.2  # 有症狀(I_sym)至住院(H)速率=1/將會住院者在有症狀(I_sym)的時間
f_s = 0.5  # 有症狀(I_sym)且不會住院者死亡機率
mu_s = 0.2  # 有症狀(I_sym)至不住院死亡(F)速率=1/不住院且會死亡者在有症狀(I_sym)的時間
gamma_s = 0.2  # 有症狀(I_sym)至不住院康復(R)速率=1/不住院且會康復者在有症狀(I_sym)的時間
f_h = 0.5  # 住院者(H)死亡機率
mu_h = 0.2  # 住院者(H)至死亡(F)速率=1/將會死亡者在住院(H)的時間
gamma_h = 0.2  # 住院者(H)至康復(R)速率=1/將會康復者在住院(H)的時間
xi = 0.5  # 康復(R)重回至易感染(S)速率=1/已康復期(R)時間

Q_sigma = 0.6  # 已隔離者中，潛伏期(Q_E)至症狀出現前(但有傳染力)(Q_I_pre)速率=1/潛伏期(Q_E)時間
Q_lambda_ = 0.1  # 已隔離者中，症狀出現前(但有傳染力)(Q_I_pre)至有症狀(Q_I_sym)或無症狀(Q_I_asym)速率=1/症狀出現前(但有傳染力)(Q_I_pre)時間
Q_gamma_a = 0.2  # 已隔離者中，無症狀(Q_I_asym)速率至康復(Q_R)速率=1/無症狀(Q_I_asym)時間
Q_eta = 0.2  # 已隔離者中，有症狀(Q_I_sym)至住院(H)速率=1/將會住院者在有症狀(Q_I_sym)的時間
Q_mu_s = 0.2  # 已隔離者中，有症狀(Q_I_sym)至不住院死亡(F)速率=1/不住院且會死亡者在有症狀(Q_I_sym)的時間
Q_gamma_s = 0.2  # 已隔離者中，有症狀(I_sym)至不住院康復(Q_R)速率=1/不住院且會康復者在有症狀(Q_I_sym)的時間
rho = 1  # 已隔離的康復者(Q_R)解放回康復者(R)的速率=1/康復者隔離(Q_R)時間

theta_E = 0.2
psi_E = 0.2
theta_I_pre = 0.2
psi_I_pre = 0.2
theta_I_asym = 0.2
psi_I_asym = 0.2
theta_I_sym = 0.2
psi_I_sym = 0.2

S = 150
E = 0
I_pre = 0
I_asym = 0
I_sym = 135
H = 0
F = 0
R = 15
Q_E = 0
Q_I_pre = 0
Q_I_asym = 0
Q_I_sym = 0
Q_R = 0

day = 100
interval = 0.2


def runSEIRS(S, E, I_pre, I_asym, I_sym, H, F, R, Q_E, Q_I_pre, Q_I_asym, Q_I_sym, Q_R, N, dt, t, params):
    locals().update(params)  # set values of the simulation parameters to custom values


    d_S2E_dt = beta * S / N * (I_pre + I_asym + I_sym)
    d_E2I_pre_dt = sigma * E
    d_I_pre2I_asym_dt = a * lambda_ * I_pre
    d_I_pre2I_sym_dt = (1 - a) * lambda_ * I_pre
    d_I_asym2R_dt = gamma_a * I_asym
    d_I_sym2R_dt = (1 - h) * (1 - f_s) * gamma_s * I_sym
    d_I_sym2F_dt = (1 - h) * f_s * mu_s * I_sym
    d_I_sym2H_dt = h * eta * I_sym
    d_H2F_dt = f_h * mu_h * H
    d_H2R_dt = (1 - f_h) * gamma_h * H
    d_E2Q_E = theta_E * psi_E * E
    d_I_pre2Q_I_pre = theta_I_pre * psi_I_pre * I_pre
    d_I_asym2Q_I_asym = theta_I_asym * psi_I_asym * I_asym
    d_I_sym2Q_I_sym = theta_I_sym * psi_I_sym * I_sym
    d_Q_E2Q_I_pre = Q_sigma * Q_E
    d_Q_I_pre2Q_I_asym_dt = a * Q_lambda_ * Q_I_pre
    d_Q_I_pre2Q_I_sym_dt = (1 - a) * Q_lambda_ * Q_I_pre
    d_Q_I_asym2Q_R_dt = Q_gamma_a * Q_I_asym
    d_Q_I_sym2Q_R_dt = (1 - h) * (1 - f_s) * Q_gamma_s * Q_I_sym
    d_Q_I_sym2F_dt = (1 - h) * (f_s) * Q_mu_s * Q_I_sym
    d_Q_I_sym2H_dt = h * Q_eta * Q_I_sym
    d_Q_R2R_dt = rho * Q_R
    d_R2S_dt = xi * R

    S += (-d_S2E_dt + d_R2S_dt) * dt
    E += (d_S2E_dt - d_E2I_pre_dt - d_E2Q_E) * dt
    I_pre += (d_E2I_pre_dt - d_I_pre2I_asym_dt - d_I_pre2I_sym_dt - d_I_pre2Q_I_pre) * dt
    I_asym += (d_I_pre2I_asym_dt - d_I_asym2R_dt - d_I_asym2Q_I_asym) * dt
    I_sym += (d_I_pre2I_sym_dt - d_I_sym2F_dt - d_I_sym2R_dt - d_I_sym2H_dt - d_I_sym2Q_I_sym) * dt
    Q_E += (d_E2Q_E - d_Q_E2Q_I_pre) * dt
    Q_I_pre += (d_I_pre2Q_I_pre + d_Q_E2Q_I_pre - d_Q_I_pre2Q_I_asym_dt - d_Q_I_pre2Q_I_sym_dt) * dt
    Q_I_asym += (d_I_asym2Q_I_asym + d_Q_I_pre2Q_I_asym_dt - d_Q_I_asym2Q_R_dt) * dt
    Q_I_sym += (d_I_sym2Q_I_sym + d_Q_I_pre2Q_I_sym_dt - d_Q_I_sym2Q_R_dt - d_Q_I_sym2F_dt - d_Q_I_sym2H_dt) * dt
    H += (d_I_sym2H_dt + d_Q_I_sym2H_dt - d_H2F_dt - d_H2R_dt) * dt
    F += (d_H2F_dt + d_I_sym2F_dt + d_Q_I_sym2F_dt) * dt
    Q_R += (d_Q_I_asym2Q_R_dt + d_Q_I_sym2Q_R_dt - d_Q_R2R_dt) * dt
    R += (d_I_asym2R_dt + d_I_sym2R_dt + d_H2R_dt + d_Q_R2R_dt - d_R2S_dt) * dt

    def fix(n):
        if n < 0:
            return 0
        return n

    S, E, I_pre, I_asym, I_sym, H, F, R, Q_E, Q_I_pre, Q_I_asym, Q_I_sym, Q_R, N = map(fix,
                                                                                       [S, E, I_pre, I_asym, I_sym, H,
                                                                                        F, R, Q_E, Q_I_pre, Q_I_asym,
                                                                                        Q_I_sym, Q_R, N])

    N = sum((S, E, I_pre, I_asym, I_sym, H, F, R, Q_E, Q_I_pre, Q_I_asym, Q_I_sym, Q_R))
    return S, E, I_pre, I_asym, I_sym, H, F, R, Q_E, Q_I_pre, Q_I_asym, Q_I_sym, Q_R, N


def integrateSEIRS(SEIRS, total_time, interval, params):
    t = 0
    for i in range(int(total_time / interval)):
        new_SEIRS = np.array(runSEIRS(*SEIRS[-1], interval, t, params))
        # print(new_SEIRS)
        # print(SEIRS)
        SEIRS = np.vstack((SEIRS, new_SEIRS))

        t += interval

    return SEIRS


def simulateSEIRS(S=S, E=E, I_pre=I_pre, I_asym=I_asym, I_sym=I_sym, H=H, F=F, R=R, Q_E=Q_E, Q_I_pre=Q_I_pre,
                  Q_I_asym=Q_I_asym, Q_I_sym=Q_I_sym, Q_R=Q_R, day=day, interval=interval, draw_graph=False, **params):
    now = time.time()

    '''運行SIR模擬'''
    init_SEIRS = {
        'S': S,
        'E': E,
        'I_pre': I_pre,
        'I_asym': I_asym,
        'I_sym': I_sym,
        'H': H,
        'F': F,
        'R': R,
        'Q_E': Q_E,
        'Q_I_pre': Q_I_pre,
        'Q_I_asym': Q_I_asym,
        'Q_I_sym': Q_I_sym,
        'Q_R': Q_R
    }
    init_SEIRS.update({'total': sum(list(init_SEIRS.values()))})  # add 'total' as a key
    print(init_SEIRS)
    SEIRS = np.array([list(init_SEIRS.values())], dtype=float)

    SEIRS = integrateSEIRS(SEIRS, day, interval, params)
    print(SEIRS[-1])
    print('execution time: %fs' % (time.time() - now))

    if draw_graph:
        '''畫圖'''
        length = np.size(SEIRS, axis=0)
        types_amount = np.size(SEIRS, axis=1)
        t = arange(0, length, 1) * (interval)  # [0,0+interval,0+interval*2,...,day]，即每筆資料對應之時刻

        # 運用內差法
        # 繪製曲線
        dt = 0.05
        tnew = arange(0, int(day / dt) + 1, 1) * dt  # [0,0+dt,0+dt*2,...,day]

        lines = [
            plot(tnew, interpolate.InterpolatedUnivariateSpline(t, SEIRS[:, i])(tnew), label=list(init_SEIRS.keys())[i])[0]
            for i in range(types_amount)]
        print(lines)

        # 圖例
        legend(handles=lines,
               shadow=True, loc=(0.85, 0.4))  # handle

        # # 標示
        # text(16.5, 290, 'Beta=%g' % (Beta))
        # text(16.5, 270, 'Gamma=%g' % (Gamma))
        # text(16.5, 250, 'Mu=%g' % (Mu))
        # text(16.5, 230, 'Alpha=%g' % (Alpha))
        # text(16.5, 210, 'Birth Rate=%g' % (birth))
        # text(16.5, 190, 'Death Rate=%g' % (death))
        # text(16.5, 170, 'Infected_0:%d' % (I0))
        # text(16.5, 150,
        #      'Infected_MAX:%d' % (np.max(SEIRS[:, [i for i, val in enumerate(init_SEIRS.keys()) if val == "I"][0]])))
        # text(16.5, 130,
        #      'Exposed_MAX:%d' % (np.max(SEIRS[:, [i for i, val in enumerate(init_SEIRS.keys()) if val == "E"][0]])))

        v = [0, day, 0, np.max(SEIRS[:, -1]) * 1.2]  # [x刻度start,x刻度end,y刻度start,y刻度end]
        axis(v)
        xlabel('time (days)')
        ylabel('Population(people)')
        title('SEIRS Model')
        grid(True)
        show()
    return SEIRS
simulateSEIRS(draw_graph=True)