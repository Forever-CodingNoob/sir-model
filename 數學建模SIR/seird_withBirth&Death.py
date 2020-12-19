from pylab import *
from scipy import interpolate
import numpy as np
import time
import math

# Parameters
Beta = 0.6  # beta
Gamma = 0.1  # gamma
Mu = 0.1  # 疾病致死率
Alpha = 1 / 7  # 1/潛伏期
N = 300
day = 100
interval = 0.2
I0 = 135
E0 = 0
R0 = 15
D0 = 0
S0 = N - I0 - E0 - R0 - D0
birth = 0.5 # 某年出生率 = 隔年新生人口數 / 該年總人數 = (隔年總人數-該年總人數) / 該年總人數 (若不考慮死亡率)
death = 0.5 # 某年死亡率 = 隔年死亡人口數 / 該年總人數 = (該年總人數-隔年總人數) / 該年總人數 (若不考慮出生率)


init_SEIRD = {'S':S0,'E':E0,'I':I0,'R':R0,'D':D0,'total':N}

def runSEIRD(S, E, I, R, D, N, dt, t):
    dBorndt = lambda population: (1/365)*population*((birth+1)**0)*math.log(birth+1)
    dDeaddt = lambda population: -(1/365)*population*((1-death)**0)*math.log(1-death)
    d_S2E_dt = Beta * S * I / N
    d_I2R_dt = Gamma * I
    d_I2D_dt = Mu * I
    d_E2I_dt = Alpha*E

    dSdt = dBorndt(S+E+I+R) - dDeaddt(S) - d_S2E_dt
    dEdt = d_S2E_dt - d_E2I_dt - dDeaddt(E)
    dIdt = d_E2I_dt - d_I2R_dt - d_I2D_dt - dDeaddt(I)
    dRdt = d_I2R_dt - dDeaddt(R)
    dDdt = d_I2D_dt + dDeaddt(S+E+I+R)

    S += dSdt * dt
    E += dEdt * dt
    I += dIdt * dt
    R += dRdt * dt
    D += dDdt * dt

    # S,I,R,D=map(int,[S,I,R,D]) # turn S,I,R,D into integer

    def fix(n):
        if n<0:
            return 0
        return n
    S,E,I,R,D=map(fix,[S,E,I,R,D])

    if I<=0:
        print('disease gone!')

    N = S+E+I+R+D
    return S, E, I, R, D, N


def integrateSEIRD(SEIRD, total_time, interval):
    t = 0
    for i in range(int(total_time / interval)):
        new_SEIRD = np.array(runSEIRD(*SEIRD[-1], interval, t))
        # print(new_SEIRD)
        # print(SEIRD)
        SEIRD = np.vstack((SEIRD, new_SEIRD))

        t += interval

    return SEIRD


if __name__ == '__main__':
    now = time.time()

    '''運行SIR模擬'''
    SEIRD = np.array([list(init_SEIRD.values())], dtype=float)

    SEIRD = integrateSEIRD(SEIRD, day, interval)
    print(SEIRD[-1], N)

    '''畫圖'''
    length = np.size(SEIRD, axis=0)
    types_amount = np.size(SEIRD, axis=1)
    t = arange(0, length, 1) * (interval)  # [0,0+interval,0+interval*2,...,day]，即每筆資料對應之時刻
    # print(id(LS))
    # print(t)
    # print(LR.size)
    # print(t.size)

    # 運用內差法
    # 繪製曲線
    dt = 0.05
    tnew = arange(0, int(day / dt) + 1, 1) * dt  # [0,0+dt,0+dt*2,...,day]

    lines = [plot(tnew, interpolate.InterpolatedUnivariateSpline(t, SEIRD[:,i])(tnew), label=list(init_SEIRD.keys())[i])[0] for i in range(types_amount)]
    print(lines)

    # 圖例
    legend(handles=lines,
           shadow=True, loc=(0.85, 0.4))  # handle

    # 標示
    text(16.5, 290, 'Beta=%g' % (Beta))
    text(16.5, 270, 'Gamma=%g' % (Gamma))
    text(16.5, 250, 'Mu=%g' % (Mu))
    text(16.5, 230, 'Alpha=%g' % (Alpha))
    text(16.5, 210, 'Infected:%d' % (I0))

    v = [0, day, 0, np.max(SEIRD[:,-1])*1.2]  # [x刻度start,x刻度end,y刻度start,y刻度end]
    axis(v)
    xlabel('time (days)')
    ylabel('Population(people)')
    title('SEIRD Model')
    grid(True)
    print('execution time: %fs' % (time.time() - now))
    show()