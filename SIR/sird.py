from pylab import *
from scipy import interpolate
import numpy as np
import time

# Parameters
Beta = 0.6
Gamma = 0.1
Mu = 0.1
N = 300
day = 30
interval = 0.2
I0 = 135
R0 = 15
D0 = 0


def runSIRD(S, I, R, D, dt):
    dS = - Beta * S * I / N
    dR = Gamma * I
    dD = Mu * I
    dI = - dS - dR - dD
    S += dS * dt
    I += dI * dt
    R += dR * dt
    D += dD * dt

    if S < 0:
        S = 0
    if I > N:
        I = N
    if R > N:
        R = N
    if D > N:
        D = N

    return S, I, R, D


def integrateSIRD(SIRD, total_time, interval):
    t = 0
    for i in range(int(total_time / interval)):
        new_SIRD = np.array([[i] for i in runSIRD(*SIRD[:, -1], interval)])
        # print(new_SIRD)
        # print(SIRD)
        SIRD = np.hstack((SIRD, new_SIRD))

        t += interval

    return SIRD


if __name__ == '__main__':
    now = time.time()

    '''運行SIR模擬'''
    I = I0
    S = N - I0 - R0
    R = R0
    D = D0

    SIRD = np.array([[S], [I], [R], [D]], dtype=float)

    SIRD = integrateSIRD(SIRD, day, interval)

    '''畫圖'''
    length = np.size(SIRD, axis=1)
    t = arange(0, length, 1) * (interval)  # [0,0+interval,0+interval*2,...,day]，即每筆資料對應之時刻
    # print(id(LS))
    # print(t)
    # print(LR.size)
    # print(t.size)

    # 運用內差法
    ls = interpolate.InterpolatedUnivariateSpline(t, SIRD[0])
    li = interpolate.InterpolatedUnivariateSpline(t, SIRD[1])
    lr = interpolate.InterpolatedUnivariateSpline(t, SIRD[2])
    ld = interpolate.InterpolatedUnivariateSpline(t, SIRD[3])

    # 繪製曲線
    dt = 0.05
    tnew = arange(0, int(day / dt) + 1, 1) * dt  # [0,0+dt,0+dt*2,...,day]
    lsnew = ls(tnew)
    linew = li(tnew)
    lrnew = lr(tnew)
    ldnew = ld(tnew)
    line1, = plot(tnew, lsnew, label='S')
    line2, = plot(tnew, linew, label='I')
    line3, = plot(tnew, lrnew, label='R')
    line4, = plot(tnew, ldnew, label='D')

    # 圖例
    legend(handles=[line1, line2, line3, line4],
           shadow=True, loc=(0.85, 0.4))  # handle

    # 標示
    text(16.5, 290, 'Beta=%g' % (Beta))
    text(16.5, 240, 'Gamma=%g' % (Gamma))
    text(16.5, 190, 'Mu=%g' % (Mu))
    text(16.5, 140, 'Infected:%d' % (I0))

    v = [0, day, 0, N]  # [x刻度start,x刻度end,y刻度start,y刻度end]
    axis(v)
    xlabel('time (days)')
    ylabel('Population(people)')
    title('SIR Model')
    grid(True)
    print('execution time: %fs' % (time.time() - now))
    show()