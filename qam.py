import wave
from math import sin, cos, pi, sqrt
import struct
import sys
import pylab as pl
from biquad_module import Biquad


class iir:
    def __init__(self):
        self.xv = [0] * 5
        self.yv = [0] * 5

    def push(self, val):
        self.xv[0] = self.xv[1]; self.xv[1] = self.xv[2]; self.xv[2] = self.xv[3]; self.xv[3] = self.xv[4]
        self.xv[4] = val / 9.794817390e+01
        self.yv[0] = self.yv[1]; self.yv[1] = self.yv[2]; self.yv[2] = self.yv[3]; self.yv[3] = self.yv[4]
        self.yv[4] = (self.xv[0] + self.xv[4]) + 4 * (self.xv[1] + self.xv[3]) + 6 * self.xv[2] \
                     + ( -0.1203895999 * self.yv[0]) + (  0.7244708295 * self.yv[1]) \
                     + ( -1.7358607092 * self.yv[2]) + (  1.9684277869 * self.yv[3])
        return self.yv[4]

    def value(self):
        return self.yv[4]


class pid:
    def __init__(self, kp, ki, kd, sample_rate, target=0):
        self.i = 0
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.target = target
        self.prev_err = 0
        self.sample_rate = sample_rate

    def set(self, target):
        self.target = target

    def update(self, measure):
        err = self.target - measure
        p = err * self.kp
        self.i += err / self.sample_rate
        i = self.i * self.ki
        d = (err - self.prev_err) * self.sample_rate * self.kd
        self.prev_err = err
        return p + i + d

def phase_detector(d0, d1, d2):
    return (d2 - d0) * d1

class qam:
    def __init__(self, filename, samplerate=48000, mode='r'):

        self.samplerate = samplerate
        self.carrier = 1450*11
        self.symrate = 960*2
        self.symlen = self.samplerate/self.symrate
        w = wave.open(filename, mode + "b")
        if mode == 'w':
            w.setnchannels(1)
            w.setframerate(samplerate)
            w.setsampwidth(2)
        elif mode == 'r':
            assert(w.getnchannels() == 1)
            assert(w.getframerate() == samplerate)
            assert(w.getsampwidth() == 2)
            invsqr2 = 1.0 / sqrt(2.0)
            cutoff = self.symrate
            self.iir_i = Biquad(Biquad.LOWPASS, cutoff*2, samplerate, invsqr2)
            self.iir_q = Biquad(Biquad.LOWPASS, cutoff*2, samplerate, invsqr2)
            self.curr_phase = 0
            self.phase_max = self.symlen
        else:
            assert(0)

        self.w = w


        self.t = 0
        self.bit_per_symbol = 3
        symbols = 1 << self.bit_per_symbol
        self.amps = [(1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1), (0,-1), (1,-1)]
        self.symbols = symbols
        self.carry = 0
        self.carry_len = 0

    def modulate_symbol(self, l):
        assert(l < self.symbols)

        for j in range(self.symlen):
            i = self.amps[l][0] * sin(2 * pi * self.carrier * self.t / self.samplerate)
            q = self.amps[l][1] * cos(2 * pi * self.carrier * self.t / self.samplerate)
            self.t += 1

            c = (i + q) / 2
            c = int((c + 1) / 2 * 65535) - 32768
            c = struct.pack("<h", c)
            self.w.writeframes(c)

    def modulate(self, data):
        for c in data:
            for i in range(8):
                self.carry >>= 1
                self.carry |= ((ord(c) >> i) & 1) << (self.bit_per_symbol - 1)
                self.carry_len += 1
                if self.carry_len == self.bit_per_symbol:
                    self.modulate_symbol(self.carry)
                    self.carry_len = 0
                    self.carry = 0


    def demodulate(self):
        di = []
        dq = []
        fi = []
        ph_ck = []
        si = []
        sq = []
        ph_err = 0

        last_clock = 0

        while 1:
            if ph_err > 0:
                d = 0
                ph_err -= 1
            else:
                d = self.w.readframes(1)
                if not d:
                    break
                d = struct.unpack("<h", d)[0]

            d = 2 * float(d + 32768) / 65535 - 1
            i = d * sin(2 * pi * self.carrier * self.t / self.samplerate)
            q = d * cos(2 * pi * self.carrier * self.t / self.samplerate)
            self.t += 1
            li = self.iir_i(i)
            lq = self.iir_q(q)

            fi.append(li)

            self.curr_phase += 1
            si.append(li)
            sq.append(lq)

            if self.curr_phase >= self.phase_max:
                self.curr_phase %= self.phase_max
                sample_cnt = self.t - last_clock
                last_clock = self.t
                i_min = len(si)/2-8
                i_max = len(si)/2+8
                assert(i_min >= 0)
                assert(i_max < len(si))

                l = i_max - i_min
                ii = sum(si[i_min:i_max])/l
                qq = sum(sq[i_min:i_max])/l
                di += [ii]*sample_cnt
                dq += [qq]*sample_cnt
                ph_ck.append(0.4)
                ph_ck += [0] * (sample_cnt-1)
                g = len(si)/4
                ei = phase_detector(si[g]-ii, si[len(si)/2]-ii, si[-g]-ii)
                eq = phase_detector(sq[g]-qq, sq[len(sq)/2]-qq, sq[-g]-qq)
                e = ei + eq
                if e < 0:
                    self.curr_phase += 1
                if e > 0:
                    self.curr_phase -= 1

                si = []
                sq = []


        if 1:
            pl.plot(di, dq, label='IQ', marker='o', color='b', ls='')
        else:
            pl.plot(di)
            pl.plot(fi)
            pl.plot(ph_ck)
        pl.grid(True)
        pl.show()

if sys.argv[1] == 'w':
    q = qam("qam.wav", mode='w')
    f = open(sys.argv[2])
    for c in f:
        q.modulate(c)
else:
    q = qam("qam.wav")
    q.demodulate()

