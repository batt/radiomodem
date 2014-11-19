import wave
from math import sin, cos, pi, sqrt
import struct
import sys
import pylab as pl
from biquad_module import Biquad
import random


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

def time_error_detector(d0, d1, d2):
    return (d2 - d0) * d1

class Sma:
    def __init__(self, size):
        self.len = size
        self.reset()

    def add(self, p):
        p2 = p*p
        self.sum += p
        self.sum2 += p2
        self.points.append(p)
        self.points2.append(p2)
        if len(self.points) == self.len + 1:
            self.sum -= self.points.pop(0)
            self.sum2 -= self.points2.pop(0)

        return self.sum / len(self.points)

    def avg(self):
        return self.sum / len(self.points)

    def variance(self):
        l = len(self.points)
        return self.sum2/l - (self.sum/l)**2

    def reset(self):
        self.points = []
        self.points2 = []
        self.sum = 0.0
        self.sum2 = 0.0

class qam:
    def __init__(self, filename, samplerate=48000, mode='r'):

        self.samplerate = samplerate
        self.carrier = 1450*11
        self.symrate = 960*2
        self.symlen = self.samplerate/self.symrate
        w = wave.open(filename, mode + "b")
        if mode == 'w':
            snr = 10 #dB
            self.att = 10**(snr/20.0)

            w.setnchannels(1)
            w.setframerate(samplerate)
            w.setsampwidth(2)
        elif mode == 'r':
            assert(w.getnchannels() == 1)
            assert(w.getframerate() == samplerate)
            assert(w.getsampwidth() == 2)
            invsqr2 = 1.0 / sqrt(2.0)
            cutoff = self.symrate
            self.iir_i = Biquad(Biquad.LOWPASS, cutoff, samplerate, invsqr2)
            self.iir_q = Biquad(Biquad.LOWPASS, cutoff, samplerate, invsqr2)
            self.curr_phase = 0
            self.phase_max = self.symlen
            self.curr_p = 0
            self.p = []
            self.p.append({'i':Sma(16), 'q':Sma(16)})
            self.p.append({'i':Sma(16), 'q':Sma(16)})
            self.preamble_sync = 0
            self.a = 0
            self.b = 0
            self.c = 0
            self.d = 0

        else:
            assert(0)

        self.w = w


        self.t = 0
        self.bit_per_symbol = 3
        symbols = 1 << self.bit_per_symbol
        self.amps = [(1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1), (0,-1), (1,-1)]
        self.ref0 = self.amps[1]
        self.ref1 = self.amps[4]
        #00110000 11000011 00001100
        self.preamble = "\x30\xC3\x0C" * 8 + "\xff\xff\xff"
        self.symbols = symbols
        self.carry = 0
        self.carry_len = 0

    def modulate_symbol(self, l):
        assert(l < self.symbols)
        #sys.stdout.write("%d " % l)
        for j in range(self.symlen):
            i = self.amps[l][0] * sin(2 * pi * self.carrier * self.t / self.samplerate)
            q = self.amps[l][1] * cos(2 * pi * self.carrier * self.t / self.samplerate)
            i += random.randrange(-1,1)/self.att
            q += random.randrange(-1,1)/self.att
            i = max(min(i,1), -1)
            q = max(min(q,1), -1)

            self.t += 1

            c = (i + q) / 2
            c = int((c + 1) / 2 * 65535) - 32768
            c = struct.pack("<h", c)
            self.w.writeframes(c)

    def modulate(self, data):
        data = self.preamble + data
        for c in data:
            for i in range(8):
                self.carry <<= 1
                self.carry |= ((ord(c) >> (7-i)) & 1)
                self.carry_len += 1
                if self.carry_len == self.bit_per_symbol:
                    self.modulate_symbol(self.carry)
                    self.carry_len = 0
                    self.carry = 0


    def get_iq(self, d):
        i = d * sin(2 * pi * self.carrier * self.t / self.samplerate)
        q = d * cos(2 * pi * self.carrier * self.t / self.samplerate)
        self.t += 1
        i = self.iir_i(i)
        q = self.iir_q(q)
        return i, q

    def compute_coeff(self, i0, i1, q0, q1, r0, r1):
        b = (i0*r1 - r0*i1) / (q1*i0 - q0*i1)
        a = (r0 - b*q0) / i0
        return a, b

    def read_data(self, i, q, lock):
        self.p[self.curr_p]['i'].add(i)
        self.p[self.curr_p]['q'].add(q)
        self.curr_p = 1 - self.curr_p
        if self.curr_p == 0:
            var = self.p[0]['i'].variance() + self.p[0]['q'].variance() + \
                  self.p[1]['i'].variance() + self.p[1]['q'].variance()

            self.preamble_sync = var < 0.01

        if lock and self.preamble_sync:
            i0 = self.p[0]['i'].avg()
            i1 = self.p[1]['i'].avg()
            q0 = self.p[0]['q'].avg()
            q1 = self.p[1]['q'].avg()
            #check for correct preamble points (one segment is greater than the other)
            l0 = i0*i0+q0*q0
            l1 = i1*i1+q1*q1
            if l1 > l0:
                i = i1
                i1 = i0
                i0 = i

                q = q1
                q1 = q0
                q0 = q

            self.a, self.b = self.compute_coeff(i0, i1, q0, q1, self.ref0[0], self.ref1[0])
            self.c, self.d = self.compute_coeff(i0, i1, q0, q1, self.ref0[1], self.ref1[1])



    def demodulate(self):
        di = []
        dq = []
        fi = []
        fq = []
        ph_ck = []
        si = []
        sq = []
        lock = []
        sync = []
        err_sma = Sma(32)
        ph_err = 1111+random.randint(0, 50)

        last_clock = 0

        while 1:
            if ph_err > 0:
                d = random.randint(-32768, 32767)
                ph_err -= 1
            else:
                d = self.w.readframes(1)
                if not d:
                    break
                d = struct.unpack("<h", d)[0]

            d = 2 * float(d + 32768) / 65535 - 1
            i, q = self.get_iq(d)

            fi.append(i)
            fq.append(q)

            self.curr_phase += 1
            si.append(i)
            sq.append(q)

            if self.curr_phase >= self.phase_max:
                self.curr_phase %= self.phase_max
                sample_cnt = self.t - last_clock
                last_clock = self.t
                idx2 = self.symlen/2
                idx4 = self.symlen/4

                ii = sum(si[idx2-4:idx2+4]) / 8
                qq = sum(sq[idx2-4:idx2+4]) / 8

                ph_ck.append(0.4)
                ph_ck += [0] * (sample_cnt-1)

                ei = time_error_detector(si[idx4]-ii, si[idx2]-ii, si[-idx4]-ii)
                eq = time_error_detector(sq[idx4]-qq, sq[idx2]-qq, sq[-idx4]-qq)
                e = ei + eq
                if e < 0:
                    self.curr_phase += 1
                if e > 0:
                    self.curr_phase -= 1

                err_sma.add(abs(e))
                ll = err_sma.variance() < 5e-6
                lock += [ll] * sample_cnt
                self.read_data(ii, qq, ll)
                sync += [self.preamble_sync] * sample_cnt
                if self.a:
                    iin = self.a*ii + self.b*qq
                    qqn = self.c*ii + self.d*qq
                else:
                    iin = qqn = 0

                di += [iin]*sample_cnt
                dq += [qqn]*sample_cnt

                si = []
                sq = []


        if 1:
            pl.plot(di, dq, label='IQ', marker='o', color='b', ls='')
        else:
            pl.plot(di)
            #pl.plot(dq)
            pl.plot(fi)
            #pl.plot(fq)
            pl.plot(ph_ck)
            pl.plot(lock)
            pl.plot(sync)
        pl.grid(True)
        pl.show()

if sys.argv[1] == 'w':
    q = qam("qam.wav", mode='w')
    f = open(sys.argv[2]).read()
    q.modulate(f)
else:
    q = qam("qam.wav")
    q.demodulate()

