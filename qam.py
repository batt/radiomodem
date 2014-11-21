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

class IqPoint:
    def __init__(self, i=0, q=0, avg_len=16):
        self.i = Sma(avg_len)
        self.q = Sma(avg_len)

    def avg(self):
        return self.i.avg(), self.q.avg()

    def add(self, i, q):
        self.i.add(i)
        self.q.add(q)
        return self.avg()

    def variance(self):
        return self.i.variance() + self.q.variance()

class Lock:
    def __init__(self, size=32):
        self.sma = Sma(size)
        self.cnt = 0
        self.curr_lock = 0
        self.out = 0

    def curr_lock(self):
        return self.curr_lock

    def locked(self):
        if self.cnt >= 20:
            self.out = 1
        elif self.cnt < 1:
            self.out = 0
        return self.out

    def add(self, e):
        self.sma.add(e)
        self.curr_lock = (abs(self.sma.avg()) < 1e-3) and (self.sma.variance() < 1e-4)
        if self.curr_lock:
            self.cnt += 1
        else:
            self.cnt -= 1

        self.cnt = min(max(self.cnt, 0), 64)

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
        else:
            assert(0)

        self.w = w


        self.t = 0
        self.bit_per_symbol = 3
        symbols = 1 << self.bit_per_symbol
        self.amps = [(1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1), (0,-1), (1,-1), (0,0)]
        self.curr_p = 0
        self.p = [IqPoint() for i in range(symbols)]
        self.ref = [1, 4]

        #00110000 11000011 00001100
        self.preamble = "\x30\xC3\x0C" * 7 + "\x30\xC3\x0F"
        self.preamble_sync = 0
        self.a = 0
        self.b = 0
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

    def compute_coeff(self, i, q, ref_sym):
        ir = self.amps[ref_sym][0]
        qr = self.amps[ref_sym][1]
        a = (q*qr + i*ir) / (i*i + q*q)
        b = (ir - a*i) / q
        return a, b

    def find_symbol(self, i, q):
        i = round(i, 0)
        q = round(q, 0)
        i = int(min(max(i, -1), 1))
        q = int(min(max(q, -1), 1))

        return self.amps.index((i, q))

    def read_data(self, i, q, lock):
        r0 = self.ref[0]
        r1 = self.ref[1]
        cp = self.ref[self.curr_p]
        self.p[cp].add(i, q)

        self.curr_p = 1 - self.curr_p
        if self.curr_p == 0:
            var = self.p[r0].variance() + self.p[r1].variance()
            self.preamble_sync = var < 0.01

        if lock:
            if self.preamble_sync:
                i0, q0 = self.p[r0].avg()
                i1, q1 = self.p[r1].avg()
                #check for correct preamble points (one segment is greater than the other)
                l0 = i0*i0+q0*q0
                l1 = i1*i1+q1*q1
                if l1 > l0:
                    r = r1
                    r1 = r0
                    r0 = r

                a0, b0 = self.compute_coeff(i0, q0, r0)
                a1, b1 = self.compute_coeff(i1, q1, r1)
                if self.a == 0:
                    sys.stdout.write("Preamble lock a0:%f, b0:%f, a1:%f, b1:%f\n" % (a0,b0,a1,b1))
                self.a = (a0+a1)/2
                self.b = (b0+b1)/2
        else:
            if self.a != 0:
                sys.stdout.write("Sync lost!\n")
            self.a = 0



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
        data_present = []
        edge_lock = Lock()
        ph_err = 11111+random.randint(0, 50)
        trail = 11111
        last_clock = 0

        while 1:
            if ph_err > 0:
                d = random.randint(-32768, 32767)
                ph_err -= 1
                data_present.append(0)
            else:
                d = self.w.readframes(1)
                if d:
                    d = struct.unpack("<h", d)[0]
                    data_present.append(1.5)
                else:
                    if trail:
                        d = random.randint(-32768, 32767)
                        trail -=1
                        data_present.append(0)
                    else:
                        break
                

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

                edge_lock.add(e)
                ll = edge_lock.locked()
                lock += [ll] * sample_cnt
                self.read_data(ii, qq, ll)
                sync += [self.preamble_sync] * sample_cnt
                if self.a:
                    iin = self.a*ii + self.b*qq
                    qqn = -self.b*ii + self.a*qq
                else:
                    iin = qqn = 0


                di += [iin]*sample_cnt
                dq += [qqn]*sample_cnt

                si = []
                sq = []


        if 1:
            pl.plot(di, dq, label='IQ', marker='o', color='b', ls='')
            pl.axis((-2,2,-2,2))
        else:
            pl.plot(di)
            #pl.plot(dq)
            pl.plot(fi)
            #pl.plot(fq)
            #pl.plot(ph_ck)
            pl.plot(data_present)
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

