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
     
class qam:
    def __init__(self, filename, samplerate=48000, mode='r'):

        self.samplerate = samplerate
        self.carrier = 1450*11
        self.symrate = 960
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
            self.iir_i = Biquad(Biquad.LOWPASS, cutoff, samplerate, invsqr2)
            self.iir_q = Biquad(Biquad.LOWPASS, cutoff, samplerate, invsqr2)
            self.iir_i_ph = Biquad(Biquad.HIGHPASS, cutoff, samplerate, invsqr2)
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
        ph_d = []
        ph_i = []
        ph_ck = []
        si = []
        sq = []
        prev_ph_i = 0
        prev_ph_q = 0
        sample_cnt = 0
        ph_err = 0
        prev_li = 0
        prev_lq = 0

        last_clock = 0
        err_ph = 0
        i_pid = pid(0.47, 8.5e-3, 8.5e-7, self.samplerate, 0)

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
            sample_cnt += 1
            li = self.iir_i(i)
            lq = self.iir_q(q)
            
            fi.append(li)
            ph_filter_i = li - prev_li
            ph_filter_q = lq - prev_lq
            prev_li = li
            prev_lq = lq

            ph_d.append(ph_filter_i)

            phi = abs(ph_filter_i) > (0.013 - prev_ph_i * 0.003)
            phq = abs(ph_filter_q) > (0.013 - prev_ph_q * 0.003)
            ph_i.append(phi)
            if (prev_ph_i == 0 and phi == 1) and (prev_ph_q == 0 and phq == 1):
                delta_ph = self.t - last_clock
                err_ph = i_pid.update(delta_ph)
                #sys.stdout.write("%f\n" % err_ph)


            self.curr_phase += 1
            prev_ph_i = phi
            prev_ph_q = phq
            si.append(li)
            sq.append(lq)
            #err_ph = 0
            if self.curr_phase + err_ph >= self.phase_max:
                last_clock = self.t
                self.curr_phase %= self.phase_max
                i_min = len(si)/2-5
                i_max = len(si)/2+5
                l = i_max - i_min
                ii = sum(si[i_min:i_max])/l
                qq = sum(sq[i_min:i_max])/l
                #sys.stdout.write("%.4f, %.4f\n" % (si, sq))
                di += [ii]*sample_cnt
                dq += [qq]*sample_cnt
                ph_ck.append(1)
                ph_ck += [0] * (sample_cnt-1)
                si = []
                sq = []
                sample_cnt = 0

        if 1:
            pl.plot(di, dq, label='IQ', marker='o', color='b', ls='')
        else:
            pl.plot(ph_d)
            pl.plot(di)
            pl.plot(fi)
            pl.plot(ph_i)
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

