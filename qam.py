import wave
from math import sin
from math import cos
from math import pi
import struct
import sys



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
     
     
class qam:
    def __init__(self, filename, samplerate=48000, mode='r'):
        w = wave.open(filename, mode + "b")
        if mode == 'w':
            w.setnchannels(1)
            w.setframerate(samplerate)
            w.setsampwidth(2)
        elif mode == 'r':
            assert(w.getnchannels() == 1)
            assert(w.getframerate() == samplerate)
            assert(w.getsampwidth() == 2)
            self.iir_i = iir()
            self.iir_q = iir()
        else:
            assert(0)

        self.w = w

        self.samplerate = samplerate
        self.carrier = 14500
        self.symrate = 9600
        self.symlen = self.samplerate/self.symrate
        self.t = 0
        self.bit_per_symbol = 3
        symbols = 1 << self.bit_per_symbol;
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
        while 1:
            d = self.w.readframes(1)
            if not d:
                return
            d = struct.unpack("<h", d)[0]
            d = 2 * float(d + 32768) / 65535 - 1
            i = d * sin(2 * pi * self.carrier * self.t / self.samplerate)
            q = d * cos(2 * pi * self.carrier * self.t / self.samplerate)
            self.t += 1
            self.iir_i.push(i)
            self.iir_q.push(q)
            sys.stdout.write("%.4f, %.4f\n" % (self.iir_i.value(), self.iir_q.value()))

if sys.argv[1] == 'w':
    q = qam("qam.wav", mode='w')
    f = open(sys.argv[2])
    for c in f:
        q.modulate(c)
else:
    q = qam("qam.wav")
    q.demodulate()

