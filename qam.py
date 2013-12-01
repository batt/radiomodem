import wave
from math import sin
from math import cos
from math import pi
import struct
import sys

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
            i=0
            q=0
            for j in range(self.symlen):
                d = self.w.readframes(1)
                if not d:
                    return
                d = struct.unpack("<h", d)[0]
                d = 2 * float(d + 32768) / 65535 - 1
                i += d * sin(2 * pi * self.carrier * self.t / self.samplerate)
                q += d * cos(2 * pi * self.carrier * self.t / self.samplerate)
                self.t += 1
            sys.stdout.write("%.4f, %.4f\n" % (i, q))

if sys.argv[1] == 'w':
    q = qam("qam.wav", mode='w')
    f = open(sys.argv[2])
    for c in f:
        q.modulate(c)
else:
    q = qam("qam.wav")
    q.demodulate()

