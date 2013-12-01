import wave
from math import sin
from math import cos
from math import pi
import struct

class qam:
    def __init__(self, filename, samplerate=44100):
        w = wave.open(filename, "wb")
        w.setnchannels(1)
        w.setframerate(samplerate)
        w.setsampwidth(2)
        self.w = w

        self.samplerate = samplerate
        self.carrier = 21000
        self.symrate = 14700
        self.symlen = self.samplerate/self.symrate
        self.t = 0
        self.bit_per_level = 3
        levels = 1 << self.bit_per_level;
        self.amps = [(1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1), (0,-1), (1,-1)]
        self.levels = levels
        self.carry = 0
        self.carry_len = 0

    def _mod(self, l):
        assert(l < self.levels)

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
                self.carry <<= 1
                self.carry |= (ord(c) >> i) & 1
                self.carry_len += 1
                if self.carry_len == self.bit_per_level:
                    self._mod(self.carry)
                    self.carry_len = 0
                    self.carry = 0

q = qam("qam.wav")
f = open("gpl-3.0.txt")
for c in f:
    q.modulate(c)

print "done."
