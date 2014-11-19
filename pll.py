#!/usr/bin/env python
# -*- coding: utf-8 -*-

from biquad_module import Biquad

from pylab import *

from math import *

import random, re
import wave
import sys
import struct

def ntrp(x,xa,xb,ya,yb):
  return (x-xa) * (yb-ya) / (xb-xa) + ya

sample_rate = 48000.0 # sampling frequency


pll_integral = 0
old_ref = 0
pll_cf = 1450
pll_loop_gain = 0.00003
ref_sig = 0

invsqr2 = 1.0 / sqrt(2.0)

cutoff = .06 # Units Hz

loop_lowpass = Biquad(Biquad.LOWPASS,cutoff,sample_rate,invsqr2)
lock_lowpass = Biquad(Biquad.LOWPASS,cutoff,sample_rate,invsqr2)

cutoff_iq = 960
i_lowpass = Biquad(Biquad.LOWPASS,cutoff_iq,sample_rate,invsqr2)
q_lowpass = Biquad(Biquad.LOWPASS,cutoff_iq,sample_rate,invsqr2)

ta = []
da = []
db = []
di = []
dq = []

noise_level = 0 # +40 db


w = wave.open(sys.argv[1], "rb")
n = 0
while 1:
  t = n / sample_rate

  # BEGIN test signal block
  d = w.readframes(1)
  if not d:
    break
  d = struct.unpack("<h", d)[0]
  test_sig = 2 * float(d + 32768) / 65535 - 1
  #noise = (random.random() * 2 - 1) * noise_level
  #test_sig += noise
  # END test signal block
  
  # BEGIN PLL block
  pll_loop_control = test_sig * ref_sig * pll_loop_gain
  pll_loop_control = loop_lowpass(pll_loop_control)
  pll_integral += pll_loop_control / sample_rate
  ref_sig = sin(2 * pi * pll_cf * (t + pll_integral))
  ref_sig_q = cos(2 * pi * pll_cf * (t + pll_integral))

  quad_ref = (ref_sig-old_ref) * sample_rate / (2 * pi * pll_cf)
  old_ref = ref_sig
  pll_lock = lock_lowpass(-quad_ref * test_sig)
  # END PLL block

  i = test_sig * ref_sig
  q = test_sig * ref_sig_q
  i = i_lowpass(i)
  q = q_lowpass(q)

  di.append(i)
  dq.append(q)
  da.append(ref_sig)
  ta.append(t)
  db.append(pll_lock * 2)

  n += 1


plot(ta, db, label='PLL lock')
#plot(di, dq, label='IQ', marker='o', color='r', ls='')

grid(True)

legend(loc='lower right')
setp(gca().get_legend().get_texts(),fontsize=9)
locs, labels = xticks()
setp(labels,fontsize=8)
locs, labels = yticks()
setp(labels,fontsize=8)

#gcf().set_size_inches(5,3.75)

name = re.sub('.*?(\w+).*','\\1',sys.argv[0])
savefig(name+'.png')

show()
