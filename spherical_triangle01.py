#! /usr/bin/env python 
""" 
This program calculate the distance between two astro
    object in degrees.

    object position in spherical coordinates.
        e.x. a pair of R.A.(Right Ascension) & Decl.(Declination)
    R.A. == sign(R.A.), h(hour), m(minutes), s(second)
    Decl.== sign(Decl.), d(degree), m(arc-minute), s(arc-second)
    object1 : R.A. & Decl.
    object2 : R.A.& Decl.

    result: distance in arc-degrees between two astro objects.
            and Spherical astronomy position angle

"""

import numpy as np

pi = np.pi

def sign(x):
   if x > 0:
      return +1
   elif x < 0:
      return -1
   else:
      return 0


class RA(object):
    def __init__(self, name, sgn, h, m, s):
        self.name = name
        self.sign = sign(sgn)
        self.h    = h
        self.m    = m
        self.s    = s
        self.deg  = self.sign * 15.0 *(self.h + self.m / 60.0 + self.s / 3600.0)
        self.rad  = self.deg * 2.0 * 3.14159 / 360

class DECL(object):
    def __init__(self, name, sgn, d, m, s):
        self.name = name
        self.sign = sign(sgn)
        self.d    = d
        self.m    = m
        self.s    = s
        self.deg  = self.sign * (self.d + self.m / 60.0 + self.s / 3600.0)
        self.rad  = self.deg * 2.0 * 3.14159 / 360

if __name__ == '__main__':

    betelguse_ra = RA('betelguse_ra', +1, 5, 52, 30)
    betelguse_decl = DECL('betelguse_decl',+1, 7, 24, 00)
 
    sirius_ra = RA('sirius_ra', +1, 6, 42, 54)
    sirius_decl = DECL('sirius_decl', -1, 16, 39, 00)

    betel_cos_decl = np.cos(betelguse_decl.rad)
    betel_sin_decl = np.sin(betelguse_decl.rad)

    sirius_cos_decl = np.cos(sirius_decl.rad)
    sirius_sin_decl = np.sin(sirius_decl.rad)

    cos_sra_bra = np.cos(sirius_ra.rad - betelguse_ra.rad)
    sin_sra_bra = np.sin(sirius_ra.rad - betelguse_ra.rad)

    sin_d_sin_th = np.cos(sirius_decl.rad) * sin_sra_bra 
    sin_d_cos_th = np.cos(betelguse_decl.rad) * np.sin(sirius_decl.rad) - betel_sin_decl * sirius_cos_decl * cos_sra_bra
    cos_d = betel_sin_decl * sirius_sin_decl + betel_cos_decl * sirius_cos_decl * cos_sra_bra

    check = ((sin_d_sin_th)**2 +(sin_d_cos_th )**2 +(cos_d)** 2)
    print "((sin_d_sin_th)**2 +(sin_d_cos_th )**2 +(cos_d)** 2)=%f" % check

    tan_th = sin_d_sin_th / sin_d_cos_th
    th = np.arctan(tan_th)
    th_deg = th * 360 / 2 / pi
    if th_deg < 0:
      th_deg = th_deg + 180

    cos_th = np.cos(th)
    sin_d = np.abs(sin_d_cos_th / cos_th )
    tan_d = sin_d / cos_d 
    diff_rad = np.arctan(tan_d)
    diff_deg = diff_rad * 360 / 2.0 / pi

    print "sin_d_sin_th=%f" % sin_d_sin_th
    print "sin_d_cos_th=%f" % sin_d_cos_th
    print "cos_d       =%f" % cos_d
    print "tan_th      =%f" % tan_th
    print "th          =%f" % th
    print "th_deg      =%f" % th_deg
    print "cos_th      =%f" % cos_th
    print "sin_d       =%f" % sin_d
    print "tan_d       =%f" % tan_d

    print "diff_rad    =%f" % diff_rad
    print "diff_deg    =%f" % diff_deg
