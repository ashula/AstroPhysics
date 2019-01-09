#! /usr/bin/env python 
""" 
In private repository.

This program calculate the distance between two astro
    object in degrees.

    object position in astro-spherical coordinates.
        e.x. a pair of R.A.(Right Ascension) & Decl.(Declination)
    R.A. == sign(R.A.), h(hour), m(minutes), s(second)
    Decl.== sign(Decl.), d(degree), m(arc-minute), s(arc-second)
    AstroObjectCoordinate == name, RA, DECL

    func. arc_distance(AstroObject1, AstroObject2)
        AstroObject1 : name, R.A. & Decl.
        AstroObject2 : name, R.A.& Decl.

    result: distance in arc-degrees between two astro objects.
            and Spherical astronomy position angle(polar angle).

    Ver-0.02:
        modify class RA, DECL data structure
        add    class AstroObjectCoordinate
        Add    func.arc_distance()
        correct star name spels.
    Ver-0.01:
        Initial implementation to calculate arc_distance and polar angular between
            Betelgeuse and Sirius.


"""

import numpy as np

pi = np.pi
debug = 0

def sign(x):
   if x > 0:
      return +1
   elif x < 0:
      return -1
   else:
      return 0

class RA(object):
    def __init__(self, h, m, s):
        self.sign = sign(h)
        self.h    = abs(h)
        self.m    = m
        self.s    = s
        self.deg  = self.sign * 15.0 *(self.h + self.m / 60.0 + self.s / 3600.0)
        self.rad  = self.deg * 2.0 * 3.14159 / 360

class DECL(object):
    def __init__(self, d, m, s):
        self.sign = sign(d)
        self.d    = abs(d)
        self.m    = m
        self.s    = s
        self.deg  = self.sign * (self.d + self.m / 60.0 + self.s / 3600.0)
        self.rad  = self.deg * 2.0 * 3.14159 / 360

class AstroObjectCoordinate(object):
    def __init__(self, name, RA, DECL):
        self.name = name
        self.ra = RA
        self.decl = DECL


"""
function arc_distance(Star1, Star2)
    Star1, Star2 : AstroObjectCoordinate;
    
    return  
        arc_distance and plar angle in arc digree and radian.

"""


def arc_distance(Star1, Star2):
    star1_cos_decl = np.cos(Star1.decl.rad)
    star1_sin_decl = np.sin(Star1.decl.rad)

    star2_cos_decl = np.cos(Star2.decl.rad)
    star2_sin_decl = np.sin(Star2.decl.rad)

    cos_sra_bra = np.cos(Star2.ra.rad - Star1.ra.rad)
    sin_sra_bra = np.sin(Star2.ra.rad - Star1.ra.rad)

    sin_d_sin_th = np.cos(Star2.decl.rad) * sin_sra_bra
    sin_d_cos_th = np.cos(Star1.decl.rad) * np.sin(Star2.decl.rad) - star1_sin_decl * star2_cos_decl * cos_sra_bra
    cos_d = star1_sin_decl * star2_sin_decl + star1_cos_decl * star2_cos_decl * cos_sra_bra

    check = ((sin_d_sin_th)**2 +(sin_d_cos_th )**2 +(cos_d)** 2)

    if (debug!=0):
        print "((sin_d_sin_th)**2 +(sin_d_cos_th )**2 +(cos_d)** 2)=%f" % check

    tan_th = sin_d_sin_th / sin_d_cos_th
    th_rad = np.arctan(tan_th)
    th_deg = th_rad * 360 / 2 / pi
    if (th_deg < 0):
        th_deg = th_deg + 180
        th_rad = th_rad + pi

    cos_th = np.cos(th_rad)
    sin_d = np.abs(sin_d_cos_th / cos_th )
    tan_d = sin_d / cos_d
    diff_rad = np.arctan(tan_d)
    diff_deg = diff_rad * 360 / 2.0 / pi

    if (debug!=0):
        print "func.arc_distance(): sin_d_sin_th=%f" % sin_d_sin_th
        print "func.arc_distance(): sin_d_cos_th=%f" % sin_d_cos_th
        print "func.arc_distance(): cos_d       =%f" % cos_d
        print "func.arc_distance(): tan_th      =%f" % tan_th
        print "func.arc_distance(): th_rad      =%f" % th_rad
        print "func.arc_distance(): th_deg      =%f" % th_deg
        print "func.arc_distance(): cos_th      =%f" % cos_th
        print "func.arc_distance(): sin_d       =%f" % sin_d
        print "func.arc_distance(): tan_d       =%f" % tan_d
        print "func.arc_distance(): diff_rad    =%f" % diff_rad
        print "func.arc_distance(): diff_deg    =%f" % diff_deg

    return diff_deg, diff_rad, th_deg, th_rad


if __name__ == '__main__':

    betelgeuse_ra = RA(+5, 52, 30)
    betelgeuse_decl = DECL(+7, 24, 00)
 
    sirius_ra = RA(+6, 42, 54)
    sirius_decl = DECL(-16, 39, 00)

    betel_cos_decl = np.cos(betelgeuse_decl.rad)
    betel_sin_decl = np.sin(betelgeuse_decl.rad)

    sirius_cos_decl = np.cos(sirius_decl.rad)
    sirius_sin_decl = np.sin(sirius_decl.rad)

    cos_sra_bra = np.cos(sirius_ra.rad - betelgeuse_ra.rad)
    sin_sra_bra = np.sin(sirius_ra.rad - betelgeuse_ra.rad)

    sin_d_sin_th = np.cos(sirius_decl.rad) * sin_sra_bra 
    sin_d_cos_th = np.cos(betelgeuse_decl.rad) * np.sin(sirius_decl.rad) - betel_sin_decl * sirius_cos_decl * cos_sra_bra
    cos_d = betel_sin_decl * sirius_sin_decl + betel_cos_decl * sirius_cos_decl * cos_sra_bra

    check = ((sin_d_sin_th)**2 +(sin_d_cos_th )**2 +(cos_d)** 2)
    print "((sin_d_sin_th)**2 +(sin_d_cos_th )**2 +(cos_d)** 2)=%f" % check

    tan_th = sin_d_sin_th / sin_d_cos_th
    th_rad = np.arctan(tan_th)
    th_deg = th_rad * 360 / 2 / pi
    if th_deg < 0:
        th_deg = th_deg + 180
        th_rad = th_rad + pi

    cos_th = np.cos(th_rad)
    sin_d = np.abs(sin_d_cos_th / cos_th )
    tan_d = sin_d / cos_d 
    diff_rad = np.arctan(tan_d)
    diff_deg = diff_rad * 360 / 2.0 / pi

    print "sin_d_sin_th=%f" % sin_d_sin_th
    print "sin_d_cos_th=%f" % sin_d_cos_th
    print "cos_d       =%f" % cos_d
    print "tan_th      =%f" % tan_th
    print "th_ra       =%f" % th_rad
    print "th_deg      =%f" % th_deg
    print "cos_th      =%f" % cos_th
    print "sin_d       =%f" % sin_d
    print "tan_d       =%f" % tan_d

    print "diff_rad    =%f" % diff_rad
    print "diff_deg    =%f" % diff_deg

    #   Test001 function arc_distance().
    print "*** Test001 ****************"
    betelgeuse = AstroObjectCoordinate('betelgeuse',betelgeuse_ra,betelgeuse_decl)
    sirius  = AstroObjectCoordinate('sirius', sirius_ra, sirius_decl)

    (distance_degree, distance_radian, th_degree, th_radian) = arc_distance(betelgeuse, sirius)

    print "arc distance[rad] =%f" % distance_radian
    print "arc distance[deg] =%f" % distance_degree
    print "polar angle[rad]  =%f" % th_radian
    print "polar angle[deg]  =%f" % th_degree

    #   Test002 function arc_distance().
    print "*** Test002 ****************"
    betelgeuse = AstroObjectCoordinate('Betelgeuse', RA(+5, 52, 30), DECL(+7, 24, 00))
    sirius    = AstroObjectCoordinate('Sirius', RA(+6, 42, 54), DECL(-16, 39, 00))

    (distance_degree, distance_radian, th_degree, th_radian) = arc_distance(betelgeuse, sirius)

    print "arc distance[rad] =%f" % distance_radian
    print "arc distance[deg] =%f" % distance_degree
    print "polar angle[rad]  =%f" % th_radian
    print "polar angle[deg]  =%f" % th_degree

