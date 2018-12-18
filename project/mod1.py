# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 10:58:03 2018

@author: aleksandra

"""

import numpy as np


kB = 1.380648*1e-23  # Дж/K = кг*м^2/с^2/К
aem = 1.660539040*1e-27  # кг
m = 1.00784*aem  # кг
ns = 0.5*0.14*1e6  # м^-3
Ts = 6500  # K
G = 6.67384*1e-11  # м^3*кг^-1*с^-2
Msol = 1.9885*1e30  # кг
V = 2.64*1e4  # V*, м/c
L = 149597870700.0  # L*, м
norm = ns*(m/(2*np.pi*kB*Ts))**(3/2)  # нормировочный множитель ф-и Максвелла
e = 1.6021766208*1e-19  # Кл
me = 5.4857990907*1e-4*aem  # кг
ligsp = 299792458  # м/с
fos = 0.416  # сила осциллятора

# в системе сгс
ligsp_cgs = ligsp*1e2  # см/c
e_cgs = 4.80320427*1e-10  # Фр = гр^(1/2)*см^(3/2)*с^(-1)
me_cgs = me*1e3  # гр
sigtot = fos*np.pi*e_cgs**(2)/me_cgs/ligsp_cgs  # см^2/c


def matrix3(r, v, alph):
    r_abs = np.linalg.norm(r)
    ez = np.cross(r, v)
    ez_abs = np.linalg.norm(ez)
    ez /= ez_abs
    vrv = np.cross(v, np.cross(r, v))
    anr = alph*r/r_abs
    ex = vrv - anr
    ex_abs = np.linalg.norm(ex)
    ex /= ex_abs
    ey = np.cross(ez, ex)
    ey_abs = np.linalg.norm(ey)
    ey /= ey_abs
    matr = np.stack([ex, ey, ez])
    return matr


def distrfunc2(r, v, Vs, mu, beta, Rs):
    """
    Calculate the distribution function of neutral hydrogen within the 
    heliosphere
    r&v are position and velocity vectors in Cartesian coordinate system with 
    Sun at the center
    Vs is velocity vector of local ISM
    mu is gravitational force to radiative pressure ratio
    
    Rs is the radius of the heliosphere    
    """
    alph = G*Msol*(1-mu)/(V)**2/L
    matr = matrix3(r, v, alph)
    r = matr@r
    v = matr@v
    Vs = matr@Vs
    rabs = np.sqrt(np.inner(r, r))
    E = np.inner(v, v)/2 - alph/rabs
    M = np.cross(r, v)
    e = np.sqrt(1 + 2*E*np.inner(M, M)/alph**2)
    a = abs(alph/2/E)
    vabs = np.linalg.norm(v)  # скорость частицы на сфере
    if (vabs**2/2 - alph/rabs + alph/Rs) <= 0:
        return 0
    if E < 0:
        vs = np.sqrt(2*(vabs**2/2 - alph/rabs + alph/Rs))
        ksis = - np.arccos((1 - Rs/a)/e)
        ksir = np.sign(r[1])*np.arccos(r[0]/a+e)
        tggam = - np.sqrt(1 - e**2)/np.tan(ksis)
        us3 = - np.sign(ksis)*vs/np.sqrt(1 + tggam**2)
        vs3 = - np.sign(ksis)*tggam*vs/np.sqrt(1 + tggam**2)
        ws3 = 0
        v = np.array([us3, vs3, ws3])
#        par = (Vs[0] - v[0])**2 + (Vs[1] - v[1])**2 + (Vs[2] - v[2])**2
        par = np.inner(Vs-v, Vs-v)
        fm = norm*np.exp(-m*par*V**2/2/kB/Ts)  # значение ф-и на границе
        xs = a*(np.cos(ksis) - e)
        ys = a*np.sqrt(1 - e**2)*np.sin(ksis)
        rs = np.array([xs, ys])
        xr = a*(np.cos(ksir) - e)
        yr = a*np.sqrt(1 - e**2)*np.sin(ksir)
    else:
        vs = np.sqrt(2*(vabs**2/2 - alph/rabs + alph/Rs))
        ksis = - np.arccosh(((Rs + a*np.sign(alph))/a/e))
        ksir = np.arcsinh(r[1]/a/np.sqrt(e**2 - 1))
        tggam = -np.sign(alph)*np.sqrt(e**2 - 1)*np.cosh(ksis)/np.sinh(ksis)
        us3 = np.sign(tggam)*vs/np.sqrt(1 + tggam**2)
        vs3 = abs(tggam)*vs/np.sqrt(1 + tggam**2)
        ws3 = 0
        v = np.array([us3, vs3, ws3])
        par = np.inner(Vs-v, Vs-v)
        fm = norm*np.exp(-m*par*V**2/2/kB/Ts)
        xs = a*(e - np.sign(alph)*np.cosh(ksis))
        ys = a*np.sqrt(e**2 - 1)*np.sinh(ksis)
        rs = np.array([xs, ys])
        xr = a*(e - np.sign(alph)*np.cosh(ksir))
        yr = a*np.sqrt(e**2 - 1)*np.sinh(ksir)
    r = np.array([xr, yr])
    tets = polarangle(rs)
    tetr = polarangle(r)
    M0 = np.sqrt(np.inner(M, M))
    if tetr > np.pi:
        extin = abs(tetr - tets)/M0
    else:
        extin = (2*np.pi - abs(tetr - tets))/M0
    return fm*np.exp(-beta*L/V*extin)


def polarangle(r):
    """
    Calculate the polar angle for given Cartesian coordinates
    r is position vector in Cartesian coordinate system
    """
    alph = np.arctan2(r[1], r[0])
    if alph < 0:
        alph += 2*np.pi
    return alph
