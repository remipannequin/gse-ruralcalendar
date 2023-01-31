/* extension.js
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * SPDX-License-Identifier: GPL-2.0-or-later
 */

/* exported init */

const {Clutter, GObject, St, GLib, Gio, GdkPixbuf, Cogl, GnomeDesktop} = imports.gi;
const Gettext = imports.gettext;
const PanelMenu = imports.ui.panelMenu;
const PopupMenu = imports.ui.popupMenu;
const Main = imports.ui.main;
const Mainloop = imports.mainloop;
const ExtensionUtils = imports.misc.extensionUtils;
const Me = ExtensionUtils.getCurrentExtension();
const Domain = Gettext.domain(Me.metadata.uuid);
const _ = Domain.gettext;



// This is useful to get the file with the right language, and generate the url to wikipedia
// translator should make sure that such file exists...
let lang;

// Code from extension FRC@jcdubacq.dubacq.fr

// Astro code lifted from https://www.fourmilab.ch/documents/calendar/
// The copyright there states this code belongs to the public domain.
// Everything between BEGIN BORROWED CODE and up to END BORROWED CODE is
// therefore put in the public domain too. Modifications were made to
// put it in a more object-oriented form.

// BEGIN BORROWED CODE

/**
 *
 */
function Astro() {
    'use strict';
    this._init();
}

Astro.prototype = {
    _init() {
        'use strict';
        this.cache = {};
        this.cachestamp = {};
    },

    J2000: 2451545.0,              // Julian day of J2000 epoch
    JulianCentury: 36525.0,                // Days in Julian century
    JulianMillennium: 36250,   // Days in Julian millennium
    AstronomicalUnit: 149597870.0,            // Astronomical unit in kilometres
    TropicalYear: 365.24219878,           // Mean solar tropical year

    rtd(r) {
        'use strict';
        return (r * 180.0) / Math.PI;
    },
    dtr(d) {
        'use strict';
        return (d * Math.PI) / 180.0;
    },
    dcos(d) {
        'use strict';
        return Math.cos(this.dtr(d));
    },
    dsin(d) {
        'use strict';
        return Math.sin(this.dtr(d));
    },
    fixangle(a) {
        'use strict';
        return a - 360.0 * Math.floor(a / 360.0);
    },
    fixangr(a) {
        'use strict';
        return a - (2 * Math.PI) * Math.floor(a / (2 * Math.PI));
    },
    mod(a, b) {
        return a - (b * Math.floor(a / b));
    },

    /*  DELTAT  --  Determine the difference, in seconds, between
        Dynamical time and Universal time.  */
    /*  Table of observed Delta T values at the beginning of
        even numbered years from 1620 through 2002.  */
    deltaTtab: [
        121, 112, 103, 95, 88, 82, 77, 72, 68, 63, 60, 56, 53, 51, 48, 46,
        44, 42, 40, 38, 35, 33, 31, 29, 26, 24, 22, 20, 18, 16, 14, 12,
        11, 10, 9, 8, 7, 7, 7, 7, 7, 7, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10,
        10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13,
        13, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16,
        16, 16, 15, 15, 14, 13, 13.1, 12.5, 12.2, 12, 12, 12, 12, 12, 12,
        11.9, 11.6, 11, 10.2, 9.2, 8.2, 7.1, 6.2, 5.6, 5.4, 5.3, 5.4, 5.6,
        5.9, 6.2, 6.5, 6.8, 7.1, 7.3, 7.5, 7.6, 7.7, 7.3, 6.2, 5.2, 2.7,
        1.4, -1.2, -2.8, -3.8, -4.8, -5.5, -5.3, -5.6, -5.7, -5.9, -6,
        -6.3, -6.5, -6.2, -4.7, -2.8, -0.1, 2.6, 5.3, 7.7, 10.4, 13.3, 16,
        18.2, 20.2, 21.1, 22.4, 23.5, 23.8, 24.3, 24, 23.9, 23.9, 23.7,
        24, 24.3, 25.3, 26.2, 27.3, 28.2, 29.1, 30, 30.7, 31.4, 32.2,
        33.1, 34, 35, 36.5, 38.3, 40.2, 42.2, 44.5, 46.5, 48.5, 50.5,
        52.2, 53.8, 54.9, 55.8, 56.9, 58.3, 60, 61.6, 63, 65, 66.6,
    ],
    deltat(year) {
        var dt, f, i, t;
        if ((year >= 1620) && (year <= 2000)) {
            i = Math.floor((year - 1620) / 2);
            f = ((year - 1620) / 2) - i;  /* Fractional part of year */
            dt = deltaTtab[i] + ((deltaTtab[i + 1] - deltaTtab[i]) * f);
        } else {
            t = (year - 2000) / 100;
            if (year < 948) {
                dt = 2177 + (497 * t) + (44.1 * t * t);
            } else {
                dt = 102 + (102 * t) + (25.3 * t * t);
                if ((year > 2000) && (year < 2100))
                    dt += 0.37 * (year - 2100);
            }
        }
        return dt;
    },

    EquinoxpTerms: [
        485, 324.96,   1934.136,
        203, 337.23,  32964.467,
        199, 342.08,     20.186,
        182,  27.85, 445267.112,
        156,  73.14,  45036.886,
        136, 171.52,  22518.443,
        77, 222.54,  65928.934,
        74, 296.72,   3034.906,
        70, 243.58,   9037.513,
        58, 119.81,  33718.147,
        52, 297.17,    150.678,
        50,  21.02,   2281.226,
        45, 247.54,  29929.562,
        44, 325.15,  31555.956,
        29,  60.93,   4443.417,
        18, 155.12,  67555.328,
        17, 288.79,   4562.452,
        16, 198.04,  62894.029,
        14, 199.76,  31436.921,
        12,  95.39,  14577.848,
        12, 287.11,  31931.756,
        12, 320.81,  34777.259,
        9, 227.73,   1222.114,
        8,  15.45,  16859.074,
    ],
    JDE0tab1000: [
        [1721139.29189, 365242.13740,  0.06134,  0.00111, -0.00071],
        [1721233.25401, 365241.72562, -0.05323,  0.00907,  0.00025],
        [1721325.70455, 365242.49558, -0.11677, -0.00297,  0.00074],
        [1721414.39987, 365242.88257, -0.00769, -0.00933, -0.00006],
    ],
    JDE0tab2000: [
        [2451623.80984, 365242.37404,  0.05169, -0.00411, -0.00057],
        [2451716.56767, 365241.62603,  0.00325,  0.00888, -0.00030],
        [2451810.21715, 365242.01767, -0.11575,  0.00337,  0.00078],
        [2451900.05952, 365242.74049, -0.06223, -0.00823,  0.00032],
    ],
    equinox(year, which) {
        'use strict';
        var deltaL, i, j, JDE0, JDE, JDE0tab, S, T, W, Y;
        /*  Initialise terms for mean equinox and solstices.  We
            have two sets: one for years prior to 1000 and a second
            for subsequent years.  */
        if (year < 1000) {
            JDE0tab = this.JDE0tab1000;
            Y = year / 1000;
        } else {
            JDE0tab = this.JDE0tab2000;
            Y = (year - 2000) / 1000;
        }
        JDE0 =  JDE0tab[which][0] +
            (JDE0tab[which][1] * Y) +
            (JDE0tab[which][2] * Y * Y) +
            (JDE0tab[which][3] * Y * Y * Y) +
            (JDE0tab[which][4] * Y * Y * Y * Y);
        T = (JDE0 - 2451545.0) / 36525;
        W = (35999.373 * T) - 2.47;
        deltaL = 1 + (0.0334 * this.dcos(W)) + (0.0007 * this.dcos(2 * W));
        //  Sum the periodic terms for time T
        S = 0;
        for (i = j = 0; i < 24; i++) {
            S += this.EquinoxpTerms[j] * this.dcos(this.EquinoxpTerms[j + 1] + (this.EquinoxpTerms[j + 2] * T));
            j += 3;
        }
        JDE = JDE0 + ((S * 0.00001) / deltaL);
        return JDE;
    },

    /*  SUNPOS  --  Position of the Sun.  Please see the comments
        on the return statement at the end of this function
        which describe the array it returns.  We return
        intermediate values because they are useful in a
        variety of other contexts.  */
    sunpos(jd) {
        var T, T2, L0, M, e, C, sunLong, sunAnomaly, sunR,
            Omega, Lambda, epsilon, epsilon0, Alpha, Delta,
            AlphaApp, DeltaApp;
        T = (jd - this.J2000) / this.JulianCentury;
        T2 = T * T;
        L0 = 280.46646 + (36000.76983 * T) + (0.0003032 * T2);
        L0 = this.fixangle(L0);
        M = 357.52911 + (35999.05029 * T) + (-0.0001537 * T2);
        M = this.fixangle(M);
        e = 0.016708634 + (-0.000042037 * T) + (-0.0000001267 * T2);
        C = ((1.914602 + (-0.004817 * T) + (-0.000014 * T2)) * this.dsin(M)) +
            ((0.019993 - (0.000101 * T)) * this.dsin(2 * M)) +
            (0.000289 * this.dsin(3 * M));
        sunLong = L0 + C;
        sunAnomaly = M + C;
        sunR = (1.000001018 * (1 - (e * e))) / (1 + (e * this.dcos(sunAnomaly)));
        Omega = 125.04 - (1934.136 * T);
        Lambda = sunLong + -0.00569 + (-0.00478 * this.dsin(Omega));
        epsilon0 = this.obliqeq(jd);
        epsilon = epsilon0 + (0.00256 * this.dcos(Omega));
        Alpha = this.rtd(Math.atan2(this.dcos(epsilon0) * this.dsin(sunLong), this.dcos(sunLong)));
        Alpha = this.fixangle(Alpha);
        Delta = this.rtd(Math.asin(this.dsin(epsilon0) * this.dsin(sunLong)));
        AlphaApp = this.rtd(Math.atan2(this.dcos(epsilon) * this.dsin(Lambda), this.dcos(Lambda)));
        AlphaApp = this.fixangle(AlphaApp);
        DeltaApp = this.rtd(Math.asin(this.dsin(epsilon) * this.dsin(Lambda)));
        return [
            L0,                           //  [0] Geometric mean longitude of the Sun
            M,                            //  [1] Mean anomaly of the Sun
            e,                            //  [2] Eccentricity of the Earth's orbit
            C,                            //  [3] Sun's equation of the Centre
            sunLong,                      //  [4] Sun's true longitude
            sunAnomaly,                   //  [5] Sun's true anomaly
            sunR,                         //  [6] Sun's radius vector in AU
            Lambda,                       //  [7] Sun's apparent longitude at true equinox of the date
            Alpha,                        //  [8] Sun's true right ascension
            Delta,                        //  [9] Sun's true declination
            AlphaApp,                     // [10] Sun's apparent right ascension
            DeltaApp,                      // [11] Sun's apparent declination
        ];
    },

    /*  OBLIQEQ  --  Calculate the obliquity of the ecliptic for a given
        Julian date.  This uses Laskar's tenth-degree
        polynomial fit (J. Laskar, Astronomy and
        Astrophysics, Vol. 157, page 68 [1986]) which is
        accurate to within 0.01 arc second between AD 1000
        and AD 3000, and within a few seconds of arc for
        +/-10000 years around AD 2000.  If we're outside the
        range in which this fit is valid (deep time) we
        simply return the J2000 value of the obliquity, which
        happens to be almost precisely the mean.  */

    oterms: [-4680.93, -1.55, 1999.25, -51.38, -249.67, -39.05, 7.12, 27.87, 5.79, 2.45],
    obliqeq(jd) {
        var eps, u, v, i;
        v = u = (jd - this.J2000) / (this.JulianCentury * 100);
        eps = 23 + (26 / 60.0) + (21.448 / 3600.0);
        if (Math.abs(u) < 1.0) {
            for (i = 0; i < 10; i++) {
                eps += (this.oterms[i] / 3600.0) * v;
                v *= u;
            }
        }
        return eps;
    },

    /*  NUTATION  --  Calculate the nutation in longitude, deltaPsi, and
        obliquity, deltaEpsilon for a given Julian date
        jd.  Results are returned as a two element Array
        giving (deltaPsi, deltaEpsilon) in degrees.  */
    /* Periodic terms for nutation in longiude (delta \Psi) and
       obliquity (delta \Epsilon) as given in table 21.A of
       Meeus, "Astronomical Algorithms", first edition. */
    nutArgMult: [0, 0, 0, 0, 1, -2, 0, 0, 2, 2, 0, 0, 0, 2, 2, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, -2, 1, 0, 2, 2, 0, 0, 0, 2, 1, 0, 0, 1, 2, 2, -2, -1, 0, 2, 2, -2, 0, 1, 0, 0, -2, 0, 0, 2, 1, 0, 0, -1, 2, 2, 2, 0, 0, 0, 0, 0, 0, 1, 0, 1, 2, 0, -1, 2, 2, 0, 0, -1, 0, 1, 0, 0, 1, 2, 1, -2, 0, 2, 0, 0, 0, 0, -2, 2, 1, 2, 0, 0, 2, 2, 0, 0, 2, 2, 2, 0, 0, 2, 0, 0, -2, 0, 1, 2, 2, 0, 0, 0, 2, 0, -2, 0, 0, 2, 0, 0, 0, -1, 2, 1, 0, 2, 0, 0, 0, 2, 0, -1, 0, 1, -2, 2, 0, 2, 2, 0, 1, 0, 0, 1, -2, 0, 1, 0, 1, 0, -1, 0, 0, 1, 0, 0, 2, -2, 0, 2, 0, -1, 2, 1, 2, 0, 1, 2, 2, 0, 1, 0, 2, 2, -2, 1, 1, 0, 0, 0, -1, 0, 2, 2, 2, 0, 0, 2, 1, 2, 0, 1, 0, 0, -2, 0, 2, 2, 2, -2, 0, 1, 2, 1, 2, 0, -2, 0, 1, 2, 0, 0, 0, 1, 0, -1, 1, 0, 0, -2, -1, 0, 2, 1, -2, 0, 0, 0, 1, 0, 0, 2, 2, 1, -2, 0, 2, 0, 1, -2, 1, 0, 2, 1, 0, 0, 1, -2, 0, -1, 0, 1, 0, 0, -2, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 2, 0, -1, -1, 1, 0, 0, 0, 1, 1, 0, 0, 0, -1, 1, 2, 2, 2, -1, -1, 2, 2, 0, 0, -2, 2, 2, 0, 0, 3, 2, 2, 2, -1, 0, 2, 2],
    nutArgCoeff: [-171996, -1742, 92095, 89, -13187, -16, 5736, -31, -2274, -2, 977, -5, 2062, 2, -895, 5, 1426, -34, 54, -1, 712, 1, -7, 0, -517, 12, 224, -6, -386, -4, 200, 0, -301, 0, 129, -1, 217, -5, -95, 3, -158, 0, 0, 0, 129, 1, -70, 0, 123, 0, -53, 0, 63, 0, 0, 0, 63, 1, -33, 0, -59, 0, 26, 0, -58, -1, 32, 0, -51, 0, 27, 0, 48, 0, 0, 0, 46, 0, -24, 0, -38, 0, 16, 0, -31, 0, 13, 0, 29, 0, 0, 0, 29, 0, -12, 0, 26, 0, 0, 0, -22, 0, 0, 0, 21, 0, -10, 0, 17, -1, 0, 0, 16, 0, -8, 0, -16, 1, 7, 0, -15, 0, 9, 0, -13, 0, 7, 0, -12, 0, 6, 0, 11, 0, 0, 0, -10, 0, 5, 0, -8, 0, 3, 0, 7, 0, -3, 0, -7, 0, 0, 0, -7, 0, 3, 0, -7, 0, 3, 0, 6, 0, 0, 0, 6, 0, -3, 0, 6, 0, -3, 0, -6, 0, 3, 0, -6, 0, 3, 0, 5, 0, 0, 0, -5, 0, 3, 0, -5, 0, 3, 0, -5, 0, 3, 0, 4, 0, 0, 0, 4, 0, 0, 0, 4, 0, 0, 0, -4, 0, 0, 0, -4, 0, 0, 0, -4, 0, 0, 0, 3, 0, 0, 0, -3, 0, 0, 0, -3, 0, 0, 0, -3, 0, 0, 0, -3, 0, 0, 0, -3, 0, 0, 0, -3, 0, 0, 0, -3, 0, 0, 0],
    nutation(jd) {
        var deltaPsi, deltaEpsilon,
            i, j,
            t = (jd - 2451545.0) / 36525.0, t2, t3, to10,
            ta = [],
            dp = 0, de = 0, ang;
        t3 = t * (t2 = t * t);
        /* Calculate angles.  The correspondence between the elements
           of our array and the terms cited in Meeus are:
           ta[0] = D  ta[0] = M  ta[2] = M'  ta[3] = F  ta[4] = \Omega
        */
        ta[0] = this.dtr(297.850363 + 445267.11148 * t - 0.0019142 * t2 +
                         t3 / 189474.0);
        ta[1] = this.dtr(357.52772 + 35999.05034 * t - 0.0001603 * t2 -
                         t3 / 300000.0);
        ta[2] = this.dtr(134.96298 + 477198.867398 * t + 0.0086972 * t2 +
                         t3 / 56250.0);
        ta[3] = this.dtr(93.27191 + 483202.017538 * t - 0.0036825 * t2 +
                         t3 / 327270);
        ta[4] = this.dtr(125.04452 - 1934.136261 * t + 0.0020708 * t2 +
                         t3 / 450000.0);
        /* Range reduce the angles in case the sine and cosine functions
           don't do it as accurately or quickly. */
        for (i = 0; i < 5; i++)
            ta[i] = this.fixangr(ta[i]);

        to10 = t / 10.0;
        for (i = 0; i < 63; i++) {
            ang = 0;
            for (j = 0; j < 5; j++) {
                if (this.nutArgMult[(i * 5) + j] !== 0)
                    ang += this.nutArgMult[(i * 5) + j] * ta[j];
            }
            dp += (this.nutArgCoeff[(i * 4) + 0] + this.nutArgCoeff[(i * 4) + 1] * to10) * Math.sin(ang);
            de += (this.nutArgCoeff[(i * 4) + 2] + this.nutArgCoeff[(i * 4) + 3] * to10) * Math.cos(ang);
        }
        /* Return the result, converting from ten thousandths of arc
           seconds to radians in the process. */
        deltaPsi = dp / (3600.0 * 10000.0);
        deltaEpsilon = de / (3600.0 * 10000.0);
        return [deltaPsi, deltaEpsilon];
    },

    /*  EQUATIONOFTIME  --  Compute equation of time for a given moment.
        Returns the equation of time as a fraction of
        a day.  */
    equationOfTime(jd) {
        var alpha, deltaPsi, E, epsilon, L0, tau;
        tau = (jd - this.J2000) / this.JulianMillennium;
        L0 = 280.4664567 + (360007.6982779 * tau) +
            (0.03032028 * tau * tau) +
            ((tau * tau * tau) / 49931) +
            -((tau * tau * tau * tau) / 15300) +
            -((tau * tau * tau * tau * tau) / 2000000);
        L0 = this.fixangle(L0);
        alpha = this.sunpos(jd)[10];
        deltaPsi = this.nutation(jd)[0];
        epsilon = this.obliqeq(jd) + this.nutation(jd)[1];
        E = L0 + -0.0057183 + -alpha + (deltaPsi * this.dcos(epsilon));
        E -= 20.0 * Math.floor(E / 20.0);
        E /= 24 * 60;
        return E;
    },

    //  GREGORIAN_TO_JD  --  Determine Julian day number from Gregorian calendar date
    gregorianEpoch: 1721425.5,
    gregorian_to_jd(year, month, day) {
        return (this.gregorianEpoch - 1) +
            (365 * (year - 1)) +
            Math.floor((year - 1) / 4) +
            -Math.floor((year - 1) / 100) +
            Math.floor((year - 1) / 400) +
            Math.floor((((367 * month) - 362) / 12) +
                       (month <= 2 ? 0
                           : this.leap_gregorian(year) ? -1 : -2
                       ) +
                       day);
    },

    //  LEAP_GREGORIAN  --  Is a given year in the Gregorian calendar a leap year ?
    leap_gregorian(year) {
        return ((year % 4) === 0) &&
            !(((year % 100) === 0) && ((year % 400) !== 0));
    },

    //  JD_TO_GREGORIAN  --  Calculate Gregorian calendar date from Julian day
    jd_to_gregorian(jd) {
        var wjd, depoch, quadricent, dqc, cent, dcent, quad, dquad, yindex, year, yearday, leapadj, month, day;

        wjd = Math.floor(jd - 0.5) + 0.5;
        depoch = wjd - this.gregorianEpoch;
        quadricent = Math.floor(depoch / 146097);
        dqc = this.mod(depoch, 146097);
        cent = Math.floor(dqc / 36524);
        dcent = this.mod(dqc, 36524);
        quad = Math.floor(dcent / 1461);
        dquad = this.mod(dcent, 1461);
        yindex = Math.floor(dquad / 365);
        year = (quadricent * 400) + (cent * 100) + (quad * 4) + yindex;
        if (!((cent === 4) || (yindex === 4)))
            year++;

        yearday = wjd - this.gregorian_to_jd(year, 1, 1);
        leapadj = wjd < this.gregorian_to_jd(year, 3, 1) ? 0
            :                   this.leap_gregorian(year) ? 1 : 2;
        month = Math.floor((((yearday + leapadj) * 12) + 373) / 367);
        day = (wjd - this.gregorian_to_jd(year, month, 1)) + 1;

        return [year, month, day];
    },

    /*  EQUINOXE_A_PARIS  --  Determine Julian day and fraction of the
        September equinox at the Paris meridian in
        a given Gregorian year.  */
    equinoxe_a_paris(year) {
        'use strict';
        var equJED, equJD, equAPP, equParis, dtParis;
        //  September equinox in dynamical time
        equJED = this.equinox(year, 2);
        //  Correct for delta T to obtain Universal time
        equJD = equJED - (this.deltat(year) / (24 * 60 * 60));
        //  Apply the equation of time to yield the apparent time at Greenwich
        equAPP = equJD + this.equationOfTime(equJED);
        /*  Finally, we must correct for the constant difference between
            the Greenwich meridian and that of Paris, 2°20'15" to the
            East.  */
        dtParis = (2 + (20 / 60.0) + (15 / (60 * 60.0))) / 360;
        equParis = equAPP + dtParis;
        return equParis;
    },

    /*  PARIS_EQUINOXE_JD  --  Calculate Julian day during which the
        September equinox, reckoned from the Paris
        meridian, occurred for a given Gregorian
        year.  */
    paris_equinoxe_jd(year) {
        var ep, epg;
        ep = this.equinoxe_a_paris(year);
        epg = Math.floor(ep - 0.5) + 0.5;
        return epg;
    },

    frenchRevolutionaryEpoch: 2375839.5,
    anneeDeLaRevolution(jd) {
        var guess = this.jd_to_gregorian(jd)[0] - 2,
            lasteq, nexteq, adr;
        lasteq = this.paris_equinoxe_jd(guess);
        while (lasteq > jd) {
            guess--;
            lasteq = this.paris_equinoxe_jd(guess);
        }
        nexteq = lasteq - 1;
        while (!((lasteq <= jd) && (jd < nexteq))) {
            lasteq = nexteq;
            guess++;
            nexteq = this.paris_equinoxe_jd(guess);
        }
        adr = Math.round((lasteq - this.frenchRevolutionaryEpoch) / this.TropicalYear) + 1;
        return [adr, lasteq];
    },

    /*  JD_TO_FRENCH_REVOLUTIONARY  --  Calculate date in the French Revolutionary
        calendar from Julian day.  The five or six
        "sansculottides" are considered a thirteenth
        month in the results of this function.  */
    jd_to_french_revolutionary(jd) {
        var an, mois, decade, jour, adr, equinoxe;
        jd = Math.floor(jd) + 0.5;
        if (this.cachestamp['FRC'] === jd)
            return this.cache['FRC'];

        adr = this.anneeDeLaRevolution(jd);
        an = adr[0];
        equinoxe = adr[1];
        mois = Math.floor((jd - equinoxe) / 30) + 1;
        jour = (jd - equinoxe) % 30;
        decade = Math.floor(jour / 10) + 1;
        jour = (jour % 10) + 1;
        this.cachestamp['FRC'] = jd;
        this.cache['FRC'] = [an, mois, decade, jour];

        return this.cache['FRC'];
    },

};

class RuralCalendar {
    constructor(data) {
        this.monthNames = [
            _('Vendémiaire'), _('Brumaire'), _('Frimaire'),
            _('Nivôse'), _('Pluviôse'), _('Ventôse'),
            _('Germinal'), _('Floréal'), _('Prairial'),
            _('Messidor'), _('Thermidor'), _('Fructidor'),
            _('Sans-culottides'),
        ];
        this.sansculottidesNames = [
            _('jour de la vertu'), _('jour du génie'),
            _('jour du travail'), _('jour de l´opinion'),
            _('jour des récompenses'), _('jour de la révolution'),
        ];
        this.dayNames = [
            'Primidi', 'Duodi', 'Tridi', 'Quartidi', 'Quintidi',
            'Sextidi', 'Septidi', 'Octidi', 'Nonidi', 'Décadi',
        ];
        this.decadeNames = [_('first'), _('second'), _('third')];
        this.saintsNames = data;
        // set the date to today,
        this.update(new Date());
    }

    romanNumeral(n) {
        var val, s = '', limit = 3999, i = 0;
        var v = [1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1];
        var r = ['M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I'];
        if (n < 1 || n > limit)
            return '';
        while (i < 13) {
            val = v[i];
            while (n >= val) {
                n -= val;
                s += r[i];
            }
            if (n === 0)
                return s;
            ++i;
        }
        return '';
    }


    daymonth() {
        let daymonth = (this.jrr[2] - 1) * 10 + this.jrr[3];
        if (this.jrr[1] !== 13) {
            if (daymonth === 1)
                daymonth += '<sup>er</sup>';
        } else if (daymonth === 1) {
            daymonth = `${daymonth}<sup>er</sup>`;
        } else {
            daymonth = `${daymonth}<sup>e</sup>`;
        }
        return daymonth;
    }

    // END BORROWED CODE

    /**
     * Update the date.
     *
     * @param {Date} time the current date
     */
    update(time) {
        //  Update Julian day
        let j = astro.gregorian_to_jd(time.getFullYear(), time.getMonth() + 1, time.getDate());
        this.jrr = astro.jd_to_french_revolutionary(j);
    }

    /**
     *
     * @returns the date string
     */
    getDate() {
        let dayName = this.daymonth();
        let monthName = this.monthNames[this.jrr[1] - 1];
        return `${dayName} ${monthName}`;
    }

    /**
     * Get the name of the object.
     *
     * @returns the name of the plant, animal (every 5 days), tool (every 10 days)
     */
    getTitle() {
        if (this.jrr[1] === 13)
            return this.sansculottidesNames[this.jrr[1] - 1];

        let i = (this.jrr[1] - 1) * 30 + (this.jrr[2] - 1) * 10 + (this.jrr[3] - 1);
        return this.saintsNames[i][0];
    }

    /**
     * Get the wikipedia URL for the summary page.
     *
     * @returns the url of the summary page.
     */
    getSummaryUrl() {
        if (this.jrr[1] === 13)
            return null;

        let i = (this.jrr[1] - 1) * 30 + (this.jrr[2] - 1) * 10 + (this.jrr[3] - 1);
        let pageTitle = this.saintsNames[i][1];
        return _(`https://${lang}.wikipedia.org/api/rest_v1/page/summary/${pageTitle}?redirect=true`);
    }

    getWikiUrl() {
        if (this.jrr[1] === 13)
            return null;

        let i = (this.jrr[1] - 1) * 30 + (this.jrr[2] - 1) * 10 + (this.jrr[3] - 1);
        let pageTitle = this.saintsNames[i][1];
        return _(`https://${lang}.wikipedia.org/wiki/${pageTitle}`);
    }
}

/**
 * Filter html markup to keep only those that are in pango.
 *
 * @param {str} text html text to convert
 * @returns the filtered string
 */
function htmlToPango(text) {
    let result = '';
    var enabledTags = ['i', 'b'];
    const regexp = /<[^>]*>/g;
    const tagValue = /<[/]?(.*)>/;
    const matches = text.matchAll(regexp);
    var i = 0;
    for (const match of matches) {
        // log(`Found ${match[0]} start=${match.index} end=${match.index + match[0].length}.`);
        // log(i, match.index, text.substring(i, match.index));

        // copy what is before the tag
        result = result.concat(text.substring(i, match.index));
        i = match.index + match[0].length;
        // copy the tag if enabled
        var value = match[0].match(tagValue);
        if (enabledTags.includes(value[1]))
            result = result.concat(value[0]);
    }
    result += result.substring(i, text.length);

    return result;
}

const RuralCalendarTopMenu = GObject.registerClass(
class RuralCalendarTopMenu extends PanelMenu.Button {
    _init(calendar) {
        this.calendar = calendar;
        this.detail_url = null;
        this.timeout = null;
        this.toptext = 'Rural Calendar';
        super._init(0.5, 'Rural Calendar');
        this.kill = false;
        this.toplabel = new St.Label({
            text: this.toptext,
            y_expand: true,
            y_align: Clutter.ActorAlign.CENTER,
        });
        this.toplabel.clutter_text.set_use_markup(true);
        // this.toplabel.clutter_text.set_width(200);
        this.add_child(this.toplabel);

        // Add popup menu items
        let box = new St.BoxLayout({style_class: 'detail-panel', vertical: true});
        this.detail_menu_item = new PopupMenu.PopupMenuItem('');
        this.detail_menu_item.setOrnament(PopupMenu.Ornament.HIDDEN);
        this.detail_menu_item.add_child(box);
        this.detail_menu_item.connect('button-press-event', this.open_wiki.bind(this));

        // title
        this.detail_title_label = new St.Label();
        box.add_child(this.detail_title_label);
        // description
        this.detail_desc_label = new St.Label();
        box.add_child(this.detail_desc_label);
        // image
        this.detail_icon = new Clutter.Actor();
        this.detail_icon.content_gravity = Clutter.ContentGravity.CENTER;
        box.add_child(this.detail_icon);
        // Extract
        this.detail_extract_label = new St.Label();
        box.add_child(this.detail_extract_label);
        // Add item to menu
        this.menu.addMenuItem(this.detail_menu_item);

        this.timeout = Mainloop.timeout_add_seconds(15, this.update.bind(this));
        this.update();
    }

    /**
     * Open the wikipedia page.
     */
    open_wiki() {
        let url = this.calendar.getWikiUrl();
        log(url);
        Gio.AppInfo.launch_default_for_uri_async(url, null, null, null);
    }

    /**
     * Update teh widget, fetching detail data if required.
     *
     * @returns true
     */
    update() {
        if (this.kill)
            return false;

        try {
            let time = new Date();
            this.calendar.update(time);
            this.date = this.calendar.getDate();
            this.title = this.calendar.getTitle();
            this.toplabel.clutter_text.set_markup(`${this.date}, <i>${this.title}</i>`);

            let link = this.calendar.getSummaryUrl();
            if (!link) {
                // No detail for this day, hide details
                this.detail_menu_item.hide();
            } else {
                this.detail_menu_item.show();
                if (this.detail_url !== link) {
                    this.detail_url = link;
                    this.fetch_detail();
                }
            }
            return true;
        } catch {
            // we may be on an old call of update, stoping timeout
            log('stopping update loop');
            return false;
        }
    }

    /**
     * Fetch the details from wikipedia.
     */
    fetch_detail() {
        if (!this.detail_url)
            log("Invalid URL, won't fetch details.");


        let detailFile = Gio.File.new_for_uri(this.detail_url);
        detailFile.read(null);
        const [, contents] = detailFile.load_contents(null);
        const decoder = new TextDecoder('utf-8');
        const contentsString = decoder.decode(contents);
        let rspData = JSON.parse(contentsString);

        if (this.title !== rspData.title)
            this.detail_title_label.clutter_text.set_markup(`<big>${this.title} (${rspData.title})</big>`);
        else
            this.detail_title_label.clutter_text.set_markup(`<big>${this.title}</big>`);

        let w;
        if (rspData.thumbnail) {
            let thumbSrcUrl = rspData.thumbnail.source;
            w = rspData.thumbnail.width;
            let h = rspData.thumbnail.height;

            let f = Gio.File.new_for_uri(thumbSrcUrl);
            let pixbuf = GdkPixbuf.Pixbuf.new_from_stream(f.read(null), null);
            let img = new Clutter.Image();
            img.set_data(pixbuf.get_pixels(),
                pixbuf.get_has_alpha()
                    ? Cogl.PixelFormat.RGBA_8888
                    : Cogl.PixelFormat.RGB_888,
                w,
                h,
                pixbuf.get_rowstride());

            this.detail_icon.content = img;
            this.detail_icon.width = w;
            this.detail_icon.height = h;
            this.detail_icon.show();
        } else {
            w = 200;
            this.detail_icon.hide();
        }
        if (rspData.description) {
            this.detail_desc_label.show();
            this.detail_desc_label.text = rspData.description;
            this.detail_desc_label.clutter_text.set_width(w);
            this.detail_desc_label.clutter_text.set_line_wrap(true);
            // this.detail_desc_label.clutter_text.set_line_wrap_mode(Pango.WrapMode.CHAR);
        } else {
            this.detail_desc_label.hide();
        }

        if (rspData.extract) {
            this.detail_extract_label.show();
            this.detail_extract_label.clutter_text.set_markup(htmlToPango(rspData.extract_html));
            this.detail_extract_label.clutter_text.set_width(w);
            this.detail_extract_label.clutter_text.set_line_wrap(true);
            this.detail_extract_label.clutter_text.set_justify(true);
            // this.detail_extract_label.clutter_text.set_line_wrap_mode(Pango.WrapMode.CHAR);
        } else {
            this.detail_extract_label.hide();
        }
    }

    removeTimeout() {
        if (this.timeout) {
            Mainloop.source_remove(this.timeout);
            this.timeout = null;
        }
        this.widget.kill = true;
    }
});



let astro = new Astro();


class Extension {
    constructor() {
        // Read names and links
        const file = Gio.File.new_for_path(`${Me.path}/french-republican-calendar_${lang}.json`);
        // Synchronous, blocking method
        const [, contents] = file.load_contents(null);
        const decoder = new TextDecoder('utf-8');
        const contentsString = decoder.decode(contents);
        let data = JSON.parse(contentsString);
        this.calendar = new RuralCalendar(data);
    }

    enable() {
        this.widget = new RuralCalendarTopMenu(this.calendar);

        let pos = 1;
        if ('apps-menu' in Main.panel.statusArea)
            pos = 2;
        Main.panel.addToStatusArea('ruralcalendar-menu', this.widget, pos, 'center');
    }

    disable() {
        if (this.widget) {
            this.widget.removeTimeout();
            this.widget.destroy();
            this.widget = null;
        }
    }
}

/**
 *
 */
function getSupportedLang() {
    let [, langCode, , ,] = GnomeDesktop.parse_locale(GLib.getenv('LANG'));
    let f = Gio.File.new_for_path(`${Me.path}/french-republican-calendar_${langCode}.json`);
    if (f.query_exists(null)) {
        return langCode;
    } else {
        log(`unable to find data file for language ${langCode}.`);
        // Default language
        return 'en';
    }
}

/**
 *
 */
function init() {
    ExtensionUtils.initTranslations(Me.metadata.uuid);
    lang = getSupportedLang();
    return new Extension();
}
