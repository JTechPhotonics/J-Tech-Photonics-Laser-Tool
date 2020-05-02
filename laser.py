#!/usr/bin/env python
"""
Modified by Jay Johnson 2015, J Tech Photonics, Inc., jtechphotonics.com
modified by Adam Polak 2014, polakiumengineering.org

based on Copyright (C) 2009 Nick Drobchenko, nick@cnc-club.ru
based on gcode.py (C) 2007 hugomatic...
based on addnodes.py (C) 2005,2007 Aaron Spike, aaron@ekips.org
based on dots.py (C) 2005 Aaron Spike, aaron@ekips.org
based on interp.py (C) 2005 Aaron Spike, aaron@ekips.org
based on bezmisc.py (C) 2005 Aaron Spike, aaron@ekips.org
based on cubicsuperpath.py (C) 2005 Aaron Spike, aaron@ekips.org

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""
import inkex
import simplestyle
import cubicsuperpath
import simpletransform

import os
import math
import bezmisc
import re
import sys
import time
import numpy
import gettext

_ = gettext.gettext

# Check if inkex has error messages. (0.46 version does not have one) Could be removed later.
if "errormsg" not in dir(inkex):
    inkex.errormsg = lambda msg: sys.stderr.write((unicode(msg) + "\n").encode("UTF-8"))


def bezierslopeatt(((bx0, by0), (bx1, by1), (bx2, by2), (bx3, by3)), t):
    ax, ay, bx, by, cx, cy, x0, y0 = bezmisc.bezierparameterize(((bx0, by0), (bx1, by1), (bx2, by2), (bx3, by3)))
    dx = 3 * ax * (t ** 2) + 2 * bx * t + cx
    dy = 3 * ay * (t ** 2) + 2 * by * t + cy
    if dx == dy == 0:
        dx = 6 * ax * t + 2 * bx
        dy = 6 * ay * t + 2 * by
        if dx == dy == 0:
            dx = 6 * ax
            dy = 6 * ay
            if dx == dy == 0:
                print_("Slope error x = %s*t^3+%s*t^2+%s*t+%s, y = %s*t^3+%s*t^2+%s*t+%s,  t = %s, dx==dy==0" % (
                    ax, bx, cx, dx, ay, by, cy, dy, t))
                print_(((bx0, by0), (bx1, by1), (bx2, by2), (bx3, by3)))
                dx, dy = 1, 1

    return dx, dy


bezmisc.bezierslopeatt = bezierslopeatt

################################################################################
#
#        Styles and additional parameters
#
################################################################################

math.pi2 = math.pi * 2
straight_tolerance = 0.0001
straight_distance_tolerance = 0.0001
engraving_tolerance = 0.0001
loft_lengths_tolerance = 0.0000001
options = {}
defaults = {
    'header': """
G90
""",
    'footer': """G1 X0 Y0

"""
}

intersection_recursion_depth = 10
intersection_tolerance = 0.00001

styles = {
    "loft_style": {
        'main curve': simplestyle.formatStyle(
            {'stroke': '#88f', 'fill': 'none', 'stroke-width': '1', 'marker-end': 'url(#Arrow2Mend)'}),
    },
    "biarc_style": {
        'biarc0': simplestyle.formatStyle(
            {'stroke': '#88f', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '1'}),
        'biarc1': simplestyle.formatStyle(
            {'stroke': '#8f8', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '1'}),
        'line': simplestyle.formatStyle(
            {'stroke': '#f88', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '1'}),
        'area': simplestyle.formatStyle(
            {'stroke': '#777', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '0.1'}),
    },
    "biarc_style_dark": {
        'biarc0': simplestyle.formatStyle(
            {'stroke': '#33a', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '1'}),
        'biarc1': simplestyle.formatStyle(
            {'stroke': '#3a3', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '1'}),
        'line': simplestyle.formatStyle(
            {'stroke': '#a33', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '1'}),
        'area': simplestyle.formatStyle(
            {'stroke': '#222', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '0.3'}),
    },
    "biarc_style_dark_area": {
        'biarc0': simplestyle.formatStyle(
            {'stroke': '#33a', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '0.1'}),
        'biarc1': simplestyle.formatStyle(
            {'stroke': '#3a3', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '0.1'}),
        'line': simplestyle.formatStyle(
            {'stroke': '#a33', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '0.1'}),
        'area': simplestyle.formatStyle(
            {'stroke': '#222', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '0.3'}),
    },
    "biarc_style_i": {
        'biarc0': simplestyle.formatStyle(
            {'stroke': '#880', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '1'}),
        'biarc1': simplestyle.formatStyle(
            {'stroke': '#808', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '1'}),
        'line': simplestyle.formatStyle(
            {'stroke': '#088', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '1'}),
        'area': simplestyle.formatStyle(
            {'stroke': '#999', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '0.3'}),
    },
    "biarc_style_dark_i": {
        'biarc0': simplestyle.formatStyle(
            {'stroke': '#dd5', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '1'}),
        'biarc1': simplestyle.formatStyle(
            {'stroke': '#d5d', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '1'}),
        'line': simplestyle.formatStyle(
            {'stroke': '#5dd', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '1'}),
        'area': simplestyle.formatStyle(
            {'stroke': '#aaa', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '0.3'}),
    },
    "biarc_style_lathe_feed": {
        'biarc0': simplestyle.formatStyle(
            {'stroke': '#07f', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '.4'}),
        'biarc1': simplestyle.formatStyle(
            {'stroke': '#0f7', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '.4'}),
        'line': simplestyle.formatStyle(
            {'stroke': '#f44', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '.4'}),
        'area': simplestyle.formatStyle(
            {'stroke': '#aaa', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '0.3'}),
    },
    "biarc_style_lathe_passing feed": {
        'biarc0': simplestyle.formatStyle(
            {'stroke': '#07f', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '.4'}),
        'biarc1': simplestyle.formatStyle(
            {'stroke': '#0f7', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '.4'}),
        'line': simplestyle.formatStyle(
            {'stroke': '#f44', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '.4'}),
        'area': simplestyle.formatStyle(
            {'stroke': '#aaa', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '0.3'}),
    },
    "biarc_style_lathe_fine feed": {
        'biarc0': simplestyle.formatStyle(
            {'stroke': '#7f0', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '.4'}),
        'biarc1': simplestyle.formatStyle(
            {'stroke': '#f70', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '.4'}),
        'line': simplestyle.formatStyle(
            {'stroke': '#744', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '.4'}),
        'area': simplestyle.formatStyle(
            {'stroke': '#aaa', 'fill': 'none', "marker-end": "url(#DrawCurveMarker)", 'stroke-width': '0.3'}),
    },
    "area artefact": simplestyle.formatStyle({'stroke': '#ff0000', 'fill': '#ffff00', 'stroke-width': '1'}),
    "area artefact arrow": simplestyle.formatStyle({'stroke': '#ff0000', 'fill': '#ffff00', 'stroke-width': '1'}),
    "dxf_points": simplestyle.formatStyle({"stroke": "#ff0000", "fill": "#ff0000"}),

}


################################################################################
#        Cubic Super Path additional functions
################################################################################

def csp_segment_to_bez(sp1, sp2):
    return sp1[1:] + sp2[:2]


def csp_split(sp1, sp2, t=.5):
    [x1, y1], [x2, y2], [x3, y3], [x4, y4] = sp1[1], sp1[2], sp2[0], sp2[1]
    x12 = x1 + (x2 - x1) * t
    y12 = y1 + (y2 - y1) * t
    x23 = x2 + (x3 - x2) * t
    y23 = y2 + (y3 - y2) * t
    x34 = x3 + (x4 - x3) * t
    y34 = y3 + (y4 - y3) * t
    x1223 = x12 + (x23 - x12) * t
    y1223 = y12 + (y23 - y12) * t
    x2334 = x23 + (x34 - x23) * t
    y2334 = y23 + (y34 - y23) * t
    x = x1223 + (x2334 - x1223) * t
    y = y1223 + (y2334 - y1223) * t
    return [sp1[0], sp1[1], [x12, y12]], [[x1223, y1223], [x, y], [x2334, y2334]], [[x34, y34], sp2[1], sp2[2]]


def csp_curvature_at_t(sp1, sp2, t, depth=3):
    ax, ay, bx, by, cx, cy, dx, dy = bezmisc.bezierparameterize(csp_segment_to_bez(sp1, sp2))

    # curvature = (x'y''-y'x'') / (x'^2+y'^2)^1.5

    f1x = 3 * ax * t ** 2 + 2 * bx * t + cx
    f1y = 3 * ay * t ** 2 + 2 * by * t + cy
    f2x = 6 * ax * t + 2 * bx
    f2y = 6 * ay * t + 2 * by
    d = (f1x ** 2 + f1y ** 2) ** 1.5
    if d != 0:
        return (f1x * f2y - f1y * f2x) / d
    else:
        t1 = f1x * f2y - f1y * f2x
        if t1 > 0: return 1e100
        if t1 < 0: return -1e100
        # Use the Lapitals rule to solve 0/0 problem for 2 times...
        t1 = 2 * (bx * ay - ax * by) * t + (ay * cx - ax * cy)
        if t1 > 0: return 1e100
        if t1 < 0: return -1e100
        t1 = bx * ay - ax * by
        if t1 > 0: return 1e100
        if t1 < 0: return -1e100
        if depth > 0:
            # little hack ;^) hope it wont influence anything...
            return csp_curvature_at_t(sp1, sp2, t * 1.004, depth - 1)
        return 1e100


def csp_at_t(sp1, sp2, t):
    ax, bx, cx, dx = sp1[1][0], sp1[2][0], sp2[0][0], sp2[1][0]
    ay, by, cy, dy = sp1[1][1], sp1[2][1], sp2[0][1], sp2[1][1]

    x1, y1 = ax + (bx - ax) * t, ay + (by - ay) * t
    x2, y2 = bx + (cx - bx) * t, by + (cy - by) * t
    x3, y3 = cx + (dx - cx) * t, cy + (dy - cy) * t

    x4, y4 = x1 + (x2 - x1) * t, y1 + (y2 - y1) * t
    x5, y5 = x2 + (x3 - x2) * t, y2 + (y3 - y2) * t

    x, y = x4 + (x5 - x4) * t, y4 + (y5 - y4) * t
    return [x, y]


def cspseglength(sp1, sp2, tolerance=0.001):
    bez = (sp1[1][:], sp1[2][:], sp2[0][:], sp2[1][:])
    return bezmisc.bezierlength(bez, tolerance)


#        Distance calculation from point to arc
def point_to_arc_distance(p, arc):
    P0, P2, c, a = arc
    dist = None
    p = P(p)
    r = (P0 - c).mag()
    if r > 0:
        i = c + (p - c).unit() * r
        alpha = ((i - c).angle() - (P0 - c).angle())
        if a * alpha < 0:
            if alpha > 0:
                alpha = alpha - math.pi2
            else:
                alpha = math.pi2 + alpha
        if between(alpha, 0, a) or min(abs(alpha), abs(alpha - a)) < straight_tolerance:
            return (p - i).mag(), [i.x, i.y]
        else:
            d1, d2 = (p - P0).mag(), (p - P2).mag()
            if d1 < d2:
                return (d1, [P0.x, P0.y])
            else:
                return (d2, [P2.x, P2.y])


def csp_to_arc_distance(sp1, sp2, arc1, arc2, tolerance=0.01):  # arc = [start,end,center,alpha]
    n, i = 10, 0
    d, d1, dl = (0, (0, 0)), (0, (0, 0)), 0
    while i < 1 or (abs(d1[0] - dl[0]) > tolerance and i < 4):
        i += 1
        dl = d1 * 1
        for j in range(n + 1):
            t = float(j) / n
            p = csp_at_t(sp1, sp2, t)
            d = min(point_to_arc_distance(p, arc1), point_to_arc_distance(p, arc2))
            d1 = max(d1, d)
        n = n * 2
    return d1[0]


################################################################################
#    Common functions
################################################################################

def atan2(*arg):
    if len(arg) == 1 and (type(arg[0]) == type([0., 0.]) or type(arg[0]) == type((0., 0.))):
        return (math.pi / 2 - math.atan2(arg[0][0], arg[0][1])) % math.pi2
    elif len(arg) == 2:

        return (math.pi / 2 - math.atan2(arg[0], arg[1])) % math.pi2
    else:
        raise ValueError, "Bad argumets for atan! (%s)" % arg


def between(c, x, y):
    return x - straight_tolerance <= c <= y + straight_tolerance or y - straight_tolerance <= c <= x + straight_tolerance


# Print arguments into specified log file
def print_(*arg):
    f = open(options.log_filename, "a")
    for s in arg:
        s = str(unicode(s).encode('unicode_escape')) + " "
        f.write(s)
    f.write("\n")
    f.close()


################################################################################
#        Point (x,y) operations
################################################################################

class P:
    def __init__(self, x, y=None):
        if not y == None:
            self.x, self.y = float(x), float(y)
        else:
            self.x, self.y = float(x[0]), float(x[1])

    def __add__(self, other):
        return P(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return P(self.x - other.x, self.y - other.y)

    def __neg__(self):
        return P(-self.x, -self.y)

    def __mul__(self, other):
        if isinstance(other, P):
            return self.x * other.x + self.y * other.y
        return P(self.x * other, self.y * other)

    __rmul__ = __mul__

    def __div__(self, other):
        return P(self.x / other, self.y / other)

    def mag(self):
        return math.hypot(self.x, self.y)

    def unit(self):
        h = self.mag()
        if h:
            return self / h
        else:
            return P(0, 0)

    def angle(self):
        return math.atan2(self.y, self.x)

    def __repr__(self):
        return '%f,%f' % (self.x, self.y)

    def l2(self):
        return self.x * self.x + self.y * self.y


################################################################################
#
#        Biarc function
#
#        Calculates biarc approximation of cubic super path segment
#        splits segment if needed or approximates it with straight line
#
################################################################################
def biarc(sp1, sp2, z1, z2, depth=0):
    def biarc_split(sp1, sp2, z1, z2, depth):
        if depth < options.biarc_max_split_depth:
            sp1, sp2, sp3 = csp_split(sp1, sp2)
            l1, l2 = cspseglength(sp1, sp2), cspseglength(sp2, sp3)
            if l1 + l2 == 0:
                zm = z1
            else:
                zm = z1 + (z2 - z1) * l1 / (l1 + l2)
            return biarc(sp1, sp2, z1, zm, depth + 1) + biarc(sp2, sp3, zm, z2, depth + 1)
        else:
            return [[sp1[1], 'line', 0, 0, sp2[1], [z1, z2]]]

    P0, P4 = P(sp1[1]), P(sp2[1])
    TS, TE, v = (P(sp1[2]) - P0), -(P(sp2[0]) - P4), P0 - P4
    tsa, tea, va = TS.angle(), TE.angle(), v.angle()
    if TE.mag() < straight_distance_tolerance and TS.mag() < straight_distance_tolerance:
        # Both tangents are zerro - line straight
        return [[sp1[1], 'line', 0, 0, sp2[1], [z1, z2]]]
    if TE.mag() < straight_distance_tolerance:
        TE = -(TS + v).unit()
        r = TS.mag() / v.mag() * 2
    elif TS.mag() < straight_distance_tolerance:
        TS = -(TE + v).unit()
        r = 1 / (TE.mag() / v.mag() * 2)
    else:
        r = TS.mag() / TE.mag()
    TS, TE = TS.unit(), TE.unit()
    tang_are_parallel = (
            (tsa - tea) % math.pi < straight_tolerance or math.pi - (tsa - tea) % math.pi < straight_tolerance)
    if (tang_are_parallel and
            ((
                     v.mag() < straight_distance_tolerance or TE.mag() < straight_distance_tolerance or TS.mag() < straight_distance_tolerance) or
             1 - abs(TS * v / (TS.mag() * v.mag())) < straight_tolerance)):
        # Both tangents are parallel and start and end are the same - line straight
        # or one of tangents still smaller then tollerance

        # Both tangents and v are parallel - line straight
        return [[sp1[1], 'line', 0, 0, sp2[1], [z1, z2]]]

    c, b, a = v * v, 2 * v * (r * TS + TE), 2 * r * (TS * TE - 1)
    if v.mag() == 0:
        return biarc_split(sp1, sp2, z1, z2, depth)
    asmall, bsmall, csmall = abs(a) < 10 ** -10, abs(b) < 10 ** -10, abs(c) < 10 ** -10
    if asmall and b != 0:
        beta = -c / b
    elif csmall and a != 0:
        beta = -b / a
    elif not asmall:
        discr = b * b - 4 * a * c
        if discr < 0:    raise ValueError, (a, b, c, discr)
        disq = discr ** .5
        beta1 = (-b - disq) / 2 / a
        beta2 = (-b + disq) / 2 / a
        if beta1 * beta2 > 0:    raise ValueError, (a, b, c, disq, beta1, beta2)
        beta = max(beta1, beta2)
    elif asmall and bsmall:
        return biarc_split(sp1, sp2, z1, z2, depth)
    alpha = beta * r
    ab = alpha + beta
    P1 = P0 + alpha * TS
    P3 = P4 - beta * TE
    P2 = (beta / ab) * P1 + (alpha / ab) * P3


    def calculate_arc_params(P0, P1, P2):
        D = (P0 + P2) / 2
        if (D - P1).mag() == 0: return None, None
        R = D - ((D - P0).mag() ** 2 / (D - P1).mag()) * (P1 - D).unit()
        p0a, p1a, p2a = (P0 - R).angle() % (2 * math.pi), (P1 - R).angle() % (2 * math.pi), (P2 - R).angle() % (
                2 * math.pi)
        alpha = (p2a - p0a) % (2 * math.pi)
        if (p0a < p2a and (p1a < p0a or p2a < p1a)) or (p2a < p1a < p0a):
            alpha = -2 * math.pi + alpha
        if abs(R.x) > 1000000 or abs(R.y) > 1000000 or (R - P0).mag < .1:
            return None, None
        else:
            return R, alpha

    R1, a1 = calculate_arc_params(P0, P1, P2)
    R2, a2 = calculate_arc_params(P2, P3, P4)
    if R1 == None or R2 == None or (R1 - P0).mag() < straight_tolerance or (
            R2 - P2).mag() < straight_tolerance: return [[sp1[1], 'line', 0, 0, sp2[1], [z1, z2]]]

    d = csp_to_arc_distance(sp1, sp2, [P0, P2, R1, a1], [P2, P4, R2, a2])
    if d > 1 and depth < options.biarc_max_split_depth:
        return biarc_split(sp1, sp2, z1, z2, depth)
    else:
        if R2.mag() * a2 == 0:
            zm = z2
        else:
            zm = z1 + (z2 - z1) * (abs(R1.mag() * a1)) / (abs(R2.mag() * a2) + abs(R1.mag() * a1))
        return [[sp1[1], 'arc', [R1.x, R1.y], a1, [P2.x, P2.y], [z1, zm]],
                [[P2.x, P2.y], 'arc', [R2.x, R2.y], a2, [P4.x, P4.y], [zm, z2]]]


################################################################################
#        Polygon class
################################################################################
class Polygon:
    def __init__(self, polygon=None):
        self.polygon = [] if polygon == None else polygon[:]

    def add(self, add):
        if type(add) == type([]):
            self.polygon += add[:]
        else:
            self.polygon += add.polygon[:]


class ArrangementGenetic:
    # gene = [fittness, order, rotation, xposition]
    # spieces = [gene]*shapes count
    # population = [spieces]
    def __init__(self, polygons, material_width):
        self.population = []
        self.genes_count = len(polygons)
        self.polygons = polygons
        self.width = material_width
        self.mutation_factor = 0.1
        self.order_mutate_factor = 1.
        self.move_mutate_factor = 1.


################################################################################
###
###        Gcodetools class
###
################################################################################

class LaserGcode(inkex.Effect):

    def export_gcode(self, gcode):
        gcode_pass = gcode
        for x in range(1, self.options.passes):
            gcode += "G91\nG1 Z-" + self.options.pass_depth + "\nG90\n" + gcode_pass
        f = open(self.options.directory + self.options.file, "w")
        f.write(
            self.options.laser_off_command + " S0" + "\n" + self.header +
            "G1 F" + self.options.travel_speed + "\n" + gcode + self.footer)
        f.close()

    def __init__(self):
        inkex.Effect.__init__(self)
        add_option = self.OptionParser.add_option

        add_option("-d", "--directory", action="store", type="string", dest="directory", default="",
                   help="Output directory")

        add_option("-f", "--filename", action="store", type="string", dest="file",
                   default="output.gcode", help="File name")

        add_option("", "--add-numeric-suffix-to-filename", action="store", type="inkbool",
                   dest="add_numeric_suffix_to_filename", default=False,
                   help="Add numeric suffix to file name")

        add_option("", "--laser-command", action="store", type="string", dest="laser_command",
                   default="M03", help="Laser gcode command")

        add_option("", "--laser-off-command", action="store", type="string", dest="laser_off_command",
                   default="M05", help="Laser gcode end command")

        add_option("", "--laser-speed", action="store", type="int", dest="laser_speed", default="750",
                   help="Laser speed (mm/min)")

        add_option("", "--travel-speed", action="store", type="string", dest="travel_speed",
                   default="3000", help="Travel speed (mm/min)")

        add_option("", "--laser-power", action="store", type="int", dest="laser_power", default="255",
                   help="S# is 256 or 10000 for full power")

        add_option("", "--passes", action="store", type="int", dest="passes", default="1",
                   help="Quantity of passes")

        add_option("", "--pass-depth", action="store", type="string", dest="pass_depth", default="1",
                   help="Depth of laser cut")

        add_option("", "--power-delay", action="store", type="string", dest="power_delay",
                   default="0", help="Laser power-on delay (ms)")

        add_option("", "--suppress-all-messages", action="store", type="inkbool",
                   dest="suppress_all_messages", default=True,
                   help="Hide messages during g-code generation")

        add_option("", "--create-log", action="store", type="inkbool", dest="log_create_log",
                   default=False, help="Create log files")

        add_option("", "--log-filename", action="store", type="string", dest="log_filename",
                   default='', help="Create log files")

        add_option("", "--engraving-draw-calculation-paths", action="store", type="inkbool",
                   dest="engraving_draw_calculation_paths", default=False,
                   help="Draw additional graphics to debug engraving path")

        add_option("", "--unit", action="store", type="string", dest="unit",
                   default="G21 (All units in mm)", help="Units either mm or inches")

        add_option("", "--active-tab", action="store", type="string", dest="active_tab", default="",
                   help="Defines which tab is active")

        add_option("", "--biarc-max-split-depth", action="store", type="int",
                   dest="biarc_max_split_depth", default="4",
                   help="Defines maximum depth of splitting while approximating using biarcs.")

    def parse_curve(self, p, layer, w=None, f=None):
        c = []
        if len(p) == 0:
            return []
        p = self.transform_csp(p, layer)

        # Sort to reduce Rapid distance
        k = range(1, len(p))
        keys = [0]
        while len(k) > 0:
            end = p[keys[-1]][-1][1]
            dist = None
            for i in range(len(k)):
                start = p[k[i]][0][1]
                dist = max((-((end[0] - start[0]) ** 2 + (end[1] - start[1]) ** 2), i), dist)
            keys += [k[dist[1]]]
            del k[dist[1]]
        for k in keys:
            subpath = p[k]
            c += [[[subpath[0][1][0], subpath[0][1][1]], 'move', 0, 0]]
            for i in range(1, len(subpath)):
                sp1 = [[subpath[i - 1][j][0], subpath[i - 1][j][1]] for j in range(3)]
                sp2 = [[subpath[i][j][0], subpath[i][j][1]] for j in range(3)]
                c += biarc(sp1, sp2, 0, 0) if w == None else biarc(sp1, sp2, -f(w[k][i - 1]), -f(w[k][i]))
            #                    l1 = biarc(sp1,sp2,0,0) if w==None else biarc(sp1,sp2,-f(w[k][i-1]),-f(w[k][i]))
            #                    print_((-f(w[k][i-1]),-f(w[k][i]), [i1[5] for i1 in l1]) )
            c += [[[subpath[-1][1][0], subpath[-1][1][1]], 'end', 0, 0]]
            print_("Curve: " + str(c))
        return c

    def draw_curve(self, curve, layer, group=None, style=styles["biarc_style"]):

        self.get_defs()
        # Add marker to defs if it does not exist
        if "DrawCurveMarker" not in self.defs:
            defs = inkex.etree.SubElement(self.document.getroot(), inkex.addNS("defs", "svg"))
            marker = inkex.etree.SubElement(defs, inkex.addNS("marker", "svg"),
                                            {"id": "DrawCurveMarker", "orient": "auto", "refX": "-8",
                                             "refY": "-2.41063", "style": "overflow:visible"})
            inkex.etree.SubElement(marker, inkex.addNS("path", "svg"),
                                   {
                                       "d": "m -6.55552,-2.41063 0,0 L -13.11104,0 c 1.0473,-1.42323 1.04126,-3.37047 0,-4.82126",
                                       "style": "fill:#000044; fill-rule:evenodd;stroke-width:0.62500000;stroke-linejoin:round;"}
                                   )
        if "DrawCurveMarker_r" not in self.defs:
            defs = inkex.etree.SubElement(self.document.getroot(), inkex.addNS("defs", "svg"))
            marker = inkex.etree.SubElement(defs, inkex.addNS("marker", "svg"),
                                            {"id": "DrawCurveMarker_r", "orient": "auto", "refX": "8",
                                             "refY": "-2.41063", "style": "overflow:visible"})
            inkex.etree.SubElement(marker, inkex.addNS("path", "svg"),
                                   {
                                       "d": "m 6.55552,-2.41063 0,0 L 13.11104,0 c -1.0473,-1.42323 -1.04126,-3.37047 0,-4.82126",
                                       "style": "fill:#000044; fill-rule:evenodd;stroke-width:0.62500000;stroke-linejoin:round;"}
                                   )
        for i in [0, 1]:
            style['biarc%s_r' % i] = simplestyle.parseStyle(style['biarc%s' % i])
            style['biarc%s_r' % i]["marker-start"] = "url(#DrawCurveMarker_r)"
            del (style['biarc%s_r' % i]["marker-end"])
            style['biarc%s_r' % i] = simplestyle.formatStyle(style['biarc%s_r' % i])

        if group is None:
            group = inkex.etree.SubElement(self.layers[min(1, len(self.layers) - 1)], inkex.addNS('g', 'svg'),
                                           {"gcodetools": "Preview group"})
        s, arcn = '', 0

        a, b, c = [0., 0.], [1., 0.], [0., 1.]
        k = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])
        a, b, c = self.transform(a, layer, True), self.transform(b, layer, True), self.transform(c, layer, True)
        if ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) * k > 0:
            reverse_angle = 1
        else:
            reverse_angle = -1
        for sk in curve:
            si = sk[:]
            si[0], si[2] = self.transform(si[0], layer, True), (
                self.transform(si[2], layer, True) if type(si[2]) == type([]) and len(si[2]) == 2 else si[2])

            if s != '':
                if s[1] == 'line':
                    inkex.etree.SubElement(group, inkex.addNS('path', 'svg'),
                                           {
                                               'style': style['line'],
                                               'd': 'M %s,%s L %s,%s' % (s[0][0], s[0][1], si[0][0], si[0][1]),
                                               "gcodetools": "Preview",
                                           }
                                           )
                elif s[1] == 'arc':
                    arcn += 1
                    sp = s[0]
                    c = s[2]
                    s[3] = s[3] * reverse_angle

                    a = ((P(si[0]) - P(c)).angle() - (P(s[0]) - P(c)).angle()) % math.pi2  # s[3]
                    if s[3] * a < 0:
                        if a > 0:
                            a = a - math.pi2
                        else:
                            a = math.pi2 + a
                    r = math.sqrt((sp[0] - c[0]) ** 2 + (sp[1] - c[1]) ** 2)
                    a_st = (math.atan2(sp[0] - c[0], - (sp[1] - c[1])) - math.pi / 2) % (math.pi * 2)
                    st = style['biarc%s' % (arcn % 2)][:]
                    if a > 0:
                        a_end = a_st + a
                        st = style['biarc%s' % (arcn % 2)]
                    else:
                        a_end = a_st * 1
                        a_st = a_st + a
                        st = style['biarc%s_r' % (arcn % 2)]
                    inkex.etree.SubElement(group, inkex.addNS('path', 'svg'),
                                           {
                                               'style': st,
                                               inkex.addNS('cx', 'sodipodi'): str(c[0]),
                                               inkex.addNS('cy', 'sodipodi'): str(c[1]),
                                               inkex.addNS('rx', 'sodipodi'): str(r),
                                               inkex.addNS('ry', 'sodipodi'): str(r),
                                               inkex.addNS('start', 'sodipodi'): str(a_st),
                                               inkex.addNS('end', 'sodipodi'): str(a_end),
                                               inkex.addNS('open', 'sodipodi'): 'true',
                                               inkex.addNS('type', 'sodipodi'): 'arc',
                                               "gcodetools": "Preview",
                                           })
            s = si


    def check_dir(self):
        if self.options.directory[-1] not in ["/", "\\"]:
            if "\\" in self.options.directory:
                self.options.directory += "\\"
            else:
                self.options.directory += "/"
        print_("Checking direcrory: '%s'" % self.options.directory)
        if (os.path.isdir(self.options.directory)):
            if (os.path.isfile(self.options.directory + 'header')):
                f = open(self.options.directory + 'header', 'r')
                self.header = f.read()
                f.close()
            else:
                self.header = defaults['header']
            if (os.path.isfile(self.options.directory + 'footer')):
                f = open(self.options.directory + 'footer', 'r')
                self.footer = f.read()
                f.close()
            else:
                self.footer = defaults['footer']

            if self.options.unit == "G21 (All units in mm)":
                self.header += "G21\n"
            elif self.options.unit == "G20 (All units in inches)":
                self.header += "G20\n"
        else:
            self.error(_("Directory does not exist! Please specify existing directory at options tab!"), "error")
            return False

        if self.options.add_numeric_suffix_to_filename:
            dir_list = os.listdir(self.options.directory)
            if "." in self.options.file:
                r = re.match(r"^(.*)(\..*)$", self.options.file)
                ext = r.group(2)
                name = r.group(1)
            else:
                ext = ""
                name = self.options.file
            max_n = 0
            for s in dir_list:
                r = re.match(r"^%s_0*(\d+)%s$" % (re.escape(name), re.escape(ext)), s)
                if r:
                    max_n = max(max_n, int(r.group(1)))
            filename = name + "_" + ("0" * (4 - len(str(max_n + 1))) + str(max_n + 1)) + ext
            self.options.file = filename

        print_("Testing writing rights on '%s'" % (self.options.directory + self.options.file))
        try:
            f = open(self.options.directory + self.options.file, "w")
            f.close()
        except:
            self.error(_("Can not write to specified file!\n%s" % (self.options.directory + self.options.file)),
                       "error")
            return False
        return True


    ################################################################################
    #
    #       Generate Gcode
    #
    #       Curve definition
    #       [start point, type = {'arc','line','move','end'}, arc center, arc angle, end point, [zstart, zend]]
    #
    ################################################################################

    def generate_gcode(self, curve, layer, depth):
        tool = self.tools
        print_("Tool in g-code generator: " + str(tool))

        def c(c):
            c = [c[i] if i < len(c) else None for i in range(6)]
            if c[5] == 0: c[5] = None
            s = [" X", " Y", " Z", " I", " J", " K"]
            r = ''
            for i in range(6):
                if c[i] != None:
                    r += s[i] + ("%f" % (round(c[i], 4))).rstrip('0')
            return r


        if len(curve) == 0: return ""

        try:
            self.last_used_tool == None
        except:
            self.last_used_tool = None
        print_("working on curve")
        print_("Curve: " + str(curve))
        g = ""

        lg, f = 'G00', "F%f" % tool['penetration feed']
        penetration_feed = "F%s" % tool['penetration feed']
        current_a = 0
        for i in range(1, len(curve)):
            #    Creating Gcode for curve between s=curve[i-1] and si=curve[i] start at s[0] end at s[4]=si[0]
            s, si = curve[i - 1], curve[i]
            feed = f if lg not in ['G01', 'G02', 'G03'] else ''
            if s[1] == 'move':
                g += "G1 " + c(si[0]) + "\n" + tool['gcode before path'] + "\n"
                lg = 'G00'
            elif s[1] == 'end':
                g += tool['gcode after path'] + "\n"
                lg = 'G00'
            elif s[1] == 'line':
                if lg == "G00": g += "G1 " + feed + "\n"
                g += "G1 " + c(si[0]) + "\n"
                lg = 'G01'
            elif s[1] == 'arc':
                r = [(s[2][0] - s[0][0]), (s[2][1] - s[0][1])]
                if lg == "G00": g += "G1 " + feed + "\n"
                if (r[0] ** 2 + r[1] ** 2) > .1:
                    r1, r2 = (P(s[0]) - P(s[2])), (P(si[0]) - P(s[2]))
                    if abs(r1.mag() - r2.mag()) < 0.001:
                        g += ("G2" if s[3] < 0 else "G3") + c(
                            si[0] + [None, (s[2][0] - s[0][0]), (s[2][1] - s[0][1])]) + "\n"
                    else:
                        r = (r1.mag() + r2.mag()) / 2
                        g += ("G2" if s[3] < 0 else "G3") + c(si[0]) + " R%f" % (r) + "\n"
                    lg = 'G02'
                else:
                    g += "G1 " + c(si[0]) + " " + feed + "\n"
                    lg = 'G01'
        if si[1] == 'end':
            g += tool['gcode after path'] + "\n"
        return g


    def get_transforms(self, g):
        root = self.document.getroot()
        trans = []
        while (g != root):
            if 'transform' in g.keys():
                t = g.get('transform')
                t = simpletransform.parseTransform(t)
                trans = simpletransform.composeTransform(t, trans) if trans != [] else t
                print_(trans)
            g = g.getparent()
        return trans


    def apply_transforms(self, g, csp):
        trans = self.get_transforms(g)
        if trans != []:
            simpletransform.applyTransformToPath(trans, csp)
        return csp


    def transform(self, source_point, layer, reverse=False):
        if layer == None:
            layer = self.current_layer if self.current_layer is not None else self.document.getroot()
        if layer not in self.transform_matrix:
            for i in range(self.layers.index(layer), -1, -1):
                if self.layers[i] in self.orientation_points:
                    break

            print_(str(self.layers))
            print_(str("I: " + str(i)))
            print_("Transform: " + str(self.layers[i]))
            if self.layers[i] not in self.orientation_points:
                self.error(_(
                    "Orientation points for '%s' layer have not been found! Please add orientation points using Orientation tab!") % layer.get(
                    inkex.addNS('label', 'inkscape')), "no_orientation_points")
            elif self.layers[i] in self.transform_matrix:
                self.transform_matrix[layer] = self.transform_matrix[self.layers[i]]
            else:
                orientation_layer = self.layers[i]
                if len(self.orientation_points[orientation_layer]) > 1:
                    self.error(
                        _("There are more than one orientation point groups in '%s' layer") % orientation_layer.get(
                            inkex.addNS('label', 'inkscape')), "more_than_one_orientation_point_groups")
                points = self.orientation_points[orientation_layer][0]
                if len(points) == 2:
                    points += [[[(points[1][0][1] - points[0][0][1]) + points[0][0][0],
                                 -(points[1][0][0] - points[0][0][0]) + points[0][0][1]],
                                [-(points[1][1][1] - points[0][1][1]) + points[0][1][0],
                                 points[1][1][0] - points[0][1][0] + points[0][1][1]]]]
                if len(points) == 3:
                    print_("Layer '%s' Orientation points: " % orientation_layer.get(inkex.addNS('label', 'inkscape')))
                    for point in points:
                        print_(point)
                    #    Zcoordinates definition taken from Orientatnion point 1 and 2
                    self.Zcoordinates[layer] = [max(points[0][1][2], points[1][1][2]),
                                                min(points[0][1][2], points[1][1][2])]
                    matrix = numpy.array([
                        [points[0][0][0], points[0][0][1], 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, points[0][0][0], points[0][0][1], 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, points[0][0][0], points[0][0][1], 1],
                        [points[1][0][0], points[1][0][1], 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, points[1][0][0], points[1][0][1], 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, points[1][0][0], points[1][0][1], 1],
                        [points[2][0][0], points[2][0][1], 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, points[2][0][0], points[2][0][1], 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, points[2][0][0], points[2][0][1], 1]
                    ])

                    if numpy.linalg.det(matrix) != 0:
                        m = numpy.linalg.solve(matrix,
                                               numpy.array(
                                                   [[points[0][1][0]], [points[0][1][1]], [1], [points[1][1][0]],
                                                    [points[1][1][1]], [1], [points[2][1][0]], [points[2][1][1]], [1]]
                                               )
                                               ).tolist()
                        self.transform_matrix[layer] = [[m[j * 3 + i][0] for i in range(3)] for j in range(3)]

                    else:
                        self.error(_(
                            "Orientation points are wrong! (if there are two orientation points they sould not be the same. If there are three orientation points they should not be in a straight line.)"),
                            "wrong_orientation_points")
                else:
                    self.error(_(
                        "Orientation points are wrong! (if there are two orientation points they sould not be the same. If there are three orientation points they should not be in a straight line.)"),
                        "wrong_orientation_points")

            self.transform_matrix_reverse[layer] = numpy.linalg.inv(self.transform_matrix[layer]).tolist()
            print_("\n Layer '%s' transformation matrixes:" % layer.get(inkex.addNS('label', 'inkscape')))
            print_(self.transform_matrix)
            print_(self.transform_matrix_reverse)

            # Zautoscale is absolute
            self.Zauto_scale[layer] = 1
            print_("Z automatic scale = %s (computed according orientation points)" % self.Zauto_scale[layer])

        x, y = source_point[0], source_point[1]
        if not reverse:
            t = self.transform_matrix[layer]
        else:
            t = self.transform_matrix_reverse[layer]
        return [t[0][0] * x + t[0][1] * y + t[0][2], t[1][0] * x + t[1][1] * y + t[1][2]]


    def transform_csp(self, csp_, layer, reverse=False):
        csp = [[[csp_[i][j][0][:], csp_[i][j][1][:], csp_[i][j][2][:]] for j in range(len(csp_[i]))] for i in
               range(len(csp_))]
        for i in xrange(len(csp)):
            for j in xrange(len(csp[i])):
                for k in xrange(len(csp[i][j])):
                    csp[i][j][k] = self.transform(csp[i][j][k], layer, reverse)
        return csp

    ################################################################################
    #        Errors handling function, notes are just printed into Logfile,
    #        warnings are printed into log file and warning message is displayed but
    #        extension continues working, errors causes log and execution is halted
    #        Notes, warnings adn errors could be assigned to space or comma or dot
    #        sepparated strings (case is ignoreg).
    ################################################################################
    def error(self, s, type_="Warning"):
        notes = "Note "
        warnings = """
                        Warning tools_warning
                        bad_orientation_points_in_some_layers
                        more_than_one_orientation_point_groups
                        more_than_one_tool
                        orientation_have_not_been_defined
                        tool_have_not_been_defined
                        selection_does_not_contain_paths
                        selection_does_not_contain_paths_will_take_all
                        selection_is_empty_will_comupe_drawing
                        selection_contains_objects_that_are_not_paths
                        """
        errors = """
                        Error
                        wrong_orientation_points
                        area_tools_diameter_error
                        no_tool_error
                        active_layer_already_has_tool
                        active_layer_already_has_orientation_points
                    """
        if type_.lower() in re.split("[\s\n,\.]+", errors.lower()):
            print_(s)
            inkex.errormsg(s + "\n")
            sys.exit()
        elif type_.lower() in re.split("[\s\n,\.]+", warnings.lower()):
            print_(s)
            if not self.options.suppress_all_messages:
                inkex.errormsg(s + "\n")
        elif type_.lower() in re.split("[\s\n,\.]+", notes.lower()):
            print_(s)
        else:
            print_(s)
            inkex.errormsg(s)
            sys.exit()

    ################################################################################
    #        Get defs from svg
    ################################################################################
    def get_defs(self):
        self.defs = {}

        def recursive(g):
            for i in g:
                if i.tag == inkex.addNS("defs", "svg"):
                    for j in i:
                        self.defs[j.get("id")] = i
                if i.tag == inkex.addNS("g", 'svg'):
                    recursive(i)

        recursive(self.document.getroot())


    ################################################################################
    #
    #        Get Gcodetools info from the svg
    #
    ################################################################################
    def get_info(self):
        self.selected_paths = {}
        self.paths = {}
        self.orientation_points = {}
        self.layers = [self.document.getroot()]
        self.Zcoordinates = {}
        self.transform_matrix = {}
        self.transform_matrix_reverse = {}
        self.Zauto_scale = {}

        def recursive_search(g, layer, selected=False):
            items = g.getchildren()
            items.reverse()
            for i in items:
                if selected:
                    self.selected[i.get("id")] = i
                if i.tag == inkex.addNS("g", 'svg') and i.get(inkex.addNS('groupmode', 'inkscape')) == 'layer':
                    self.layers += [i]
                    recursive_search(i, i)
                elif i.get('gcodetools') == "Gcodetools orientation group":
                    points = self.get_orientation_points(i)
                    if points != None:
                        self.orientation_points[layer] = self.orientation_points[layer] + [
                            points[:]] if layer in self.orientation_points else [points[:]]
                        print_("Found orientation points in '%s' layer: %s" % (
                            layer.get(inkex.addNS('label', 'inkscape')), points))
                    else:
                        self.error(_(
                            "Warning! Found bad orientation points in '%s' layer. Resulting Gcode could be corrupt!") % layer.get(
                            inkex.addNS('label', 'inkscape')), "bad_orientation_points_in_some_layers")
                elif i.tag == inkex.addNS('path', 'svg'):
                    if "gcodetools" not in i.keys():
                        self.paths[layer] = self.paths[layer] + [i] if layer in self.paths else [i]
                        if i.get("id") in self.selected:
                            self.selected_paths[layer] = self.selected_paths[layer] + [
                                i] if layer in self.selected_paths else [i]
                elif i.tag == inkex.addNS("g", 'svg'):
                    recursive_search(i, layer, (i.get("id") in self.selected))
                elif i.get("id") in self.selected:
                    self.error(_(
                        "This extension works with Paths and Dynamic Offsets and groups of them only! All other objects will be ignored!\nSolution 1: press Path->Object to path or Shift+Ctrl+C.\nSolution 2: Path->Dynamic offset or Ctrl+J.\nSolution 3: export all contours to PostScript level 2 (File->Save As->.ps) and File->Import this file."),
                        "selection_contains_objects_that_are_not_paths")


        recursive_search(self.document.getroot(), self.document.getroot())

    def get_orientation_points(self, g):
        items = g.getchildren()
        items.reverse()
        p2, p3 = [], []
        p = None
        for i in items:
            if i.tag == inkex.addNS("g", 'svg') and i.get("gcodetools") == "Gcodetools orientation point (2 points)":
                p2 += [i]
            if i.tag == inkex.addNS("g", 'svg') and i.get("gcodetools") == "Gcodetools orientation point (3 points)":
                p3 += [i]
        if len(p2) == 2:
            p = p2
        elif len(p3) == 3:
            p = p3
        if p == None: return None
        points = []
        for i in p:
            point = [[], []]
            for node in i:
                if node.get('gcodetools') == "Gcodetools orientation point arrow":
                    point[0] = self.apply_transforms(node, cubicsuperpath.parsePath(node.get("d")))[0][0][1]
                if node.get('gcodetools') == "Gcodetools orientation point text":
                    r = re.match(
                        r'(?i)\s*\(\s*(-?\s*\d*(?:,|\.)*\d*)\s*;\s*(-?\s*\d*(?:,|\.)*\d*)\s*;\s*(-?\s*\d*(?:,|\.)*\d*)\s*\)\s*',
                        node.text)
                    point[1] = [float(r.group(1)), float(r.group(2)), float(r.group(3))]
            if point[0] != [] and point[1] != []:    points += [point]
        if len(points) == len(p2) == 2 or len(points) == len(p3) == 3:
            return points
        else:
            return None

    ################################################################################
    #
    #        dxfpoints
    #
    ################################################################################
    def dxfpoints(self):
        if self.selected_paths == {}:
            self.error(_(
                "Noting is selected. Please select something to convert to drill point (dxfpoint) or clear point sign."),
                "warning")
        for layer in self.layers:
            if layer in self.selected_paths:
                for path in self.selected_paths[layer]:
                    if self.options.dxfpoints_action == 'replace':
                        path.set("dxfpoint", "1")
                        r = re.match("^\s*.\s*(\S+)", path.get("d"))
                        if r != None:
                            print_(("got path=", r.group(1)))
                            path.set("d",
                                     "m %s 2.9375,-6.343750000001 0.8125,1.90625 6.843748640396,-6.84374864039 0,0 0.6875,0.6875 -6.84375,6.84375 1.90625,0.812500000001 z" % r.group(
                                         1))
                            path.set("style", styles["dxf_points"])

                    if self.options.dxfpoints_action == 'save':
                        path.set("dxfpoint", "1")

                    if self.options.dxfpoints_action == 'clear' and path.get("dxfpoint") == "1":
                        path.set("dxfpoint", "0")

    ################################################################################
    #
    #        Laser
    #
    ################################################################################
    def laser(self):

        def get_boundaries(points):
            minx, miny, maxx, maxy = None, None, None, None
            out = [[], [], [], []]
            for p in points:
                if minx == p[0]:
                    out[0] += [p]
                if minx == None or p[0] < minx:
                    minx = p[0]
                    out[0] = [p]

                if miny == p[1]:
                    out[1] += [p]
                if miny == None or p[1] < miny:
                    miny = p[1]
                    out[1] = [p]

                if maxx == p[0]:
                    out[2] += [p]
                if maxx == None or p[0] > maxx:
                    maxx = p[0]
                    out[2] = [p]

                if maxy == p[1]:
                    out[3] += [p]
                if maxy == None or p[1] > maxy:
                    maxy = p[1]
                    out[3] = [p]
            return out


        def remove_duplicates(points):
            i = 0
            out = []
            for p in points:
                for j in xrange(i, len(points)):
                    if p == points[j]: points[j] = [None, None]
                if p != [None, None]: out += [p]
            i += 1
            return (out)


        def get_way_len(points):
            l = 0
            for i in xrange(1, len(points)):
                l += math.sqrt((points[i][0] - points[i - 1][0]) ** 2 + (points[i][1] - points[i - 1][1]) ** 2)
            return l


        def sort_dxfpoints(points):
            points = remove_duplicates(points)

            ways = [
                # l=0, d=1, r=2, u=3
                [3, 0],  # ul
                [3, 2],  # ur
                [1, 0],  # dl
                [1, 2],  # dr
                [0, 3],  # lu
                [0, 1],  # ld
                [2, 3],  # ru
                [2, 1],  # rd
            ]

            minimal_way = []
            minimal_len = None
            minimal_way_type = None
            for w in ways:
                tpoints = points[:]
                cw = []
                for j in xrange(0, len(points)):
                    p = get_boundaries(get_boundaries(tpoints)[w[0]])[w[1]]
                    tpoints.remove(p[0])
                    cw += p
                curlen = get_way_len(cw)
                if minimal_len == None or curlen < minimal_len:
                    minimal_len = curlen
                    minimal_way = cw
                    minimal_way_type = w

            return minimal_way

        if self.selected_paths == {}:
            paths = self.paths
            self.error(_("No paths are selected! Trying to work on all available paths."), "warning")
        else:
            paths = self.selected_paths

        self.check_dir()
        gcode = ""

        biarc_group = inkex.etree.SubElement(
            self.selected_paths.keys()[0] if len(self.selected_paths.keys()) > 0 else self.layers[0],
            inkex.addNS('g', 'svg'))
        print_(("self.layers=", self.layers))
        print_(("paths=", paths))
        for layer in self.layers:
            if layer in paths:
                print_(("layer", layer))
                p = []
                dxfpoints = []
                for path in paths[layer]:
                    print_(str(layer))
                    if "d" not in path.keys():
                        self.error(_(
                            "Warning: One or more paths dont have 'd' parameter, try to Ungroup (Ctrl+Shift+G) and Object to Path (Ctrl+Shift+C)!"),
                            "selection_contains_objects_that_are_not_paths")
                        continue
                    csp = cubicsuperpath.parsePath(path.get("d"))
                    csp = self.apply_transforms(path, csp)
                    if path.get("dxfpoint") == "1":
                        tmp_curve = self.transform_csp(csp, layer)
                        x = tmp_curve[0][0][0][0]
                        y = tmp_curve[0][0][0][1]
                        print_("got dxfpoint (scaled) at (%f,%f)" % (x, y))
                        dxfpoints += [[x, y]]
                    else:
                        p += csp
                dxfpoints = sort_dxfpoints(dxfpoints)
                curve = self.parse_curve(p, layer)
                self.draw_curve(curve, layer, biarc_group)
                gcode += self.generate_gcode(curve, layer, 0)

        self.export_gcode(gcode)

    ################################################################################
    #
    #        Orientation
    #
    ################################################################################
    def orientation(self, layer=None):
        print_("entering orientations")
        if layer == None:
            layer = self.current_layer if self.current_layer is not None else self.document.getroot()
        if layer in self.orientation_points:
            self.error(_("Active layer already has orientation points! Remove them or select another layer!"),
                       "active_layer_already_has_orientation_points")

        orientation_group = inkex.etree.SubElement(layer, inkex.addNS('g', 'svg'),
                                                   {"gcodetools": "Gcodetools orientation group"})

        # translate == ['0', '-917.7043']
        if layer.get("transform") != None:
            translate = layer.get("transform").replace("translate(", "").replace(")", "").split(",")
        else:
            translate = [0, 0]

        # doc height in pixels (38 mm == 143.62204724px)
        doc_height = self.unittouu(self.document.getroot().xpath('@height', namespaces=inkex.NSS)[0])

        if self.document.getroot().get('height') == "100%":
            doc_height = 1052.3622047
            print_("Overruding height from 100 percents to %s" % doc_height)

        print_("Document height: " + str(doc_height));

        if self.options.unit == "G21 (All units in mm)":
            points = [[0., 0., 0.], [100., 0., 0.], [0., 100., 0.]]
            orientation_scale = 1
            print_("orientation_scale < 0 ===> switching to mm units=%0.10f" % orientation_scale)
        elif self.options.unit == "G20 (All units in inches)":
            points = [[0., 0., 0.], [5., 0., 0.], [0., 5., 0.]]
            orientation_scale = 90
            print_("orientation_scale < 0 ===> switching to inches units=%0.10f" % orientation_scale)

        points = points[:2]

        print_(("using orientation scale", orientation_scale, "i=", points))
        for i in points:
            # X == Correct!
            # si == x,y coordinate in px
            # si have correct coordinates
            # if layer have any tranform it will be in translate so lets add that
            si = [i[0] * orientation_scale, (i[1] * orientation_scale) + float(translate[1])]
            g = inkex.etree.SubElement(orientation_group, inkex.addNS('g', 'svg'),
                                       {'gcodetools': "Gcodetools orientation point (2 points)"})
            inkex.etree.SubElement(g, inkex.addNS('path', 'svg'),
                                   {
                                       'style': "stroke:none;fill:#000000;",
                                       'd': 'm %s,%s 2.9375,-6.343750000001 0.8125,1.90625 6.843748640396,-6.84374864039 0,0 0.6875,0.6875 -6.84375,6.84375 1.90625,0.812500000001 z z' % (
                                           si[0], -si[1] + doc_height),
                                       'gcodetools': "Gcodetools orientation point arrow"
                                   })
            t = inkex.etree.SubElement(g, inkex.addNS('text', 'svg'),
                                       {
                                           'style': "font-size:10px;font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;fill:#000000;fill-opacity:1;stroke:none;",
                                           inkex.addNS("space", "xml"): "preserve",
                                           'x': str(si[0] + 10),
                                           'y': str(-si[1] - 10 + doc_height),
                                           'gcodetools': "Gcodetools orientation point text"
                                       })
            t.text = "(%s; %s; %s)" % (i[0], i[1], i[2])


    ################################################################################
    #
    #        Effect
    #
    #        Main function of Gcodetools class
    #
    ################################################################################
    def effect(self):
        global options
        options = self.options
        options.self = self
        options.doc_root = self.document.getroot()
        # define print_ function
        global print_
        if self.options.log_create_log:
            try:
                if os.path.isfile(self.options.log_filename): os.remove(self.options.log_filename)
                f = open(self.options.log_filename, "a")
                f.write("Gcodetools log file.\nStarted at %s.\n%s\n" % (
                    time.strftime("%d.%m.%Y %H:%M:%S"), options.log_filename))
                f.write("%s tab is active.\n" % self.options.active_tab)
                f.close()
            except:
                print_ = lambda *x: None
        else:
            print_ = lambda *x: None
        self.get_info()
        if self.orientation_points == {}:
            self.error(_(
                "Orientation points have not been defined! A default set of orientation points has been automatically added."),
                "warning")
            self.orientation(self.layers[min(0, len(self.layers) - 1)])
            self.get_info()

        self.tools = {
            "name": "Laser Engraver",
            "id": "Laser Engraver",
            "penetration feed": self.options.laser_speed,
            "feed": self.options.laser_speed,
            "gcode before path": ("G4 P0 \n" + self.options.laser_command + " S" + str(
                int(self.options.laser_power)) + "\nG4 P" + self.options.power_delay),
            "gcode after path": (
                    "G4 P0 \n" + self.options.laser_off_command + " S0" + "\n" + "G1 F" + self.options.travel_speed),
        }

        self.get_info()
        self.laser()


e = LaserGcode()
e.affect()
