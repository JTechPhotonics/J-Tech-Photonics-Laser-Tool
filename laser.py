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
import inkex, simplestyle, simplepath
import cubicsuperpath, simpletransform, bezmisc

import os
import math
import bezmisc
import re
import copy
import sys
import time
import cmath
import numpy
import codecs
import random
import gettext

_ = gettext.gettext

### Check if inkex has errormsg (0.46 version doesnot have one.) Could be removed later.
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


def ireplace(self, old, new, count=0):
    pattern = re.compile(re.escape(old), re.I)
    return re.sub(pattern, new, self, count)


################################################################################
###
###        Styles and additional parameters
###
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
###        Cubic Super Path additional functions
################################################################################

def csp_simple_bound(csp):
    minx, miny, maxx, maxy = None, None, None, None
    for subpath in csp:
        for sp in subpath:
            for p in sp:
                minx = min(minx, p[0]) if minx != None else p[0]
                miny = min(miny, p[1]) if miny != None else p[1]
                maxx = max(maxx, p[0]) if maxx != None else p[0]
                maxy = max(maxy, p[1]) if maxy != None else p[1]
    return minx, miny, maxx, maxy


def csp_segment_to_bez(sp1, sp2):
    return sp1[1:] + sp2[:2]


def csp_to_point_distance(csp, p, dist_bounds=[0, 1e100], tolerance=.01):
    min_dist = [1e100, 0, 0, 0]
    for j in range(len(csp)):
        for i in range(1, len(csp[j])):
            d = csp_seg_to_point_distance(csp[j][i - 1], csp[j][i], p, sample_points=5, tolerance=.01)
            if d[0] < dist_bounds[0]:
                #                draw_pointer( list(csp_at_t(subpath[dist[2]-1],subpath[dist[2]],dist[3]))
                #                    +list(csp_at_t(csp[dist[4]][dist[5]-1],csp[dist[4]][dist[5]],dist[6])),"red","line", comment = math.sqrt(dist[0]))
                return [d[0], j, i, d[1]]
            else:
                if d[0] < min_dist[0]: min_dist = [d[0], j, i, d[1]]
    return min_dist


def csp_seg_to_point_distance(sp1, sp2, p, sample_points=5, tolerance=.01):
    ax, ay, bx, by, cx, cy, dx, dy = csp_parameterize(sp1, sp2)
    dx, dy = dx - p[0], dy - p[1]
    if sample_points < 2: sample_points = 2
    d = min([(p[0] - sp1[1][0]) ** 2 + (p[1] - sp1[1][1]) ** 2, 0.],
            [(p[0] - sp2[1][0]) ** 2 + (p[1] - sp2[1][1]) ** 2, 1.])
    for k in range(sample_points):
        t = float(k) / (sample_points - 1)
        i = 0
        while i == 0 or abs(f) > 0.000001 and i < 20:
            t2, t3 = t ** 2, t ** 3
            f = (ax * t3 + bx * t2 + cx * t + dx) * (3 * ax * t2 + 2 * bx * t + cx) + (
                        ay * t3 + by * t2 + cy * t + dy) * (3 * ay * t2 + 2 * by * t + cy)
            df = (6 * ax * t + 2 * bx) * (ax * t3 + bx * t2 + cx * t + dx) + (3 * ax * t2 + 2 * bx * t + cx) ** 2 + (
                        6 * ay * t + 2 * by) * (ay * t3 + by * t2 + cy * t + dy) + (3 * ay * t2 + 2 * by * t + cy) ** 2
            if df != 0:
                t = t - f / df
            else:
                break
            i += 1
        if 0 <= t <= 1:
            p1 = csp_at_t(sp1, sp2, t)
            d1 = (p1[0] - p[0]) ** 2 + (p1[1] - p[1]) ** 2
            if d1 < d[0]:
                d = [d1, t]
    return d


def csp_seg_to_csp_seg_distance(sp1, sp2, sp3, sp4, dist_bounds=[0, 1e100], sample_points=5, tolerance=.01):
    # check the ending points first
    dist = csp_seg_to_point_distance(sp1, sp2, sp3[1], sample_points, tolerance)
    dist += [0.]
    if dist[0] <= dist_bounds[0]: return dist
    d = csp_seg_to_point_distance(sp1, sp2, sp4[1], sample_points, tolerance)
    if d[0] < dist[0]:
        dist = d + [1.]
        if dist[0] <= dist_bounds[0]: return dist
    d = csp_seg_to_point_distance(sp3, sp4, sp1[1], sample_points, tolerance)
    if d[0] < dist[0]:
        dist = [d[0], 0., d[1]]
        if dist[0] <= dist_bounds[0]: return dist
    d = csp_seg_to_point_distance(sp3, sp4, sp2[1], sample_points, tolerance)
    if d[0] < dist[0]:
        dist = [d[0], 1., d[1]]
        if dist[0] <= dist_bounds[0]: return dist
    sample_points -= 2
    if sample_points < 1: sample_points = 1
    ax1, ay1, bx1, by1, cx1, cy1, dx1, dy1 = csp_parameterize(sp1, sp2)
    ax2, ay2, bx2, by2, cx2, cy2, dx2, dy2 = csp_parameterize(sp3, sp4)
    #    try to find closes points using Newtons method
    for k in range(sample_points):
        for j in range(sample_points):
            t1, t2 = float(k + 1) / (sample_points + 1), float(j) / (sample_points + 1)
            t12, t13, t22, t23 = t1 * t1, t1 * t1 * t1, t2 * t2, t2 * t2 * t2
            i = 0
            F1, F2, F = [0, 0], [[0, 0], [0, 0]], 1e100
            x, y = ax1 * t13 + bx1 * t12 + cx1 * t1 + dx1 - (
                        ax2 * t23 + bx2 * t22 + cx2 * t2 + dx2), ay1 * t13 + by1 * t12 + cy1 * t1 + dy1 - (
                               ay2 * t23 + by2 * t22 + cy2 * t2 + dy2)
            while i < 2 or abs(F - Flast) > tolerance and i < 30:
                # draw_pointer(csp_at_t(sp1,sp2,t1))
                f1x = 3 * ax1 * t12 + 2 * bx1 * t1 + cx1
                f1y = 3 * ay1 * t12 + 2 * by1 * t1 + cy1
                f2x = 3 * ax2 * t22 + 2 * bx2 * t2 + cx2
                f2y = 3 * ay2 * t22 + 2 * by2 * t2 + cy2
                F1[0] = 2 * f1x * x + 2 * f1y * y
                F1[1] = -2 * f2x * x - 2 * f2y * y
                F2[0][0] = 2 * (6 * ax1 * t1 + 2 * bx1) * x + 2 * f1x * f1x + 2 * (
                            6 * ay1 * t1 + 2 * by1) * y + 2 * f1y * f1y
                F2[0][1] = -2 * f1x * f2x - 2 * f1y * f2y
                F2[1][0] = -2 * f2x * f1x - 2 * f2y * f1y
                F2[1][1] = -2 * (6 * ax2 * t2 + 2 * bx2) * x + 2 * f2x * f2x - 2 * (
                            6 * ay2 * t2 + 2 * by2) * y + 2 * f2y * f2y
                F2 = inv_2x2(F2)
                if F2 != None:
                    t1 -= (F2[0][0] * F1[0] + F2[0][1] * F1[1])
                    t2 -= (F2[1][0] * F1[0] + F2[1][1] * F1[1])
                    t12, t13, t22, t23 = t1 * t1, t1 * t1 * t1, t2 * t2, t2 * t2 * t2
                    x, y = ax1 * t13 + bx1 * t12 + cx1 * t1 + dx1 - (
                                ax2 * t23 + bx2 * t22 + cx2 * t2 + dx2), ay1 * t13 + by1 * t12 + cy1 * t1 + dy1 - (
                                       ay2 * t23 + by2 * t22 + cy2 * t2 + dy2)
                    Flast = F
                    F = x * x + y * y
                else:
                    break
                i += 1
            if F < dist[0] and 0 <= t1 <= 1 and 0 <= t2 <= 1:
                dist = [F, t1, t2]
                if dist[0] <= dist_bounds[0]:
                    return dist
    return dist


def csp_to_csp_distance(csp1, csp2, dist_bounds=[0, 1e100], tolerance=.01):
    dist = [1e100, 0, 0, 0, 0, 0, 0]
    for i1 in range(len(csp1)):
        for j1 in range(1, len(csp1[i1])):
            for i2 in range(len(csp2)):
                for j2 in range(1, len(csp2[i2])):
                    d = csp_seg_bound_to_csp_seg_bound_max_min_distance(csp1[i1][j1 - 1], csp1[i1][j1],
                                                                        csp2[i2][j2 - 1], csp2[i2][j2])
                    if d[0] >= dist_bounds[1]: continue
                    if d[1] < dist_bounds[0]: return [d[1], i1, j1, 1, i2, j2, 1]
                    d = csp_seg_to_csp_seg_distance(csp1[i1][j1 - 1], csp1[i1][j1], csp2[i2][j2 - 1], csp2[i2][j2],
                                                    dist_bounds, tolerance=tolerance)
                    if d[0] < dist[0]:
                        dist = [d[0], i1, j1, d[1], i2, j2, d[2]]
                    if dist[0] <= dist_bounds[0]:
                        return dist
            if dist[0] >= dist_bounds[1]:
                return dist
    return dist


#    draw_pointer( list(csp_at_t(csp1[dist[1]][dist[2]-1],csp1[dist[1]][dist[2]],dist[3]))
#                + list(csp_at_t(csp2[dist[4]][dist[5]-1],csp2[dist[4]][dist[5]],dist[6])), "#507","line")


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


def csp_true_bounds(csp):
    # Finds minx,miny,maxx,maxy of the csp and return their (x,y,i,j,t)
    minx = [float("inf"), 0, 0, 0]
    maxx = [float("-inf"), 0, 0, 0]
    miny = [float("inf"), 0, 0, 0]
    maxy = [float("-inf"), 0, 0, 0]
    for i in range(len(csp)):
        for j in range(1, len(csp[i])):
            ax, ay, bx, by, cx, cy, x0, y0 = bezmisc.bezierparameterize(
                (csp[i][j - 1][1], csp[i][j - 1][2], csp[i][j][0], csp[i][j][1]))
            roots = cubic_solver(0, 3 * ax, 2 * bx, cx) + [0, 1]
            for root in roots:
                if type(root) is complex and abs(root.imag) < 1e-10:
                    root = root.real
                if type(root) is not complex and 0 <= root <= 1:
                    y = ay * (root ** 3) + by * (root ** 2) + cy * root + y0
                    x = ax * (root ** 3) + bx * (root ** 2) + cx * root + x0
                    maxx = max([x, y, i, j, root], maxx)
                    minx = min([x, y, i, j, root], minx)

            roots = cubic_solver(0, 3 * ay, 2 * by, cy) + [0, 1]
            for root in roots:
                if type(root) is complex and root.imag == 0:
                    root = root.real
                if type(root) is not complex and 0 <= root <= 1:
                    y = ay * (root ** 3) + by * (root ** 2) + cy * root + y0
                    x = ax * (root ** 3) + bx * (root ** 2) + cx * root + x0
                    maxy = max([y, x, i, j, root], maxy)
                    miny = min([y, x, i, j, root], miny)
    maxy[0], maxy[1] = maxy[1], maxy[0]
    miny[0], miny[1] = miny[1], miny[0]

    return minx, miny, maxx, maxy


############################################################################
### csp_segments_intersection(sp1,sp2,sp3,sp4)
###
### Returns array containig all intersections between two segmets of cubic
### super path. Results are [ta,tb], or [ta0, ta1, tb0, tb1, "Overlap"]
### where ta, tb are values of t for the intersection point.
############################################################################
def csp_segments_intersection(sp1, sp2, sp3, sp4):
    a, b = csp_segment_to_bez(sp1, sp2), csp_segment_to_bez(sp3, sp4)

    def polish_intersection(a, b, ta, tb, tolerance=intersection_tolerance):
        ax, ay, bx, by, cx, cy, dx, dy = bezmisc.bezierparameterize(a)
        ax1, ay1, bx1, by1, cx1, cy1, dx1, dy1 = bezmisc.bezierparameterize(b)
        i = 0
        F, F1 = [.0, .0], [[.0, .0], [.0, .0]]
        while i == 0 or (abs(F[0]) ** 2 + abs(F[1]) ** 2 > tolerance and i < 10):
            ta3, ta2, tb3, tb2 = ta ** 3, ta ** 2, tb ** 3, tb ** 2
            F[0] = ax * ta3 + bx * ta2 + cx * ta + dx - ax1 * tb3 - bx1 * tb2 - cx1 * tb - dx1
            F[1] = ay * ta3 + by * ta2 + cy * ta + dy - ay1 * tb3 - by1 * tb2 - cy1 * tb - dy1
            F1[0][0] = 3 * ax * ta2 + 2 * bx * ta + cx
            F1[0][1] = -3 * ax1 * tb2 - 2 * bx1 * tb - cx1
            F1[1][0] = 3 * ay * ta2 + 2 * by * ta + cy
            F1[1][1] = -3 * ay1 * tb2 - 2 * by1 * tb - cy1
            det = F1[0][0] * F1[1][1] - F1[0][1] * F1[1][0]
            if det != 0:
                F1 = [[F1[1][1] / det, -F1[0][1] / det], [-F1[1][0] / det, F1[0][0] / det]]
                ta = ta - (F1[0][0] * F[0] + F1[0][1] * F[1])
                tb = tb - (F1[1][0] * F[0] + F1[1][1] * F[1])
            else:
                break
            i += 1

        return ta, tb


    def recursion(a, b, ta0, ta1, tb0, tb1, depth_a, depth_b):
        global bezier_intersection_recursive_result
        if a == b:
            bezier_intersection_recursive_result += [[ta0, tb0, ta1, tb1, "Overlap"]]
            return
        tam, tbm = (ta0 + ta1) / 2, (tb0 + tb1) / 2
        if depth_a > 0 and depth_b > 0:
            a1, a2 = bez_split(a, 0.5)
            b1, b2 = bez_split(b, 0.5)
            if bez_bounds_intersect(a1, b1): recursion(a1, b1, ta0, tam, tb0, tbm, depth_a - 1, depth_b - 1)
            if bez_bounds_intersect(a2, b1): recursion(a2, b1, tam, ta1, tb0, tbm, depth_a - 1, depth_b - 1)
            if bez_bounds_intersect(a1, b2): recursion(a1, b2, ta0, tam, tbm, tb1, depth_a - 1, depth_b - 1)
            if bez_bounds_intersect(a2, b2): recursion(a2, b2, tam, ta1, tbm, tb1, depth_a - 1, depth_b - 1)
        elif depth_a > 0:
            a1, a2 = bez_split(a, 0.5)
            if bez_bounds_intersect(a1, b): recursion(a1, b, ta0, tam, tb0, tb1, depth_a - 1, depth_b)
            if bez_bounds_intersect(a2, b): recursion(a2, b, tam, ta1, tb0, tb1, depth_a - 1, depth_b)
        elif depth_b > 0:
            b1, b2 = bez_split(b, 0.5)
            if bez_bounds_intersect(a, b1): recursion(a, b1, ta0, ta1, tb0, tbm, depth_a, depth_b - 1)
            if bez_bounds_intersect(a, b2): recursion(a, b2, ta0, ta1, tbm, tb1, depth_a, depth_b - 1)
        else:  # Both segments have been subdevided enougth. Let's get some intersections :).
            intersection, t1, t2 = straight_segments_intersection([a[0]] + [a[3]], [b[0]] + [b[3]])
            if intersection:
                if intersection == "Overlap":
                    t1 = (max(0, min(1, t1[0])) + max(0, min(1, t1[1]))) / 2
                    t2 = (max(0, min(1, t2[0])) + max(0, min(1, t2[1]))) / 2
                bezier_intersection_recursive_result += [[ta0 + t1 * (ta1 - ta0), tb0 + t2 * (tb1 - tb0)]]

    global bezier_intersection_recursive_result
    bezier_intersection_recursive_result = []
    recursion(a, b, 0., 1., 0., 1., intersection_recursion_depth, intersection_recursion_depth)
    intersections = bezier_intersection_recursive_result
    for i in range(len(intersections)):
        if len(intersections[i]) < 5 or intersections[i][4] != "Overlap":
            intersections[i] = polish_intersection(a, b, intersections[i][0], intersections[i][1])
    return intersections


def csp_segments_true_intersection(sp1, sp2, sp3, sp4):
    intersections = csp_segments_intersection(sp1, sp2, sp3, sp4)
    res = []
    for intersection in intersections:
        if (
                (len(intersection) == 5 and intersection[4] == "Overlap" and (
                        0 <= intersection[0] <= 1 or 0 <= intersection[1] <= 1) and (
                         0 <= intersection[2] <= 1 or 0 <= intersection[3] <= 1))
                or (0 <= intersection[0] <= 1 and 0 <= intersection[1] <= 1)
        ):
            res += [intersection]
    return res


def csp_get_t_at_curvature(sp1, sp2, c, sample_points=16):
    # returns a list containning [t1,t2,t3,...,tn],  0<=ti<=1...
    if sample_points < 2: sample_points = 2
    tolerance = .0000000001
    res = []
    ax, ay, bx, by, cx, cy, dx, dy = csp_parameterize(sp1, sp2)
    for k in range(sample_points):
        t = float(k) / (sample_points - 1)
        i, F = 0, 1e100
        while i < 2 or abs(F) > tolerance and i < 17:
            try:  # some numerical calculation could exceed the limits
                t2 = t * t
                # slopes...
                f1x = 3 * ax * t2 + 2 * bx * t + cx
                f1y = 3 * ay * t2 + 2 * by * t + cy
                f2x = 6 * ax * t + 2 * bx
                f2y = 6 * ay * t + 2 * by
                f3x = 6 * ax
                f3y = 6 * ay
                d = (f1x ** 2 + f1y ** 2) ** 1.5
                F1 = (
                        ((f1x * f3y - f3x * f1y) * d - (f1x * f2y - f2x * f1y) * 3. * (f2x * f1x + f2y * f1y) * (
                                    (f1x ** 2 + f1y ** 2) ** .5)) /
                        ((f1x ** 2 + f1y ** 2) ** 3)
                )
                F = (f1x * f2y - f1y * f2x) / d - c
                t -= F / F1
            except:
                break
            i += 1
        if 0 <= t <= 1 and F <= tolerance:
            if len(res) == 0:
                res.append(t)
            for i in res:
                if abs(t - i) <= 0.001:
                    break
            if not abs(t - i) <= 0.001:
                res.append(t)
    return res


def csp_max_curvature(sp1, sp2):
    ax, ay, bx, by, cx, cy, dx, dy = csp_parameterize(sp1, sp2)
    tolerance = .0001
    F = 0.
    i = 0
    while i < 2 or F - Flast < tolerance and i < 10:
        t = .5
        f1x = 3 * ax * t ** 2 + 2 * bx * t + cx
        f1y = 3 * ay * t ** 2 + 2 * by * t + cy
        f2x = 6 * ax * t + 2 * bx
        f2y = 6 * ay * t + 2 * by
        f3x = 6 * ax
        f3y = 6 * ay
        d = pow(f1x ** 2 + f1y ** 2, 1.5)
        if d != 0:
            Flast = F
            F = (f1x * f2y - f1y * f2x) / d
            F1 = (
                    (d * (f1x * f3y - f3x * f1y) - (f1x * f2y - f2x * f1y) * 3. * (f2x * f1x + f2y * f1y) * pow(
                        f1x ** 2 + f1y ** 2, .5)) /
                    (f1x ** 2 + f1y ** 2) ** 3
            )
            i += 1
            if F1 != 0:
                t -= F / F1
            else:
                break
        else:
            break
    return t


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


def csp_curvature_radius_at_t(sp1, sp2, t):
    c = csp_curvature_at_t(sp1, sp2, t)
    if c == 0:
        return 1e100
    else:
        return 1 / c


def csp_special_points(sp1, sp2):
    # special points = curvature == 0
    ax, ay, bx, by, cx, cy, dx, dy = bezmisc.bezierparameterize((sp1[1], sp1[2], sp2[0], sp2[1]))
    a = 3 * ax * by - 3 * ay * bx
    b = 3 * ax * cy - 3 * cx * ay
    c = bx * cy - cx * by
    roots = cubic_solver(0, a, b, c)
    res = []
    for i in roots:
        if type(i) is complex and i.imag == 0:
            i = i.real
        if type(i) is not complex and 0 <= i <= 1:
            res.append(i)
    return res


def csp_subpath_ccw(subpath):
    # Remove all zerro length segments
    s = 0
    # subpath = subpath[:]
    if (P(subpath[-1][1]) - P(subpath[0][1])).l2() > 1e-10:
        subpath[-1][2] = subpath[-1][1]
        subpath[0][0] = subpath[0][1]
        subpath += [[subpath[0][1], subpath[0][1], subpath[0][1]]]
    pl = subpath[-1][2]
    for sp1 in subpath:
        for p in sp1:
            s += (p[0] - pl[0]) * (p[1] + pl[1])
            pl = p
    return s < 0


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


def csp_splitatlength(sp1, sp2, l=0.5, tolerance=0.01):
    bez = (sp1[1][:], sp1[2][:], sp2[0][:], sp2[1][:])
    t = bezmisc.beziertatlength(bez, l, tolerance)
    return csp_split(sp1, sp2, t)


def cspseglength(sp1, sp2, tolerance=0.001):
    bez = (sp1[1][:], sp1[2][:], sp2[0][:], sp2[1][:])
    return bezmisc.bezierlength(bez, tolerance)


def csplength(csp):
    total = 0
    lengths = []
    for sp in csp:
        for i in xrange(1, len(sp)):
            l = cspseglength(sp[i - 1], sp[i])
            lengths.append(l)
            total += l
    return lengths, total


def csp_segments(csp):
    l, seg = 0, [0]
    for sp in csp:
        for i in xrange(1, len(sp)):
            l += cspseglength(sp[i - 1], sp[i])
            seg += [l]

    if l > 0:
        seg = [seg[i] / l for i in xrange(len(seg))]
    return seg, l


def rebuild_csp(csp, segs, s=None):
    # rebuild_csp() adds to csp control points making it's segments looks like segs
    if s == None: s, l = csp_segments(csp)

    if len(s) > len(segs): return None
    segs = segs[:]
    segs.sort()
    for i in xrange(len(s)):
        d = None
        for j in xrange(len(segs)):
            d = min([abs(s[i] - segs[j]), j], d) if d != None else [abs(s[i] - segs[j]), j]
        del segs[d[1]]
    for i in xrange(len(segs)):
        for j in xrange(0, len(s)):
            if segs[i] < s[j]: break
        if s[j] - s[j - 1] != 0:
            t = (segs[i] - s[j - 1]) / (s[j] - s[j - 1])
            sp1, sp2, sp3 = csp_split(csp[j - 1], csp[j], t)
            csp = csp[:j - 1] + [sp1, sp2, sp3] + csp[j + 1:]
            s = s[:j] + [s[j - 1] * (1 - t) + s[j] * t] + s[j:]
    return csp, s


def csp_slope(sp1, sp2, t):
    bez = (sp1[1][:], sp1[2][:], sp2[0][:], sp2[1][:])
    return bezmisc.bezierslopeatt(bez, t)


def csp_line_intersection(l1, l2, sp1, sp2):
    dd = l1[0]
    cc = l2[0] - l1[0]
    bb = l1[1]
    aa = l2[1] - l1[1]
    if aa == cc == 0: return []
    if aa:
        coef1 = cc / aa
        coef2 = 1
    else:
        coef1 = 1
        coef2 = aa / cc
    bez = (sp1[1][:], sp1[2][:], sp2[0][:], sp2[1][:])
    ax, ay, bx, by, cx, cy, x0, y0 = bezmisc.bezierparameterize(bez)
    a = coef1 * ay - coef2 * ax
    b = coef1 * by - coef2 * bx
    c = coef1 * cy - coef2 * cx
    d = coef1 * (y0 - bb) - coef2 * (x0 - dd)
    roots = cubic_solver(a, b, c, d)
    retval = []
    for i in roots:
        if type(i) is complex and abs(i.imag) < 1e-7:
            i = i.real
        if type(i) is not complex and -1e-10 <= i <= 1. + 1e-10:
            retval.append(i)
    return retval


def csp_from_arc(start, end, center, r, slope_st):
    # Creates csp that approximise specified arc
    r = abs(r)
    alpha = (atan2(end[0] - center[0], end[1] - center[1]) - atan2(start[0] - center[0],
                                                                   start[1] - center[1])) % math.pi2

    sectors = int(abs(alpha) * 2 / math.pi) + 1
    alpha_start = atan2(start[0] - center[0], start[1] - center[1])
    cos_, sin_ = math.cos(alpha_start), math.sin(alpha_start)
    k = (4. * math.tan(alpha / sectors / 4.) / 3.)
    if dot(slope_st, [- sin_ * k * r, cos_ * k * r]) < 0:
        if alpha > 0:
            alpha -= math.pi2
        else:
            alpha += math.pi2
    if abs(alpha * r) < 0.001:
        return []

    sectors = int(abs(alpha) * 2 / math.pi) + 1
    k = (4. * math.tan(alpha / sectors / 4.) / 3.)
    result = []
    for i in range(sectors + 1):
        cos_, sin_ = math.cos(alpha_start + alpha * i / sectors), math.sin(alpha_start + alpha * i / sectors)
        sp = [[], [center[0] + cos_ * r, center[1] + sin_ * r], []]
        sp[0] = [sp[1][0] + sin_ * k * r, sp[1][1] - cos_ * k * r]
        sp[2] = [sp[1][0] - sin_ * k * r, sp[1][1] + cos_ * k * r]
        result += [sp]
    result[0][0] = result[0][1][:]
    result[-1][2] = result[-1][1]

    return result


def point_to_arc_distance(p, arc):
    ###        Distance calculattion from point to arc
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


def csp_simple_bound_to_point_distance(p, csp):
    minx, miny, maxx, maxy = None, None, None, None
    for subpath in csp:
        for sp in subpath:
            for p_ in sp:
                minx = min(minx, p_[0]) if minx != None else p_[0]
                miny = min(miny, p_[1]) if miny != None else p_[1]
                maxx = max(maxx, p_[0]) if maxx != None else p_[0]
                maxy = max(maxy, p_[1]) if maxy != None else p_[1]
    return math.sqrt(max(minx - p[0], p[0] - maxx, 0) ** 2 + max(miny - p[1], p[1] - maxy, 0) ** 2)


def csp_point_inside_bound(sp1, sp2, p):
    bez = [sp1[1], sp1[2], sp2[0], sp2[1]]
    x, y = p
    c = 0
    for i in range(4):
        [x0, y0], [x1, y1] = bez[i - 1], bez[i]
        if x0 - x1 != 0 and (y - y0) * (x1 - x0) >= (x - x0) * (y1 - y0) and x > min(x0, x1) and x <= max(x0, x1):
            c += 1
    return c % 2 == 0


def csp_bound_to_point_distance(sp1, sp2, p):
    if csp_point_inside_bound(sp1, sp2, p):
        return 0.
    bez = csp_segment_to_bez(sp1, sp2)
    min_dist = 1e100
    for i in range(0, 4):
        d = point_to_line_segment_distance_2(p, bez[i - 1], bez[i])
        if d <= min_dist: min_dist = d
    return min_dist


def line_line_intersect(p1, p2, p3, p4):  # Return only true intersection.
    if (p1[0] == p2[0] and p1[1] == p2[1]) or (p3[0] == p4[0] and p3[1] == p4[1]): return False
    x = (p2[0] - p1[0]) * (p4[1] - p3[1]) - (p2[1] - p1[1]) * (p4[0] - p3[0])
    if x == 0:  # Lines are parallel
        if (p3[0] - p1[0]) * (p2[1] - p1[1]) == (p3[1] - p1[1]) * (p2[0] - p1[0]):
            if p3[0] != p4[0]:
                t11 = (p1[0] - p3[0]) / (p4[0] - p3[0])
                t12 = (p2[0] - p3[0]) / (p4[0] - p3[0])
                t21 = (p3[0] - p1[0]) / (p2[0] - p1[0])
                t22 = (p4[0] - p1[0]) / (p2[0] - p1[0])
            else:
                t11 = (p1[1] - p3[1]) / (p4[1] - p3[1])
                t12 = (p2[1] - p3[1]) / (p4[1] - p3[1])
                t21 = (p3[1] - p1[1]) / (p2[1] - p1[1])
                t22 = (p4[1] - p1[1]) / (p2[1] - p1[1])
            return ("Overlap" if (0 <= t11 <= 1 or 0 <= t12 <= 1) and (0 <= t21 <= 1 or 0 <= t22 <= 1) else False)
        else:
            return False
    else:
        return (
                0 <= ((p4[0] - p3[0]) * (p1[1] - p3[1]) - (p4[1] - p3[1]) * (p1[0] - p3[0])) / x <= 1 and
                0 <= ((p2[0] - p1[0]) * (p1[1] - p3[1]) - (p2[1] - p1[1]) * (p1[0] - p3[0])) / x <= 1)


def line_line_intersection_points(p1, p2, p3, p4):  # Return only points [ (x,y) ]
    if (p1[0] == p2[0] and p1[1] == p2[1]) or (p3[0] == p4[0] and p3[1] == p4[1]): return []
    x = (p2[0] - p1[0]) * (p4[1] - p3[1]) - (p2[1] - p1[1]) * (p4[0] - p3[0])
    if x == 0:  # Lines are parallel
        if (p3[0] - p1[0]) * (p2[1] - p1[1]) == (p3[1] - p1[1]) * (p2[0] - p1[0]):
            if p3[0] != p4[0]:
                t11 = (p1[0] - p3[0]) / (p4[0] - p3[0])
                t12 = (p2[0] - p3[0]) / (p4[0] - p3[0])
                t21 = (p3[0] - p1[0]) / (p2[0] - p1[0])
                t22 = (p4[0] - p1[0]) / (p2[0] - p1[0])
            else:
                t11 = (p1[1] - p3[1]) / (p4[1] - p3[1])
                t12 = (p2[1] - p3[1]) / (p4[1] - p3[1])
                t21 = (p3[1] - p1[1]) / (p2[1] - p1[1])
                t22 = (p4[1] - p1[1]) / (p2[1] - p1[1])
            res = []
            if (0 <= t11 <= 1 or 0 <= t12 <= 1) and (0 <= t21 <= 1 or 0 <= t22 <= 1):
                if 0 <= t11 <= 1: res += [p1]
                if 0 <= t12 <= 1: res += [p2]
                if 0 <= t21 <= 1: res += [p3]
                if 0 <= t22 <= 1: res += [p4]
            return res
        else:
            return []
    else:
        t1 = ((p4[0] - p3[0]) * (p1[1] - p3[1]) - (p4[1] - p3[1]) * (p1[0] - p3[0])) / x
        t2 = ((p2[0] - p1[0]) * (p1[1] - p3[1]) - (p2[1] - p1[1]) * (p1[0] - p3[0])) / x
        if 0 <= t1 <= 1 and 0 <= t2 <= 1:
            return [[p1[0] * (1 - t1) + p2[0] * t1, p1[1] * (1 - t1) + p2[1] * t1]]
        else:
            return []


def point_to_point_d2(a, b):
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2


def point_to_point_d(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)


def point_to_line_segment_distance_2(p1, p2, p3):
    # p1 - point, p2,p3 - line segment
    # draw_pointer(p1)
    w0 = [p1[0] - p2[0], p1[1] - p2[1]]
    v = [p3[0] - p2[0], p3[1] - p2[1]]
    c1 = w0[0] * v[0] + w0[1] * v[1]
    if c1 <= 0:
        return w0[0] * w0[0] + w0[1] * w0[1]
    c2 = v[0] * v[0] + v[1] * v[1]
    if c2 <= c1:
        return (p1[0] - p3[0]) ** 2 + (p1[1] - p3[1]) ** 2
    return (p1[0] - p2[0] - v[0] * c1 / c2) ** 2 + (p1[1] - p2[1] - v[1] * c1 / c2)


def line_to_line_distance_2(p1, p2, p3, p4):
    if line_line_intersect(p1, p2, p3, p4): return 0
    return min(
        point_to_line_segment_distance_2(p1, p3, p4),
        point_to_line_segment_distance_2(p2, p3, p4),
        point_to_line_segment_distance_2(p3, p1, p2),
        point_to_line_segment_distance_2(p4, p1, p2))


def csp_seg_bound_to_csp_seg_bound_max_min_distance(sp1, sp2, sp3, sp4):
    bez1 = csp_segment_to_bez(sp1, sp2)
    bez2 = csp_segment_to_bez(sp3, sp4)
    min_dist = 1e100
    max_dist = 0.
    for i in range(4):
        if csp_point_inside_bound(sp1, sp2, bez2[i]) or csp_point_inside_bound(sp3, sp4, bez1[i]):
            min_dist = 0.
            break
    for i in range(4):
        for j in range(4):
            d = line_to_line_distance_2(bez1[i - 1], bez1[i], bez2[j - 1], bez2[j])
            if d < min_dist: min_dist = d
            d = (bez2[j][0] - bez1[i][0]) ** 2 + (bez2[j][1] - bez1[i][1]) ** 2
            if max_dist < d: max_dist = d
    return min_dist, max_dist


def csp_reverse(csp):
    for i in range(len(csp)):
        n = []
        for j in csp[i]:
            n = [[j[2][:], j[1][:], j[0][:]]] + n
        csp[i] = n[:]
    return csp


def csp_normalized_slope(sp1, sp2, t):
    ax, ay, bx, by, cx, cy, dx, dy = bezmisc.bezierparameterize((sp1[1][:], sp1[2][:], sp2[0][:], sp2[1][:]))
    if sp1[1] == sp2[1] == sp1[2] == sp2[0]: return [1., 0.]
    f1x = 3 * ax * t * t + 2 * bx * t + cx
    f1y = 3 * ay * t * t + 2 * by * t + cy
    if abs(f1x * f1x + f1y * f1y) > 1e-20:
        l = math.sqrt(f1x * f1x + f1y * f1y)
        return [f1x / l, f1y / l]

    if t == 0:
        f1x = sp2[0][0] - sp1[1][0]
        f1y = sp2[0][1] - sp1[1][1]
        if abs(f1x * f1x + f1y * f1y) > 1e-20:
            l = math.sqrt(f1x * f1x + f1y * f1y)
            return [f1x / l, f1y / l]
        else:
            f1x = sp2[1][0] - sp1[1][0]
            f1y = sp2[1][1] - sp1[1][1]
            if f1x * f1x + f1y * f1y != 0:
                l = math.sqrt(f1x * f1x + f1y * f1y)
                return [f1x / l, f1y / l]
    elif t == 1:
        f1x = sp2[1][0] - sp1[2][0]
        f1y = sp2[1][1] - sp1[2][1]
        if abs(f1x * f1x + f1y * f1y) > 1e-20:
            l = math.sqrt(f1x * f1x + f1y * f1y)
            return [f1x / l, f1y / l]
        else:
            f1x = sp2[1][0] - sp1[1][0]
            f1y = sp2[1][1] - sp1[1][1]
            if f1x * f1x + f1y * f1y != 0:
                l = math.sqrt(f1x * f1x + f1y * f1y)
                return [f1x / l, f1y / l]
    else:
        return [1., 0.]


def csp_normalized_normal(sp1, sp2, t):
    nx, ny = csp_normalized_slope(sp1, sp2, t)
    return [-ny, nx]


def csp_parameterize(sp1, sp2):
    return bezmisc.bezierparameterize(csp_segment_to_bez(sp1, sp2))


def csp_concat_subpaths(*s):
    def concat(s1, s2):
        if s1 == []: return s2
        if s2 == []: return s1
        if (s1[-1][1][0] - s2[0][1][0]) ** 2 + (s1[-1][1][1] - s2[0][1][1]) ** 2 > 0.00001:
            return s1[:-1] + [[s1[-1][0], s1[-1][1], s1[-1][1]], [s2[0][1], s2[0][1], s2[0][2]]] + s2[1:]
        else:
            return s1[:-1] + [[s1[-1][0], s2[0][1], s2[0][2]]] + s2[1:]

    if len(s) == 0: return []
    if len(s) == 1: return s[0]
    result = s[0]
    for s1 in s[1:]:
        result = concat(result, s1)
    return result


def csp_draw(csp, color="#05f", group=None, style="fill:none;", width=.1, comment=""):
    if csp != [] and csp != [[]]:
        if group == None: group = options.doc_root
        style += "stroke:" + color + ";" + "stroke-width:%0.4fpx;" % width
        args = {"d": cubicsuperpath.formatPath(csp), "style": style}
        if comment != "": args["comment"] = str(comment)
        inkex.etree.SubElement(group, inkex.addNS('path', 'svg'), args)


def csp_subpaths_end_to_start_distance2(s1, s2):
    return (s1[-1][1][0] - s2[0][1][0]) ** 2 + (s1[-1][1][1] - s2[0][1][1]) ** 2


def csp_subpath_line_to(subpath, points):
    # Appends subpath with line or polyline.
    if len(points) > 0:
        if len(subpath) > 0:
            subpath[-1][2] = subpath[-1][1][:]
        if type(points[0]) == type([1, 1]):
            for p in points:
                subpath += [[p[:], p[:], p[:]]]
        else:
            subpath += [[points, points, points]]
    return subpath


def csp_join_subpaths(csp):
    result = csp[:]
    done_smf = True
    joined_result = []
    while done_smf:
        done_smf = False
        while len(result) > 0:
            s1 = result[-1][:]
            del (result[-1])
            j = 0
            joined_smf = False
            while j < len(joined_result):
                if csp_subpaths_end_to_start_distance2(joined_result[j], s1) < 0.000001:
                    joined_result[j] = csp_concat_subpaths(joined_result[j], s1)
                    done_smf = True
                    joined_smf = True
                    break
                if csp_subpaths_end_to_start_distance2(s1, joined_result[j]) < 0.000001:
                    joined_result[j] = csp_concat_subpaths(s1, joined_result[j])
                    done_smf = True
                    joined_smf = True
                    break
                j += 1
            if not joined_smf: joined_result += [s1[:]]
        if done_smf:
            result = joined_result[:]
            joined_result = []
    return joined_result


def triangle_cross(a, b, c):
    return (a[0] - b[0]) * (c[1] - b[1]) - (c[0] - b[0]) * (a[1] - b[1])


def csp_segment_convex_hull(sp1, sp2):
    a, b, c, d = sp1[1][:], sp1[2][:], sp2[0][:], sp2[1][:]

    abc = triangle_cross(a, b, c)
    abd = triangle_cross(a, b, d)
    bcd = triangle_cross(b, c, d)
    cad = triangle_cross(c, a, d)
    if abc == 0 and abd == 0: return [min(a, b, c, d), max(a, b, c, d)]
    if abc == 0: return [d, min(a, b, c), max(a, b, c)]
    if abd == 0: return [c, min(a, b, d), max(a, b, d)]
    if bcd == 0: return [a, min(b, c, d), max(b, c, d)]
    if cad == 0: return [b, min(c, a, d), max(c, a, d)]

    m1, m2, m3 = abc * abd > 0, abc * bcd > 0, abc * cad > 0
    if m1 and m2 and m3: return [a, b, c]
    if m1 and m2 and not m3: return [a, b, c, d]
    if m1 and not m2 and m3: return [a, b, d, c]
    if not m1 and m2 and m3: return [a, d, b, c]
    if m1 and not (m2 and m3): return [a, b, d]
    if not (m1 and m2) and m3: return [c, a, d]
    if not (m1 and m3) and m2: return [b, c, d]

    raise ValueError, "csp_segment_convex_hull happend something that shouldnot happen!"


################################################################################
###        Bezier additional functions
################################################################################

def bez_bounds_intersect(bez1, bez2):
    return bounds_intersect(bez_bound(bez2), bez_bound(bez1))


def bez_bound(bez):
    return [
        min(bez[0][0], bez[1][0], bez[2][0], bez[3][0]),
        min(bez[0][1], bez[1][1], bez[2][1], bez[3][1]),
        max(bez[0][0], bez[1][0], bez[2][0], bez[3][0]),
        max(bez[0][1], bez[1][1], bez[2][1], bez[3][1]),
    ]


def bounds_intersect(a, b):
    return not ((a[0] > b[2]) or (b[0] > a[2]) or (a[1] > b[3]) or (b[1] > a[3]))


def tpoint((x1, y1), (x2, y2), t):
    return [x1 + t * (x2 - x1), y1 + t * (y2 - y1)]


def bez_to_csp_segment(bez):
    return [bez[0], bez[0], bez[1]], [bez[2], bez[3], bez[3]]


def bez_split(a, t=0.5):
    a1 = tpoint(a[0], a[1], t)
    at = tpoint(a[1], a[2], t)
    b2 = tpoint(a[2], a[3], t)
    a2 = tpoint(a1, at, t)
    b1 = tpoint(b2, at, t)
    a3 = tpoint(a2, b1, t)
    return [a[0], a1, a2, a3], [a3, b1, b2, a[3]]


def bez_at_t(bez, t):
    return csp_at_t([bez[0], bez[0], bez[1]], [bez[2], bez[3], bez[3]], t)


def bez_to_point_distance(bez, p, needed_dist=[0., 1e100]):
    # returns [d^2,t]
    return csp_seg_to_point_distance(bez_to_csp_segment(bez), p, needed_dist)


def bez_normalized_slope(bez, t):
    return csp_normalized_slope([bez[0], bez[0], bez[1]], [bez[2], bez[3], bez[3]], t)


################################################################################
###    Some vector functions
################################################################################

def normalize((x, y)):
    l = math.sqrt(x ** 2 + y ** 2)
    if l == 0:
        return [0., 0.]
    else:
        return [x / l, y / l]


def cross(a, b):
    return a[1] * b[0] - a[0] * b[1]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def rotate_ccw(d):
    return [-d[1], d[0]]


def vectors_ccw(a, b):
    return a[0] * b[1] - b[0] * a[1] < 0


def vector_from_to_length(a, b):
    return math.sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]))


################################################################################
###    Common functions
################################################################################

def matrix_mul(a, b):
    return [[sum([a[i][k] * b[k][j] for k in range(len(a[0]))]) for j in range(len(b[0]))] for i in range(len(a))]
    try:
        return [[sum([a[i][k] * b[k][j] for k in range(len(a[0]))]) for j in range(len(b[0]))] for i in range(len(a))]
    except:
        return None


def transpose(a):
    try:
        return [[a[i][j] for i in range(len(a))] for j in range(len(a[0]))]
    except:
        return None


def det_3x3(a):
    return float(
        a[0][0] * a[1][1] * a[2][2] + a[0][1] * a[1][2] * a[2][0] + a[1][0] * a[2][1] * a[0][2]
        - a[0][2] * a[1][1] * a[2][0] - a[0][0] * a[2][1] * a[1][2] - a[0][1] * a[2][2] * a[1][0]
    )


def inv_3x3(a):  # invert matrix 3x3
    det = det_3x3(a)
    if det == 0: return None
    return [
        [(a[1][1] * a[2][2] - a[2][1] * a[1][2]) / det, -(a[0][1] * a[2][2] - a[2][1] * a[0][2]) / det,
         (a[0][1] * a[1][2] - a[1][1] * a[0][2]) / det],
        [-(a[1][0] * a[2][2] - a[2][0] * a[1][2]) / det, (a[0][0] * a[2][2] - a[2][0] * a[0][2]) / det,
         -(a[0][0] * a[1][2] - a[1][0] * a[0][2]) / det],
        [(a[1][0] * a[2][1] - a[2][0] * a[1][1]) / det, -(a[0][0] * a[2][1] - a[2][0] * a[0][1]) / det,
         (a[0][0] * a[1][1] - a[1][0] * a[0][1]) / det]
    ]


def inv_2x2(a):  # invert matrix 2x2
    det = a[0][0] * a[1][1] - a[1][0] * a[0][1]
    if det == 0: return None
    return [
        [a[1][1] / det, -a[0][1] / det],
        [-a[1][0] / det, a[0][0] / det]
    ]


def small(a):
    global small_tolerance
    return abs(a) < small_tolerance


def atan2(*arg):
    if len(arg) == 1 and (type(arg[0]) == type([0., 0.]) or type(arg[0]) == type((0., 0.))):
        return (math.pi / 2 - math.atan2(arg[0][0], arg[0][1])) % math.pi2
    elif len(arg) == 2:

        return (math.pi / 2 - math.atan2(arg[0], arg[1])) % math.pi2
    else:
        raise ValueError, "Bad argumets for atan! (%s)" % arg


def draw_text(text, x, y, style=None, font_size=20):
    if style == None:
        style = "font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;fill:#000000;fill-opacity:1;stroke:none;"
    style += "font-size:%fpx;" % font_size
    t = inkex.etree.SubElement(options.doc_root, inkex.addNS('text', 'svg'), {
        'x': str(x),
        inkex.addNS("space", "xml"): "preserve",
        'y': str(y)
    })
    text = str(text).split("\n")
    for s in text:
        span = inkex.etree.SubElement(t, inkex.addNS('tspan', 'svg'),
                                      {
                                          'x': str(x),
                                          'y': str(+y),
                                          inkex.addNS("role", "sodipodi"): "line",
                                      })
        y += font_size
        span.text = s


def draw_pointer(x, color="#f00", figure="cross", comment="", width=.1):
    if figure == "line":
        s = ""
        for i in range(1, len(x) / 2):
            s += " %s, %s " % (x[i * 2], x[i * 2 + 1])
        inkex.etree.SubElement(options.doc_root, inkex.addNS('path', 'svg'), {"d": "M %s,%s L %s" % (x[0], x[1], s),
                                                                              "style": "fill:none;stroke:%s;stroke-width:%f;" % (
                                                                              color, width), "comment": str(comment)})
    else:
        inkex.etree.SubElement(options.doc_root, inkex.addNS('path', 'svg'),
                               {"d": "m %s,%s l 10,10 -20,-20 10,10 -10,10, 20,-20" % (x[0], x[1]),
                                "style": "fill:none;stroke:%s;stroke-width:%f;" % (color, width),
                                "comment": str(comment)})


def straight_segments_intersection(a, b,
                                   true_intersection=True):  # (True intersection means check ta and tb are in [0,1])
    ax, bx, cx, dx, ay, by, cy, dy = a[0][0], a[1][0], b[0][0], b[1][0], a[0][1], a[1][1], b[0][1], b[1][1]
    if (ax == bx and ay == by) or (cx == dx and cy == dy): return False, 0, 0
    if (bx - ax) * (dy - cy) - (by - ay) * (dx - cx) == 0:  # Lines are parallel
        ta = (ax - cx) / (dx - cx) if cx != dx else (ay - cy) / (dy - cy)
        tb = (bx - cx) / (dx - cx) if cx != dx else (by - cy) / (dy - cy)
        tc = (cx - ax) / (bx - ax) if ax != bx else (cy - ay) / (by - ay)
        td = (dx - ax) / (bx - ax) if ax != bx else (dy - ay) / (by - ay)
        return (
                   "Overlap" if 0 <= ta <= 1 or 0 <= tb <= 1 or 0 <= tc <= 1 or 0 <= td <= 1 or not true_intersection else False), (
               ta, tb), (tc, td)
    else:
        ta = ((ay - cy) * (dx - cx) - (ax - cx) * (dy - cy)) / ((bx - ax) * (dy - cy) - (by - ay) * (dx - cx))
        tb = (ax - cx + ta * (bx - ax)) / (dx - cx) if dx != cx else (ay - cy + ta * (by - ay)) / (dy - cy)
        return (0 <= ta <= 1 and 0 <= tb <= 1 or not true_intersection), ta, tb


def isnan(x): return type(x) is float and x != x


def isinf(x): inf = 1e5000; return x == inf or x == -inf


def between(c, x, y):
    return x - straight_tolerance <= c <= y + straight_tolerance or y - straight_tolerance <= c <= x + straight_tolerance


def cubic_solver(a, b, c, d):
    if a != 0:
        #    Monics formula see http://en.wikipedia.org/wiki/Cubic_function#Monic_formula_of_roots
        a, b, c = (b / a, c / a, d / a)
        m = 2 * a ** 3 - 9 * a * b + 27 * c
        k = a ** 2 - 3 * b
        n = m ** 2 - 4 * k ** 3
        w1 = -.5 + .5 * cmath.sqrt(3) * 1j
        w2 = -.5 - .5 * cmath.sqrt(3) * 1j
        if n >= 0:
            t = m + math.sqrt(n)
            m1 = pow(t / 2, 1. / 3) if t >= 0 else -pow(-t / 2, 1. / 3)
            t = m - math.sqrt(n)
            n1 = pow(t / 2, 1. / 3) if t >= 0 else -pow(-t / 2, 1. / 3)
        else:
            m1 = pow(complex((m + cmath.sqrt(n)) / 2), 1. / 3)
            n1 = pow(complex((m - cmath.sqrt(n)) / 2), 1. / 3)
        x1 = -1. / 3 * (a + m1 + n1)
        x2 = -1. / 3 * (a + w1 * m1 + w2 * n1)
        x3 = -1. / 3 * (a + w2 * m1 + w1 * n1)
        return [x1, x2, x3]
    elif b != 0:
        det = c ** 2 - 4 * b * d
        if det > 0:
            return [(-c + math.sqrt(det)) / (2 * b), (-c - math.sqrt(det)) / (2 * b)]
        elif d == 0:
            return [-c / (b * b)]
        else:
            return [(-c + cmath.sqrt(det)) / (2 * b), (-c - cmath.sqrt(det)) / (2 * b)]
    elif c != 0:
        return [-d / c]
    else:
        return []


################################################################################
###        print_ prints any arguments into specified log file
################################################################################

def print_(*arg):
    f = open(options.log_filename, "a")
    for s in arg:
        s = str(unicode(s).encode('unicode_escape')) + " "
        f.write(s)
    f.write("\n")
    f.close()


################################################################################
###        Point (x,y) operations
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

    def dot(self, other):
        return self.x * other.x + self.y * other.y

    def rot(self, theta):
        c = math.cos(theta)
        s = math.sin(theta)
        return P(self.x * c - self.y * s, self.x * s + self.y * c)

    def angle(self):
        return math.atan2(self.y, self.x)

    def __repr__(self):
        return '%f,%f' % (self.x, self.y)

    def pr(self):
        return "%.2f,%.2f" % (self.x, self.y)

    def to_list(self):
        return [self.x, self.y]

    def ccw(self):
        return P(-self.y, self.x)

    def l2(self):
        return self.x * self.x + self.y * self.y


################################################################################
###
### Offset function
###
### This function offsets given cubic super path.
### It's based on src/livarot/PathOutline.cpp from Inkscape's source code.
###
###
################################################################################
def csp_offset(csp, r):
    offset_tolerance = 0.05
    offset_subdivision_depth = 10
    time_ = time.time()
    time_start = time_
    print_("Offset start at %s" % time_)
    print_("Offset radius %s" % r)


    def create_offset_segment(sp1, sp2, r):
        # See    Gernot Hoffmann "Bezier Curves"  p.34 -> 7.1 Bezier Offset Curves
        p0, p1, p2, p3 = P(sp1[1]), P(sp1[2]), P(sp2[0]), P(sp2[1])
        s0, s1, s3 = p1 - p0, p2 - p1, p3 - p2
        n0 = s0.ccw().unit() if s0.l2() != 0 else P(csp_normalized_normal(sp1, sp2, 0))
        n3 = s3.ccw().unit() if s3.l2() != 0 else P(csp_normalized_normal(sp1, sp2, 1))
        n1 = s1.ccw().unit() if s1.l2() != 0 else (n0.unit() + n3.unit()).unit()

        q0, q3 = p0 + r * n0, p3 + r * n3
        c = csp_curvature_at_t(sp1, sp2, 0)
        q1 = q0 + (p1 - p0) * (1 - (r * c if abs(c) < 100 else 0))
        c = csp_curvature_at_t(sp1, sp2, 1)
        q2 = q3 + (p2 - p3) * (1 - (r * c if abs(c) < 100 else 0))

        return [[q0.to_list(), q0.to_list(), q1.to_list()], [q2.to_list(), q3.to_list(), q3.to_list()]]


    def csp_get_subapths_last_first_intersection(s1, s2):
        _break = False
        for i in range(1, len(s1)):
            sp11, sp12 = s1[-i - 1], s1[-i]
            for j in range(1, len(s2)):
                sp21, sp22 = s2[j - 1], s2[j]
                intersection = csp_segments_true_intersection(sp11, sp12, sp21, sp22)
                if intersection != []:
                    _break = True
                    break
            if _break: break
        if _break:
            intersection = max(intersection)
            return [len(s1) - i, intersection[0], j, intersection[1]]
        else:
            return []


    def csp_join_offsets(prev, next, sp1, sp2, sp1_l, sp2_l, r):
        if len(next) > 1:
            if (P(prev[-1][1]) - P(next[0][1])).l2() < 0.001:
                return prev, [], next
            intersection = csp_get_subapths_last_first_intersection(prev, next)
            if intersection != []:
                i, t1, j, t2 = intersection
                sp1_, sp2_, sp3_ = csp_split(prev[i - 1], prev[i], t1)
                sp3_, sp4_, sp5_ = csp_split(next[j - 1], next[j], t2)
                return prev[:i - 1] + [sp1_, sp2_], [], [sp4_, sp5_] + next[j + 1:]

        # Offsets do not intersect... will add an arc...
        start = (P(csp_at_t(sp1_l, sp2_l, 1.)) + r * P(csp_normalized_normal(sp1_l, sp2_l, 1.))).to_list()
        end = (P(csp_at_t(sp1, sp2, 0.)) + r * P(csp_normalized_normal(sp1, sp2, 0.))).to_list()
        arc = csp_from_arc(start, end, sp1[1], r, csp_normalized_slope(sp1_l, sp2_l, 1.))
        if arc == []:
            return prev, [], next
        else:
            # Clip prev by arc
            if csp_subpaths_end_to_start_distance2(prev, arc) > 0.00001:
                intersection = csp_get_subapths_last_first_intersection(prev, arc)
                if intersection != []:
                    i, t1, j, t2 = intersection
                    sp1_, sp2_, sp3_ = csp_split(prev[i - 1], prev[i], t1)
                    sp3_, sp4_, sp5_ = csp_split(arc[j - 1], arc[j], t2)
                    prev = prev[:i - 1] + [sp1_, sp2_]
                    arc = [sp4_, sp5_] + arc[j + 1:]
                # else : raise ValueError, "Offset curvature clipping error"
            # Clip next by arc
            if next == []:
                return prev, [], arc
            if csp_subpaths_end_to_start_distance2(arc, next) > 0.00001:
                intersection = csp_get_subapths_last_first_intersection(arc, next)
                if intersection != []:
                    i, t1, j, t2 = intersection
                    sp1_, sp2_, sp3_ = csp_split(arc[i - 1], arc[i], t1)
                    sp3_, sp4_, sp5_ = csp_split(next[j - 1], next[j], t2)
                    arc = arc[:i - 1] + [sp1_, sp2_]
                    next = [sp4_, sp5_] + next[j + 1:]
                # else : raise ValueError, "Offset curvature clipping error"

            return prev, arc, next


    def offset_segment_recursion(sp1, sp2, r, depth, tolerance):
        sp1_r, sp2_r = create_offset_segment(sp1, sp2, r)
        err = max(
            csp_seg_to_point_distance(sp1_r, sp2_r, (
                        P(csp_at_t(sp1, sp2, .25)) + P(csp_normalized_normal(sp1, sp2, .25)) * r).to_list())[0],
            csp_seg_to_point_distance(sp1_r, sp2_r, (
                        P(csp_at_t(sp1, sp2, .50)) + P(csp_normalized_normal(sp1, sp2, .50)) * r).to_list())[0],
            csp_seg_to_point_distance(sp1_r, sp2_r, (
                        P(csp_at_t(sp1, sp2, .75)) + P(csp_normalized_normal(sp1, sp2, .75)) * r).to_list())[0],
        )

        if err > tolerance ** 2 and depth > 0:
            # print_(csp_seg_to_point_distance(sp1_r,sp2_r, (P(csp_at_t(sp1,sp2,.25)) + P(csp_normalized_normal(sp1,sp2,.25))*r).to_list())[0], tolerance)
            if depth > offset_subdivision_depth - 2:
                t = csp_max_curvature(sp1, sp2)
                t = max(.1, min(.9, t))
            else:
                t = .5
            sp3, sp4, sp5 = csp_split(sp1, sp2, t)
            r1 = offset_segment_recursion(sp3, sp4, r, depth - 1, tolerance)
            r2 = offset_segment_recursion(sp4, sp5, r, depth - 1, tolerance)
            return r1[:-1] + [[r1[-1][0], r1[-1][1], r2[0][2]]] + r2[1:]
        else:
            # csp_draw([[sp1_r,sp2_r]])
            # draw_pointer(sp1[1]+sp1_r[1], "#057", "line")
            # draw_pointer(sp2[1]+sp2_r[1], "#705", "line")
            return [sp1_r, sp2_r]


    ############################################################################
    # Some small definitions
    ############################################################################
    csp_len = len(csp)

    ############################################################################
    # Prepare the path
    ############################################################################
    # Remove all small segments (segment length < 0.001)

    for i in xrange(len(csp)):
        for j in xrange(len(csp[i])):
            sp = csp[i][j]
            if (P(sp[1]) - P(sp[0])).mag() < 0.001:
                csp[i][j][0] = sp[1]
            if (P(sp[2]) - P(sp[0])).mag() < 0.001:
                csp[i][j][2] = sp[1]
    for i in xrange(len(csp)):
        for j in xrange(1, len(csp[i])):
            if cspseglength(csp[i][j - 1], csp[i][j]) < 0.001:
                csp[i] = csp[i][:j] + csp[i][j + 1:]
        if cspseglength(csp[i][-1], csp[i][0]) > 0.001:
            csp[i][-1][2] = csp[i][-1][1]
            csp[i] += [[csp[i][0][1], csp[i][0][1], csp[i][0][1]]]

    # TODO Get rid of self intersections.

    original_csp = csp[:]
    # Clip segments which has curvature>1/r. Because their offset will be selfintersecting and very nasty.

    print_("Offset prepared the path in %s" % (time.time() - time_))
    print_("Path length = %s" % sum([len(i) for i in csp]))
    time_ = time.time()

    ############################################################################
    # Offset
    ############################################################################
    # Create offsets for all segments in the path. And join them together inside each subpath.
    unclipped_offset = [[] for i in xrange(csp_len)]
    offsets_original = [[] for i in xrange(csp_len)]
    join_points = [[] for i in xrange(csp_len)]
    intersection = [[] for i in xrange(csp_len)]
    for i in xrange(csp_len):
        subpath = csp[i]
        subpath_offset = []
        last_offset_len = 0
        for sp1, sp2 in zip(subpath, subpath[1:]):
            segment_offset = csp_offset_segment(sp1, sp2, r)
            if subpath_offset == []:
                subpath_offset = segment_offset

                prev_l = len(subpath_offset)
            else:
                prev, arc, next = csp_join_offsets(subpath_offset[-prev_l:], segment_offset, sp1, sp2, sp1_l, sp2_l, r)
                # csp_draw([prev],"Blue")
                # csp_draw([arc],"Magenta")
                subpath_offset = csp_concat_subpaths(subpath_offset[:-prev_l + 1], prev, arc, next)
                prev_l = len(next)
            sp1_l, sp2_l = sp1[:], sp2[:]

        # Join last and first offsets togother to close the curve

        prev, arc, next = csp_join_offsets(subpath_offset[-prev_l:], subpath_offset[:2], subpath[0], subpath[1], sp1_l,
                                           sp2_l, r)
        subpath_offset[:2] = next[:]
        subpath_offset = csp_concat_subpaths(subpath_offset[:-prev_l + 1], prev, arc)
        # csp_draw([prev],"Blue")
        # csp_draw([arc],"Red")
        # csp_draw([next],"Red")

        # Collect subpath's offset and save it to unclipped offset list.
        unclipped_offset[i] = subpath_offset[:]

        # for k,t in intersection[i]:
        #    draw_pointer(csp_at_t(subpath_offset[k-1], subpath_offset[k], t))

    # inkex.etree.SubElement( options.doc_root, inkex.addNS('path','svg'), {"d": cubicsuperpath.formatPath(unclipped_offset), "style":"fill:none;stroke:#0f0;"} )
    print_("Offsetted path in %s" % (time.time() - time_))
    time_ = time.time()

    # for i in range(len(unclipped_offset)):
    #    csp_draw([unclipped_offset[i]], color = ["Green","Red","Blue"][i%3], width = .1)
    # return []
    ############################################################################
    # Now to the clipping.
    ############################################################################
    # First of all find all intersection's between all segments of all offseted subpaths, including self intersections.

    # TODO define offset tolerance here
    global small_tolerance
    small_tolerance = 0.01
    summ = 0
    summ1 = 0
    for subpath_i in xrange(csp_len):
        for subpath_j in xrange(subpath_i, csp_len):
            subpath = unclipped_offset[subpath_i]
            subpath1 = unclipped_offset[subpath_j]
            for i in xrange(1, len(subpath)):
                # If subpath_i==subpath_j we are looking for self intersections, so
                # we'll need search intersections only for xrange(i,len(subpath1))
                for j in (xrange(i, len(subpath1)) if subpath_i == subpath_j else xrange(len(subpath1))):
                    if subpath_i == subpath_j and j == i:
                        # Find self intersections of a segment
                        sp1, sp2, sp3 = csp_split(subpath[i - 1], subpath[i], .5)
                        intersections = csp_segments_intersection(sp1, sp2, sp2, sp3)
                        summ += 1
                        for t in intersections:
                            summ1 += 1
                            if not (small(t[0] - 1) and small(t[1])) and 0 <= t[0] <= 1 and 0 <= t[1] <= 1:
                                intersection[subpath_i] += [[i, t[0] / 2], [j, t[1] / 2 + .5]]
                    else:
                        intersections = csp_segments_intersection(subpath[i - 1], subpath[i], subpath1[j - 1],
                                                                  subpath1[j])
                        summ += 1
                        for t in intersections:
                            summ1 += 1
                            # TODO tolerance dependence to cpsp_length(t)
                            if len(t) == 2 and 0 <= t[0] <= 1 and 0 <= t[1] <= 1 and not (
                                    subpath_i == subpath_j and (
                                    (j - i - 1) % (len(subpath) - 1) == 0 and small(t[0] - 1) and small(t[1]) or
                                    (i - j - 1) % (len(subpath) - 1) == 0 and small(t[1] - 1) and small(t[0]))):
                                intersection[subpath_i] += [[i, t[0]]]
                                intersection[subpath_j] += [[j, t[1]]]
                                # draw_pointer(csp_at_t(subpath[i-1],subpath[i],t[0]),"#f00")
                                # print_(t)
                                # print_(i,j)
                            elif len(t) == 5 and t[4] == "Overlap":
                                intersection[subpath_i] += [[i, t[0]], [i, t[1]]]
                                intersection[subpath_j] += [[j, t[1]], [j, t[3]]]

    print_("Intersections found in %s" % (time.time() - time_))
    print_("Examined %s segments" % (summ))
    print_("found %s intersections" % (summ1))
    time_ = time.time()

    ########################################################################
    # Split unclipped offset by intersection points into splitted_offset
    ########################################################################
    splitted_offset = []
    for i in xrange(csp_len):
        subpath = unclipped_offset[i]
        if len(intersection[i]) > 0:
            parts = csp_subpath_split_by_points(subpath, intersection[i])
            # Close    parts list to close path (The first and the last parts are joined together)
            if [1, 0.] not in intersection[i]:
                parts[0][0][0] = parts[-1][-1][0]
                parts[0] = csp_concat_subpaths(parts[-1], parts[0])
                splitted_offset += parts[:-1]
            else:
                splitted_offset += parts[:]
        else:
            splitted_offset += [subpath[:]]

    # for i in range(len(splitted_offset)):
    #    csp_draw([splitted_offset[i]], color = ["Green","Red","Blue"][i%3])
    print_("Splitted in %s" % (time.time() - time_))
    time_ = time.time()

    ########################################################################
    # Clipping
    ########################################################################
    result = []
    for subpath_i in range(len(splitted_offset)):
        clip = False
        s1 = splitted_offset[subpath_i]
        for subpath_j in range(len(splitted_offset)):
            s2 = splitted_offset[subpath_j]
            if (P(s1[0][1]) - P(s2[-1][1])).l2() < 0.0001 and ((subpath_i + 1) % len(splitted_offset) != subpath_j):
                if dot(csp_normalized_normal(s2[-2], s2[-1], 1.), csp_normalized_slope(s1[0], s1[1], 0.)) * r < -0.0001:
                    clip = True
                    break
            if (P(s2[0][1]) - P(s1[-1][1])).l2() < 0.0001 and ((subpath_j + 1) % len(splitted_offset) != subpath_i):
                if dot(csp_normalized_normal(s2[0], s2[1], 0.), csp_normalized_slope(s1[-2], s1[-1], 1.)) * r > 0.0001:
                    clip = True
                    break

        if not clip:
            result += [s1[:]]
        elif options.offset_draw_clippend_path:
            csp_draw([s1], color="Red", width=.1)
            draw_pointer(csp_at_t(s2[-2], s2[-1], 1.) +
                         (P(csp_at_t(s2[-2], s2[-1], 1.)) + P(
                             csp_normalized_normal(s2[-2], s2[-1], 1.)) * 10).to_list(), "Green", "line")
            draw_pointer(csp_at_t(s1[0], s1[1], 0.) +
                         (P(csp_at_t(s1[0], s1[1], 0.)) + P(csp_normalized_slope(s1[0], s1[1], 0.)) * 10).to_list(),
                         "Red", "line")

    # Now join all together and check closure and orientation of result
    joined_result = csp_join_subpaths(result)
    # Check if each subpath from joined_result is closed
    # csp_draw(joined_result,color="Green",width=1)


    for s in joined_result[:]:
        if csp_subpaths_end_to_start_distance2(s, s) > 0.001:
            # Remove open parts
            if options.offset_draw_clippend_path:
                csp_draw([s], color="Orange", width=1)
                draw_pointer(s[0][1], comment=csp_subpaths_end_to_start_distance2(s, s))
                draw_pointer(s[-1][1], comment=csp_subpaths_end_to_start_distance2(s, s))
            joined_result.remove(s)
        else:
            # Remove small parts
            minx, miny, maxx, maxy = csp_true_bounds([s])
            if (minx[0] - maxx[0]) ** 2 + (miny[1] - maxy[1]) ** 2 < 0.1:
                joined_result.remove(s)
    print_("Clipped and joined path in %s" % (time.time() - time_))
    time_ = time.time()

    ########################################################################
    # Now to the Dummy cliping: remove parts from splitted offset if their
    # centers are  closer to the original path than offset radius.
    ########################################################################

    r1, r2 = ((0.99 * r) ** 2, (1.01 * r) ** 2) if abs(r * .01) < 1 else ((abs(r) - 1) ** 2, (abs(r) + 1) ** 2)
    for s in joined_result[:]:
        dist = csp_to_point_distance(original_csp, s[int(len(s) / 2)][1], dist_bounds=[r1, r2], tolerance=.000001)
        if not r1 < dist[0] < r2:
            joined_result.remove(s)
            if options.offset_draw_clippend_path:
                csp_draw([s], comment=math.sqrt(dist[0]))
                draw_pointer(
                    csp_at_t(csp[dist[1]][dist[2] - 1], csp[dist[1]][dist[2]], dist[3]) + s[int(len(s) / 2)][1], "blue",
                    "line", comment=[math.sqrt(dist[0]), i, j, sp])

    print_("-----------------------------")
    print_("Total offset time %s" % (time.time() - time_start))
    print_()
    return joined_result


################################################################################
###
###        Biarc function
###
###        Calculates biarc approximation of cubic super path segment
###        splits segment if needed or approximates it with straight line
###
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


def biarc_curve_segment_length(seg):
    if seg[1] == "arc":
        return math.sqrt((seg[0][0] - seg[2][0]) ** 2 + (seg[0][1] - seg[2][1]) ** 2) * seg[3]
    elif seg[1] == "line":
        return math.sqrt((seg[0][0] - seg[4][0]) ** 2 + (seg[0][1] - seg[4][1]) ** 2)
    else:
        return 0


################################################################################
###        Polygon class
################################################################################
class Polygon:
    def __init__(self, polygon=None):
        self.polygon = [] if polygon == None else polygon[:]


    def move(self, x, y):
        for i in range(len(self.polygon)):
            for j in range(len(self.polygon[i])):
                self.polygon[i][j][0] += x
                self.polygon[i][j][1] += y


    def bounds(self):
        minx, miny, maxx, maxy = 1e400, 1e400, -1e400, -1e400
        for poly in self.polygon:
            for p in poly:
                if minx > p[0]: minx = p[0]
                if miny > p[1]: miny = p[1]
                if maxx < p[0]: maxx = p[0]
                if maxy < p[1]: maxy = p[1]
        return minx * 1, miny * 1, maxx * 1, maxy * 1


    def width(self):
        b = self.bounds()
        return b[2] - b[0]


    def rotate_(self, sin, cos):
        for i in range(len(self.polygon)):
            for j in range(len(self.polygon[i])):
                x, y = self.polygon[i][j][0], self.polygon[i][j][1]
                self.polygon[i][j][0] = x * cos - y * sin
                self.polygon[i][j][1] = x * sin + y * cos


    def rotate(self, a):
        cos, sin = math.cos(a), math.sin(a)
        self.rotate_(sin, cos)


    def drop_into_direction(self, direction, surface):
        # Polygon is a list of simple polygons
        # Surface is a polygon + line y = 0
        # Direction is [dx,dy]
        if len(self.polygon) == 0 or len(self.polygon[0]) == 0: return
        if direction[0] ** 2 + direction[1] ** 2 < 1e-10: return
        direction = normalize(direction)
        sin, cos = direction[0], -direction[1]
        self.rotate_(-sin, cos)
        surface.rotate_(-sin, cos)
        self.drop_down(surface, zerro_plane=False)
        self.rotate_(sin, cos)
        surface.rotate_(sin, cos)


    def centroid(self):
        centroids = []
        sa = 0
        for poly in self.polygon:
            cx, cy, a = 0, 0, 0
            for i in range(len(poly)):
                [x1, y1], [x2, y2] = poly[i - 1], poly[i]
                cx += (x1 + x2) * (x1 * y2 - x2 * y1)
                cy += (y1 + y2) * (x1 * y2 - x2 * y1)
                a += (x1 * y2 - x2 * y1)
            a *= 3.
            if abs(a) > 0:
                cx /= a
                cy /= a
                sa += abs(a)
                centroids += [[cx, cy, a]]
        if sa == 0: return [0., 0.]
        cx, cy = 0., 0.
        for c in centroids:
            cx += c[0] * c[2]
            cy += c[1] * c[2]
        cx /= sa
        cy /= sa
        return [cx, cy]


    def drop_down(self, surface, zerro_plane=True):
        # Polygon is a list of simple polygons
        # Surface is a polygon + line y = 0
        # Down means min y (0,-1)
        if len(self.polygon) == 0 or len(self.polygon[0]) == 0: return
        # Get surface top point
        top = surface.bounds()[3]
        if zerro_plane: top = max(0, top)
        # Get polygon bottom point
        bottom = self.bounds()[1]
        self.move(0, top - bottom + 10)
        # Now get shortest distance from surface to polygon in positive x=0 direction
        # Such distance = min(distance(vertex, edge)...)  where edge from surface and
        # vertex from polygon and vice versa...
        dist = 1e300
        for poly in surface.polygon:
            for i in range(len(poly)):
                for poly1 in self.polygon:
                    for i1 in range(len(poly1)):
                        st, end = poly[i - 1], poly[i]
                        vertex = poly1[i1]
                        if st[0] <= vertex[0] <= end[0] or end[0] <= vertex[0] <= st[0]:
                            if st[0] == end[0]:
                                d = min(vertex[1] - st[1], vertex[1] - end[1])
                            else:
                                d = vertex[1] - st[1] - (end[1] - st[1]) * (vertex[0] - st[0]) / (end[0] - st[0])
                            if dist > d: dist = d
                        # and vice versa just change the sign because vertex now under the edge
                        st, end = poly1[i1 - 1], poly1[i1]
                        vertex = poly[i]
                        if st[0] <= vertex[0] <= end[0] or end[0] <= vertex[0] <= st[0]:
                            if st[0] == end[0]:
                                d = min(- vertex[1] + st[1], -vertex[1] + end[1])
                            else:
                                d = - vertex[1] + st[1] + (end[1] - st[1]) * (vertex[0] - st[0]) / (end[0] - st[0])
                            if dist > d: dist = d

        if zerro_plane and dist > 10 + top: dist = 10 + top
        # print_(dist, top, bottom)
        # self.draw()
        self.move(0, -dist)


    def draw(self, color="#075", width=.1):
        for poly in self.polygon:
            csp_draw([csp_subpath_line_to([], poly + [poly[0]])], color=color, width=width)


    def add(self, add):
        if type(add) == type([]):
            self.polygon += add[:]
        else:
            self.polygon += add.polygon[:]


    def point_inside(self, p):
        inside = False
        for poly in self.polygon:
            for i in range(len(poly)):
                st, end = poly[i - 1], poly[i]
                if p == st or p == end: return True  # point is a vertex = point is on the edge
                if st[0] > end[0]: st, end = end, st  # This will be needed to check that edge if open only at rigth end
                c = (p[1] - st[1]) * (end[0] - st[0]) - (end[1] - st[1]) * (p[0] - st[0])
                # print_(c)
                if st[0] <= p[0] < end[0]:
                    if c < 0:
                        inside = not inside
                    elif c == 0:
                        return True  # point is on the edge
                elif st[0] == end[0] == p[0] and (
                        st[1] <= p[1] <= end[1] or end[1] <= p[1] <= st[1]):  # point is on the edge
                    return True
        return inside


    def hull(self):
        # Add vertices at all self intersection points.
        hull = []
        for i1 in range(len(self.polygon)):
            poly1 = self.polygon[i1]
            poly_ = []
            for j1 in range(len(poly1)):
                s, e = poly1[j1 - 1], poly1[j1]
                poly_ += [s]

                # Check self intersections
                for j2 in range(j1 + 1, len(poly1)):
                    s1, e1 = poly1[j2 - 1], poly1[j2]
                    int_ = line_line_intersection_points(s, e, s1, e1)
                    for p in int_:
                        if point_to_point_d2(p, s) > 0.000001 and point_to_point_d2(p, e) > 0.000001:
                            poly_ += [p]
                # Check self intersections with other polys
                for i2 in range(len(self.polygon)):
                    if i1 == i2: continue
                    poly2 = self.polygon[i2]
                    for j2 in range(len(poly2)):
                        s1, e1 = poly2[j2 - 1], poly2[j2]
                        int_ = line_line_intersection_points(s, e, s1, e1)
                        for p in int_:
                            if point_to_point_d2(p, s) > 0.000001 and point_to_point_d2(p, e) > 0.000001:
                                poly_ += [p]
            hull += [poly_]
        # Create the dictionary containing all edges in both directions
        edges = {}
        for poly in self.polygon:
            for i in range(len(poly)):
                s, e = tuple(poly[i - 1]), tuple(poly[i])
                if (point_to_point_d2(e, s) < 0.000001): continue
                break_s, break_e = False, False
                for p in edges:
                    if point_to_point_d2(p, s) < 0.000001:
                        break_s = True
                        s = p
                    if point_to_point_d2(p, e) < 0.000001:
                        break_e = True
                        e = p
                    if break_s and break_e: break
                l = point_to_point_d(s, e)
                if not break_s and not break_e:
                    edges[s] = [[s, e, l]]
                    edges[e] = [[e, s, l]]
                    # draw_pointer(s+e,"red","line")
                    # draw_pointer(s+e,"red","line")
                else:
                    if e in edges:
                        for edge in edges[e]:
                            if point_to_point_d2(edge[1], s) < 0.000001:
                                break
                        if point_to_point_d2(edge[1], s) > 0.000001:
                            edges[e] += [[e, s, l]]
                            # draw_pointer(s+e,"red","line")

                    else:
                        edges[e] = [[e, s, l]]
                        # draw_pointer(s+e,"green","line")
                    if s in edges:
                        for edge in edges[s]:
                            if point_to_point_d2(edge[1], e) < 0.000001:
                                break
                        if point_to_point_d2(edge[1], e) > 0.000001:
                            edges[s] += [[s, e, l]]
                            # draw_pointer(s+e,"red","line")
                    else:
                        edges[s] = [[s, e, l]]
                        # draw_pointer(s+e,"green","line")


        def angle_quadrant(sin, cos):
            # quadrants are (0,pi/2], (pi/2,pi], (pi,3*pi/2], (3*pi/2, 2*pi], i.e. 0 is in the 4-th quadrant
            if sin > 0 and cos >= 0: return 1
            if sin >= 0 and cos < 0: return 2
            if sin < 0 and cos <= 0: return 3
            if sin <= 0 and cos > 0: return 4


        def angle_is_less(sin, cos, sin1, cos1):
            # 0 = 2*pi is the largest angle
            if [sin1, cos1] == [0, 1]: return True
            if [sin, cos] == [0, 1]: return False
            if angle_quadrant(sin, cos) > angle_quadrant(sin1, cos1):
                return False
            if angle_quadrant(sin, cos) < angle_quadrant(sin1, cos1):
                return True
            if sin >= 0 and cos > 0: return sin < sin1
            if sin > 0 and cos <= 0: return sin > sin1
            if sin <= 0 and cos < 0: return sin > sin1
            if sin < 0 and cos >= 0: return sin < sin1


        def get_closes_edge_by_angle(edges, last):
            # Last edge is normalized vector of the last edge.
            min_angle = [0, 1]
            next = last
            last_edge = [(last[0][0] - last[1][0]) / last[2], (last[0][1] - last[1][1]) / last[2]]
            for p in edges:
                # draw_pointer(list(p[0])+[p[0][0]+last_edge[0]*40,p[0][1]+last_edge[1]*40], "Red", "line", width=1)
                # print_("len(edges)=",len(edges))
                cur = [(p[1][0] - p[0][0]) / p[2], (p[1][1] - p[0][1]) / p[2]]
                cos, sin = dot(cur, last_edge), cross(cur, last_edge)
                # draw_pointer(list(p[0])+[p[0][0]+cur[0]*40,p[0][1]+cur[1]*40], "Orange", "line", width=1, comment = [sin,cos])
                # print_("cos, sin=",cos,sin)
                # print_("min_angle_before=",min_angle)

                if angle_is_less(sin, cos, min_angle[0], min_angle[1]):
                    min_angle = [sin, cos]
                    next = p
                # print_("min_angle=",min_angle)

            return next

        # Join edges together into new polygon cutting the vertexes inside new polygon
        self.polygon = []
        len_edges = sum([len(edges[p]) for p in edges])
        loops = 0

        while len(edges) > 0:
            poly = []
            if loops > len_edges: raise ValueError, "Hull error"
            loops += 1
            # Find left most vertex.
            start = (1e100, 1)
            for edge in edges:
                start = min(start, min(edges[edge]))
            last = [(start[0][0] - 1, start[0][1]), start[0], 1]
            first_run = True
            loops1 = 0
            while (last[1] != start[0] or first_run):
                first_run = False
                if loops1 > len_edges: raise ValueError, "Hull error"
                loops1 += 1
                next = get_closes_edge_by_angle(edges[last[1]], last)
                # draw_pointer(next[0]+next[1],"Green","line", comment=i, width= 1)
                # print_(next[0],"-",next[1])

                last = next
                poly += [list(last[0])]
            self.polygon += [poly]
            # Remove all edges that are intersects new poly (any vertex inside new poly)
            poly_ = Polygon([poly])
            for p in edges.keys()[:]:
                if poly_.point_inside(list(p)): del edges[p]
        self.draw(color="Green", width=1)


class Arangement_Genetic:
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


    def add_random_species(self, count):
        for i in range(count):
            specimen = []
            order = range(self.genes_count)
            random.shuffle(order)
            for j in order:
                specimen += [[j, random.random(), random.random()]]
            self.population += [[None, specimen]]


    def species_distance2(self, sp1, sp2):
        # retun distance, each component is normalized
        s = 0
        for j in range(self.genes_count):
            s += ((sp1[j][0] - sp2[j][0]) / self.genes_count) ** 2 + ((sp1[j][1] - sp2[j][1])) ** 2 + (
            (sp1[j][2] - sp2[j][2])) ** 2
        return s

    def leave_top_species(self, count):
        self.population.sort()
        res = [copy.deepcopy(self.population[0])]
        del self.population[0]
        for i in range(count - 1):
            t = []
            for j in range(20):
                i1 = random.randint(0, len(self.population) - 1)
                t += [[self.population[i1][0], i1]]
            t.sort()
            res += [copy.deepcopy(self.population[t[0][1]])]
            del self.population[t[0][1]]
        self.population = res
        # del self.population[0]
        # for c in range(count-1) :
        #    rank = []
        #    for i in range(len(self.population)) :
        #        sim = self.similarity(self.population[i][1],res)
        #        rank += [ [self.population[i][0] / sim if sim>0 else 1e100,i] ]
        #    rank.sort()
        #    res += [  copy.deepcopy(self.population[rank[0][1]]) ]
        #    print_(rank[0],self.population[rank[0][1]][0])
        #    print_(res[-1])
        #    del self.population[rank[0][1]]

        self.population = res


    def populate_species(self, count, parent_count):
        self.population.sort()
        self.inc = 0
        for c in range(count):
            parent1 = random.randint(0, parent_count - 1)
            parent2 = random.randint(0, parent_count - 1)
            if parent1 == parent2: parent2 = (parent2 + 1) % parent_count
            parent1, parent2 = self.population[parent1][1], self.population[parent2][1]
            i1, i2 = 0, 0
            genes_order = []
            specimen = [[0, 0., 0.] for i in range(self.genes_count)]

            self.incest_mutation_multiplyer = 1.
            self.incest_mutation_count_multiplyer = 1.

            if self.species_distance2(parent1, parent2) <= .01 / self.genes_count:
                # OMG it's a incest :O!!!
                # Damn you bastards!
                self.inc += 1
                self.incest_mutation_multiplyer = 2.
                self.incest_mutation_count_multiplyer = 2.
            else:
                if random.random() < .01: print_(self.species_distance2(parent1, parent2))
            start_gene = random.randint(0, self.genes_count)
            end_gene = (max(1, random.randint(0, self.genes_count),
                            int(self.genes_count / 4)) + start_gene) % self.genes_count
            if end_gene < start_gene:
                end_gene, start_gene = start_gene, end_gene
                parent1, parent2 = parent2, parent1
            for i in range(start_gene, end_gene):
                # rotation_mutate_param = random.random()/100
                # xposition_mutate_param = random.random()/100
                tr = 1.  # - rotation_mutate_param
                tp = 1.  # - xposition_mutate_param
                specimen[i] = [parent1[i][0], parent1[i][1] * tr + parent2[i][1] * (1 - tr),
                               parent1[i][2] * tp + parent2[i][2] * (1 - tp)]
                genes_order += [parent1[i][0]]

            for i in range(0, start_gene) + range(end_gene, self.genes_count):
                tr = 0.  # rotation_mutate_param
                tp = 0.  # xposition_mutate_param
                j = i
                while parent2[j][0] in genes_order:
                    j = (j + 1) % self.genes_count
                specimen[i] = [parent2[j][0], parent1[i][1] * tr + parent2[i][1] * (1 - tr),
                               parent1[i][2] * tp + parent2[i][2] * (1 - tp)]
                genes_order += [parent2[j][0]]

            for i in range(random.randint(self.mutation_genes_count[0],
                                          self.mutation_genes_count[0] * self.incest_mutation_count_multiplyer)):
                if random.random() < self.order_mutate_factor * self.incest_mutation_multiplyer:
                    i1, i2 = random.randint(0, self.genes_count - 1), random.randint(0, self.genes_count - 1)
                    specimen[i1][0], specimen[i2][0] = specimen[i2][0], specimen[i1][0]
                if random.random() < self.move_mutation_factor * self.incest_mutation_multiplyer:
                    i1 = random.randint(0, self.genes_count - 1)
                    specimen[i1][1] = (specimen[i1][
                                           1] + random.random() * math.pi2 * self.move_mutation_multiplier) % 1.
                    specimen[i1][2] = (specimen[i1][2] + random.random() * self.move_mutation_multiplier) % 1.
            self.population += [[None, specimen]]


    def test_spiece_drop_down(self, spiece):
        surface = Polygon()
        for p in spiece:
            time_ = time.time()
            poly = Polygon(copy.deepcopy(self.polygons[p[0]].polygon))
            poly.rotate(p[1] * math.pi2)
            w = poly.width()
            left = poly.bounds()[0]
            poly.move(-left + (self.width - w) * p[2], 0)
            poly.drop_down(surface)
            surface.add(poly)
        return surface


    def test(self, test_function):
        for i in range(len(self.population)):
            if self.population[i][0] == None:
                surface = test_function(self.population[i][1])
                b = surface.bounds()
                self.population[i][0] = (b[3] - b[1]) * (b[2] - b[0])
        self.population.sort()


    def test_spiece_centroid(self, spiece):
        poly = Polygon(copy.deepcopy(self.polygons[spiece[0][0]].polygon))
        poly.rotate(spiece[0][2] * math.pi2)
        surface = Polygon(poly.polygon)
        i = 0
        for p in spiece[1:]:
            i += 1
            poly = Polygon(copy.deepcopy(self.polygons[p[0]].polygon))
            poly.rotate(p[2] * math.pi2)
            c = surface.centroid()
            c1 = poly.centroid()
            direction = [math.cos(p[1] * math.pi2), -math.sin(p[1] * math.pi2)]
            poly.move(c[0] - c1[0] - direction[0] * 100, c[1] - c1[1] - direction[1] * 100)
            poly.drop_into_direction(direction, surface)
            surface.add(poly)
        return surface


        # surface.draw()


################################################################################
###
###        Gcodetools class
###
################################################################################

class laser_gcode(inkex.Effect):

    def export_gcode(self, gcode):
        gcode_pass = gcode
        for x in range(1, self.options.passes):
            gcode += "G91\nG1 Z-" + self.options.pass_depth + "\nG90\n" + gcode_pass
        f = open(self.options.directory + self.options.file, "w")
        f.write(
            self.options.laser_off_command + " S0" + "\n" + self.header + "G1 F" + self.options.travel_speed + "\n" + gcode + self.footer)
        f.close()

    def __init__(self):
        inkex.Effect.__init__(self)
        self.OptionParser.add_option("-d", "--directory", action="store", type="string", dest="directory", default="",
                                     help="Output directory")
        self.OptionParser.add_option("-f", "--filename", action="store", type="string", dest="file",
                                     default="output.gcode", help="File name")
        self.OptionParser.add_option("", "--add-numeric-suffix-to-filename", action="store", type="inkbool",
                                     dest="add_numeric_suffix_to_filename", default=False,
                                     help="Add numeric suffix to file name")
        self.OptionParser.add_option("", "--laser-command", action="store", type="string", dest="laser_command",
                                     default="M03", help="Laser gcode command")
        self.OptionParser.add_option("", "--laser-off-command", action="store", type="string", dest="laser_off_command",
                                     default="M05", help="Laser gcode end command")
        self.OptionParser.add_option("", "--laser-speed", action="store", type="int", dest="laser_speed", default="750",
                                     help="Laser speed (mm/min)")
        self.OptionParser.add_option("", "--travel-speed", action="store", type="string", dest="travel_speed",
                                     default="3000", help="Travel speed (mm/min)")
        self.OptionParser.add_option("", "--laser-power", action="store", type="int", dest="laser_power", default="255",
                                     help="S# is 256 or 10000 for full power")
        self.OptionParser.add_option("", "--passes", action="store", type="int", dest="passes", default="1",
                                     help="Quantity of passes")
        self.OptionParser.add_option("", "--pass-depth", action="store", type="string", dest="pass_depth", default="1",
                                     help="Depth of laser cut")
        self.OptionParser.add_option("", "--power-delay", action="store", type="string", dest="power_delay",
                                     default="0", help="Laser power-on delay (ms)")
        self.OptionParser.add_option("", "--suppress-all-messages", action="store", type="inkbool",
                                     dest="suppress_all_messages", default=True,
                                     help="Hide messages during g-code generation")
        self.OptionParser.add_option("", "--create-log", action="store", type="inkbool", dest="log_create_log",
                                     default=False, help="Create log files")
        self.OptionParser.add_option("", "--log-filename", action="store", type="string", dest="log_filename",
                                     default='', help="Create log files")
        self.OptionParser.add_option("", "--engraving-draw-calculation-paths", action="store", type="inkbool",
                                     dest="engraving_draw_calculation_paths", default=False,
                                     help="Draw additional graphics to debug engraving path")
        self.OptionParser.add_option("", "--unit", action="store", type="string", dest="unit",
                                     default="G21 (All units in mm)", help="Units either mm or inches")
        self.OptionParser.add_option("", "--active-tab", action="store", type="string", dest="active_tab", default="",
                                     help="Defines which tab is active")
        self.OptionParser.add_option("", "--biarc-max-split-depth", action="store", type="int",
                                     dest="biarc_max_split_depth", default="4",
                                     help="Defines maximum depth of splitting while approximating using biarcs.")

    def parse_curve(self, p, layer, w=None, f=None):
        c = []
        if len(p) == 0:
            return []
        p = self.transform_csp(p, layer)

        ### Sort to reduce Rapid distance
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
        # Add marker to defs if it doesnot exists
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

        if group == None:
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
    ###
    ###        Generate Gcode
    ###        Generates Gcode on given curve.
    ###
    ###        Crve defenitnion [start point, type = {'arc','line','move','end'}, arc center, arc angle, end point, [zstart, zend]]
    ###
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

        def calculate_angle(a, current_a):
            return min(
                [abs(a - current_a % math.pi2 + math.pi2), a + current_a - current_a % math.pi2 + math.pi2],
                [abs(a - current_a % math.pi2 - math.pi2), a + current_a - current_a % math.pi2 - math.pi2],
                [abs(a - current_a % math.pi2), a + current_a - current_a % math.pi2])[1]

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

            ###self.Zauto_scale[layer]  = math.sqrt( (self.transform_matrix[layer][0][0]**2 + self.transform_matrix[layer][1][1]**2)/2 )
            ### Zautoscale is absolete
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
    ###        Errors handling function, notes are just printed into Logfile,
    ###        warnings are printed into log file and warning message is displayed but
    ###        extension continues working, errors causes log and execution is halted
    ###        Notes, warnings adn errors could be assigned to space or comma or dot
    ###        sepparated strings (case is ignoreg).
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
    ###        Get defs from svg
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
    ###
    ###        Get Gcodetools info from the svg
    ###
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
    ###
    ###        dxfpoints
    ###
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
    ###
    ###        Laser
    ###
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
    ###
    ###        Orientation
    ###
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
    ###
    ###        Effect
    ###
    ###        Main function of Gcodetools class
    ###
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


e = laser_gcode()
e.affect()
