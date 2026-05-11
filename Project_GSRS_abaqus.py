# -*- coding: utf-8 -*-
"""Geodesic lattice shell generator for Abaqus (短程线网壳).

Builds a spherical-cap wire frame consisting of:
  - An upper icosahedral cap subdivided at frequency f1
  - A lower antiprism-like band subdivided at frequency f2
All nodes are projected onto a sphere sized from the given span and rise.
"""
from __future__ import division, print_function
import math

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup


PENTAGON = 5
TWO_PI = 2.0 * math.pi


# --- vector helpers -------------------------------------------------------

def vec_add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])

def vec_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

def vec_scale(a, k):
    return (a[0] * k, a[1] * k, a[2] * k)

def distance(a, b):
    dx, dy, dz = a[0] - b[0], a[1] - b[1], a[2] - b[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)

def project_on_sphere(p, radius):
    """Scale p onto the sphere of given radius (equivalent to r * p / |p|)."""
    length = math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2])
    if length == 0:
        return (0.0, 0.0, 0.0)
    k = radius / length
    return (p[0] * k, p[1] * k, p[2] * k)


def _interpolate(start, end, steps):
    """Return `steps` points start, start+step, ..., end-step (open at end)."""
    step = vec_scale(vec_sub(end, start), 1.0 / steps)
    return [vec_add(start, vec_scale(step, j)) for j in range(steps)]


def _variance(values):
    n = len(values)
    if n == 0:
        return 0.0
    mean = sum(values) / n
    return sum([(v - mean) ** 2 for v in values]) / n


# --- geometry -------------------------------------------------------------

class ShellGeometry(object):
    """Pure geometry: builds node coordinates and bar lengths, no Abaqus calls."""

    def __init__(self, f1, f2, span, rise, adjustment):
        self.f1 = f1
        self.f2 = f2
        self.radius = rise / 2.0 + (span * span / 8.0) / rise
        self.z_angle = math.acos((self.radius - rise) / self.radius)
        base_split = f1 / float(f1 + f2)
        self.z_angle_upper = (base_split + adjustment) * self.z_angle
        self.upper_layers = self._build_upper()   # upper_layers[m][q] : m in 0..f1
        self.lower_layers = self._build_lower()   # lower_layers[i][q] : i in 0..f2

    def _build_upper(self):
        """Upper cap: icosahedral subdivision, ring m has PENTAGON*m nodes."""
        f1, radius = self.f1, self.radius
        apex = (0.0, 0.0, 1.0)

        # Five rays from the apex down to the base ring of the cap.
        base_ring = [
            (math.sin(self.z_angle_upper) * math.cos(TWO_PI / PENTAGON * k),
             math.sin(self.z_angle_upper) * math.sin(TWO_PI / PENTAGON * k),
             math.cos(self.z_angle_upper))
            for k in range(PENTAGON)
        ]
        # Points along each edge from apex to base_ring[k], f1+1 points each.
        edges = [
            [vec_add(apex, vec_scale(vec_sub(base, apex), m / float(f1)))
             for m in range(f1 + 1)]
            for base in base_ring
        ]

        layers = [[(0.0, 0.0, radius)]]  # ring 0 is just the apex
        for m in range(1, f1 + 1):
            ring = []
            for k in range(PENTAGON):
                a = edges[k][m]
                b = edges[(k + 1) % PENTAGON][m]
                ring.extend(_interpolate(a, b, m))
            layers.append(ring)

        for m in range(1, f1 + 1):
            layers[m] = [project_on_sphere(p, radius) for p in layers[m]]
        return layers

    def _build_lower(self):
        """Lower band between the cap's base ring and the equator-side ring."""
        f1, f2, radius = self.f1, self.f2, self.radius
        pentagon_step = TWO_PI / PENTAGON
        offset = f2 / float(f1) * math.pi / PENTAGON
        aux_offset = (f1 - f2) / float(f1) * pentagon_step

        # ring_top: 5 vertices of the cap's base (unit sphere), one per pentagon slot.
        ring_top = [
            vec_scale(self.upper_layers[f1][i], 1.0 / radius)
            for i in range(PENTAGON * f1) if i % f1 == 0
        ]
        # ring_bot: main vertices on the lower ring.
        ring_bot = [
            (math.sin(self.z_angle) * math.cos(offset + i * pentagon_step),
             math.sin(self.z_angle) * math.sin(offset + i * pentagon_step),
             math.cos(self.z_angle))
            for i in range(PENTAGON)
        ]
        # ring_aux: auxiliary vertices used for odd seed columns.
        ring_aux = [
            (math.sin(self.z_angle) * math.cos(offset + aux_offset + i * pentagon_step),
             math.sin(self.z_angle) * math.sin(offset + aux_offset + i * pentagon_step),
             math.cos(self.z_angle))
            for i in range(PENTAGON)
        ]

        # seeds: 2*PENTAGON columns from ring_top down to ring_bot (f2+1 points each).
        # Even-indexed columns use ring_bot; odd-indexed columns use ring_aux.
        seeds = []
        m = n = 0
        while n < PENTAGON:
            bottom = ring_bot[n]
            top = ring_top[m]
            step = vec_scale(vec_sub(bottom, top), 1.0 / f2)
            seeds.append([vec_add(top, vec_scale(step, j)) for j in range(f2 + 1)])
            if abs(m - n) > 2:
                break
            if m > n:
                n += 1
                continue
            if m == n:
                m += 1
                ring_bot[n] = ring_aux[n]     # next seed uses the auxiliary vertex
            if m == PENTAGON:
                m = 0

        # Fill each layer i in [0, f2] by sweeping between adjacent seed columns.
        layers = []
        total = 2 * PENTAGON
        for i in range(f2 + 1):
            ring = []
            for m in range(total):
                next_m = (m + 1) % total
                if m % 2 == 0 and i != f1:
                    ring.extend(_interpolate(seeds[m][i], seeds[next_m][i], f1 - i))
                elif m % 2 != 0 and i != 0:
                    ring.extend(_interpolate(seeds[m][i], seeds[next_m][i], i))
            layers.append(ring)

        for i in range(f2 + 1):
            for q in range(PENTAGON * f1):
                layers[i][q] = project_on_sphere(layers[i][q], radius)
        # Stitch top ring of lower band to the cap's base ring (shared nodes).
        layers[0] = list(self.upper_layers[f1])
        return layers

    def bar_lengths(self):
        """Return the length of every bar (for variance-based adjustment search)."""
        lengths = []
        self._collect_upper_bars(lengths, connect_only=False)
        self._collect_lower_bars(lengths, connect_only=False)
        return lengths

    def _collect_upper_bars(self, out, connect_only, connect=None):
        """Emit bars for the upper cap. If connect_only, call connect(p1, p2) instead."""
        f1 = self.f1
        layers = self.upper_layers
        for m in range(1, f1 + 1):
            ring_len = PENTAGON * m
            for q in range(ring_len):
                self._emit(out, connect_only, connect,
                           layers[m][q], layers[m][(q + 1) % ring_len])

        # Radial bars between layer m and m-1.
        for m in range(f1, 0, -1):
            ring_len = PENTAGON * m
            prev_len = PENTAGON * (m - 1)
            non_corner = []
            for q in range(ring_len):
                if q % m == 0:
                    # corner of the pentagonal patch: connects to the previous-ring corner
                    prev_idx = (q // m) * (m - 1) if prev_len > 0 else 0
                    self._emit(out, connect_only, connect,
                               layers[m][q], layers[m - 1][prev_idx])
                else:
                    non_corner.append(q)
            # Remaining nodes in this ring connect to two neighbours on ring m-1.
            for idx, q in enumerate(non_corner):
                p = layers[m][q]
                a = layers[m - 1][idx] if prev_len > 0 else layers[m - 1][0]
                b = layers[m - 1][(idx + 1) % prev_len] if prev_len > 0 else layers[m - 1][0]
                self._emit(out, connect_only, connect, p, a)
                self._emit(out, connect_only, connect, p, b)

    def _collect_lower_bars(self, out, connect_only, connect=None):
        f1, f2 = self.f1, self.f2
        layers = self.lower_layers
        ring_len = PENTAGON * f1
        # Horizontal bars per ring (skip the top ring; already emitted as the cap base).
        for m in range(1, f2 + 1):
            for q in range(ring_len):
                self._emit(out, connect_only, connect,
                           layers[m][q], layers[m][(q + 1) % ring_len])
        # Diagonal bars between adjacent rings.
        for m in range(f2):
            for q in range(ring_len):
                q_prev = q - 1 if q > 0 else ring_len - 1
                self._emit(out, connect_only, connect,
                           layers[m][q], layers[m + 1][q])
                self._emit(out, connect_only, connect,
                           layers[m][q], layers[m + 1][q_prev])

    @staticmethod
    def _emit(out, connect_only, connect, a, b):
        if connect_only:
            connect(a, b)
        else:
            out.append(distance(a, b))


# --- adjustment search ----------------------------------------------------

def search_best_adjustment(f1, f2, span, rise,
                           candidates=None):
    """Pick the adjustment coefficient that minimises bar-length variance."""
    if candidates is None:
        candidates = [i / 100.0 for i in range(-2, 13)]  # -0.02 .. 0.12 step 0.01
    best_adj = 0.0
    best_var = float("inf")
    for adj in candidates:
        geom = ShellGeometry(f1, f2, span, rise, adj)
        var = _variance(geom.bar_lengths())
        if var < best_var:
            best_var = var
            best_adj = adj
    return best_adj


# --- Abaqus part builder --------------------------------------------------

class GeodesicShellPart(object):
    def __init__(self, part_name, f1, f2, span, rise, model_name='Model-1'):
        adjustment = search_best_adjustment(f1, f2, span, rise)
        self.geometry = ShellGeometry(f1, f2, span, rise, adjustment)
        mdb.models[model_name].Part(
            name=part_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        self.part = mdb.models[model_name].parts[part_name]

    def build(self):
        geom = self.geometry
        connect = lambda a, b: self.part.WirePolyLine(points=(a, b), meshable=ON)
        geom._collect_upper_bars(None, connect_only=True, connect=connect)
        geom._collect_lower_bars(None, connect_only=True, connect=connect)


def main(part_name, f1, f2, span, rise):
    GeodesicShellPart(part_name, f1, f2, span, rise).build()
    print("Completed!")


if __name__ == "__main__":
    main(part_name='duanchengxian_project01', f1=5, f2=2, span=20, rise=9)
