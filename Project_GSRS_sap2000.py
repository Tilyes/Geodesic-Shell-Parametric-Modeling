# -*- coding: utf-8 -*-
"""Geodesic lattice shell generator for SAP2000 OAPI.

Mirrors the Abaqus version geometry-for-geometry; differences are only in
how frames/nodes are materialised (SAP2000 FrameObj.AddByCoord vs. Abaqus
WirePolyLine) and in SAP2000 application bootstrap.

Units: N, mm, C  (SAP2000 eUnits = 9)
"""
from __future__ import division, print_function
import math
import os
import sys

import comtypes.client


PENTAGON = 5
TWO_PI = 2.0 * math.pi
UNITS_N_MM_C = 9


# --- vector helpers -------------------------------------------------------

def vec_add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])

def vec_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

def vec_scale(a, k):
    return (a[0] * k, a[1] * k, a[2] * k)

def project_on_sphere(p, radius):
    length = math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2])
    if length == 0:
        return (0.0, 0.0, 0.0)
    k = radius / length
    return (p[0] * k, p[1] * k, p[2] * k)

def _interpolate(start, end, steps):
    step = vec_scale(vec_sub(end, start), 1.0 / steps)
    return [vec_add(start, vec_scale(step, j)) for j in range(steps)]

def _variance(values):
    n = len(values)
    if n == 0:
        return 0.0
    mean = sum(values) / n
    return sum([(v - mean) ** 2 for v in values]) / n

# --- geometry (pure math, identical to the Abaqus version) ---------------

class ShellGeometry(object):
    def __init__(self, f1, f2, span, rise, adjustment):
        self.f1 = f1
        self.f2 = f2
        self.radius = rise / 2.0 + (span * span / 8.0) / rise
        self.z_angle = math.acos((self.radius - rise) / self.radius)
        base_split = f1 / float(f1 + f2)
        self.z_angle_upper = (base_split + adjustment) * self.z_angle
        self.upper_layers = self._build_upper()
        self.lower_layers = self._build_lower()

    def _build_upper(self):
        f1, radius = self.f1, self.radius
        apex = (0.0, 0.0, 1.0)
        base_ring = [
            (math.sin(self.z_angle_upper) * math.cos(TWO_PI / PENTAGON * k),
             math.sin(self.z_angle_upper) * math.sin(TWO_PI / PENTAGON * k),
             math.cos(self.z_angle_upper))
            for k in range(PENTAGON)
        ]
        edges = [
            [vec_add(apex, vec_scale(vec_sub(base, apex), m / float(f1)))
             for m in range(f1 + 1)]
            for base in base_ring
        ]
        layers = [[(0.0, 0.0, radius)]]
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
        f1, f2, radius = self.f1, self.f2, self.radius
        pentagon_step = TWO_PI / PENTAGON
        offset = f2 / float(f1) * math.pi / PENTAGON
        aux_offset = (f1 - f2) / float(f1) * pentagon_step

        ring_top = [
            vec_scale(self.upper_layers[f1][i], 1.0 / radius)
            for i in range(PENTAGON * f1) if i % f1 == 0
        ]
        ring_bot = [
            (math.sin(self.z_angle) * math.cos(offset + i * pentagon_step),
             math.sin(self.z_angle) * math.sin(offset + i * pentagon_step),
             math.cos(self.z_angle))
            for i in range(PENTAGON)
        ]
        ring_aux = [
            (math.sin(self.z_angle) * math.cos(offset + aux_offset + i * pentagon_step),
             math.sin(self.z_angle) * math.sin(offset + aux_offset + i * pentagon_step),
             math.cos(self.z_angle))
            for i in range(PENTAGON)
        ]

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
                ring_bot[n] = ring_aux[n]
            if m == PENTAGON:
                m = 0

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
        layers[0] = list(self.upper_layers[f1])
        return layers

    def iter_bars(self):
        """Yield (p1, p2) for every bar — upper cap then lower band."""
        f1, f2 = self.f1, self.f2
        upper = self.upper_layers
        for m in range(1, f1 + 1):
            ring_len = PENTAGON * m
            for q in range(ring_len):
                yield upper[m][q], upper[m][(q + 1) % ring_len]
        for m in range(f1, 0, -1):
            ring_len = PENTAGON * m
            prev_len = PENTAGON * (m - 1)
            non_corner = []
            for q in range(ring_len):
                if q % m == 0:
                    prev_idx = (q // m) * (m - 1) if prev_len > 0 else 0
                    yield upper[m][q], upper[m - 1][prev_idx]
                else:
                    non_corner.append(q)
            for idx, q in enumerate(non_corner):
                p = upper[m][q]
                if prev_len > 0:
                    yield p, upper[m - 1][idx]
                    yield p, upper[m - 1][(idx + 1) % prev_len]
                else:
                    yield p, upper[m - 1][0]
                    yield p, upper[m - 1][0]

        lower = self.lower_layers
        ring_len = PENTAGON * f1
        for m in range(1, f2 + 1):
            for q in range(ring_len):
                yield lower[m][q], lower[m][(q + 1) % ring_len]
        for m in range(f2):
            for q in range(ring_len):
                q_prev = q - 1 if q > 0 else ring_len - 1
                yield lower[m][q], lower[m + 1][q]
                yield lower[m][q], lower[m + 1][q_prev]

    def bar_lengths(self):
        return [math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)
                for a, b in self.iter_bars()]


def search_best_adjustment(f1, f2, span, rise, candidates=None):
    if candidates is None:
        candidates = [i / 100.0 for i in range(-2, 13)]
    best_adj, best_var = 0.0, float("inf")
    for adj in candidates:
        var = _variance(ShellGeometry(f1, f2, span, rise, adj).bar_lengths())
        if var < best_var:
            best_var, best_adj = var, adj
    return best_adj


# --- SAP2000 bootstrap ----------------------------------------------------

def start_sap2000(visible=True):
    """Start a new SAP2000 instance via the OAPI helper and return (SapObject, SapModel)."""
    helper = comtypes.client.CreateObject('SAP2000v1.Helper')
    helper = helper.QueryInterface(comtypes.gen.SAP2000v1.cHelper)
    sap_object = helper.CreateObjectProgID('CSI.SAP2000.API.SapObject')
    sap_object.ApplicationStart()
    if visible:
        sap_object.Hide()   # SAP starts hidden by default for OAPI; force show below
        sap_object.Unhide()
    sap_model = sap_object.SapModel
    sap_model.InitializeNewModel(UNITS_N_MM_C)
    sap_model.File.NewBlank()
    return sap_object, sap_model


# --- builder --------------------------------------------------------------

class GeodesicShellPart(object):
    """Build the geodesic shell as SAP2000 frame elements."""

    def __init__(self, sap_model, f1, f2, span, rise, frame_prefix='B'):
        adjustment = search_best_adjustment(f1, f2, span, rise)
        self.geometry = ShellGeometry(f1, f2, span, rise, adjustment)
        self.sap_model = sap_model
        self.frame_prefix = frame_prefix

    def build(self):
        frame = self.sap_model.FrameObj
        created = 0
        for p1, p2 in self.geometry.iter_bars():
            # AddByCoord signature: (x1,y1,z1,x2,y2,z2, Name, PropName, UserName, CSys)
            ret = frame.AddByCoord(
                p1[0], p1[1], p1[2],
                p2[0], p2[1], p2[2],
                '',          # Name — let SAP auto-assign
                'Default',   # Section property
                '',          # UserName
                'Global')
            # AddByCoord returns [Name, ret] via comtypes by-ref outputs; older
            # builds may return just the int ret code.
            code = ret[-1] if isinstance(ret, (tuple, list)) else ret
            if code != 0:
                raise RuntimeError('FrameObj.AddByCoord failed (code={0})'.format(code))
            created += 1
        return created


def main(f1=5, f2=2, span=20.0, rise=9.0,
         save_path=r'D:\Personal\Show\GSRS\geodesic_shell.sdb'):
    sap_object, sap_model = start_sap2000(visible=True)
    print('SAP2000 OAPI version: {0}'.format(sap_object.GetOAPIVersionNumber()))
    try:
        builder = GeodesicShellPart(sap_model, f1, f2, span, rise)
        n = builder.build()
        print('Created {0} frame elements.'.format(n))
        if save_path:
            parent = os.path.dirname(save_path)
            if parent and not os.path.isdir(parent):
                os.makedirs(parent)
            ret = sap_model.File.Save(save_path)
            if ret != 0:
                raise RuntimeError('File.Save failed (code={0})'.format(ret))
            print('Saved to {0}'.format(save_path))
        print('Completed! SAP2000 window left open.')
    except Exception:
        # Leave SAP open so the partial state is inspectable.
        raise


if __name__ == '__main__':
    main(f1=5, f2=2, span=20.0, rise=9.0)
