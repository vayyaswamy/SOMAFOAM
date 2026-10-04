"""Readers for SOMAFOAM (foam-extend) cases written in ASCII.

Used by plot_1d.py, plot_2d.py and plot_electrodes.py; it can also be imported
directly:

    import somafoam_io as sio
    case = sio.Case("path/to/case")
    t = case.latest_time()
    mesh = case.mesh(t)
    ne = case.field("electron", t)          # numpy array, one value per cell

Only numpy is required. Binary cases are not supported: set
"writeFormat ascii;" in system/controlDict.
"""

import os
import re

import numpy as np

_NUMBER = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")
_LIST_START = re.compile(r"\n\s*(\d+)\s*\n?\s*\(")
_INTERNAL = re.compile(
    r"internalField\s+(uniform\s+(\([^)]*\)|[^;]+);|nonuniform\s+List<(\w+)>\s*(\d+)\s*\()"
)


def _read(path):
    with open(path, "rb") as f:
        raw = f.read()
    text = raw.decode("latin-1")
    if re.search(r"format\s+binary", text[:2000]):
        raise ValueError(
            f"{path} is binary; set 'writeFormat ascii;' in system/controlDict"
        )
    return text


def _numbers(text, count):
    """First `count` numbers of a parenthesis-free or parenthesised block."""
    return np.array(text.replace("(", " ").replace(")", " ").split()[:count],
                    dtype=float)


def read_internal_field(path, n_cells=None):
    """Internal field of a vol field file as an array of shape (n,) or (n, 3).

    A uniform field is expanded to n_cells values if n_cells is given, and
    returned as a single value (or 3-vector) otherwise.
    """
    text = _read(path)
    m = _INTERNAL.search(text)
    if not m:
        raise ValueError(f"no internalField in {path}")

    if m.group(1).startswith("uniform"):
        value = m.group(2).strip()
        if value.startswith("("):
            v = np.array(value.strip("()").split(), dtype=float)
            return np.tile(v, (n_cells, 1)) if n_cells else v
        v = float(value)
        return np.full(n_cells, v) if n_cells else v

    kind, n = m.group(3), int(m.group(4))
    body = text[m.end():]
    if kind == "scalar":
        return np.array(body[:body.index(")")].split()[:n], dtype=float)

    width = {"vector": 3, "symmTensor": 6, "tensor": 9}.get(kind)
    if width is None:
        raise ValueError(f"unsupported field type {kind} in {path}")
    end = body.index("boundaryField") if "boundaryField" in body else len(body)
    return _numbers(body[:end], n*width).reshape(n, width)


def read_scalar_list(path, dtype=float):
    """A plain list file (owner, neighbour, cellLevel, ...)."""
    text = _read(path)
    m = _LIST_START.search(text)
    if not m:
        raise ValueError(f"no list in {path}")
    n = int(m.group(1))
    body = text[m.end():]
    return np.array(body[:body.index(")")].split()[:n], dtype=dtype)


def read_points(path):
    text = _read(path)
    m = _LIST_START.search(text)
    n = int(m.group(1))
    return _numbers(text[m.end():], 3*n).reshape(n, 3)


def read_faces(path):
    """List of integer arrays, one per face."""
    text = _read(path)
    m = _LIST_START.search(text)
    n = int(m.group(1))
    faces = [
        np.array(points.split(), dtype=int)
        for points in re.findall(r"\d+\s*\(([^()]*)\)", text[m.end():])
    ]
    return faces[:n]


class Mesh:
    """Polyhedral mesh of one case (or one processor) at one time."""

    def __init__(self, poly_mesh_dir):
        self.directory = poly_mesh_dir
        self.points = read_points(os.path.join(poly_mesh_dir, "points"))
        self.faces = read_faces(os.path.join(poly_mesh_dir, "faces"))
        self.owner = read_scalar_list(
            os.path.join(poly_mesh_dir, "owner"), int)
        self.neighbour = read_scalar_list(
            os.path.join(poly_mesh_dir, "neighbour"), int)
        self.n_cells = int(self.owner.max()) + 1
        self._cell_points = None
        self._centres = None

    @property
    def cell_points(self):
        """Point labels of each cell (list of sorted integer arrays)."""
        if self._cell_points is None:
            sets = [set() for _ in range(self.n_cells)]
            for face, cell in zip(self.faces, self.owner):
                sets[cell].update(face.tolist())
            for face, cell in zip(self.faces, self.neighbour):
                sets[cell].update(face.tolist())
            self._cell_points = [np.array(sorted(s)) for s in sets]
        return self._cell_points

    @property
    def centres(self):
        """Cell centres, as the mean of the cell's points; shape (n, 3)."""
        if self._centres is None:
            self._centres = np.array(
                [self.points[p].mean(axis=0) for p in self.cell_points])
        return self._centres

    def extent(self):
        return self.points.max(axis=0) - self.points.min(axis=0)

    def axes_by_extent(self):
        """Coordinate axes ordered from the longest to the shortest extent
        of the cell centres (the one-cell-thick directions come last)."""
        spread = self.centres.max(axis=0) - self.centres.min(axis=0)
        return list(np.argsort(-spread))

    def polygons(self, axes):
        """Outline of each cell in the plane of the two given axes, for a 2D
        (one cell thick) mesh: list of (k, 2) arrays, ordered around the
        cell. Cells with hanging nodes give polygons with more than four
        corners."""
        normal = ({0, 1, 2} - set(axes)).pop()
        z = self.points[:, normal]
        front = z <= z.min() + 1e-6*max(z.max() - z.min(), 1e-300)
        polygons = []
        for labels, centre in zip(self.cell_points, self.centres):
            p = self.points[labels[front[labels]]][:, axes]
            c = centre[axes]
            angle = np.arctan2(p[:, 1] - c[1], p[:, 0] - c[0])
            polygons.append(p[np.argsort(angle)])
        return polygons


class Case:
    """A case directory, or one processor directory of a parallel case."""

    def __init__(self, path):
        self.path = os.path.abspath(path)
        if not os.path.isdir(os.path.join(self.path, "constant")):
            raise ValueError(f"{path} is not a case directory")
        self._meshes = {}

    # -- times ---------------------------------------------------------------

    def times(self):
        """Time directory names, sorted by value."""
        names = [
            d for d in os.listdir(self.path)
            if _NUMBER.match(d) and os.path.isdir(os.path.join(self.path, d))
        ]
        return sorted(names, key=float)

    def latest_time(self):
        return self.times()[-1]

    def time_name(self, value):
        """Name of the time directory closest to a value (or 'latest')."""
        names = self.times()
        if not names:
            raise ValueError(f"no time directories in {self.path}")
        if str(value) == "latest":
            return names[-1]
        if str(value) in names:
            return str(value)
        target = float(value)
        return min(names, key=lambda n: abs(float(n) - target))

    def times_with(self, field):
        return [t for t in self.times()
                if os.path.isfile(os.path.join(self.path, t, field))]

    # -- mesh and fields -----------------------------------------------------

    def poly_mesh_dir(self, time=None):
        """polyMesh directory in use at a time: that of the latest time not
        after it that has one (moving or refined meshes), else constant."""
        if time is not None:
            for name in reversed(self.times()):
                if float(name) <= float(time)*(1 + 1e-12):
                    d = os.path.join(self.path, name, "polyMesh")
                    if os.path.isfile(os.path.join(d, "points")):
                        return d
        return os.path.join(self.path, "constant", "polyMesh")

    def mesh(self, time=None):
        d = self.poly_mesh_dir(time)
        if d not in self._meshes:
            self._meshes[d] = Mesh(d)
        return self._meshes[d]

    def has_field(self, name, time):
        return os.path.isfile(os.path.join(self.path, str(time), name))

    def field(self, name, time):
        """Internal field at a time. 'cellLevel' is read from the mesh
        directory of an adaptively refined case."""
        mesh = self.mesh(time)
        if name == "cellLevel":
            path = os.path.join(self.poly_mesh_dir(time), "cellLevel")
            if not os.path.isfile(path):
                return np.zeros(mesh.n_cells)
            return read_scalar_list(path)
        return read_internal_field(
            os.path.join(self.path, str(time), name), mesh.n_cells)

    def scalar(self, name, time, component=None):
        """A field as one value per cell: a component (0, 1, 2) or, by
        default, the magnitude of a vector field."""
        f = self.field(name, time)
        if f.ndim == 1:
            return f
        if component is None:
            return np.linalg.norm(f, axis=1)
        return f[:, component]

    # -- parallel cases ------------------------------------------------------

    def processors(self):
        """Cases of the processor directories, if any."""
        cases = []
        i = 0
        while os.path.isdir(os.path.join(self.path, f"processor{i}")):
            cases.append(Case(os.path.join(self.path, f"processor{i}")))
            i += 1
        return cases

    def parts(self, field, time):
        """The cases that hold a field at a time: this case if it has been
        reconstructed (or was run in serial), else its processors."""
        here = os.path.isdir(os.path.join(self.path, str(time)))
        if self.has_field(field, time) or (field == "cellLevel" and here):
            return [self]
        procs = [
            p for p in self.processors()
            if p.has_field(field, time) or (
                field == "cellLevel"
                and os.path.isdir(os.path.join(p.path, str(time))))
        ]
        if not procs:
            raise FileNotFoundError(
                f"{field} not found at time {time} in {self.path} "
                "or its processor directories")
        return procs


_UNITS = {
    "Te": "K", "T": "K", "Tion": "K", "Phi": "V", "E": "V/m", "p": "Pa",
    "cellLevel": "",
}


def label(field, component=None):
    """Axis label with units for the usual SOMAFOAM fields. Files named
    after a species hold its number density."""
    base = field[:-4] if field.endswith("Mean") else field
    if base.startswith("N_"):
        base = base[2:]
    if base in _UNITS:
        unit = _UNITS[base]
    elif base.startswith(("J_", "F_", "U_", "Sy_", "D_", "mu_")) or \
            base in ("Jnet", "Jtot", "refinementIndicator"):
        unit = None
    else:
        unit = "m$^{-3}$"
    text = field
    if component is not None:
        text += "_" + "xyz"[component]
    if field.endswith("Mean"):
        text += " (average)"
    return f"{text} ({unit})" if unit else text


def all_times(case):
    """Time names of a case and of its processor directories, sorted."""
    names = set(case.times())
    for p in case.processors():
        names.update(p.times())
    return sorted(names, key=float)


def nearest_time(names, value):
    """Entry of a sorted list of time names for a value or 'latest'."""
    if not names:
        raise ValueError("no time directories found")
    if str(value) == "latest":
        return names[-1]
    if str(value) in names:
        return str(value)
    return min(names, key=lambda n: abs(float(n) - float(value)))


def read_electrode_file(path):
    """Columns of an electrodeVoltageCurrent output file as a dict of
    arrays: t, V, I_conduction, I_displacement, I_total, j_total."""
    data = np.loadtxt(path, comments="#", ndmin=2)
    names = ["t", "V", "I_conduction", "I_displacement", "I_total", "j_total"]
    return {n: data[:, i] for i, n in enumerate(names)}


def read_electrode(case_path, patch, name="electrodes"):
    """Voltage and current history of a patch, joined over all the start
    times of the run (a restart writes into a new directory)."""
    root = os.path.join(case_path, name)
    if not os.path.isdir(root):
        raise FileNotFoundError(
            f"{root} not found: add the electrodeVoltageCurrent function "
            "object to system/controlDict")
    starts = sorted((d for d in os.listdir(root) if _NUMBER.match(d)),
                    key=float)
    pieces = []
    for i, start in enumerate(starts):
        f = os.path.join(root, start, patch + ".dat")
        if not os.path.isfile(f):
            continue
        d = read_electrode_file(f)
        if i + 1 < len(starts):
            keep = d["t"] < float(starts[i + 1])
            d = {k: v[keep] for k, v in d.items()}
        pieces.append(d)
    if not pieces:
        raise FileNotFoundError(f"no {patch}.dat under {root}")
    return {k: np.concatenate([p[k] for p in pieces]) for k in pieces[0]}
