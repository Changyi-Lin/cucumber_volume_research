#!/usr/bin/env python3
"""Print connected-component size, volume, and bounds for CGAL triangle PLY files."""

from __future__ import annotations

import argparse
from collections import defaultdict, deque
from pathlib import Path
import struct
import numpy as np


def read_cgal_ply(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open("rb") as handle:
        header = []
        while True:
            line = handle.readline().decode("ascii").strip()
            header.append(line)
            if line == "end_header":
                break
        if "format binary_little_endian 1.0" not in header:
            raise ValueError("Only CGAL binary_little_endian PLY is supported")
        vertex_count = int(next(x.split()[2] for x in header if x.startswith("element vertex ")))
        face_count = int(next(x.split()[2] for x in header if x.startswith("element face ")))
        vertices = np.fromfile(handle, dtype="<f8", count=vertex_count * 3).reshape(-1, 3)
        faces = np.empty((face_count, 3), dtype=np.int32)
        for i in range(face_count):
            count = struct.unpack("<B", handle.read(1))[0]
            if count != 3:
                raise ValueError(f"Non-triangle face with {count} vertices")
            faces[i] = struct.unpack("<iii", handle.read(12))
    return vertices, faces


def components(faces: np.ndarray) -> list[list[int]]:
    edges: dict[tuple[int, int], list[int]] = defaultdict(list)
    for fi, (a, b, c) in enumerate(faces):
        for u, v in ((a, b), (b, c), (c, a)):
            edges[min(u, v), max(u, v)].append(fi)
    adjacency = [[] for _ in faces]
    for linked in edges.values():
        for a in linked:
            adjacency[a].extend(b for b in linked if b != a)
    seen = np.zeros(len(faces), dtype=bool)
    output = []
    for start in range(len(faces)):
        if seen[start]:
            continue
        seen[start] = True
        queue = deque([start])
        component = []
        while queue:
            face = queue.popleft()
            component.append(face)
            for neighbor in adjacency[face]:
                if not seen[neighbor]:
                    seen[neighbor] = True
                    queue.append(neighbor)
        output.append(component)
    return sorted(output, key=len, reverse=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("meshes", nargs="+", type=Path)
    args = parser.parse_args()
    for path in args.meshes:
        vertices, faces = read_cgal_ply(path)
        groups = components(faces)
        print(f"\n{path}: vertices={len(vertices)} faces={len(faces)} components={len(groups)}")
        for rank, group in enumerate(groups, 1):
            tri = faces[np.asarray(group)]
            used = np.unique(tri)
            a, b, c = vertices[tri[:, 0]], vertices[tri[:, 1]], vertices[tri[:, 2]]
            signed_volume = np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0
            lo, hi = vertices[used].min(axis=0), vertices[used].max(axis=0)
            print(f"  #{rank}: faces={len(group):6d} ({len(group)/len(faces):7.3%}) "
                  f"vertices={len(used):6d} signed_volume_ml={signed_volume/1000:10.4f} "
                  f"AABB=({', '.join(f'{x:.2f}' for x in (hi-lo))})")


if __name__ == "__main__":
    main()

