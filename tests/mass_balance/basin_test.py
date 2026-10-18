#!/usr/bin/env python3
"""Mass-balance test for EF5's kinematic wave routing.

Builds a 100 x 100 cell basin of 0.01 degree cells near 35.5 N with a single
outlet: hillslopes drain sideways into one trunk channel. It rains 10 mm in
an hour, EF5 (CREST + kinematic wave) runs until the basin has drained, and
the outlet discharge volume must match the rain within 1%. The cases cover
subsurface flow into the channels (all rain becomes interflow) and overland
flow (fully impervious), north-south and east-west trunks, and 2-, 15- and
60-minute steps.

Usage: basin_test.py <ef5 binary> <work dir>
"""
import csv
import math
import os
import subprocess
import sys

N = 100  # cells per side
MID = N // 2  # trunk row or column
RAIN_MM = 10.0

CASES = [
    # name, trunk, pathway, time step, end of run
    ("interflow N-S 2 min", "ns", "interflow", "2u", "202308010000"),
    ("interflow E-W 60 min", "ew", "interflow", "60u", "202308010000"),
    ("overland N-S 15 min", "ns", "overland", "15u", "202307010000"),
    ("overland E-W 60 min", "ew", "overland", "60u", "202307010000"),
]


def write_asc(path, rows, fmt):
    with open(path, "w") as f:
        f.write(f"ncols {N}\nnrows {N}\nxllcorner -98.0\nyllcorner 35.0\n"
                f"cellsize 0.01\nNODATA_value -9999\n")
        for row in rows:
            f.write(" ".join(fmt % v for v in row) + "\n")


def build_basin(folder, trunk):
    """Fishbone basin: the trunk is column MID flowing south (ns) or row MID
    flowing east (ew), and every other cell flows straight towards it."""
    os.makedirs(f"{folder}/precip", exist_ok=True)
    os.makedirs(f"{folder}/pet", exist_ok=True)
    ddm, dem, fam = [], [], []
    for y in range(N):
        ddm_row, dem_row, fam_row = [], [], []
        for x in range(N):
            across, along = (x, y) if trunk == "ns" else (y, x)
            if across < MID:
                d, a = (1 if trunk == "ns" else 4), across + 1  # E or S
            elif across > MID:
                d, a = (16 if trunk == "ns" else 64), N - across  # W or N
            else:
                d, a = (4 if trunk == "ns" else 1), (along + 1) * N  # S or E
            ddm_row.append(d)
            fam_row.append(a)
            dem_row.append(10 + 2 * abs(across - MID) + 2 * (N - along))
        ddm.append(ddm_row)
        dem.append(dem_row)
        fam.append(fam_row)
    write_asc(f"{folder}/dem.asc", dem, "%d")
    write_asc(f"{folder}/ddm.asc", ddm, "%d")
    write_asc(f"{folder}/fam.asc", fam, "%d")
    write_asc(f"{folder}/precip/P_2023060101.asc",
              [[RAIN_MM] * N for _ in range(N)], "%g")
    for month in (6, 7, 8):
        write_asc(f"{folder}/pet/PET_{month:02d}.asc",
                  [[0] * N for _ in range(N)], "%d")
    return (MID, N - 1) if trunk == "ns" else (N - 1, MID)


def basin_area_m2():
    """Total cell area as EF5 computes it: GeographicProjection uses 110574 m
    per degree, and the latitude of each row's northern edge."""
    side = 110574 * 0.01
    return sum(side * side * math.cos((36.0 - 0.01 * y) / 180 * 3.1415926)
               for y in range(N)) * N


def run_case(ef5, work, name, trunk, pathway, step, end):
    folder = os.path.join(work, name.replace(" ", "_"))
    outlet_x, outlet_y = build_basin(folder, trunk)
    # Interflow: saturated thin soil with a huge Ksat sends all rain to
    # interflow. Overland: all impervious, with quick overland flow.
    im = 100 if pathway == "overland" else 0
    alpha0 = 0.1 if pathway == "overland" else 1.5
    control = f"""[Basic]
DEM={folder}/dem.asc
DDM={folder}/ddm.asc
FAM={folder}/fam.asc
PROJ=geographic
ESRIDDM=true
SELFFAM=true
[PrecipForcing P]
TYPE=ASC
UNIT=mm/h
FREQ=1h
LOC={folder}/precip
NAME=P_YYYYMMDDHH.asc
[PETForcing PET]
TYPE=ASC
UNIT=mm/h
FREQ=1m
LOC={folder}/pet
NAME=PET_MM.asc
[Gauge out]
cellx={outlet_x}
celly={outlet_y}
outputts=true
[Basin B]
gauge=out
[CrestParamSet PC]
gauge=out
wm=10
b=0.5
im={im}
ke=0
fc=1000
iwu=100
[KWParamSet R]
gauge=out
under=10
leaki=0.1
th=150
isu=0
alpha=2.0
beta=0.6
alpha0={alpha0}
[Task T]
STYLE=SIMU
MODEL=CREST
ROUTING=KW
BASIN=B
PRECIP=P
PET=PET
OUTPUT={folder}/
PARAM_SET=PC
ROUTING_PARAM_SET=R
TIMESTEP={step}
TIME_BEGIN=202306010000
TIME_END={end}
[Execute]
TASK=T
"""
    with open(f"{folder}/control.txt", "w") as f:
        f.write(control)
    with open(f"{folder}/ef5.log", "w") as log:
        result = subprocess.run([ef5, f"{folder}/control.txt"], stdout=log,
                                stderr=subprocess.STDOUT)
    with open(f"{folder}/ef5.log", errors="replace") as log:
        errors = [line for line in log if "ERROR:" in line]
    if result.returncode != 0 or errors:
        print(f"FAIL {name}: EF5 exited with {result.returncode}")
        print("".join(errors[:10]))
        return False

    seconds = int(step[:-1]) * 60
    with open(f"{folder}/ts.out.crest.csv") as f:
        rows = list(csv.reader(f))[1:]
    discharge = [float(row[1]) for row in rows]
    volume_out = sum(discharge) * seconds
    volume_in = RAIN_MM / 1000.0 * basin_area_m2()
    ratio = volume_out / volume_in
    ok = abs(ratio - 1.0) <= 0.01 and discharge[-1] < 0.01 * max(discharge)
    print(f"{'ok  ' if ok else 'FAIL'} {name}: outlet/rain = {ratio:.4f}, "
          f"peak {max(discharge):.1f} m3/s, final {discharge[-1]:.2f} m3/s")
    return ok


def main():
    ef5 = os.path.realpath(sys.argv[1])
    work = os.path.realpath(sys.argv[2])
    os.makedirs(work, exist_ok=True)
    failed = sum(not run_case(ef5, work, *case) for case in CASES)
    print(f"{failed} of {len(CASES)} mass-balance cases failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
