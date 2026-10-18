#!/bin/bash
# Runs EF5 (CRESTPHYS + kinematic wave) on the synthetic 1 km CONUS domain at
# 2-minute steps, writing soil moisture and unit discharge every step, once on
# 4 threads and once on 1, and fails unless:
#   - EF5 logs no errors, finds every forcing file, and writes every expected
#     grid, state and results.json
#   - every grid covers all land cells with no NaN, infinite or negative
#     values, and the final discharge state is not all zero (it is unrounded;
#     the unit discharge grids are rounded to 0.1 m^3/s/km^2)
#   - the two runs write byte-identical grids
#
# Usage: tests/run_conus.sh <ef5 binary> <work dir> [steps] [columns rows]
# Run ./configure first, since the test tools include config.h.
set -euo pipefail

EF5=$(realpath "$1")
WORK=$(mkdir -p "$2" && realpath "$2")
STEPS=${3:-15}
SIZE=${4:+$4 $5}
TESTS=$(cd "$(dirname "$0")" && pwd)
ROOT=$(dirname "$TESTS")
CXXFLAGS="-O2 -std=gnu++11 -I$ROOT/src -I$ROOT -I/usr/include/geotiff"
LIBS="-ltiff -lgeotiff -lz"

# Wall time and peak memory, when GNU time is installed
timed() { # label command...
  local label=$1; shift
  if [ -x /usr/bin/time ]; then
    /usr/bin/time -f "$label: %e s, %M KB peak" "$@"
  else
    "$@"
  fi
}

g++ $CXXFLAGS -o "$WORK/gen_conus" "$TESTS/conus/gen_conus.cpp" "$ROOT/src/TifGrid.cpp" $LIBS
g++ $CXXFLAGS -o "$WORK/tifcheck" "$TESTS/tools/tifcheck.cpp" "$ROOT/src/TifGrid.cpp" $LIBS

echo "== Generating the domain"
timed generate "$WORK/gen_conus" "$WORK/domain" "$STEPS" $SIZE 2>&1 | grep -v "legacy Deflate"

# The model clock runs from 2023-06-01 00:00 in 2-minute steps
stamp() { date -u -d "2023-06-01 00:00 UTC + $(( 2 * $1 )) minutes" +"$2"; }
END=$(stamp "$STEPS" %Y%m%d%H%M)

run() { # threads
  local out="$WORK/out_t$1"
  rm -rf "$out" && mkdir -p "$out"
  cp "$WORK/domain/control.txt" "$out/control.txt"
  cat >> "$out/control.txt" <<EOF
[Task CONUS]
STYLE=SIMU
MODEL=CRESTPHYS
ROUTING=KW
BASIN=CONUS
PRECIP=MRMS
PET=PET
OUTPUT=$out/
STATES=$out/
TIME_STATE=$END
PARAM_SET=P
ROUTING_PARAM_SET=R
DEFAULTPARAMSGAUGE=g0
TIMESTEP=2u
TIME_BEGIN=202306010000
TIME_END=$END
OUTPUT_GRIDS=soilmoisture|unitstreamflow
[Execute]
TASK=CONUS
EOF
  echo "== Running EF5 on $1 thread(s)"
  timed "ef5 on $1 thread(s)" env OMP_NUM_THREADS="$1" "$EF5" "$out/control.txt" \
    > "$out/ef5.log" 2>&1 || {
    tail -50 "$out/ef5.log"; echo "EF5 exited with an error"; exit 1; }
  grep -a "KB peak" "$out/ef5.log" || true
  grep -a "routing sec" "$out/ef5.log" | tail -3
  if grep -aq "ERROR:" "$out/ef5.log"; then
    grep -a "ERROR:" "$out/ef5.log" | head -20; echo "EF5 logged errors"; exit 1
  fi
  # EF5 runs on with zeros when a forcing file is missing
  if grep -aq "Missing .* file" "$out/ef5.log"; then
    grep -ao "Missing [^)]*)" "$out/ef5.log" | head -5
    echo "EF5 could not find forcing files"; exit 1
  fi
}

check() { # threads
  local out="$WORK/out_t$1" missing=0 step
  local expected=("$out/results.json")
  for step in $(seq 1 "$STEPS"); do
    expected+=("$out/sm.$(stamp "$step" %Y%m%d_%H%M).crestphys.tif"
               "$out/unitq.$(stamp "$step" %Y%m%d_%H%M).crestphys.tif")
  done
  local when=$(stamp "$STEPS" %Y%m%d_%H%M)
  for state in crestphys_SM crestphys_GW kwr_pCQ kwr_pOQ kwr_IR; do
    expected+=("$out/${state}_$when.tif")
  done
  for f in "${expected[@]}"; do
    [ -s "$f" ] || { echo "Missing output $f"; missing=1; }
  done
  [ $missing -eq 0 ] || exit 1
  echo "All ${#expected[@]} expected outputs written on $1 thread(s)"
  "$WORK/tifcheck" --reference "$WORK/domain/basic/dem.tif" "$out"/*.tif | grep -v "^ok"
  "$WORK/tifcheck" --nonzero "$out/kwr_pCQ_$when.tif"
}

run 4
check 4
run 1
check 1

echo "== Comparing 1 and 4 threads"
differ=0
for f in "$WORK/out_t4"/*.tif; do
  cmp -s "$f" "$WORK/out_t1/$(basename "$f")" || { echo "Differs: $(basename "$f")"; differ=1; }
done
[ $differ -eq 0 ] || { echo "Outputs depend on the thread count"; exit 1; }
echo "All $(ls "$WORK/out_t4"/*.tif | wc -l) grids identical on 1 and 4 threads"
