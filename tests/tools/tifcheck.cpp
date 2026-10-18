// Checks EF5 output grids: every cell that is not no-data must be finite and
// not negative. With --nonzero, each grid must also have a value above zero.
// With --reference, each grid must have as many cells with data as the
// reference grid; this also catches unreadable rows, which the reader turns
// into no-data.
//
// Usage: tifcheck [--nonzero] [--reference <grid.tif>] <grid.tif>...
// Exits with 1 if any grid fails or cannot be read.
#include "TifGrid.h"
#include <cmath>
#include <cstdio>
#include <cstring>

static long CountCells(FloatGrid *grid) {
  long cells = 0;
  for (long r = 0; r < grid->numRows; r++) {
    for (long c = 0; c < grid->numCols; c++) {
      cells += (grid->data[r][c] != grid->noData);
    }
  }
  return cells;
}

int main(int argc, char **argv) {
  bool wantNonzero = false;
  long wantCells = -1;
  int failed = 0, total = 0;
  for (int a = 1; a < argc; a++) {
    if (!strcmp(argv[a], "--nonzero")) {
      wantNonzero = true;
      continue;
    }
    if (!strcmp(argv[a], "--reference") && a + 1 < argc) {
      FloatGrid *reference = ReadFloatTifGrid(argv[++a]);
      if (!reference) {
        printf("FAIL %s: cannot read reference\n", argv[a]);
        return 1;
      }
      wantCells = CountCells(reference);
      delete reference;
      continue;
    }
    total++;
    FloatGrid *grid = ReadFloatTifGrid(argv[a]);
    if (!grid) {
      printf("FAIL %s: cannot read\n", argv[a]);
      failed++;
      continue;
    }
    long valid = 0, positive = 0, bad = 0;
    double maxValue = 0.0;
    for (long r = 0; r < grid->numRows; r++) {
      for (long c = 0; c < grid->numCols; c++) {
        float value = grid->data[r][c];
        if (value == grid->noData) {
          continue;
        }
        valid++;
        if (!std::isfinite(value) || value < 0.0f) {
          bad++;
        } else if (value > 0.0f) {
          positive++;
          if (value > maxValue) {
            maxValue = value;
          }
        }
      }
    }
    bool ok = valid > 0 && bad == 0 && (!wantNonzero || positive > 0) &&
              (wantCells < 0 || valid == wantCells);
    printf("%s %s: %ld cells", ok ? "ok  " : "FAIL", argv[a], valid);
    if (wantCells >= 0 && valid != wantCells) {
      printf(" (reference has %ld)", wantCells);
    }
    printf(", %ld above zero, %ld NaN/inf/negative, max %g\n", positive, bad,
           maxValue);
    failed += !ok;
    delete grid;
  }
  printf("%d of %d grids failed\n", failed, total);
  return failed ? 1 : 0;
}
