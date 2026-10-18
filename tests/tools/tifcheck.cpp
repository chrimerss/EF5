// Checks EF5 output grids: every cell that is not no-data must be finite and
// not negative. With --nonzero, each grid must also have a value above zero.
//
// Usage: tifcheck [--nonzero] <grid.tif>...
// Exits with 1 if any grid fails or cannot be read.
#include "TifGrid.h"
#include <cmath>
#include <cstdio>
#include <cstring>

int main(int argc, char **argv) {
  bool wantNonzero = false;
  int failed = 0, total = 0;
  for (int a = 1; a < argc; a++) {
    if (!strcmp(argv[a], "--nonzero")) {
      wantNonzero = true;
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
    bool ok = valid > 0 && bad == 0 && (!wantNonzero || positive > 0);
    printf("%s %s: %ld cells, %ld above zero, %ld NaN/inf/negative, max %g\n",
           ok ? "ok  " : "FAIL", argv[a], valid, positive, bad, maxValue);
    failed += !ok;
    delete grid;
  }
  printf("%d of %d grids failed\n", failed, total);
  return failed ? 1 : 0;
}
