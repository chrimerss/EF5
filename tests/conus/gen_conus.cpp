// Synthetic CONUS domain for testing EF5 at full size: the MRMS 0.01 degree
// grid (7000 x 3500, about 9.5 million land cells at 1 km), with a
// hydrologically conditioned DEM, ESRI D8 flow directions, a self-inclusive
// flow accumulation, CRESTPHYS and kinematic wave parameter grids, monthly
// PET and 2-minute precipitation rates from moving storm cells.
//
// Usage: gen_conus <output dir> <precip steps> [columns rows]
//
// Writes <output dir>/{basic,params,pet,precip}/ and <output dir>/control.txt
// with every block except [Task] and [Execute]; one gauge sits at each outlet.
#include "Grid.h"
#include "TifGrid.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <queue>
#include <random>
#include <string>
#include <sys/stat.h>
#include <vector>

static long NC = 7000, NR = 3500;
static const double CS = 0.01, LEFT = -130.0, TOP = 55.0;
static const float ND = -9999.0f;

static FloatGrid *NewGrid(float fill) {
  FloatGrid *g = new FloatGrid();
  g->numCols = NC;
  g->numRows = NR;
  g->cellSize = CS;
  g->extent.left = LEFT;
  g->extent.top = TOP;
  g->extent.right = LEFT + NC * CS;
  g->extent.bottom = TOP - NR * CS;
  g->noData = ND;
  g->geoSet = false;
  g->data = new float *[NR];
  for (long r = 0; r < NR; r++) {
    g->data[r] = new float[NC];
    for (long c = 0; c < NC; c++) {
      g->data[r][c] = fill;
    }
  }
  return g;
}

// Smooth value noise: bilinear interpolation of a random lattice
struct ValueNoise {
  int n;
  std::vector<float> v;
  ValueNoise(int n_, unsigned seed) : n(n_), v((n_ + 2) * (n_ + 2)) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> u(-1, 1);
    for (size_t i = 0; i < v.size(); i++) {
      v[i] = u(rng);
    }
  }
  float at(double fx, double fy) const { // fx, fy in [0, 1]
    double x = fx * n, y = fy * n;
    int i = (int)x, j = (int)y;
    double tx = x - i, ty = y - j;
    tx = tx * tx * (3 - 2 * tx);
    ty = ty * ty * (3 - 2 * ty);
    int w = n + 2;
    return (1 - ty) * ((1 - tx) * v[j * w + i] + tx * v[j * w + i + 1]) +
           ty * ((1 - tx) * v[(j + 1) * w + i] + tx * v[(j + 1) * w + i + 1]);
  }
};

static const int DX[8] = {1, 1, 0, -1, -1, -1, 0, 1};
static const int DY[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static const int ESRI[8] = {1, 2, 4, 8, 16, 32, 64, 128};

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "Usage: %s <output dir> <precip steps> [columns rows]\n",
            argv[0]);
    return 1;
  }
  std::string out = argv[1];
  int nPrecip = atoi(argv[2]);
  if (argc > 4) {
    NC = atol(argv[3]);
    NR = atol(argv[4]);
  }
  const char *dirs[] = {"", "/basic", "/params", "/pet", "/precip"};
  for (int d = 0; d < 5; d++) {
    mkdir((out + dirs[d]).c_str(), 0755);
  }
  long N = NC * NR;

  // Land mask: an irregular CONUS-like blob in pseudo lon/lat, so the shape
  // does not depend on the grid size
  ValueNoise coast(24, 7);
  std::vector<unsigned char> land(N, 0);
  long nLand = 0;
  for (long r = 0; r < NR; r++) {
    double lat = 55.0 - 35.0 * (r + 0.5) / NR;
    for (long c = 0; c < NC; c++) {
      double lon = -130.0 + 70.0 * (c + 0.5) / NC;
      double ex = (lon + 96.0) / 26.5, ey = (lat - 37.5) / 11.4;
      double thr = 1.0 + 0.18 * coast.at((double)c / NC, (double)r / NR);
      if (ex * ex + ey * ey < thr) {
        land[r * NC + c] = 1;
        nLand++;
      }
    }
  }

  // Raw DEM: mountain ranges, a trunk valley, tilts to the east and south
  ValueNoise n1(12, 1), n2(60, 2), n3(300, 3), n4(1500, 4);
  std::vector<float> z(N, ND);
  for (long r = 0; r < NR; r++) {
    double lat = 55.0 - 35.0 * (r + 0.5) / NR;
    for (long c = 0; c < NC; c++) {
      if (!land[r * NC + c]) {
        continue;
      }
      double lon = -130.0 + 70.0 * (c + 0.5) / NC;
      double fx = (double)c / NC, fy = (double)r / NR;
      double rockies = 2500.0 * exp(-pow((lon + 110.0) / 7.0, 2));
      double apps = 700.0 * exp(-pow((lon + 80.0) / 3.0, 2)) *
                    exp(-pow((lat - 38.0) / 8.0, 2));
      double plains = -12.0 * (lon + 130.0) * (NC / 7000.0) + 900.0;
      double gulf = 20.0 * (lat - 29.0) * (NR / 3500.0);
      double valley = -250.0 * exp(-pow((lon + 91.0) / 2.5, 2));
      double noise = 300 * n1.at(fx, fy) + 120 * n2.at(fx, fy) +
                     40 * n3.at(fx, fy) + 8 * n4.at(fx, fy);
      z[r * NC + c] =
          std::max(1.0, rockies + apps + plains + gulf + valley + noise + 500);
    }
  }

  // Priority-flood fill with epsilon (Barnes et al., 2014): every land cell
  // drains to the edge of the mask, and the flood order gives D8 directions
  typedef std::pair<float, long> Cell;
  std::priority_queue<Cell, std::vector<Cell>, std::greater<Cell> > open;
  std::vector<unsigned char> closed(N, 0);
  std::vector<signed char> dir(N, -1);
  for (long r = 0; r < NR; r++) {
    for (long c = 0; c < NC; c++) {
      long i = r * NC + c;
      if (!land[i]) {
        continue;
      }
      for (int k = 0; k < 8; k++) {
        long rr = r + DY[k], cc = c + DX[k];
        if (rr < 0 || rr >= NR || cc < 0 || cc >= NC || !land[rr * NC + cc]) {
          closed[i] = 1;
          dir[i] = k; // Outlet, drains off the mask
          open.push(Cell(z[i], i));
          break;
        }
      }
    }
  }
  while (!open.empty()) {
    Cell e = open.top();
    open.pop();
    long i = e.second, r = i / NC, c = i % NC;
    for (int k = 0; k < 8; k++) {
      long rr = r + DY[k], cc = c + DX[k];
      if (rr < 0 || rr >= NR || cc < 0 || cc >= NC) {
        continue;
      }
      long j = rr * NC + cc;
      if (!land[j] || closed[j]) {
        continue;
      }
      closed[j] = 1;
      if (z[j] <= z[i]) {
        z[j] = nextafterf(z[i], 1e9f) + 1e-3f;
      }
      dir[j] = (k + 4) % 8; // j drains to i
      open.push(Cell(z[j], j));
    }
  }

  // Flow accumulation, counting each cell itself
  std::vector<long> down(N, -1);
  std::vector<int> pending(N, 0);
  std::vector<float> fam(N, ND);
  for (long i = 0; i < N; i++) {
    if (!land[i]) {
      continue;
    }
    long rr = i / NC + DY[(int)dir[i]], cc = i % NC + DX[(int)dir[i]];
    if (rr >= 0 && rr < NR && cc >= 0 && cc < NC && land[rr * NC + cc]) {
      down[i] = rr * NC + cc;
      pending[down[i]]++;
    }
    fam[i] = 1;
  }
  std::vector<long> ready;
  for (long i = 0; i < N; i++) {
    if (land[i] && pending[i] == 0) {
      ready.push_back(i);
    }
  }
  while (!ready.empty()) {
    long i = ready.back();
    ready.pop_back();
    if (down[i] >= 0) {
      fam[down[i]] += fam[i];
      if (--pending[down[i]] == 0) {
        ready.push_back(down[i]);
      }
    }
  }

  FloatGrid *g = NewGrid(ND);
  for (long i = 0; i < N; i++) {
    g->data[i / NC][i % NC] = land[i] ? z[i] : ND;
  }
  WriteFloatTifGrid((out + "/basic/dem.tif").c_str(), g);
  for (long i = 0; i < N; i++) {
    g->data[i / NC][i % NC] = land[i] ? (float)ESRI[(int)dir[i]] : ND;
  }
  WriteFloatTifGrid((out + "/basic/ddm.tif").c_str(), g);
  for (long i = 0; i < N; i++) {
    g->data[i / NC][i % NC] = fam[i];
  }
  WriteFloatTifGrid((out + "/basic/fam.tif").c_str(), g);

  // Spatially varying parameter grids; the scalars in the control file are 1.
  // CRESTPHYS reads im_grid as a fraction (only a scalar im is in percent).
  const char *pnames[] = {"wm",    "b",    "im",     "fc",    "hmaxaq",
                          "gwc",   "gwe",  "alpha",  "beta",  "alpha0",
                          "under", "leaki", "th"};
  const float pbase[] = {150, 0.5, 0.02, 10,  1000, 0.1,  2,
                         2.0, 0.6, 1.5, 1.0, 0.05, 50.0};
  ValueNoise pn(40, 11);
  for (int p = 0; p < 13; p++) {
    for (long r = 0; r < NR; r++) {
      for (long c = 0; c < NC; c++) {
        float m = 1.0f + 0.3f * pn.at((double)c / NC + 0.01 * p, (double)r / NR);
        g->data[r][c] = land[r * NC + c] ? pbase[p] * m : ND;
      }
    }
    WriteFloatTifGrid((out + "/params/" + pnames[p] + ".tif").c_str(), g);
  }

  // Monthly PET climatology (mm/h)
  for (long i = 0; i < N; i++) {
    g->data[i / NC][i % NC] = land[i] ? 0.15f : ND;
  }
  for (int m = 1; m <= 12; m++) {
    char name[512];
    snprintf(name, sizeof(name), "%s/pet/PET_%02d.tif", out.c_str(), m);
    WriteFloatTifGrid(name, g);
  }

  // 2-minute precipitation rate (mm/h) from storm cells drifting east, named
  // like MRMS PrecipRate, starting 2023-06-01 00:02
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> ux(0, NC), uy(0, NR), ur(30, 250),
      ua(5, 80);
  struct Storm {
    double x, y, rad, amp;
  };
  std::vector<Storm> storms(60);
  for (size_t s = 0; s < storms.size(); s++) {
    storms[s].x = ux(rng);
    storms[s].y = uy(rng);
    storms[s].rad = ur(rng);
    storms[s].amp = ua(rng);
  }
  for (int t = 0; t < nPrecip; t++) {
    for (long r = 0; r < NR; r++) {
      for (long c = 0; c < NC; c++) {
        g->data[r][c] = 0.0f;
      }
    }
    for (size_t s = 0; s < storms.size(); s++) {
      Storm &st = storms[s];
      long r0 = std::max(0L, (long)(st.y - 3 * st.rad));
      long r1 = std::min(NR, (long)(st.y + 3 * st.rad));
      long c0 = std::max(0L, (long)(st.x - 3 * st.rad));
      long c1 = std::min(NC, (long)(st.x + 3 * st.rad));
      for (long r = r0; r < r1; r++) {
        for (long c = c0; c < c1; c++) {
          double d2 = ((c - st.x) * (c - st.x) + (r - st.y) * (r - st.y)) /
                      (st.rad * st.rad);
          if (d2 < 9) {
            g->data[r][c] += st.amp * exp(-d2);
          }
        }
      }
      st.x += 3; // About 90 km/h at 1 km cells and 2-minute steps
    }
    int minutes = 2 * (t + 1);
    char name[512];
    snprintf(name, sizeof(name),
             "%s/precip/PrecipRate_202306%02d-%02d%02d00.tif", out.c_str(),
             1 + minutes / 1440, (minutes / 60) % 24, minutes % 60);
    WriteFloatTifGrid(name, g);
  }

  // Control file without [Task] and [Execute]
  const char *o = out.c_str();
  FILE *fp = fopen((out + "/control.txt").c_str(), "w");
  fprintf(fp,
          "[Basic]\nDEM=%s/basic/dem.tif\nDDM=%s/basic/ddm.tif\n"
          "FAM=%s/basic/fam.tif\nPROJ=geographic\nESRIDDM=true\nSELFFAM=true\n\n",
          o, o, o);
  fprintf(fp,
          "[PrecipForcing MRMS]\nTYPE=TIF\nUNIT=mm/h\nFREQ=2u\nLOC=%s/precip\n"
          "NAME=PrecipRate_YYYYMMDD-HHUU00.tif\n\n",
          o);
  fprintf(fp,
          "[PETForcing PET]\nTYPE=TIF\nUNIT=mm/h\nFREQ=1m\nLOC=%s/pet\n"
          "NAME=PET_MM.tif\n\n",
          o);
  long numGauges = 0;
  for (long i = 0; i < N; i++) {
    if (land[i] && down[i] < 0) {
      fprintf(fp, "[Gauge g%ld]\ncellx=%ld\ncelly=%ld\noutputts=false\n",
              numGauges++, i % NC, i / NC);
    }
  }
  fprintf(fp, "\n[Basin CONUS]\n");
  for (long k = 0; k < numGauges; k++) {
    fprintf(fp, "gauge=g%ld\n", k);
  }
  fprintf(fp, "\n[CrestphysParamSet P]\ngauge=g0\n");
  for (int p = 0; p < 7; p++) {
    fprintf(fp, "%s_grid=%s/params/%s.tif\n", pnames[p], o, pnames[p]);
  }
  fprintf(fp, "wm=1.0\nb=1.0\nim=1.0\nke=1.0\nfc=1.0\nksoil=0.001\niwu=50\n"
              "igw=25\nhmaxaq=1.0\ngwc=1.0\ngwe=1.0\n\n");
  fprintf(fp, "[KWParamSet R]\ngauge=g0\n");
  for (int p = 7; p < 13; p++) {
    fprintf(fp, "%s_grid=%s/params/%s.tif\n", pnames[p], o, pnames[p]);
  }
  fprintf(fp, "under=1.0\nleaki=1.0\nth=1.0\nisu=0\nalpha=1.0\nbeta=1.0\n"
              "alpha0=1.0\n\n");
  fclose(fp);

  printf("%ld x %ld grid, %ld land cells (%.1f%%), %ld outlets, %d precip "
         "steps\n",
         NC, NR, nLand, 100.0 * nLand / N, numGauges, nPrecip);
  return 0;
}
