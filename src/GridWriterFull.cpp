#include "GridWriterFull.h"
#include "BasicConfigSection.h"
#include <climits>
#ifndef _WIN32
#include <system_error>
#endif

extern LongGrid *g_DEM;

// Grids allowed to wait for the writer before WriteGrid blocks
static const size_t maxQueuedGrids = 8;

GridWriterFull::GridWriterFull() {
  indexedNodes = NULL;
#ifndef _WIN32
  busy = 0;
  stopping = false;
#endif
}

GridWriterFull::GridWriterFull(const GridWriterFull &other)
    : GridWriterFull() {}

GridWriterFull::~GridWriterFull() {
#ifndef _WIN32
  if (writer.joinable()) {
    {
      std::lock_guard<std::mutex> guard(lock);
      stopping = true;
    }
    wake.notify_one();
    writer.join(); // The writer drains the queue before it stops
  }
#endif
}

void GridWriterFull::Initialize() {
  if (grid.data) {
    return; // Already set up
  }

  // Initialize everything in the new grid

  grid.extent.left = g_DEM->extent.left;
  grid.extent.bottom = g_DEM->extent.bottom;
  grid.extent.right = g_DEM->extent.right;
  grid.extent.top = g_DEM->extent.top;
  grid.numCols = g_DEM->numCols;
  grid.numRows = g_DEM->numRows;
  grid.cellSize = g_DEM->cellSize;
  // 2019-04: output gridded surface runoff ---------------------------------
  //grid.noData = -9999.0f;
  grid.noData = -9.0f;
  // --------------------------------- 
  grid.geoSet = g_DEM->geoSet;
  grid.modelType = g_DEM->modelType;
  grid.geographicType = g_DEM->geographicType;
  grid.geodeticDatum = g_DEM->geodeticDatum;

  grid.data = new float *[grid.numRows];
  for (long i = 0; i < grid.numRows; i++) {
    grid.data[i] = new float[grid.numCols];
  }

  // Grid is setup! Copy over no data values everywhere
  for (long row = 0; row < grid.numRows; row++) {
    for (long col = 0; col < grid.numCols; col++) {
      grid.data[row][col] = grid.noData;
    }
  }

#ifndef _WIN32
  try {
    writer = std::thread(&GridWriterFull::WriterLoop, this);
  } catch (const std::system_error &) {
    // No threads available, WriteGrid writes synchronously instead
  }
#endif
}

void GridWriterFull::WriteGrid(std::vector<GridNode> *nodes,
                               std::vector<float> *data, const char *file,
                               bool ascii) {
  Job *job = new Job;
  job->values = *data;
  job->file = file;
  job->ascii = ascii;
  Submit(nodes, job);
}

void GridWriterFull::WriteGrid(std::vector<GridNode> *nodes,
                               std::vector<double> *data, const char *file,
                               bool ascii) {
  Job *job = new Job;
  job->values.assign(data->begin(), data->end());
  job->file = file;
  job->ascii = ascii;
  Submit(nodes, job);
}

void GridWriterFull::Flush() {
#ifndef _WIN32
  std::unique_lock<std::mutex> guard(lock);
  while (busy > 0) {
    done.wait(guard);
  }
#endif
}

void GridWriterFull::Submit(std::vector<GridNode> *nodes, Job *job) {
  if (!grid.data) {
    Initialize();
  }
  if (nodes != indexedNodes) {
    Flush(); // The writer must be idle before the cell map changes
    IndexNodes(nodes);
  }

#ifndef _WIN32
  if (writer.joinable()) {
    std::unique_lock<std::mutex> guard(lock);
    while (busy >= maxQueuedGrids) {
      done.wait(guard);
    }
    jobs.push_back(job);
    busy++;
    wake.notify_one();
    return;
  }
#endif

  Write(job);
  delete job;
}

void GridWriterFull::IndexNodes(std::vector<GridNode> *nodes) {
  // The nodes are fixed for a simulation, so their raster cells are looked up
  // once
  size_t numNodes = nodes->size();
  cells.assign(numNodes, NULL);
  for (size_t i = 0; i < numNodes; i++) {
    GridNode *node = &(nodes->at(i));
    if (node->gauge) {
      cells[i] = &(grid.data[node->y][node->x]);
    }
  }
  indexedNodes = nodes;
}

void GridWriterFull::Write(Job *job) {
  size_t numNodes = cells.size();
  for (size_t i = 0; i < numNodes; i++) {
    if (cells[i]) {
      *(cells[i]) = job->values.at(i);
    }
  }

  if (job->ascii) {
    WriteFloatAscGrid(job->file.c_str(), &grid);
  } else {
    char *artist = NULL;
    char *copyright = NULL;
    char *datetime = NULL;
    if (g_basicConfig->GetArtist()[0]) {
      artist = g_basicConfig->GetArtist();
    }
    if (g_basicConfig->GetCopyright()[0]) {
      copyright = g_basicConfig->GetCopyright();
    }
    WriteFloatTifGrid(job->file.c_str(), &grid, artist, datetime, copyright);
  }
}

#ifndef _WIN32
void GridWriterFull::WriterLoop() {
  std::unique_lock<std::mutex> guard(lock);
  while (true) {
    while (jobs.empty() && !stopping) {
      wake.wait(guard);
    }
    if (jobs.empty()) {
      return; // Stopping and nothing left to write
    }
    Job *job = jobs.front();
    jobs.pop_front();
    guard.unlock();
    Write(job);
    delete job;
    guard.lock();
    busy--;
    done.notify_all();
  }
}
#endif
