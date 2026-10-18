#ifndef GRIDWRITER_H
#define GRIDWRITER_H

#include "AscGrid.h"
#include "Grid.h"
#include "GridNode.h"
#include "TifGrid.h"
#include <string>
#include <vector>
#ifndef _WIN32
#include <condition_variable>
#include <deque>
#include <mutex>
#include <thread>
#endif

// Writes per-node values as full rasters. Outside Windows the files are
// written on a background thread: WriteGrid copies the values and returns, so
// compression and disk I/O overlap the next time step. Flush waits until all
// queued grids are on disk.
class GridWriterFull {
public:
  GridWriterFull();
  // Copies start out empty and idle; std::vector<Simulator> needs this and
  // only copies writers that have not been used yet
  GridWriterFull(const GridWriterFull &other);
  GridWriterFull &operator=(const GridWriterFull &other) = delete;
  ~GridWriterFull();
  void Initialize();
  void WriteGrid(std::vector<GridNode> *nodes, std::vector<float> *data,
                 const char *file, bool ascii = true);
  void WriteGrid(std::vector<GridNode> *nodes, std::vector<double> *data,
                 const char *file, bool ascii = true);
  void Flush();

private:
  struct Job {
    std::vector<float> values; // one per node
    std::string file;
    bool ascii;
  };

  void Submit(std::vector<GridNode> *nodes, Job *job);
  void IndexNodes(std::vector<GridNode> *nodes);
  void Write(Job *job);

  FloatGrid grid;
  std::vector<GridNode> *indexedNodes;
  std::vector<float *> cells; // raster cell of each node, NULL if not written

#ifndef _WIN32
  void WriterLoop();

  std::deque<Job *> jobs;
  size_t busy; // jobs queued or being written
  bool stopping;
  std::mutex lock;
  std::condition_variable wake; // signals the writer
  std::condition_variable done; // signals threads waiting on the writer
  std::thread writer;
#endif
};

#endif
