#ifndef KW_MODEL_H
#define KW_MODEL_H

#include "ModelBase.h"

enum KW_LAYER {
  KW_LAYER_FASTFLOW,
  KW_LAYER_INTERFLOW,
  KW_LAYER_BASEFLOW,
  KW_LAYER_QTY,
};

enum STATES_KW { STATE_KW_PQ, STATE_KW_PO, STATE_KW_IR, STATE_KW_QTY };

struct KWGridNode : BasicGridNode {
  float params[PARAM_KINEMATIC_QTY];
  float states[STATE_KW_QTY];

  bool channelGridCell;
  double slopeSqrt;
  double hillSlopeSqrt;

  double
      nexTime[KW_LAYER_QTY]; // This is a by product of computing cell routing
  KWGridNode *routeCNode[2][KW_LAYER_QTY]; // This is the node we route water to
  GridNode *routeNode[2][KW_LAYER_QTY];
  double routeAmount[2][KW_LAYER_QTY];
  double incomingWater[KW_LAYER_QTY];

  // double reservoirs[KW_LAYER_QTY]; // CREST has two excess storage reservoirs
  // (overland & interflow)

  // double previousStreamflow; //cms
  // double previousOverland;
  double incomingWaterOverland, incomingWaterChannel;
};

class KWRoute : public RoutingModel {

public:
  KWRoute();
  ~KWRoute();
  float SetObsInflow(long index, float inflow);
  bool InitializeModel(std::vector<GridNode> *newNodes,
                       std::map<GaugeConfigSection *, float *> *paramSettings,
                       std::vector<FloatGrid *> *paramGrids);
  void InitializeStates(TimeVar *beginTime, char *statePath,
                        std::vector<float> *fastFlow,
                        std::vector<float> *interFlow,
                        std::vector<float> *baseFlow);
  void SaveStates(TimeVar *currentTime, char *statePath,
                  GridWriterFull *gridWriter);
  bool Route(float stepHours, std::vector<float> *fastFlow,
             std::vector<float> *interFlow, std::vector<float> *baseFlow, std::vector<float> *discharge);
  float GetMaxSpeed() { return maxSpeed; }

private:
  void RouteInt(float stepSeconds, GridNode *node, KWGridNode *cNode,
                float fastFlow, float interFlow, float baseFlow);
  void
  InitializeParameters(std::map<GaugeConfigSection *, float *> *paramSettings,
                       std::vector<FloatGrid *> *paramGrids);
  void InitializeRouting(float timeSeconds);
  void BuildSchedule();
  void GatherInflow(unsigned long index, KWGridNode *cNode);

  std::vector<GridNode> *nodes;
  std::vector<KWGridNode> kwNodes;
  float maxSpeed;
  bool initialized;

  // Parallel routing schedule. Nodes are grouped into levels so that every
  // node's upstream neighbours and interflow leak sources are in earlier
  // levels; nodes within a level are routed concurrently. Each node pulls its
  // inflow in the same order the serial sweep pushed it, so results do not
  // depend on the number of threads.
  std::vector<unsigned int> levelNodes; // node indices grouped by level
  std::vector<unsigned int> levelStart; // offset of each level in levelNodes
  std::vector<unsigned int> upStart;    // offset of each node in upNodes
  std::vector<unsigned int> upNodes;    // direct upstream nodes, decreasing
  std::vector<unsigned int> leakStart;  // offset of each node in leakSrc
  std::vector<unsigned int> leakSrc;    // 2 * source + slot, decreasing source
  std::vector<double> leakOut;          // interflow leak per (source, slot)
};

#endif
