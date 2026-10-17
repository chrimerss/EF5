#include "KinematicRoute.h"
#include "AscGrid.h"
#include "DatedName.h"
#include "Messages.h"
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>

static const char *stateStrings[] = {
    "pCQ",
    "pOQ",
    "IR",
};

// Solve D * q + alpha * q^beta = rhs for q >= 0. The left side increases with
// q, so the root is unique and lies in [0, rhs / D]; Newton steps that leave
// that bracket fall back to bisection. The tolerance is relative so small
// flows are resolved as accurately as large ones.
static float SolveKW(float D, float alpha, float beta, float rhs, float guess) {
  if (!(rhs > 0.0f)) {
    return 0.0f;
  }
  double lo = 0.0, hi = rhs / D;
  if (!(alpha > 0.0f)) {
    return hi;
  }
  double q = (guess > lo && guess < hi) ? guess : 0.5 * hi;
  for (int itr = 0; itr < 60; itr++) {
    double qBeta = pow(q, (double)beta);
    double resError = D * q + alpha * qBeta - rhs;
    if (fabs(resError) <= 1e-6 * rhs || hi - lo <= 1e-12 * hi) {
      break;
    }
    if (resError > 0.0) {
      hi = q;
    } else {
      lo = q;
    }
    double next = q - resError / (D + alpha * beta * qBeta / q);
    q = (next > lo && next < hi) ? next : 0.5 * (lo + hi);
  }
  return q;
}

KWRoute::KWRoute() {}

KWRoute::~KWRoute() {}

float KWRoute::SetObsInflow(long index, float inflow) {
  KWGridNode *cNode = &(kwNodes[index]);
  GridNode *node = &nodes->at(index);
  float prev;
  if (!node->channelGridCell) {
    prev = cNode->states[STATE_KW_PQ] * cNode->flowWidth;
    cNode->states[STATE_KW_PQ] = inflow / cNode->flowWidth;
    cNode->incomingWaterOverland = inflow;
  } else {
    prev = cNode->states[STATE_KW_PQ];
    float diff = 0.0;
    if (inflow > 1.0 && prev > 1.0) {
      diff = inflow / prev - 1.0;
    }
    GaugeConfigSection *thisGauge = node->gauge;
    float multiplier = 1.0;
    if (diff > 0.2) {
      multiplier = 0.5;
    } else if (diff < -0.2) {
      multiplier = 2.0;
    }
    printf(" Multiplier %f,%f,%f,%f ", multiplier, diff, inflow, prev);
    if (multiplier != 1.0) {
      size_t numNodes = nodes->size();
      for (size_t i = 0; i < numNodes; i++) {
        node = &nodes->at(i);
        if (node->gauge == thisGauge) {
          KWGridNode *cNode = &(kwNodes[i]);
          cNode->params[PARAM_KINEMATIC_ALPHA] *= multiplier;
          if (cNode->params[PARAM_KINEMATIC_ALPHA] < 0.01) {
            cNode->params[PARAM_KINEMATIC_ALPHA] = 0.01;
          } else if (cNode->params[PARAM_KINEMATIC_ALPHA] > 200.0) {
            cNode->params[PARAM_KINEMATIC_ALPHA] = 200.0;
          }
        }
      }
    }
    cNode->states[STATE_KW_PQ] = inflow;
    cNode->incomingWaterChannel = inflow;
  }
  return prev;
}

bool KWRoute::InitializeModel(
    std::vector<GridNode> *newNodes,
    std::map<GaugeConfigSection *, float *> *paramSettings,
    std::vector<FloatGrid *> *paramGrids) {

  nodes = newNodes;
  if (kwNodes.size() != nodes->size()) {
    kwNodes.resize(nodes->size());
  }

  // Fill in modelIndex in the gridNodes
  size_t numNodes = nodes->size();
  for (size_t i = 0; i < numNodes; i++) {
    GridNode *node = &nodes->at(i);
    node->modelIndex = i;
    KWGridNode *cNode = &(kwNodes[i]);
    cNode->slopeSqrt = pow(node->slope, 0.5f);
    cNode->hillSlopeSqrt = pow(node->slope * 0.5, 0.5f);
    cNode->flowWidth = node->area * 1000000.0 / node->horLen; // km^2 to m^2
    cNode->incomingWater[KW_LAYER_BASEFLOW]=0.0;
    cNode->incomingWater[KW_LAYER_INTERFLOW] = 0.0;
    cNode->incomingWater[KW_LAYER_FASTFLOW] = 0.0;
    cNode->incomingWaterOverland = 0.0;
    cNode->incomingWaterChannel = 0.0;
    for (int p = 0; p < STATE_KW_QTY; p++) {
      cNode->states[p] = 0.0;
    }
  }

  InitializeParameters(paramSettings, paramGrids);
  initialized = false;
  maxSpeed = 1.0;

  return true;
}

void KWRoute::InitializeStates(TimeVar *beginTime, char *statePath,
                               std::vector<float> *fastFlow,
                               std::vector<float> *interFlow,
                               std::vector<float> *baseFlow) {
  DatedName timeStr;
  timeStr.SetNameStr("YYYYMMDD_HHUU");
  timeStr.ProcessNameLoose(NULL);
  timeStr.UpdateName(beginTime->GetTM());

  char buffer[300];
  for (int p = 0; p < STATE_KW_QTY; p++) {
    sprintf(buffer, "%s/kwr_%s_%s.tif", statePath, stateStrings[p],
            timeStr.GetName());

    FloatGrid *sGrid = ReadFloatTifGrid(buffer);
    if (sGrid) {
      printf("Using Kinematic Wave Routing %s State Grid %s\n", stateStrings[p],
             buffer);
      if (g_DEM->IsSpatialMatch(sGrid)) {
        for (size_t i = 0; i < nodes->size(); i++) {
          GridNode *node = &nodes->at(i);
          KWGridNode *cNode = &(kwNodes[i]);
          if (sGrid->data[node->y][node->x] != sGrid->noData) {
            cNode->states[p] = sGrid->data[node->y][node->x];
          }
        }
      } else {
        GridLoc pt;
        for (size_t i = 0; i < nodes->size(); i++) {
          GridNode *node = &(nodes->at(i));
          KWGridNode *cNode = &(kwNodes[i]);
          if (sGrid->GetGridLoc(node->refLoc.x, node->refLoc.y, &pt) &&
              sGrid->data[pt.y][pt.x] != sGrid->noData) {
            cNode->states[p] = sGrid->data[pt.y][pt.x];
          }
        }
      }
      delete sGrid;
    } else {
      printf("Kinematic Wave Routing %s State Grid %s not found!\n",
             stateStrings[p], buffer);
    }
  }
}

void KWRoute::SaveStates(TimeVar *currentTime, char *statePath,
                         GridWriterFull *gridWriter) {
  DatedName timeStr;
  timeStr.SetNameStr("YYYYMMDD_HHUU");
  timeStr.ProcessNameLoose(NULL);
  timeStr.UpdateName(currentTime->GetTM());

  std::vector<float> dataVals;
  dataVals.resize(nodes->size());

  char buffer[300];
  for (int p = 0; p < STATE_KW_QTY; p++) {
    sprintf(buffer, "%s/kwr_%s_%s.tif", statePath, stateStrings[p],
            timeStr.GetName());
    for (size_t i = 0; i < nodes->size(); i++) {
      KWGridNode *cNode = &(kwNodes[i]);
      dataVals[i] = cNode->states[p];
    }
    gridWriter->WriteGrid(nodes, &dataVals, buffer, false);
  }
}

bool KWRoute::Route(float stepHours, std::vector<float> *fastFlow,
                    std::vector<float> *interFlow,
                    std::vector<float> *baseFlow,
                    std::vector<float> *discharge) {

  if (!initialized) {
    initialized = true;
    InitializeRouting(stepHours * 3600.0f);
    BuildSchedule();
  }

  if (levelStart.empty()) {
    return false; // BuildSchedule failed and said why
  }

  long numNodes = (long)nodes->size();
  long numLevels = (long)levelStart.size() - 1;

  // Sweep from the headwaters down, one level at a time
#if _OPENMP
#pragma omp parallel
#endif
  for (long l = 0; l < numLevels; l++) {
#if _OPENMP
#pragma omp for schedule(static)
#endif
    for (long k = levelStart[l]; k < (long)levelStart[l + 1]; k++) {
      unsigned long i = levelNodes[k];
      KWGridNode *cNode = &(kwNodes[i]);
      GatherInflow(i, cNode);
      RouteInt(stepHours * 3600.0f, &(nodes->at(i)), cNode, fastFlow->at(i),
               interFlow->at(i), baseFlow->at(i));
    }
  }

#if _OPENMP
#pragma omp parallel for
#endif
  for (long i = 0; i < numNodes; i++) {
    KWGridNode *cNode = &(kwNodes[i]);
    if (!cNode->channelGridCell) {
      // Interflow leaked into this cell during the sweep
      for (unsigned int e = leakStart[i]; e < leakStart[i + 1]; e++) {
        cNode->incomingWater[KW_LAYER_BASEFLOW] += leakOut[leakSrc[e]];
      }
    }
    interFlow->at(i) = 0.0; // cNode->incomingWater[KW_LAYER_INTERFLOW];
    baseFlow ->at(i)=  0.0; // cNode->incomingWater[KW_LAYER_BASEFLOW];
    fastFlow->at(i) = 0.0; // cNode->incomingWater[KW_LAYER_FASTFLOW];
    cNode->incomingWaterOverland = 0.0;
    cNode->incomingWaterChannel = 0.0;
    if (!cNode->channelGridCell) {
      float q = cNode->incomingWater[KW_LAYER_FASTFLOW] * cNode->flowWidth;
      q += (cNode->incomingWater[KW_LAYER_INTERFLOW] * nodes->at(i).area / 3.6);
      q += (cNode -> incomingWater[KW_LAYER_BASEFLOW]*nodes->at(i).area/3.6);
      discharge->at(i) = q; // * (stepHours * 3600.0f);
    } else {
      discharge->at(i) = cNode->incomingWater[KW_LAYER_FASTFLOW];
    }
    cNode->states[STATE_KW_IR] =
        cNode->states[STATE_KW_IR] + cNode->incomingWater[KW_LAYER_INTERFLOW] + cNode->incomingWater[KW_LAYER_BASEFLOW];
    cNode->incomingWater[KW_LAYER_INTERFLOW] =
        0.0; // Zero here so we can save states
    cNode->incomingWater[KW_LAYER_FASTFLOW] =
        0.0; // Zero here so we can save states
    cNode->incomingWater[KW_LAYER_BASEFLOW] =
        0.0; // Zero here so we can save states        
  }

  // InitializeRouting(stepHours * 3600.0f);

  return true;
}

void KWRoute::RouteInt(float stepSeconds, GridNode *node, KWGridNode *cNode,
                       float fastFlow, float interFlow, float baseFlow) {

  if (!cNode->channelGridCell) {
    // overland routing
    // printf("%f\n",baseFlow);

    float beta = 0.6;
    float alpha = cNode->params[PARAM_KINEMATIC_ALPHA0];

    fastFlow /= 1000.0;          // mm to m
    float newInWater = fastFlow; // / node->horLen;

    float A, B, C, D, E;
    float backDiffq = 0.0;
    double qIn = cNode->incomingWaterOverland / cNode->flowWidth; // m^2/s
    if (qIn + cNode->states[STATE_KW_PQ] > 0.0) {
      backDiffq = pow((qIn + cNode->states[STATE_KW_PQ]) / 2.0, beta - 1.0);
      if (!std::isfinite(backDiffq)) {
        backDiffq = 0.0;
      }
    }
    A = (stepSeconds / node->horLen) * qIn;
    B = alpha * beta * cNode->states[STATE_KW_PQ] * backDiffq;
    C = stepSeconds * newInWater;
    D = stepSeconds / node->horLen;
    E = alpha * beta * backDiffq;
    float estq = (A + B + C) / (D + E); // cms/m
    float rhs = A + alpha * pow(cNode->states[STATE_KW_PQ], beta) +
                stepSeconds * newInWater;
    estq = SolveKW(D, alpha, beta, rhs, estq);

    float newq = estq;

    cNode->states[STATE_KW_PQ] = newq;

    // The downstream node pulls this in GatherInflow
    cNode->incomingWater[KW_LAYER_FASTFLOW] = newq;

    // Add Interflow Excess Water to Reservoir
    cNode->states[STATE_KW_IR] += (interFlow+baseFlow);
    double interflowLeak =
        cNode->states[STATE_KW_IR] * cNode->params[PARAM_KINEMATIC_LEAKI];
    // printf(" %f ", interflowLeak);
    cNode->states[STATE_KW_IR] -= interflowLeak;
    if (cNode->states[STATE_KW_IR] < 0) {
      cNode->states[STATE_KW_IR] = 0;
    }

    if (cNode->routeCNode[0][KW_LAYER_INTERFLOW]) {
      double interflowLeak0 =
          interflowLeak * cNode->routeAmount[0][KW_LAYER_INTERFLOW] *
          node->area / cNode->routeNode[0][KW_LAYER_INTERFLOW]->area;
      double leakAmount = interflowLeak0;
      /*if (cNode->routeNode[0][KW_LAYER_INTERFLOW]->channelGridCell) {
      printf(" 0 got %f ",
      cNode->routeCNode[0][KW_LAYER_INTERFLOW]->incomingWater[KW_LAYER_INTERFLOW]);
      }*/
      // printf("0 got here... incoming base flow: %f, leakAmount: %f\n",cNode->routeCNode[0][KW_LAYER_INTERFLOW]->incomingWater[KW_LAYER_INTERFLOW], leakAmount);
      // The receiving node pulls this from leakOut
      leakOut[2 * node->modelIndex] = leakAmount;
    }


    if (cNode->routeCNode[1][KW_LAYER_INTERFLOW]) {
      double interflowLeak1 =
          interflowLeak * cNode->routeAmount[1][KW_LAYER_INTERFLOW] *
          node->area / cNode->routeNode[1][KW_LAYER_INTERFLOW]->area;
      double leakAmount = interflowLeak1;
      // printf(" 1 got %f ", leakAmount);
      leakOut[2 * node->modelIndex + 1] = leakAmount;
    }


  } else {

    // First do overland routing

    float beta = 0.6;
    float alpha = cNode->params[PARAM_KINEMATIC_ALPHA0];
    interFlow += cNode->incomingWater[KW_LAYER_INTERFLOW] +
                 cNode->incomingWater[KW_LAYER_BASEFLOW];
    // printf(" %f ", cNode->incomingWater[KW_LAYER_INTERFLOW]);
    fastFlow /= 1000.0; // mm to m
    baseFlow /= 1000.0;
    interFlow /= 1000.0; // mm to m
    float newInWater = (fastFlow + interFlow+ baseFlow);
    // printf("fast flow: %f, interflow: %f, baseflow: %f\n",fastFlow,interFlow,baseFlow);

    float A, B, C, D, E;
    float backDiffq = 0.0;
    double qIn = cNode->incomingWaterOverland / cNode->flowWidth; // m^2/s
    if (qIn + cNode->states[STATE_KW_PO] > 0.0) {
      backDiffq = pow((qIn + cNode->states[STATE_KW_PO]) / 2.0, beta - 1.0);
      if (!std::isfinite(backDiffq)) {
        backDiffq = 0.0;
      }
    }
    A = (stepSeconds / node->horLen) * qIn;
    B = alpha * beta * cNode->states[STATE_KW_PO] * backDiffq;
    C = stepSeconds * newInWater;
    D = stepSeconds / node->horLen;
    E = alpha * beta * backDiffq;
    float estq = (A + B + C) / (D + E); // cms/m
    float rhs = A + alpha * pow(cNode->states[STATE_KW_PO], beta) +
                stepSeconds * newInWater;
    estq = SolveKW(D, alpha, beta, rhs, estq);
    float newq = estq;

    cNode->states[STATE_KW_PO] = newq;

    // Here we compute channel routing
    beta = cNode->params[PARAM_KINEMATIC_BETA];
    alpha = cNode->params[PARAM_KINEMATIC_ALPHA];

    // Channel Flow
    // Compute Q at current grid point
    float backDiffQ = 0.0;
    if (cNode->incomingWaterChannel + cNode->states[STATE_KW_PQ] > 0.0) {
      backDiffQ =
          pow((cNode->incomingWaterChannel + cNode->states[STATE_KW_PQ]) / 2.0,
              beta - 1.0);
      if (!std::isfinite(backDiffQ)) {
        backDiffQ = 0.0;
      }
    }

    A = (stepSeconds / node->horLen) * cNode->incomingWaterChannel;
    B = alpha * beta * cNode->states[STATE_KW_PQ] * backDiffQ;
    float lateral = newq * cNode->flowWidth / node->horLen; // m^2/s
    C = stepSeconds * lateral;
    D = stepSeconds / node->horLen;
    E = alpha * beta * backDiffQ;
    float estQ = (A + B + C) / (D + E); // cms
    rhs = A + alpha * pow(cNode->states[STATE_KW_PQ], beta) + C;
    estQ = SolveKW(D, alpha, beta, rhs, estQ);
    float newWater = estQ;
    /*if (newWater != newWater) {
    printf("New water is %f (%f, %f) %f %f [%f %f %f %f %f] %f %f\n", newWater,
    cNode->incomingWaterChannel, cNode->states[STATE_KW_PQ], newq,
    cNode->incomingWaterOverland, A, B, C, D, E, alpha, 0.0);
    }*/
    cNode->states[STATE_KW_PQ] =
        newWater; // Update previous Q for further routing if "steps" > 1
    // The downstream node pulls this in GatherInflow

    cNode->incomingWater[KW_LAYER_FASTFLOW] = newWater;
    cNode->incomingWater[KW_LAYER_INTERFLOW] = 0.0;
    cNode->incomingWater[KW_LAYER_BASEFLOW] = 0.0;
  }
}

void KWRoute::InitializeParameters(
    std::map<GaugeConfigSection *, float *> *paramSettings,
    std::vector<FloatGrid *> *paramGrids) {

  // This pass distributes parameters
  size_t numNodes = nodes->size();
  size_t unused = 0;
  for (size_t i = 0; i < numNodes; i++) {
    GridNode *node = &nodes->at(i);
    KWGridNode *cNode = &(kwNodes[i]);
    if (!node->gauge) {
      unused++;
      continue;
    }
    // Copy all of the parameters over
    memcpy(cNode->params, (*paramSettings)[node->gauge],
           sizeof(float) * PARAM_KINEMATIC_QTY);

    if (!paramGrids->at(PARAM_KINEMATIC_ISU)) {
      cNode->states[STATE_KW_IR] = cNode->params[PARAM_KINEMATIC_ISU];
    }
    cNode->incomingWater[KW_LAYER_INTERFLOW] = 0.0;
    cNode->incomingWater[KW_LAYER_BASEFLOW] = 0.0;
    cNode->incomingWater[KW_LAYER_FASTFLOW] = 0.0;

    // Deal with the distributed parameters here
    GridLoc pt;
    for (size_t paramI = 0; paramI < PARAM_KINEMATIC_QTY; paramI++) {
      FloatGrid *grid = paramGrids->at(paramI);
      if (grid && g_DEM->IsSpatialMatch(grid)) {
        if (grid->data[node->y][node->x] == 0) {
          grid->data[node->y][node->x] = 0.01;
        }
        cNode->params[paramI] *= grid->data[node->y][node->x];
      } else if (grid &&
                 grid->GetGridLoc(node->refLoc.x, node->refLoc.y, &pt)) {
        if (grid->data[pt.y][pt.x] == 0) {
          grid->data[pt.y][pt.x] = 0.01;
          // printf("Using nodata value in param %s\n",
          // modelParamStrings[MODEL_CREST][paramI]);
        }
        cNode->params[paramI] *= grid->data[pt.y][pt.x];
      }
    }

    if (cNode->params[PARAM_KINEMATIC_LEAKI] < 0.0) {
      // printf("Node Leak Interflow(%f) is less than 0, setting to 0.\n",
      // cNode->params[PARAM_KINEMATIC_LEAKI]);
      cNode->params[PARAM_KINEMATIC_LEAKI] = 0.0;
    } else if (cNode->params[PARAM_KINEMATIC_LEAKI] > 1.0) {
      // printf("Node Leak Interflow(%f) is greater than 1, setting to 1.\n",
      // cNode->params[PARAM_KINEMATIC_LEAKI]);
      cNode->params[PARAM_KINEMATIC_LEAKI] = 1.0;
    }

    if (cNode->params[PARAM_KINEMATIC_ALPHA] < 0.0) {
      // printf("Node Alpha(%f) is less than 0, setting to 1.\n",
      // cNode->params[PARAM_KINEMATIC_ALPHA]);
      cNode->params[PARAM_KINEMATIC_ALPHA] = 1.0;
    }

    if (cNode->params[PARAM_KINEMATIC_ALPHA0] < 0.0) {
      // printf("Node Alpha0(%f) is less than 0, setting to 1.\n",
      // cNode->params[PARAM_KINEMATIC_ALPHA0]);
      cNode->params[PARAM_KINEMATIC_ALPHA0] = 1.0;
    }

    if (cNode->params[PARAM_KINEMATIC_BETA] < 0.0) {
      // printf("Node Beta(%f) is less than 0, setting to 0.6.\n",
      // cNode->params[PARAM_KINEMATIC_BETA]);
      cNode->params[PARAM_KINEMATIC_BETA] = 0.6;
    }

    if (node->fac > cNode->params[PARAM_KINEMATIC_TH]) {
      node->channelGridCell = true;
      cNode->channelGridCell = true;
    } else {
      node->channelGridCell = false;
      cNode->channelGridCell = false;
    }
  }
}

void KWRoute::InitializeRouting(float timeSeconds) {

  // This pass distributes parameters & calculates the time it takes for water
  // to cross the grid cell.
  size_t numNodes = nodes->size();
#if _OPENMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < numNodes; i++) {
    GridNode *node = &nodes->at(i);
    KWGridNode *cNode = &(kwNodes[i]);

    // Calculate the water speed for interflow
    float speedUnder = cNode->params[PARAM_KINEMATIC_UNDER] * cNode->slopeSqrt;

    float nexTimeUnder = node->horLen / speedUnder;
    cNode->nexTime[KW_LAYER_INTERFLOW] = nexTimeUnder;
    cNode->nexTime[KW_LAYER_BASEFLOW] = nexTimeUnder;
  }

  // This pass figures out which cell water is routed to
#if _OPENMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < numNodes; i++) {
    GridNode *currentNode, *previousNode;
    float currentSeconds, previousSeconds;
    GridNode *node = &nodes->at(i);
    KWGridNode *cNode = &(kwNodes[i]);

    // Interflow routing
    previousSeconds = 0;
    currentSeconds = 0;
    currentNode = node;
    previousNode = NULL;
    while (currentSeconds < timeSeconds && currentNode &&
           !kwNodes[currentNode->modelIndex].channelGridCell) {
      if (currentNode) {
        previousSeconds = currentSeconds;
        previousNode = currentNode;
        currentSeconds +=
            kwNodes[currentNode->modelIndex].nexTime[KW_LAYER_INTERFLOW];
        if (currentNode->downStreamNode != INVALID_DOWNSTREAM_NODE) {
          currentNode = &(nodes->at(currentNode->downStreamNode));
        } else {
          currentNode = NULL;
        }
      } else {
        if (timeSeconds > currentSeconds) {
          previousNode = NULL;
        }
        break; // We have effectively run out of nodes to transverse, this is
               // done!
      }
    }

    cNode->routeNode[0][KW_LAYER_INTERFLOW] = currentNode;
    cNode->routeCNode[0][KW_LAYER_INTERFLOW] =
        (currentNode) ? &(kwNodes[currentNode->modelIndex]) : NULL;
    cNode->routeNode[1][KW_LAYER_INTERFLOW] = previousNode;
    cNode->routeCNode[1][KW_LAYER_INTERFLOW] =
        (previousNode) ? &(kwNodes[previousNode->modelIndex]) : NULL;
    cNode->routeNode[0][KW_LAYER_BASEFLOW] = currentNode;
    cNode->routeCNode[0][KW_LAYER_BASEFLOW] =
        (currentNode) ? &(kwNodes[currentNode->modelIndex]) : NULL;
    cNode->routeNode[1][KW_LAYER_BASEFLOW] = previousNode;
    cNode->routeCNode[1][KW_LAYER_BASEFLOW] =
        (previousNode) ? &(kwNodes[previousNode->modelIndex]) : NULL;        
    if (currentNode && !kwNodes[currentNode->modelIndex].channelGridCell) {
      if ((currentSeconds - previousSeconds) > 0) {
        cNode->routeAmount[0][KW_LAYER_INTERFLOW] =
            (timeSeconds - previousSeconds) /
            (currentSeconds - previousSeconds);
        cNode->routeAmount[1][KW_LAYER_INTERFLOW] =
            1.0 - cNode->routeAmount[0][KW_LAYER_INTERFLOW];

        cNode->routeAmount[0][KW_LAYER_BASEFLOW] =
            (timeSeconds - previousSeconds) /
            (currentSeconds - previousSeconds);
        cNode->routeAmount[1][KW_LAYER_BASEFLOW] =
            1.0 - cNode->routeAmount[0][KW_LAYER_BASEFLOW];
      }      
    } else {
      cNode->routeAmount[0][KW_LAYER_INTERFLOW] = 1.0;
      cNode->routeAmount[1][KW_LAYER_INTERFLOW] = 0.0;
      cNode->routeAmount[0][KW_LAYER_BASEFLOW] = 1.0;
      cNode->routeAmount[1][KW_LAYER_BASEFLOW] = 0.0;      
    }
  }
}

void KWRoute::BuildSchedule() {

  size_t numNodes = nodes->size();
  if (numNodes >= UINT_MAX / 2) {
    ERROR_LOGF("Too many nodes (%lu) for the routing schedule!",
               (unsigned long)numNodes);
    return;
  }

  // A node's level is the longest path to it from the headwaters, so all of
  // its upstream nodes sit in lower levels
  std::vector<unsigned int> level(numNodes, 0), pending(numNodes, 0);
  std::vector<unsigned int> ready;
  for (size_t i = 0; i < numNodes; i++) {
    GridNode *node = &nodes->at(i);
    if (node->downStreamNode != INVALID_DOWNSTREAM_NODE) {
      pending[node->downStreamNode]++;
    }
  }
  for (size_t i = 0; i < numNodes; i++) {
    if (pending[i] == 0) {
      ready.push_back(i);
    }
  }
  size_t visited = 0;
  unsigned int numLevels = 0;
  while (!ready.empty()) {
    unsigned int i = ready.back();
    ready.pop_back();
    visited++;
    if (level[i] + 1 > numLevels) {
      numLevels = level[i] + 1;
    }
    unsigned long down = nodes->at(i).downStreamNode;
    if (down != INVALID_DOWNSTREAM_NODE) {
      if (level[down] < level[i] + 1) {
        level[down] = level[i] + 1;
      }
      if (--pending[down] == 0) {
        ready.push_back(down);
      }
    }
  }
  if (visited != numNodes) {
    WARNING_LOGF("%lu nodes are in a flow direction loop!",
                 (unsigned long)(numNodes - visited));
  }

  // Group the nodes by level, in increasing index within each level
  levelStart.assign(numLevels + 1, 0);
  for (size_t i = 0; i < numNodes; i++) {
    levelStart[level[i] + 1]++;
  }
  for (unsigned int l = 0; l < numLevels; l++) {
    levelStart[l + 1] += levelStart[l];
  }
  std::vector<unsigned int> fill(levelStart.begin(), levelStart.end() - 1);
  levelNodes.resize(numNodes);
  for (size_t i = 0; i < numNodes; i++) {
    levelNodes[fill[level[i]]++] = i;
  }

  // Direct upstream nodes in decreasing index, the order the serial sweep
  // pushed them downstream
  upStart.assign(numNodes + 1, 0);
  for (size_t i = 0; i < numNodes; i++) {
    unsigned long down = nodes->at(i).downStreamNode;
    if (down != INVALID_DOWNSTREAM_NODE) {
      upStart[down + 1]++;
    }
  }
  for (size_t i = 0; i < numNodes; i++) {
    upStart[i + 1] += upStart[i];
  }
  fill.assign(upStart.begin(), upStart.end() - 1);
  upNodes.resize(upStart[numNodes]);
  for (long i = numNodes - 1; i >= 0; i--) {
    unsigned long down = nodes->at(i).downStreamNode;
    if (down != INVALID_DOWNSTREAM_NODE) {
      upNodes[fill[down]++] = i;
    }
  }

  // Interflow leak sources of each node, also in decreasing source index
  leakStart.assign(numNodes + 1, 0);
  for (size_t i = 0; i < numNodes; i++) {
    KWGridNode *cNode = &(kwNodes[i]);
    for (int slot = 0; slot < 2 && !cNode->channelGridCell; slot++) {
      if (cNode->routeCNode[slot][KW_LAYER_INTERFLOW]) {
        leakStart[cNode->routeNode[slot][KW_LAYER_INTERFLOW]->modelIndex + 1]++;
      }
    }
  }
  for (size_t i = 0; i < numNodes; i++) {
    leakStart[i + 1] += leakStart[i];
  }
  fill.assign(leakStart.begin(), leakStart.end() - 1);
  leakSrc.resize(leakStart[numNodes]);
  for (long i = numNodes - 1; i >= 0; i--) {
    KWGridNode *cNode = &(kwNodes[i]);
    for (int slot = 0; slot < 2 && !cNode->channelGridCell; slot++) {
      if (cNode->routeCNode[slot][KW_LAYER_INTERFLOW]) {
        unsigned long target =
            cNode->routeNode[slot][KW_LAYER_INTERFLOW]->modelIndex;
        leakSrc[fill[target]++] = 2 * i + slot;
      }
    }
  }
  leakOut.assign(2 * numNodes, 0.0);

  INFO_LOGF("Kinematic wave routing schedule has %u levels for %lu nodes",
            numLevels, (unsigned long)numNodes);
}

void KWRoute::GatherInflow(unsigned long index, KWGridNode *cNode) {
  // Overland flow (m^3/s) from non-channel upstream nodes and channel flow
  // from channel ones, summed in the order the serial sweep added them
  for (unsigned int e = upStart[index]; e < upStart[index + 1]; e++) {
    KWGridNode *upNode = &(kwNodes[upNodes[e]]);
    if (upNode->channelGridCell) {
      cNode->incomingWaterChannel += upNode->incomingWater[KW_LAYER_FASTFLOW];
    } else {
      cNode->incomingWaterOverland +=
          upNode->incomingWater[KW_LAYER_FASTFLOW] * upNode->flowWidth;
    }
  }

  // Channel cells use the interflow leaked into them while routing; other
  // cells add it to their reservoir after the sweep
  if (cNode->channelGridCell) {
    for (unsigned int e = leakStart[index]; e < leakStart[index + 1]; e++) {
      cNode->incomingWater[KW_LAYER_BASEFLOW] += leakOut[leakSrc[e]];
    }
  }
}
