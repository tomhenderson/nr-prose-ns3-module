/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */
/*
 * NIST-developed software is provided by NIST as a public
 * service. You may use, copy and distribute copies of the software in
 * any medium, provided that you keep intact this entire notice. You
 * may improve, modify and create derivative works of the software or
 * any portion of the software, and you may copy and distribute such
 * modifications or works. Modified works should carry a notice
 * stating that you changed the software and should note the date and
 * nature of any such change. Please explicitly acknowledge the
 * National Institute of Standards and Technology as the source of the
 * software.
 *
 * NIST-developed software is expressly provided "AS IS." NIST MAKES
 * NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT OR ARISING BY
 * OPERATION OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 * NON-INFRINGEMENT AND DATA ACCURACY. NIST NEITHER REPRESENTS NOR
 * WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED
 * OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES NOT
 * WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE
 * SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE
 * CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.
 *
 * You are solely responsible for determining the appropriateness of
 * using and distributing the software and you assume all risks
 * associated with its use, including but not limited to the risks and
 * costs of program errors, compliance with applicable laws, damage to
 * or loss of data, programs or equipment, and the unavailability or
 * interruption of operation. This software is not intended to be used
 * in any situation where a failure could cause risk of injury or
 * damage to property. The software developed by NIST employees is not
 * subject to copyright protection within the United States.
 */

/**
 * \ingroup examples
 * \file nr-prose-u2u-multihop.cc
 * \brief A proof of concept of sidelink multihop UE-to-UE relay
 *
 * TODO!
 *
 * \code{.unparsed}
$ ./ns3 run "nr-prose-u2u-multihop --help"
    \endcode
 *
 */

#include "ns3/core-module.h"
#include "ns3/config-store.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/applications-module.h"
#include "ns3/mobility-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/nr-module.h"
#include "ns3/nr-prose-module.h"
#include "ns3/lte-module.h"
#include "ns3/stats-module.h"
#include "ns3/config-store-module.h"
#include "ns3/log.h"
#include "ns3/antenna-module.h"
#include <iomanip>
#include <sqlite3.h>
#include "ns3/flow-monitor-module.h"

using namespace ns3;

NS_LOG_COMPONENT_DEFINE ("NrProseU2uMultihop");

/*
 * Global methods and variable to hook trace sources from different layers of
 * the protocol stack.
 */

/**
 * \brief Method to listen to the trace SlPscchScheduling of NrUeMac, which gets
 *        triggered upon the transmission of SCI format 1-A from UE MAC.
 *
 * \param pscchStats Pointer to the \link UeMacPscchTxOutputStats \endlink class,
 *        which is responsible for writing the trace source parameters to a database.
 * \param pscchStatsParams Parameters of the trace source.
 */
void
NotifySlPscchScheduling (UeMacPscchTxOutputStats *pscchStats,
                         const SlPscchUeMacStatParameters pscchStatsParams)
{
  pscchStats->Save (pscchStatsParams);
}

/**
 * \brief Method to listen to the trace SlPsschScheduling of NrUeMac, which gets
 *        triggered upon the transmission of SCI format 2-A and data from UE MAC.
 *
 * \param psschStats Pointer to the \link UeMacPsschTxOutputStats \endlink class,
 *        which is responsible for writing the trace source parameters to a database.
 * \param psschStatsParams Parameters of the trace source.
 */
void
NotifySlPsschScheduling (UeMacPsschTxOutputStats *psschStats,
                         const SlPsschUeMacStatParameters psschStatsParams)
{
  psschStats->Save (psschStatsParams);
}

/**
 * \brief Method to listen to the trace RxPscchTraceUe of NrSpectrumPhy, which gets
 *        triggered upon the reception of SCI format 1-A.
 *
 * \param pscchStats Pointer to the \link UePhyPscchRxOutputStats \endlink class,
 *        which is responsible for writing the trace source parameters to a database.
 * \param pscchStatsParams Parameters of the trace source.
 */
void
NotifySlPscchRx (UePhyPscchRxOutputStats *pscchStats,
                 const SlRxCtrlPacketTraceParams pscchStatsParams)
{
  pscchStats->Save (pscchStatsParams);
}

/**
 * \brief Method to listen to the trace RxPsschTraceUe of NrSpectrumPhy, which gets
 *        triggered upon the reception of SCI format 2-A and data.
 *
 * \param psschStats Pointer to the \link UePhyPsschRxOutputStats \endlink class,
 *        which is responsible for writing the trace source parameters to a database.
 * \param psschStatsParams Parameters of the trace source.
 */
void
NotifySlPsschRx (UePhyPsschRxOutputStats *psschStats,
                 const SlRxDataPacketTraceParams psschStatsParams)
{
  psschStats->Save (psschStatsParams);
}

std::map<std::string, uint32_t> g_relayNasPacketCounter;
/*
 * \brief Trace sink function for logging reception of data packets in the NAS
 *        layer by UE(s) acting as relay UE
 *
 * \param stream the output stream wrapper where the trace will be written
 * \param nodeIp the IP of the relay UE
 * \param srcIp the IP of the node sending the packet
 * \param dstIp the IP of the node that would be receiving the packet
 * \param srcLink the link from which the relay UE received the packet (UL, DL, or SL)
 * \param dstLink the link towards which the relay routed the packet (UL, DL, or SL)
 * \param p the packet
 */
void
TraceSinkRelayNasRxPacketTrace (Ptr<OutputStreamWrapper> stream, std::string nodeId,
                                Ipv4Address nodeIp, Ipv4Address srcIp, Ipv4Address dstIp, uint16_t dstPort,
                                std::string srcLink, std::string dstLink, Ptr<Packet> p)
{
  *stream->GetStream () << Simulator::Now ().GetSeconds ()
                        << "\t" << nodeId
                        << "\t" << nodeIp
                        << "\t" << srcIp
                        << "\t" << dstIp
                        << "\t" << dstPort
                        << "\t" << srcLink
                        << "\t" << dstLink
                        << std::endl;
  std::ostringstream  oss;
  oss  << nodeId << "         "  << nodeIp << "         " << srcIp << "->" << dstIp << ":" << dstPort << "         " <<  srcLink << "->" << dstLink;

  std::string mapKey = oss.str ();
  auto it = g_relayNasPacketCounter.find (mapKey);
  if (it == g_relayNasPacketCounter.end ())
  {
    g_relayNasPacketCounter.insert (std::pair < std::string, uint32_t> (mapKey, 1));
  }
  else
  {
    it->second += 1;
  }
}

/* 
 * Global variables to calculate simulation statistics. 
 * Global because depending on the path finding algorithm the calculation is 
 * done either in the main or in a function scheduled for simulation time 
 */
uint16_t g_foundNpaths = 0; 
uint32_t g_minPathCount = std::numeric_limits<uint32_t>::max();
uint32_t g_maxPathCount = 0;
uint64_t g_totalPathCount = 0; 
uint32_t g_minNonZeroPathCount = std::numeric_limits<uint32_t>::max();
uint64_t g_totalNonZeroPathCount = 0;
uint32_t g_nonZeroCount = 0;
double g_meanPathCount = 0;
double g_meanNonZeroPathCount = 0;
double g_d = 0;

//Global variable to propagate schedule type
bool g_dynSch = true;

/*
 * Global variables to count TX/RX packets and bytes in the application layer.
 */
uint32_t rxByteCounter = 0;
uint32_t txByteCounter = 0;
uint32_t rxPktCounter = 0;
uint32_t txPktCounter = 0;

/**
 * \brief Method to listen to the packet sink application trace Rx.
 * \param packet The packet
 * \param The address of the transmitter
 */
void
ReceivePacket (Ptr<const Packet> packet, const Address &from [[maybe_unused]])
{
  rxByteCounter += packet->GetSize ();
  rxPktCounter++;
  //std::cout <<Simulator::Now().GetMilliSeconds() << " Rx" <<std::endl;

}

/**
 * \brief Method to listen to the transmitting application trace Tx.
 * \param packet The packet
 */
void
TransmitPacket (Ptr<const Packet> packet)
{
  txByteCounter += packet->GetSize ();
  txPktCounter++;
  //std::cout <<Simulator::Now().GetMilliSeconds() << " Tx" <<std::endl;
}


/*
* TODO
*/
Ptr<ns3::Node> GetNodeById(const ns3::NodeContainer& nodes, uint32_t nodeId) {
  for (uint32_t i = 0; i < nodes.GetN(); ++i) {
      if (nodes.Get(i)->GetId() == nodeId) {
          return nodes.Get(i);
      }
  }
  return nullptr; // Return nullptr if node is not found
}

/*
* TODO
*/
std::vector<uint32_t> SelectRandomRelays(ns3::NodeContainer ueNodes, uint32_t nRelays, Ptr<ns3::UniformRandomVariable> randomVariable) 
{
  std::vector<uint32_t> selectedNodeIds;
  std::unordered_set<uint32_t> uniqueIndices;

  // Select unique random nodes
  while (selectedNodeIds.size() < nRelays) {
      uint32_t index = randomVariable->GetInteger(0, ueNodes.GetN() - 1);
      if (uniqueIndices.insert(index).second) {
          selectedNodeIds.push_back(ueNodes.Get(index)->GetId());
      }
  }

  return selectedNodeIds;
}

/*
* TODO
*/
uint32_t GetRandomNodeId(ns3::NodeContainer& nodes, Ptr<ns3::UniformRandomVariable> randomVariable) {

  uint32_t randomIndex = randomVariable->GetInteger(0, nodes.GetN() - 1);

  return nodes.Get(randomIndex)->GetId();
}

/*
* TODO
*/
uint32_t GetRandomNodeIdExcludeRelays(ns3::NodeContainer& nodes, std::vector<uint32_t>& relayNodeIds, Ptr<ns3::UniformRandomVariable> randomVariable) {

  uint32_t randomIndex;
  uint32_t nodeId;

  do {
      randomIndex = randomVariable->GetInteger(0, nodes.GetN() - 1);
      nodeId = nodes.Get(randomIndex)->GetId();
  } while (std::find(relayNodeIds.begin(), relayNodeIds.end(), nodeId) != relayNodeIds.end());

  return nodeId;
}

/*
* TODO
*/
double CalculateDistance(Ptr<ns3::Node> a, Ptr<ns3::Node> b) {
  Ptr<ns3::MobilityModel> mobilityA = a->GetObject<ns3::MobilityModel>();
  Ptr<ns3::MobilityModel> mobilityB = b->GetObject<ns3::MobilityModel>();
  return mobilityA->GetDistanceFrom(mobilityB);
}

/*
* TODO
*/
void ExportDataForPlotting(const ns3::NodeContainer& nodes, 
                           const std::vector<uint32_t>& relayNodeIds,
                           const std::map<uint16_t, std::list<Ptr<ns3::Node>>>& paths,
                           double d) {
  std::ofstream file("nodes_and_paths.csv");

  // Write node information
  for (uint32_t i = 0; i < nodes.GetN(); ++i) {
    Ptr<ns3::Node> node = nodes.Get(i);
    Ptr<ns3::MobilityModel> mobility = node->GetObject<ns3::MobilityModel>();
    ns3::Vector pos = mobility->GetPosition();
    bool isRelay = std::find(relayNodeIds.begin(), relayNodeIds.end(), node->GetId()) != relayNodeIds.end();
    file << node->GetId() << "," << pos.x << "," << pos.y << "," << (isRelay ? "relay" : "non-relay") << "\n";
  }

  // Write path information
  for (const auto& pathEntry : paths) {
      file << "Path," << pathEntry.first;
      for (const auto& node : pathEntry.second) {
          file << "," << node->GetId();
      }
      file << "\n";
  }
  file << "D," << d << "\n";;

  file.close();
}


/*
* TODO
*/
struct RelayStats {
    uint32_t relayNodeId;
    uint32_t pathCount;
};

/*
* TODO
*/
std::vector<RelayStats> GetRelayStats(const std::vector<uint32_t>& selectedRelays, 
                                          const std::map<uint16_t, std::list<Ptr<ns3::Node>>>& nodePaths) {
  
  std::vector<RelayStats> result;

  // Iterate through each relay
  for (uint32_t relay : selectedRelays) {
      RelayStats info;
      info.relayNodeId = relay;
      info.pathCount = 0;

      // Iterate through each path
      for (const auto& path : nodePaths) {
          // Check if path contains more than two nodes (i.e., there is a relay)
          if (path.second.size() <= 2) continue;

          // Get the first and last node IDs in the path
          uint32_t firstNodeId = path.second.front()->GetId();
          uint32_t lastNodeId = path.second.back()->GetId();

          // Check if the relay is the first or last node (i.e., an end node)
          if (relay == firstNodeId || relay == lastNodeId) continue;

          // Iterate through each node in the path
          for (const Ptr<ns3::Node>& node : path.second) {
              // Check if the current node is the relay
              if (node->GetId() == relay) {
                  info.pathCount++;
                  break; // Move to next path once the relay is found in the current path
              }
          }
      }

      result.push_back(info);
  }

  return result;
}

/*
* TODO
*/
typedef std::map<uint32_t, std::map<uint32_t, double>> ConnectivityMap;

/*
* TODO
*/
ConnectivityMap GetDistanceConnectivity(NodeContainer& ueNodes) {

  ConnectivityMap connectivityMap;

  for (auto itA = ueNodes.Begin(); itA != ueNodes.End(); ++itA) {
      Ptr<Node> nodeA = *itA;
      std::map<uint32_t, double> nodeDistances;

      for (auto itB = ueNodes.Begin(); itB != ueNodes.End(); ++itB) {
          Ptr<Node> nodeB = *itB;
          if (nodeA != nodeB) {
              Ptr<MobilityModel> mobilityA = nodeA->GetObject<MobilityModel>();
              Ptr<MobilityModel> mobilityB = nodeB->GetObject<MobilityModel>();
              double distance = mobilityA->GetDistanceFrom(mobilityB);
              nodeDistances[nodeB->GetId()] = distance;
          }
      }
      connectivityMap[nodeA->GetId()] = nodeDistances;
  }


  std::cout << "Full Connectivity: \n";
  for (const auto& node : connectivityMap) {
  std::cout << "Node " << node.first << " ";
    for (const auto& connection : node.second) {
        std::cout << "[" << connection.first << ", " << connection.second << "] ";
    }
  std::cout << std::endl;
  }
  return connectivityMap;
}


enum ComparisonOp { LESS_OR_EQUAL, GREATER };
/*
* TODO
*/
ConnectivityMap FilterConnectivity(ConnectivityMap& connectivityMap, double metricConstraint, ComparisonOp op) {

  ConnectivityMap filteredMap;

  for (auto& node : connectivityMap) {
      std::map<uint32_t, double> filteredDistances;

      for (auto& distance : node.second) {
          bool condition = (op == LESS_OR_EQUAL && distance.second <= metricConstraint) ||
                            (op == GREATER && distance.second > metricConstraint);
          if (condition) {
              filteredDistances.insert(distance);
          }
      }

      if (!filteredDistances.empty()) {
          filteredMap[node.first] = filteredDistances;
      }
  }

  
  std::cout << "Filtered Connectivity (" << metricConstraint << "):\n";
  for (const auto& node : filteredMap) {
  std::cout << "Node " << node.first << " ";
    for (const auto& connection : node.second) {
        std::cout << "[" << connection.first << ", " << connection.second << "] ";
    }
  std::cout << std::endl; 
  }

  return filteredMap;
}


/*
* TODO
*/
typedef std::list<Ptr<Node>> Path;
typedef std::pair<Path, double> PathInfo;
typedef std::list<PathInfo> PathList;


/*
* TODO
*/
void DFS(uint32_t currentNodeId, uint32_t tgtNodeId, const ConnectivityMap& connectivityMap, 
         std::vector<bool>& visited, Path& currentPath, PathList& allPaths, double currentMetric, int depth, int maxDepth, 
         const std::vector<uint32_t>& relayNodeIds) {

  // Base conditions
  if (depth > maxDepth || visited[currentNodeId] || connectivityMap.find(currentNodeId) == connectivityMap.end()) {
      return;
  }

  visited[currentNodeId] = true;
  currentPath.push_back(NodeList::GetNode(currentNodeId));

/*     for (const auto& node : currentPath) {
      std::cout << node->GetId() << " -> ";
  }
  std::cout << std::endl; */

  if (currentNodeId == tgtNodeId) {
      allPaths.push_back(make_pair(currentPath, currentMetric));
      NS_LOG_DEBUG ("allPaths.size() " << allPaths.size());

  } else {
      for (const auto& neighbour : connectivityMap.at(currentNodeId)) {
          // Check if the neighbour is a relay node or the target node
          if (std::find(relayNodeIds.begin(), relayNodeIds.end(), neighbour.first) != relayNodeIds.end() ||
              neighbour.first == tgtNodeId) {
              double newMetric = currentMetric + neighbour.second;
              DFS(neighbour.first, tgtNodeId, connectivityMap, visited, currentPath, allPaths, newMetric, depth + 1, maxDepth, relayNodeIds);
          }
      }
  }

  visited[currentNodeId] = false;
  currentPath.pop_back();
}

/*
* TODO
*/
PathList FindAllPaths(const NodeContainer ueNodes, const ConnectivityMap& connectivityMap, 
                      const std::vector<uint32_t>& relayNodeIds, uint32_t srcNodeId, uint32_t tgtNodeId, int maxDepth) {
   
  PathList allPaths;
  Path currentPath;
  std::vector<bool> visited(ueNodes.GetN(), false);
  double currentMetric = 0.0;

  DFS(srcNodeId, tgtNodeId, connectivityMap, visited, currentPath, allPaths, currentMetric, 0, maxDepth, relayNodeIds);

/*     std::cout << "All paths "<< srcNodeId << "->" << tgtNodeId << ":\n";

    for (const auto& pathInfo : allPaths) {
    std::cout << "Path: ";
    for (const auto& node : pathInfo.first) {
        std::cout << node->GetId() << " -> ";
    }
    std::cout << " Total Metric: " << pathInfo.second << std::endl;
  }  */

  return allPaths;
}

/*
* TODO
*/
PathList FilterPathsByRelays(PathList& paths, uint32_t maxNumRelays) {
    
  PathList filteredPaths;

  for (auto& pathInfo : paths) {
      if (pathInfo.first.size() <= maxNumRelays + 2) { // +2 includes both src and tgt
          filteredPaths.push_back(pathInfo);
      }
  }
/*   std::cout << "Filtered paths (max " << maxNumRelays <<" relays):\n";

    for (const auto& pathInfo : filteredPaths) {
      std::cout << "Filtered Paths: ";
      for (const auto& node : pathInfo.first) {
          std::cout << node->GetId() << " -> ";
      }
      std::cout << " Total Metric: " << pathInfo.second << std::endl;
  }  */
  return filteredPaths;
}

/*
* TODO
*/
PathInfo GetSmallestAggregatedMetricPath(PathList& paths) {
   
  PathInfo shortestPath;
  double minMetric = std::numeric_limits<double>::max();

  for (auto& pathInfo : paths) {
      double currentMetric = pathInfo.second;
      if (currentMetric < minMetric ||
          (currentMetric == minMetric && pathInfo.first.size() < shortestPath.first.size())) {
          shortestPath = pathInfo;
          minMetric = currentMetric;
      }
  }
  std::cout << "Shortest Path: ";
  for (auto nodeIt = shortestPath.first.begin(); nodeIt != shortestPath.first.end(); ++nodeIt) {
      std::cout << (*nodeIt)->GetId();
      if (std::next(nodeIt) != shortestPath.first.end()) {
          std::cout << " -> ";
      }
  }
  std::cout << " Metric: " << shortestPath.second << std::endl;

  return shortestPath;
}

/*
* TODO
*/
PathInfo GetLargestAggregatedMetricPath(PathList& paths) {

  PathInfo largestPath;
  // Start with the smallest possible metric (0 if metrics are always positive)
  double maxMetric = -std::numeric_limits<double>::max();

  for (auto& pathInfo : paths) {
      double currentMetric = pathInfo.second;
      // Change comparison to find the largest metric
      if (currentMetric > maxMetric ||
          (currentMetric == maxMetric && pathInfo.first.size() < largestPath.first.size())) {
          largestPath = pathInfo;
          maxMetric = currentMetric;
      }
  }
  std::cout << "Largest Path: ";
  for (auto nodeIt = largestPath.first.begin(); nodeIt != largestPath.first.end(); ++nodeIt) {
      std::cout << (*nodeIt)->GetId();
      if (std::next(nodeIt) != largestPath.first.end()) {
          std::cout << " -> ";
      }
  }
  std::cout << " Metric: " << largestPath.second << std::endl;

  return largestPath;
}

/*
* TODO
*/
PathInfo GetRandomAggregatedMetricPath(PathList& paths, Ptr<ns3::UniformRandomVariable> randomVariable) {
   
  randomVariable->SetAttribute("Min", ns3::DoubleValue(0));
  randomVariable->SetAttribute("Max", ns3::DoubleValue(paths.size() - 1));

  // Return an empty path if the list is empty
  if (paths.empty()) {
      return PathInfo(); // Returns an empty PathInfo
  }

  // Get a random index
  int randomIndex = randomVariable->GetInteger();

  // Iterate to the random index
  auto pathIt = paths.begin();
  std::advance(pathIt, randomIndex);

  // Print the path and the aggregated metric
  std::cout << "Random Path: ";
  for (auto nodeIt = pathIt->first.begin(); nodeIt != pathIt->first.end(); ++nodeIt) {
      std::cout << (*nodeIt)->GetId();
      if (std::next(nodeIt) != pathIt->first.end()) {
          std::cout << " -> ";
      }
  }
  std::cout << " Metric: " << pathIt->second << std::endl;

  return *pathIt;
}

/*
* TODO
*/
PathInfo GetShortestNodeCountPath(PathList& paths) {

  PathInfo bestPath;
  // Initialize with a very large number, assuming no path will have this many nodes
  size_t minNodeCount = std::numeric_limits<size_t>::max();
  double maxMetric = std::numeric_limits<double>::lowest(); // To track the largest metric

  for (auto& pathInfo : paths) {
      size_t currentNodeCount = pathInfo.first.size();
      double currentMetric = pathInfo.second;

      // Check if the current path has fewer nodes or equal number of nodes but a larger metric
      if (currentNodeCount < minNodeCount || 
          (currentNodeCount == minNodeCount && currentMetric > maxMetric)) {
          bestPath = pathInfo;
          minNodeCount = currentNodeCount;
          maxMetric = currentMetric;
      }
  }

  std::cout << "Shortest Node Count Path: ";
  for (auto nodeIt = bestPath.first.begin(); nodeIt != bestPath.first.end(); ++nodeIt) {
      std::cout << (*nodeIt)->GetId();
      if (std::next(nodeIt) != bestPath.first.end()) {
          std::cout << " -> ";
      }
  }
  std::cout << " Node Count: " << bestPath.first.size() << ", Metric: " << bestPath.second << std::endl;

  return bestPath;
}

/*
* TODO
*/
std::map<uint32_t, uint32_t> nodeIdToL2IdMap;

/*
* TODO
*/
void FillNodeIdToL2IdMapping(const NodeContainer ueNodes) 
{

  for (uint32_t j = 0; j < ueNodes.GetN(); ++j) 
  {
    // Get the L2 ID
    Ptr<NrSlUeProse> prose = ueNodes.Get(j)->GetDevice (0)->GetObject<NrUeNetDevice>()->GetObject<NrSlUeProse>();
    uint32_t l2Id = prose->GetL2Id();

    // Get the Node ID
    uint32_t ueNodeId = ueNodes.Get(j)->GetId();

    // Update the mapping
    nodeIdToL2IdMap[ueNodeId] = l2Id;

  }
}

/*
* TODO
*/
uint32_t FindNodeIdFromL2Id(uint32_t l2Id) {
  for (const auto& pair : nodeIdToL2IdMap) {
      if (pair.second == l2Id) {
          return pair.first;
      }
  }
  return (uint32_t)-1; // Return an invalid ID if not found
}

/*
* TODO
*/
ConnectivityMap GetDiscoveryConnectivity(const NodeContainer ueNodes) {

  ConnectivityMap connectivityMap;
  for (uint32_t j = 0; j < ueNodes.GetN(); ++j) {
    Ptr<NrSlUeProse> prose = ueNodes.Get(j)->GetDevice (0)->GetObject<NrUeNetDevice>()->GetObject<NrSlUeProse>();
      auto rsrpMap = prose->GetRsrpMeasurementsMap();
      uint32_t currentNodeId = ueNodes.Get(j)->GetId(); // Get the Node ID of the current UE

      for (const auto& entry : rsrpMap) {
          uint32_t peerL2Id = entry.first;
          uint32_t peerNodeId = FindNodeIdFromL2Id(peerL2Id); // Find the Node ID corresponding to the peer L2 ID
          double rsrpValueDbm = entry.second.first;
          double rsrpValueMw = pow(10, rsrpValueDbm / 10);

          // Populate the ConnectivityMap
          connectivityMap[currentNodeId][peerNodeId] = rsrpValueMw;
      }
  }

  std::cout << "Full Connectivity: \n";
  for (const auto& node : connectivityMap) {
  std::cout << "Node " << node.first << " ";
    for (const auto& connection : node.second) {
        std::cout << "[" << connection.first << ", " << connection.second << "] ";
    }
  std::cout << std::endl;
  }
  return connectivityMap;
}

/*
* TODO
*/
void PrintRsrpMeasurementMap (const NodeContainer ueNodes)
{
  for (uint32_t j = 0; j < ueNodes.GetN (); ++j)
  {
    Ptr<NrSlUeProse> prose = ueNodes.Get(j)->GetDevice (0)->GetObject<NrUeNetDevice>()->GetObject<NrSlUeProse>();
    auto rsrpMap = prose->GetRsrpMeasurementsMap ();
    std::cout << "UE Node ID: " << ueNodes.Get(j)->GetId() << ", L2 ID: " << prose->GetL2Id() << std::endl;
    for (const auto& entry : rsrpMap) {
      uint32_t peerL2Id = entry.first;
      double rsrpValueDbm  = entry.second.first;
      double rsrpValueMw = pow(10, rsrpValueDbm / 10);

      std::cout << "-> Peer Node ID: " << FindNodeIdFromL2Id(peerL2Id) << ", Peer L2 ID: " << peerL2Id 
                << ", RSRP Value: " << rsrpValueDbm << " dBm / " << rsrpValueMw << " mW" << std::endl;
    }
  } 
}

/*
* TODO
*/
struct TrafficStats
{
  uint32_t nRxPkts;
  uint32_t nTxPkts;
  double lossRatio;
  double meanDelay;
  double meanJitter;
};

/*
* TODO
*/
struct FlowInfo
{
  uint16_t id;
  Ipv4Address srcIp;
  Ipv4Address tgtIp;
  uint16_t srcNodeId;
  uint16_t tgtNodeId;
  std::list<Ptr<NetDevice>> ueDevPath;
  uint16_t nHops;
  TrafficStats trafficStats; //To be filled at the end of the simulation
};

/*
* TODO
*/
void FindPathsFromDiscoveryInfo (Ptr<NrSlProseHelper> nrSlProseHelper, const NodeContainer ueNodes, const std::vector<uint32_t>& relayNodeIds,
                                std::map<uint16_t, FlowInfo> *flowsInfo, uint32_t maxNumRelays, double rsrpConstraint)

{

  FillNodeIdToL2IdMapping (ueNodes);
  PrintRsrpMeasurementMap(ueNodes);

  std::map<uint16_t, std::list<Ptr<ns3::Node>>> nodePaths; //Structure used for plotting and relay statistics

  // Step 1: Get Connectivity
  ConnectivityMap connectivityMap = GetDiscoveryConnectivity (ueNodes);
  // Step 2: Filter Connectivity by Metric Constraint
  ComparisonOp op = GREATER; // We are going to keep RSRPs that are greater than the rsrpConstraint
  ConnectivityMap filteredConnectivity = FilterConnectivity(connectivityMap, rsrpConstraint, op);
  for (auto &itFlowInfo : *flowsInfo)
  {
    std::list<Ptr<ns3::Node>> pathTest;

    uint32_t srcNodeId = itFlowInfo.second.srcNodeId; 
    uint32_t tgtNodeId = itFlowInfo.second.tgtNodeId;
    std::cout << "Path ID: "<< itFlowInfo.second.id << " srcNodeId: " << srcNodeId << " tgtNodeId: " << tgtNodeId << std::endl;

    // Step 3: Find All Paths Between src and dst
    int maxDepth = maxNumRelays + 2;  // Does this do the maxNumRelays filter already?
    PathList paths = FindAllPaths(ueNodes, filteredConnectivity, relayNodeIds, srcNodeId, tgtNodeId, maxDepth);
    // Step 4: Filter Paths by Number of Relays
    PathList filteredPaths = FilterPathsByRelays(paths, maxNumRelays);

    // Step 5: Get Path
    //PathInfo shortestPathInfo = GetLargestAggregatedMetricPath(filteredPaths);
    PathInfo shortestPathInfo = GetShortestNodeCountPath(filteredPaths);

    pathTest = shortestPathInfo.first;

    if (pathTest.size() > 0)
    {
      itFlowInfo.second.nHops = pathTest.size () - 1;
      g_foundNpaths++;
    }
    else
    {
      itFlowInfo.second.nHops = 0;
    }

    std::cout << "Path (srcIp: " << itFlowInfo.second.srcIp << ", tgtIp: " << itFlowInfo.second.tgtIp << ", Flow Id: " << itFlowInfo.second.id << ") " << itFlowInfo.second.nHops <<" hops:" << std::endl;
    std::cout << " NodeID\tIPv4Address\tposition(X,Y,Z): " <<std::endl;
    for (auto it = pathTest.begin (); it != pathTest.end (); ++it)
      {
        std::cout << " " << (*it)->GetId () << "\t" 
                  << (*it)->GetObject<Ipv4L3Protocol> ()->GetAddress (1, 0).GetLocal () << "\t" 
                  << "(" << (*it)->GetObject<MobilityModel> ()->GetPosition ().x <<","
                  << (*it)->GetObject<MobilityModel> ()->GetPosition ().y <<","
                  << (*it)->GetObject<MobilityModel> ()->GetPosition ().z <<")" 
                  <<std::endl;

        Ptr<NrUeNetDevice> netDevPtr = (*it)->GetDevice (0)->GetObject<NrUeNetDevice> ();
        
        itFlowInfo.second.ueDevPath.push_back (netDevPtr);
      }  
    nodePaths.emplace(itFlowInfo.second.id, pathTest); //Structure used for plotting and scenario 
  }
    
  ExportDataForPlotting (ueNodes, relayNodeIds, nodePaths, g_d); 

  //Same slInfo template for all paths at the moment
  SidelinkInfo slInfo;
  slInfo.m_castType = SidelinkInfo::CastType::Unicast;
  slInfo.m_harqEnabled = true;

  if (g_dynSch)
  {
    slInfo.m_dynamic = true;
  }
  else
  {
    slInfo.m_dynamic = false;
    slInfo.m_rri = MilliSeconds (20);
  }

  //Configure the paths
  std::cout << "Configuring paths... "<< std::endl;
  std::cout << "slInfo.m_harqEnabled: "<< slInfo.m_harqEnabled 
            << " slInfo.m_dynamic: "<< slInfo.m_dynamic 
            << " slInfo.m_rri: "<< slInfo.m_rri.GetMilliSeconds () << " ms" 
            << std::endl;

  for (auto &itFlowInfo : *flowsInfo)
  {    
    if (itFlowInfo.second.ueDevPath.size () > 0)
    {
      nrSlProseHelper->ConfigureU2uRelayPath (itFlowInfo.second.srcIp , itFlowInfo.second.tgtIp, itFlowInfo.first, slInfo, itFlowInfo.second.ueDevPath);
    }
  }

  std::cout << "Paths configured."<< std::endl;

  //Calculate relay stats
  std::vector<RelayStats> relayStats = GetRelayStats (relayNodeIds, nodePaths);
  std::cout << "Relay Stats:\nRelayNodeID\tPathCount\n";
  for (const auto& relayInfo : relayStats) 
  {
      std::cout << relayInfo.relayNodeId << "\t"
                << relayInfo.pathCount << std::endl;
  }

  for (const auto& relayInfo : relayStats) 
  {
    g_minPathCount = std::min(g_minPathCount, relayInfo.pathCount);
    g_maxPathCount = std::max(g_maxPathCount, relayInfo.pathCount);
    g_totalPathCount += relayInfo.pathCount;
    if (relayInfo.pathCount > 0) 
    {
          g_minNonZeroPathCount = std::min(g_minNonZeroPathCount, relayInfo.pathCount);
          g_totalNonZeroPathCount += relayInfo.pathCount;
          g_nonZeroCount++;
    }
  }

  if (!relayStats.empty()) 
  {
    g_meanPathCount = static_cast<double>(g_totalPathCount) / relayStats.size();

    std::cout << "Minimum Path Count: " << g_minPathCount << std::endl;
    std::cout << "Maximum Path Count: " << g_maxPathCount << std::endl;
    std::cout << "Mean Path Count: " << g_meanPathCount << std::endl;
    if (g_nonZeroCount > 0) 
    {
      g_meanNonZeroPathCount = static_cast<double>(g_totalNonZeroPathCount) / g_nonZeroCount;
      std::cout << "Minimum Path Count (Non-Zero): " << g_minNonZeroPathCount << std::endl;
      std::cout << "Mean Path Count (Non-Zero): " << g_meanNonZeroPathCount << std::endl;
    } else {
        std::cout << "No non-zero relay path information available." << std::endl;
    }
  } else 
  {
      std::cout << "No relay path information available." << std::endl;
  }

}

/*
* TODO
*/
void
TraceSinkDiscoveryRsrpTrace (Ptr<OutputStreamWrapper> stream,
                                   uint32_t thisL2Id, uint32_t peerL2Id, double rsrp)
{
 
  *stream->GetStream () << Simulator::Now ().GetMilliSeconds () / 1000.0 << "\t" << thisL2Id << "\t" << peerL2Id <<  "\t" << rsrp  << std::endl;

}

Time g_T5;
class NodePhyStats {
public:
  struct PhyStats {
      uint64_t nTxCtrl = 0;         // PSCCH
      uint64_t nTxData = 0;         // PSSCH
      uint64_t nRxCtrl = 0;         // Successfully decoded PSCCH
      uint64_t nRxData = 0;         // Successfully decoded PSSCH
      uint64_t nCorruptRxCtrl = 0;  // PSCCH corrupted (BLER)
      uint64_t nCorruptRxData = 0;  // PSSCH corrupted (BLER)
      uint64_t nHdRxCtrl = 0;       // PSCCH ignored due to half duplex
      uint64_t nHdRxData = 0;       // PSSCH ignored due to half duplex
      uint64_t nIgnoredData = 0;    // Not expected (PSCCH not received)
  };

  uint32_t m_nodeId;
  uint32_t m_l2Id;

  PhyStats m_stats;
};

/*
* TODO
*/
void
TraceSlPscchTx (NodePhyStats *stats, const Time duration)
{
  if (Simulator::Now() > g_T5)
  {
    stats->m_stats.nTxCtrl ++;
  }
}

/*
* TODO
*/
void
TraceSlPsschTx (NodePhyStats *stats, const Time duration)
{
  if (Simulator::Now() > g_T5)
  {
    stats->m_stats.nTxData ++;
  }
}

/*
* TODO
*/
void
TraceSlPscchRx (NodePhyStats *stats, const SlRxCtrlPacketTraceParams pscchStatsParams)
{
  
  if (Simulator::Now() > g_T5)
  {
    //Only collect stats for traffic direceted to this UE 
    if (pscchStatsParams.m_dstL2Id == stats->m_l2Id)
    {
      if (pscchStatsParams.m_corrupt)
      {
          stats->m_stats.nCorruptRxCtrl ++;
      }
      else
      {
        stats->m_stats.nRxCtrl ++;
      }
    }  
  }

}

/*
* TODO
*/
void
TraceSlPsschRx (NodePhyStats *stats, const SlRxDataPacketTraceParams psschStatsParams)
{
  if (Simulator::Now() > g_T5)
  {
    //Only collect stats for traffic direceted to this UE 
    if (psschStatsParams.m_dstL2Id == stats->m_l2Id)
    {
      if (psschStatsParams.m_corrupt)
      {
        stats->m_stats.nCorruptRxData ++;
      }
      else
      {
        stats->m_stats.nRxData ++;
      }
    }  
  }
}

/*
* TODO
*/
void
TraceSlPscchRxHd (NodePhyStats *stats, uint64_t oldValue, uint64_t newValue)
{
  if (Simulator::Now() > g_T5)
  {
    stats->m_stats.nHdRxCtrl ++;
  }
}

/*
* TODO
*/
void
TraceSlPsschRxHd (NodePhyStats *stats, uint64_t oldValue, uint64_t newValue)
{
  if (Simulator::Now() > g_T5)
  {
    stats->m_stats.nHdRxData ++;
  }
}

/*
* TODO
*/
void
TraceSlPsschRxIg (NodePhyStats *stats, uint64_t oldValue, uint64_t newValue)
{
  if (Simulator::Now() > g_T5)
  {
    stats->m_stats.nIgnoredData ++;
  }
}

int
main (int argc, char *argv[])
{
  //Topology parameters
  uint16_t nUes = 3;
  uint16_t d = 60; //meters
  std::string deployment = "Linear";
  double relayDensity = 0.5;
  uint16_t nPaths = 2;
  std::string pathFindingAlgo = "DistShortest";
  bool allowReverse = false; 
  bool allowRelayEndNode = false;
  uint32_t maxNumRelays = 4; // Example maximum number of relays
  double distanceConstraint =  60.0; //m
  double rsrpConstraint =  1e-20; //mW
  bool prioToSps = false;

  // Traffic parameters
  uint32_t udpPacketSize = 60; //bytes
  double dataRate = 24; //kbps
  bool bidirectional = false;

  // Simulation timeline parameters
  Time trafficDuration = Seconds (30.0); //Total simulation time
  Time startTrafficTime = Seconds (3.0); //Time to start the traffic in the application layer. Will be adapted if discovery is used

  // NR parameters
  uint16_t numerologyBwpSl = 0; //The numerology to be used in sidelink bandwidth part
  double centralFrequencyBandSl = 5.89e9; // band n47  TDD //Here band is analogous to channel
  double bandwidthBandSl = 40e6; //40 MHz
  double txPower = 23; //dBm
  uint32_t mcs = 5; 

  // SL parameters
  bool sensing = false;
  uint32_t maxNumTx = 1;
  bool dynSch = true;

  CommandLine cmd;
  cmd.AddValue ("nUes", "The number of UEs in the topology", nUes);
  cmd.AddValue ("deployment", "The type of deployment. Options: Linear/SqrRand/SqrGrid", deployment);
  cmd.AddValue ("d", "The side of the square when using SqrGrid or SqrRand topologies, or the distance between first and last UE when using linear topology", d);
  cmd.AddValue ("relayDensity", "The relay density", relayDensity);
  cmd.AddValue ("nPaths", "The number of path to create", nPaths);
  cmd.AddValue ("pathFindingAlgo", "The algorithm used to find the paths. Options: DistShortest/DistRandom/DiscoveryShortest", pathFindingAlgo);
  cmd.AddValue ("distanceConstraint", "When distance is use for path finding: the coverage radius of each UE", distanceConstraint);
  cmd.AddValue ("rsrpConstraint", "When discovery is use for path finding: the min RSRP to consider connectivity", rsrpConstraint);
  cmd.AddValue ("maxNumRelays", "The maximum number of relays allowed in a path", maxNumRelays);
  cmd.AddValue ("allowReverse", "If true src->tgt and tgt->src are different paths considered sperately, if false only the first one is considered", allowReverse);
  cmd.AddValue ("allowRelayEndNode", "If true relays can be src or tgt in paths", allowRelayEndNode);
  cmd.AddValue ("packetSizeBe", "packet size in bytes to be used by the traffic", udpPacketSize);
  cmd.AddValue ("dataRate", "The data rate in kilobits per second", dataRate);
  cmd.AddValue ("bidirectional","If false, only a flow in the stc->tgt is installed. If true, two traffic flows are installed, one on each direction. ", bidirectional);
  cmd.AddValue ("mu", "The numerology to be used in sidelink bandwidth part", numerologyBwpSl);
  cmd.AddValue ("sensing", "If true, it enables the sensing based resource selection for SL, otherwise, no sensing is applied", sensing);
  cmd.AddValue ("maxNumTx", "The maximum numeber of transmissions in the PSSCH", maxNumTx);
  cmd.AddValue ("dynSch", "If True, all flows use dynamic scheduling. If false, all flows use SPS", dynSch);
  cmd.AddValue ("mcs", "SL MCS", mcs);


  // Parse the command line
  cmd.Parse (argc, argv);

  g_d = d;
  g_dynSch = dynSch;
  Time stopTrafficTime = startTrafficTime + trafficDuration;

  if (! (pathFindingAlgo == "DistShortest" || pathFindingAlgo == "DistRandom" || pathFindingAlgo == "DiscoveryShortest"))
  {
     NS_FATAL_ERROR ("Path finding algorithm not recognized");
  }
  //Check if the frequency is in the allowed range.
  NS_ABORT_IF (centralFrequencyBandSl > 6e9);

  //Setup large enough buffer size to avoid overflow
  Config::SetDefault ("ns3::LteRlcUm::MaxTxBufferSize", UintegerValue (999999999));
  //Setup -t-Reordering to 0ms to avoid RLC reordering delay
  Config::SetDefault ("ns3::LteRlcUm::ReorderingTimer", TimeValue (MilliSeconds (0)));

  // Discovery Frequency
  Time discInterval = Seconds (0.02); // Interval between two discovery announcements
  Config::SetDefault ("ns3::NrSlUeProse::DiscoveryInterval",TimeValue (discInterval));

  //UE nodes creation
  NodeContainer ueNodes;
  ueNodes.Create (nUes);

  //UE nodes mobility setup
  MobilityHelper mobility;
  mobility.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
  Ptr<ListPositionAllocator> listPositionAllocator = CreateObject<ListPositionAllocator> ();

  //Linear
  if (deployment == "Linear")
  {
    double interUeDistance = d / nUes;
    for (uint16_t i = 0; i < nUes; i++)
      {
        listPositionAllocator->Add (Vector (interUeDistance * i, 0.0, 1.5));
      }
    mobility.SetPositionAllocator (listPositionAllocator);
    mobility.Install (ueNodes);
  }
  else if (deployment ==  "SqrRand")
  {
    Ptr<RandomBoxPositionAllocator> sqrRandPositionAllocator = CreateObject<RandomBoxPositionAllocator> ();
    sqrRandPositionAllocator->SetAttribute ("X", StringValue ("ns3::UniformRandomVariable[Min=0.0|Max=" + std::to_string(d) + "]"));
    sqrRandPositionAllocator->SetAttribute ("Y", StringValue ("ns3::UniformRandomVariable[Min=0.0|Max=" + std::to_string(d) + "]"));
    sqrRandPositionAllocator->SetAttribute ("Z", StringValue ("ns3::ConstantRandomVariable[Constant=1.5]"));
    sqrRandPositionAllocator->AssignStreams (70000); //Fix stream to have same results for each Run indepently of previous random variables creations

    mobility.SetPositionAllocator (sqrRandPositionAllocator);
    mobility.Install (ueNodes);
  }
  else if (deployment ==  "SqrGrid")
  {
    double interUeDistance = d / (std::ceil(std::sqrt(nUes)) - 1 );

    Ptr<GridPositionAllocator>  sqrGridPositionAllocator;
    sqrGridPositionAllocator = CreateObject<GridPositionAllocator> ();
    sqrGridPositionAllocator->SetAttribute ("GridWidth", UintegerValue (std::ceil(std::sqrt(nUes))));
    sqrGridPositionAllocator->SetAttribute ("MinX", DoubleValue (0.0));
    sqrGridPositionAllocator->SetAttribute ("MinY", DoubleValue (0.0));
    sqrGridPositionAllocator->SetAttribute ("Z", DoubleValue (1.5));
    sqrGridPositionAllocator->SetAttribute ("DeltaX", DoubleValue (interUeDistance));
    sqrGridPositionAllocator->SetAttribute ("DeltaY", DoubleValue (interUeDistance));
    sqrGridPositionAllocator->SetAttribute ("LayoutType", EnumValue (GridPositionAllocator::ROW_FIRST));
    mobility.SetPositionAllocator (sqrGridPositionAllocator);
    mobility.Install (ueNodes);
  }
  else 
  {
    NS_FATAL_ERROR ("Deployment type not recognized");
  }

  /*
   * Setup the NR module. We create the various helpers needed for the
   * NR simulation:
   * - EpcHelper, which will setup the core network entities
   * - NrHelper, which takes care of creating and connecting the various
   * part of the NR stack
   */
  Ptr<NrPointToPointEpcHelper> epcHelper = CreateObject<NrPointToPointEpcHelper> ();
  Ptr<NrHelper> nrHelper = CreateObject<NrHelper> ();

  // Put the pointers inside nrHelper
  nrHelper->SetEpcHelper (epcHelper);

  /*
   * Spectrum division. We create one operational band, containing
   * one component carrier, and a single bandwidth part
   * centered at the frequency specified by the input parameters.
   */
  BandwidthPartInfoPtrVector allBwps;
  CcBwpCreator ccBwpCreator;
  const uint8_t numCcPerBand = 1;

  //Create the configuration for the CcBwpHelper. SimpleOperationBandConfcreates a single BWP per CC
  CcBwpCreator::SimpleOperationBandConf bandConfSl (centralFrequencyBandSl, bandwidthBandSl,
  //                                                  numCcPerBand, BandwidthPartInfo::V2V_Highway);
                                                    numCcPerBand, BandwidthPartInfo::RMa_LoS);

  // By using the configuration created, it is time to make the operation bands
  OperationBandInfo bandSl = ccBwpCreator.CreateOperationBandContiguousCc (bandConfSl);

  //Configure 3GPP channel model
  Config::SetDefault ("ns3::ThreeGppChannelModel::UpdatePeriod", TimeValue (MilliSeconds (10)));
  nrHelper->SetChannelConditionModelAttribute ("UpdatePeriod", TimeValue (MilliSeconds (0)));
  nrHelper->SetPathlossAttribute ("ShadowingEnabled", BooleanValue (false));

  /*
   * Initialize channel and pathloss, plus other things inside bandSl. If needed,
   * the band configuration can be done manually, but we leave it for more
   * sophisticated examples. For the moment, this method will take care
   * of all the spectrum initialization needs.
   */
  auto bandMask = NrHelper::INIT_PROPAGATION | NrHelper::INIT_CHANNEL;
  bandMask |= NrHelper::INIT_FADING;

  nrHelper->InitializeOperationBand (&bandSl, bandMask);
  allBwps = CcBwpCreator::GetAllBwps ({bandSl});

  Packet::EnableChecking ();
  Packet::EnablePrinting ();


  /*
   * Antennas for all the UEs
   * We are not using beamforming in SL, rather we are using
   * quasi-omnidirectional transmission and reception, which is the default
   * configuration of the beams.
   */
  nrHelper->SetUeAntennaAttribute ("NumRows", UintegerValue (1));
  nrHelper->SetUeAntennaAttribute ("NumColumns", UintegerValue (2));
  nrHelper->SetUeAntennaAttribute ("AntennaElement",
                                   PointerValue (CreateObject<IsotropicAntennaModel> ()));

  nrHelper->SetUePhyAttribute ("TxPower", DoubleValue (txPower));

  //NR Sidelink attribute of UE MAC, which are would be common for all the UEs
  nrHelper->SetUeMacTypeId(NrSlUeMac::GetTypeId());
  nrHelper->SetUeMacAttribute ("EnableSensing", BooleanValue (sensing));
  nrHelper->SetUeMacAttribute ("T1", UintegerValue (2));

   //Resource selection window length (given by T1 and T2) depend on the numerology
  uint16_t t2;
  switch (numerologyBwpSl)
  {
    case 0:
      //t2 = 33; // with T1 = 2, this gives a window of 32 slots with mu = 0 => 32 ms
      t2 = 17; // with T1 = 2, this gives a window of 16 slots with mu = 0 => 16 ms

      break;
    case 1:
      t2 = 33; // with T1 = 2, this gives a window of 32 slots with mu = 1 => 16 ms

      break;
    case 2:
      //t2 = 33; // with T1 = 2, this gives a window of 32 slots with mu = 2 => 8 ms
      t2 = 65; // with T1 = 2, this gives a window of 64 slots with mu = 2 => 16 ms
      break;

    default:
      NS_FATAL_ERROR ("Numerology for sidelink not recognized");
      break;
  }
  NS_LOG_DEBUG ("T2: " << t2 );
  nrHelper->SetUeMacAttribute ("T2", UintegerValue (t2));
  nrHelper->SetUeMacAttribute ("ActivePoolId", UintegerValue (0));
  nrHelper->SetUeMacAttribute ("SlThresPsschRsrp", IntegerValue (-128));

  uint8_t bwpIdForGbrMcptt = 0;

  nrHelper->SetBwpManagerTypeId (TypeId::LookupByName ("ns3::NrSlBwpManagerUe"));
  nrHelper->SetUeBwpManagerAlgorithmAttribute ("GBR_MC_PUSH_TO_TALK",
                                               UintegerValue (bwpIdForGbrMcptt));

  std::set<uint8_t> bwpIdContainer;
  bwpIdContainer.insert (bwpIdForGbrMcptt);

  NetDeviceContainer uesNetDev = nrHelper->InstallUeDevice (ueNodes, allBwps);

  // When all the configuration is done, explicitly call UpdateConfig ()
  for (auto it = uesNetDev.Begin (); it != uesNetDev.End (); ++it)
    {
      DynamicCast<NrUeNetDevice> (*it)->UpdateConfig ();
    }

  /*Create NrSlHelper which will configure the UEs protocol stack to be ready to
   *perform Sidelink related procedures
   */
  Ptr<NrSlHelper> nrSlHelper = CreateObject<NrSlHelper> ();
  // Put the pointers inside NrSlHelper
  nrSlHelper->SetEpcHelper (epcHelper);

  /*
   * Set the SL error model and AMC
   * Error model type: ns3::NrEesmCcT1, ns3::NrEesmCcT2, ns3::NrEesmIrT1,
   *                   ns3::NrEesmIrT2, ns3::NrLteMiErrorModel
   * AMC type: NrAmc::ShannonModel or NrAmc::ErrorModel
   */
  std::string errorModel = "ns3::NrEesmIrT1";
  nrSlHelper->SetSlErrorModel (errorModel);
  nrSlHelper->SetUeSlAmcAttribute ("AmcModel", EnumValue (NrAmc::ErrorModel));

  /*
   * Set the SL scheduler attributes
   * In this example we use NrSlUeMacSchedulerDefault scheduler, which uses
   * fix MCS value and schedules logical channels by priority order first and
   * then by creation order
   */
  nrSlHelper->SetNrSlSchedulerTypeId (NrSlUeMacSchedulerFixedMcs::GetTypeId ());
  nrSlHelper->SetUeSlSchedulerAttribute("Mcs", UintegerValue(14));
  nrSlHelper->SetUeSlSchedulerAttribute("PriorityToSps", BooleanValue(prioToSps));

  /*
   * Very important method to configure UE protocol stack, i.e., it would
   * configure all the SAPs among the layers, setup callbacks, configure
   * error model, configure AMC, and configure ChunkProcessor in Interference
   * API.
   */
  nrSlHelper->PrepareUeForSidelink (uesNetDev, bwpIdContainer);

  /*
   * Start preparing for all the sub Structs/RRC Information Element (IEs)
   * of LteRrcSap::SidelinkPreconfigNr. This is the main structure, which would
   * hold all the pre-configuration related to Sidelink.
   */

  //SlResourcePoolNr IE
  LteRrcSap::SlResourcePoolNr slResourcePoolNr;
  //get it from pool factory
  Ptr<NrSlCommResourcePoolFactory> ptrFactory = Create<NrSlCommResourcePoolFactory> ();
  //Configure specific parameters of interest:
  std::vector<std::bitset<1>> slBitmap = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  ptrFactory->SetSlTimeResources (slBitmap);
  ptrFactory->SetSlSensingWindow (100); // T0 in ms
  ptrFactory->SetSlSelectionWindow (5);
  ptrFactory->SetSlFreqResourcePscch (10); // PSCCH RBs
  ptrFactory->SetSlSubchannelSize (10);
  ptrFactory->SetSlMaxNumPerReserve (3);
  std::list<uint16_t> resourceReservePeriodList = {0, 20}; // in ms
  ptrFactory->SetSlResourceReservePeriodList (resourceReservePeriodList);
  //Once parameters are configured, we can create the pool
  LteRrcSap::SlResourcePoolNr pool = ptrFactory->CreatePool ();
  slResourcePoolNr = pool;

  //Configure the SlResourcePoolConfigNr IE, which holds a pool and its id
  LteRrcSap::SlResourcePoolConfigNr slresoPoolConfigNr;
  slresoPoolConfigNr.haveSlResourcePoolConfigNr = true;
  //Pool id, ranges from 0 to 15
  uint16_t poolId = 0;
  LteRrcSap::SlResourcePoolIdNr slResourcePoolIdNr;
  slResourcePoolIdNr.id = poolId;
  slresoPoolConfigNr.slResourcePoolId = slResourcePoolIdNr;
  slresoPoolConfigNr.slResourcePool = slResourcePoolNr;

  //Configure the SlBwpPoolConfigCommonNr IE, which holds an array of pools
  LteRrcSap::SlBwpPoolConfigCommonNr slBwpPoolConfigCommonNr;
  //Array for pools, we insert the pool in the array as per its poolId
  slBwpPoolConfigCommonNr.slTxPoolSelectedNormal[slResourcePoolIdNr.id] = slresoPoolConfigNr;

  //Configure the BWP IE
  LteRrcSap::Bwp bwp;
  bwp.numerology = numerologyBwpSl;
  bwp.symbolsPerSlots = 14;
  bwp.rbPerRbg = 1;
  bwp.bandwidth = bandwidthBandSl/1000/100; // SL configuration requires BW in Multiple of 100 KHz;

  //Configure the SlBwpGeneric IE
  LteRrcSap::SlBwpGeneric slBwpGeneric;
  slBwpGeneric.bwp = bwp;
  slBwpGeneric.slLengthSymbols = LteRrcSap::GetSlLengthSymbolsEnum (14);
  slBwpGeneric.slStartSymbol = LteRrcSap::GetSlStartSymbolEnum (0);

  //Configure the SlBwpConfigCommonNr IE
  LteRrcSap::SlBwpConfigCommonNr slBwpConfigCommonNr;
  slBwpConfigCommonNr.haveSlBwpGeneric = true;
  slBwpConfigCommonNr.slBwpGeneric = slBwpGeneric;
  slBwpConfigCommonNr.haveSlBwpPoolConfigCommonNr = true;
  slBwpConfigCommonNr.slBwpPoolConfigCommonNr = slBwpPoolConfigCommonNr;

  //Configure the SlFreqConfigCommonNr IE, which holds the array to store
  //the configuration of all Sidelink BWP (s).
  LteRrcSap::SlFreqConfigCommonNr slFreConfigCommonNr;
  //Array for BWPs. Here we will iterate over the BWPs, which
  //we want to use for SL.
  for (const auto &it : bwpIdContainer)
    {
      // it is the BWP id
      slFreConfigCommonNr.slBwpList[it] = slBwpConfigCommonNr;
    }

  //Configure the TddUlDlConfigCommon IE
  LteRrcSap::TddUlDlConfigCommon tddUlDlConfigCommon;
  tddUlDlConfigCommon.tddPattern = "DL|DL|DL|F|UL|UL|UL|UL|UL|UL|";

  //Configure the SlPreconfigGeneralNr IE
  LteRrcSap::SlPreconfigGeneralNr slPreconfigGeneralNr;
  slPreconfigGeneralNr.slTddConfig = tddUlDlConfigCommon;

  //Configure the SlUeSelectedConfig IE
  LteRrcSap::SlUeSelectedConfig slUeSelectedPreConfig;
  slUeSelectedPreConfig.slProbResourceKeep = 0;
  //Configure the SlPsschTxParameters IE
  LteRrcSap::SlPsschTxParameters psschParams;
  psschParams.slMaxTxTransNumPssch = maxNumTx;
  //Configure the SlPsschTxConfigList IE
  LteRrcSap::SlPsschTxConfigList pscchTxConfigList;
  pscchTxConfigList.slPsschTxParameters[0] = psschParams;
  slUeSelectedPreConfig.slPsschTxConfigList = pscchTxConfigList;

  /*
   * Finally, configure the SidelinkPreconfigNr This is the main structure
   * that needs to be communicated to NrSlUeRrc class
   */
  LteRrcSap::SidelinkPreconfigNr slPreConfigNr;
  slPreConfigNr.slPreconfigGeneral = slPreconfigGeneralNr;
  slPreConfigNr.slUeSelectedPreConfig = slUeSelectedPreConfig;
  slPreConfigNr.slPreconfigFreqInfoList[0] = slFreConfigCommonNr;

  //Communicate the above pre-configuration to the NrSlHelper
  nrSlHelper->InstallNrSlPreConfiguration (uesNetDev, slPreConfigNr);

  /****************************** End SL Configuration ***********************/

  /*
   * Fix the random streams
   */
  int64_t stream = 1;
  nrHelper->AssignStreams (uesNetDev, stream);
  nrSlHelper->AssignStreams (uesNetDev, stream + 3000);
 
  /*
   * Configure the IPv4 stack
   */
  InternetStackHelper internet;
  internet.Install (ueNodes);
  Ipv4InterfaceContainer ueIpIface;
  ueIpIface = epcHelper->AssignUeIpv4Address (uesNetDev);

  std::vector<Ipv4Address> ipv4AddressVector (nUes);

  // set the default gateway for the UE
  Ipv4StaticRoutingHelper ipv4RoutingHelper;
  for (uint32_t u = 0; u < ueNodes.GetN (); ++u)
    {
      Ptr<Node> ueNode = ueNodes.Get (u);
      // Set the default gateway for the UE
      Ptr<Ipv4StaticRouting> ueStaticRouting =
          ipv4RoutingHelper.GetStaticRouting (ueNode->GetObject<Ipv4> ());
      ueStaticRouting->SetDefaultRoute (epcHelper->GetUeDefaultGatewayAddress (), 1);

      //Obtain local IPv4 addresses that will be used to route the unicast traffic upon setup of the direct link
      ipv4AddressVector[u] =
          ueNodes.Get (u)->GetObject<Ipv4L3Protocol> ()->GetAddress (1, 0).GetLocal ();
    }

  /*
   * Configure ProSe
   */

  //Create ProSe helper
  Ptr<NrSlProseHelper> nrSlProseHelper = CreateObject<NrSlProseHelper> ();
  // Install ProSe layer and corresponding SAPs in the UEs
  nrSlProseHelper->PrepareUesForProse (uesNetDev);

  //Selecting ramson relays, finding paths, configuring them, etc. 
  std::map<uint16_t, FlowInfo> flowsInfo;
  std::cout << "Selecting random relays... "<< std::endl;

  //Random variable for selecting relays
  Ptr<ns3::UniformRandomVariable> randomVariableRelays = CreateObject<ns3::UniformRandomVariable>();
  randomVariableRelays->SetStream (70000); //Fix stream to have same results for each Run indepently of previous random variables creations

  std::vector<uint32_t> selectedRelays = SelectRandomRelays(ueNodes, std::ceil(nUes * relayDensity), randomVariableRelays);
  std::cout << "Selected Relays Node IDs: ";
  for (const auto& id : selectedRelays) {
        std::cout << id << " ";
    }
    std::cout << std::endl;
  
  std::map<uint16_t, std::list<Ptr<ns3::Node>>> nodePaths; //Structure used for plotting and relay statistics
  std::set<std::pair<uint32_t, uint32_t>> usedNodePairs;

  std::cout << "Generating src->tgt pairs (nPaths=" << nPaths << ")... "<< std::endl;

  std::cout << "*****************************************"<< std::endl;
  ConnectivityMap connectivityMap, filteredConnectivity;
  
  if (pathFindingAlgo == "DistShortest" || pathFindingAlgo == "DistRandom")
  {
    // Step 1: Get Connectivity
    connectivityMap = GetDistanceConnectivity(ueNodes);
    // Step 2: Filter Connectivity by Metric Constraint
    ComparisonOp op = LESS_OR_EQUAL; // Example operation
    filteredConnectivity = FilterConnectivity(connectivityMap, distanceConstraint, op);
   }
  uint16_t port = 8000;

  //Random variable for getting end nodes
  Ptr<ns3::UniformRandomVariable> randomVariableNodes = CreateObject<ns3::UniformRandomVariable>();
  randomVariableNodes->SetStream (70001); //Fix stream to have same results for each Run indepently of previous random variables creations

  //Random variable for random paths
  ns3::Ptr<ns3::UniformRandomVariable> randomVariablePaths = ns3::CreateObject<ns3::UniformRandomVariable>();
  randomVariablePaths->SetStream (70002); //Fix stream to have same results for each Run indepently of previous random variables creations

  for (uint16_t i = 0; i < nPaths; i++)
  {
    uint32_t srcNodeId, tgtNodeId;
    do
    {
      if (allowRelayEndNode)
      {
        srcNodeId = GetRandomNodeId (ueNodes, randomVariableNodes);
        tgtNodeId = GetRandomNodeId (ueNodes, randomVariableNodes);
      }
      else 
      {
        srcNodeId = GetRandomNodeIdExcludeRelays(ueNodes, selectedRelays, randomVariableNodes);
        tgtNodeId = GetRandomNodeIdExcludeRelays(ueNodes, selectedRelays, randomVariableNodes);
      }
    } 
    while (srcNodeId == tgtNodeId || usedNodePairs.find(std::make_pair(srcNodeId, tgtNodeId)) != usedNodePairs.end());

    
    //Configuring flow info
    FlowInfo fInfo;
    fInfo.id = i;
    fInfo.srcNodeId = srcNodeId;
    fInfo.tgtNodeId = tgtNodeId;
    fInfo.srcIp = GetNodeById(ueNodes, srcNodeId)->GetObject<Ipv4L3Protocol> ()->GetAddress (1, 0).GetLocal ();
    fInfo.tgtIp = GetNodeById(ueNodes, tgtNodeId)->GetObject<Ipv4L3Protocol> ()->GetAddress (1, 0).GetLocal ();

    std::cout << "Path ID: "<< fInfo.id << " srcNodeId: " << srcNodeId << " tgtNodeId: " << tgtNodeId << std::endl;


    if (pathFindingAlgo == "DistShortest" || pathFindingAlgo == "DistRandom")
    {
      //Finding the path
      std::list<Ptr<ns3::Node>> pathTest;

      // Step 3: Find All Paths Between src and dst
      int maxDepth = maxNumRelays + 2 ; // Does this do the maxNumRelays filter already?
      PathList paths = FindAllPaths(ueNodes, filteredConnectivity, selectedRelays, srcNodeId, tgtNodeId, maxDepth);
      // Step 4: Filter Paths by Number of Relays
      PathList filteredPaths = FilterPathsByRelays(paths, maxNumRelays);
      
      if (pathFindingAlgo == "DistShortest")
      {
        // Step 5: Get Shortest Path
        PathInfo shortestPathInfo = GetSmallestAggregatedMetricPath(filteredPaths);
        pathTest = shortestPathInfo.first;

      }
      else if (pathFindingAlgo == "DistRandom")
      {
        PathInfo shortestPathInfo = GetRandomAggregatedMetricPath(filteredPaths, randomVariablePaths);
        pathTest = shortestPathInfo.first;
      }

      if (pathTest.size() > 0)
      {
          g_foundNpaths ++;
      }
      else
      {
        std::cout << "No path found" << std::endl;
      }
  
      if (pathTest.size() > 0)
      {
        fInfo.nHops = pathTest.size () - 1;
      }
      else
      {
        fInfo.nHops = 0;
      }

      std::cout << "Path (srcIp: " << fInfo.srcIp << ", tgtIp: " << fInfo.tgtIp << ", Flow Id: " << fInfo.id << ") " << fInfo.nHops <<" hops:" << std::endl;
      std::cout << " NodeID\tIPv4Address\tposition(X,Y,Z): " <<std::endl;
      for (auto it = pathTest.begin (); it != pathTest.end (); ++it)
        {
          std::cout << " " << (*it)->GetId () << "\t" 
                    << (*it)->GetObject<Ipv4L3Protocol> ()->GetAddress (1, 0).GetLocal () << "\t" 
                    << "(" << (*it)->GetObject<MobilityModel> ()->GetPosition ().x <<","
                    << (*it)->GetObject<MobilityModel> ()->GetPosition ().y <<","
                    << (*it)->GetObject<MobilityModel> ()->GetPosition ().z <<")" 
                    <<std::endl;

          Ptr<NrUeNetDevice> netDevPtr = (*it)->GetDevice (0)->GetObject<NrUeNetDevice> ();
          
          fInfo.ueDevPath.push_back (netDevPtr);
        }  
      nodePaths.emplace(fInfo.id, pathTest); //Structure used for plotting and scenario 
    }
    
    //Regardeless of the path finding algo used (distance of discovery) 
    //we store flow info used for routing and traffic configuration
    flowsInfo.emplace(port, fInfo); 
    port++;

    //Keep track of the src->dst pair we use so that we don't repeat a path
    usedNodePairs.insert(std::make_pair(srcNodeId, tgtNodeId));
    if (!allowReverse) 
    {
      usedNodePairs.insert(std::make_pair(tgtNodeId, srcNodeId));
    }

  }

  if (pathFindingAlgo == "DistShortest" || pathFindingAlgo == "DistRandom")
  {
    std::cout << "\nFound " << g_foundNpaths << "/" << nPaths << " paths\n\n";

    //Same slInfo template for all paths at the moment
    SidelinkInfo slInfo;
    slInfo.m_castType = SidelinkInfo::CastType::Unicast;
    slInfo.m_harqEnabled = true;
    if (dynSch)
    {
      slInfo.m_dynamic = true;
    }
    else
    {
      slInfo.m_dynamic = false;
      slInfo.m_rri= MilliSeconds(20);
    }

    //Configure the paths
    std::cout << "Configuring paths... "<< std::endl;
    std::cout << "slInfo.m_harqEnabled: "<< slInfo.m_harqEnabled 
              << " slInfo.m_dynamic: "<< slInfo.m_dynamic 
              << " slInfo.m_rri: "<< slInfo.m_rri.GetMilliSeconds () << " ms" 
              << std::endl;
    for (auto itPath = flowsInfo.begin (); itPath != flowsInfo.end (); ++itPath)
      {
        if (itPath->second.ueDevPath.size () > 0)
        {
          nrSlProseHelper->ConfigureU2uRelayPath (itPath->second.srcIp , itPath->second.tgtIp, itPath->first, slInfo, itPath->second.ueDevPath);
        }
      }
  }

  if (pathFindingAlgo == "DiscoveryShortest")
  {
    Time startMonDiscTime = Seconds (2.0); //Time to start the Prose discovery procedure in seconds
    Time stopMonDiscTime = Seconds (2.0 + uesNetDev.GetN ()*1.0 + 1.0); //Time to stop the Prose discovery procedure in seconds

  // Application Codes to announce/monitor
    std::map<Ptr<NetDevice>, std::list<uint32_t> > announcePayloads; 
    std::map<Ptr<NetDevice>, std::list<uint32_t> > monitorPayloads;

    // Destination L2 IDs for each application code
    std::map<Ptr<NetDevice>, std::list<uint32_t> > announceDstL2IdsMap;
    std::map<Ptr<NetDevice>, std::list<uint32_t> > monitorDstL2IdsMap;

    for (uint32_t i = 1; i <= uesNetDev.GetN (); ++i)
      {
        //For each UE, announce one appCode and monitor all the others appCode
        announcePayloads[uesNetDev.Get (i - 1)].push_back (i);
        announceDstL2IdsMap[uesNetDev.Get (i - 1)].push_back (100*i);

        for (uint32_t j = 1; j <= uesNetDev.GetN (); ++j)
          {
            if (i != j)
              {
                monitorPayloads[uesNetDev.Get (i - 1)].push_back (j);
                monitorDstL2IdsMap[uesNetDev.Get (i - 1)].push_back (100*j);
              }
          }

        //IMPORTANT!
        Ptr<LteUeRrc> remoteRrc = uesNetDev.Get(i - 1)->GetObject <NrUeNetDevice>()->GetRrc ();
        remoteRrc->EnableUeSlRsrpMeasurements ();     
      }

    for (uint32_t i = 0; i < uesNetDev.GetN (); ++i)
      {
        //Start announcing
        Simulator::Schedule (startMonDiscTime + Seconds(i*1.0), &NrSlProseHelper::StartDiscovery, nrSlProseHelper, uesNetDev.Get (i), 
            announcePayloads[uesNetDev.Get (i)], announceDstL2IdsMap[uesNetDev.Get (i)], NrSlUeProse::Announcing);

        //Stop announcing
        Simulator::Schedule (startMonDiscTime + Seconds(i*1.0 + 1.0) ,&NrSlProseHelper::StopDiscovery, nrSlProseHelper, uesNetDev.Get (i), 
                              announcePayloads[uesNetDev.Get (i)], NrSlUeProse::Announcing);
        //Start monitoring
        Simulator::Schedule (startMonDiscTime, &NrSlProseHelper::StartDiscovery, nrSlProseHelper, uesNetDev.Get (i), 
            monitorPayloads[uesNetDev.Get (i)], monitorDstL2IdsMap[uesNetDev.Get (i)], NrSlUeProse::Monitoring);
        //Stop monitoring
        Simulator::Schedule (stopMonDiscTime,&NrSlProseHelper::StopDiscovery, nrSlProseHelper, uesNetDev.Get (i), 
                              monitorPayloads[uesNetDev.Get (i)], NrSlUeProse::Monitoring);
      }

    Time T4 = stopMonDiscTime + Seconds(1.0);
    //Simulator::Schedule(delay, &PrintRsrpMeasurementMap, uesNetDev);
    Simulator::Schedule(T4, &FindPathsFromDiscoveryInfo, nrSlProseHelper, ueNodes, selectedRelays, &flowsInfo, maxNumRelays, rsrpConstraint);

    startTrafficTime = T4 + Seconds (1.0); //Time to start the traffic in the application layer
    g_T5 = startTrafficTime;

    stopTrafficTime = startTrafficTime + trafficDuration;

    //Configure L3 filtering. 
    LteRrcSap::SlRemoteUeConfig slRemoteConfig;
    slRemoteConfig.slReselectionConfig.slRsrpThres = rsrpConstraint;
    slRemoteConfig.slReselectionConfig.slFilterCoefficientRsrp = 0;
    slRemoteConfig.slReselectionConfig.slHystMin = 0;

    LteRrcSap::SlDiscConfigCommon slDiscConfigCommon;
    slDiscConfigCommon.slRemoteUeConfigCommon = slRemoteConfig;

    nrSlProseHelper->InstallNrSlDiscoveryConfiguration (uesNetDev, uesNetDev, slDiscConfigCommon); //Reusing this function made for l3 relay

  }


  /*********************** End ProSe configuration ***************************/

  /*
   * Configure the applications:
   * - Client app: OnOff application configured to generate CBR traffic.
   *               Installed in the initiating UE of each link by default, and
   *               in target UE if 'bidirectional' flag is True.
   * - Server app: PacketSink application to consume the received traffic.
   *               Installed in all UEs
   */
  NS_LOG_INFO ("Configuring applications...");
  // Random variable to randomize a bit start times of the client applications
  //to avoid simulation artifacts of all the TX UEs transmitting at the same time.
  Ptr<UniformRandomVariable> startTimeRnd = CreateObject<UniformRandomVariable> ();
  startTimeRnd->SetStream (70003); //Fix stream to have same results for each Run indepently of previous random variables creations

  startTimeRnd->SetAttribute ("Min", DoubleValue (0));
  startTimeRnd->SetAttribute ("Max", DoubleValue (0.10));

  std::string dataRateString = std::to_string (dataRate) + "kb/s";
  ApplicationContainer clientApps;
  ApplicationContainer serverApps;

  std::cout << "Traffic flows: " << std::endl;
  for (auto itPath = flowsInfo.begin (); itPath != flowsInfo.end (); ++itPath)
  {
    OnOffHelper sidelinkClient ("ns3::UdpSocketFactory",
                                InetSocketAddress (itPath->second.tgtIp, itPath->first)); //Towards target UE and port
    sidelinkClient.SetAttribute ("EnableSeqTsSizeHeader", BooleanValue (true));
    sidelinkClient.SetConstantRate (DataRate (dataRateString), udpPacketSize);
    ApplicationContainer app = sidelinkClient.Install (GetNodeById(ueNodes, itPath->second.srcNodeId)); // Installed in source UE 
    Time appStart = startTrafficTime + Seconds (startTimeRnd->GetValue ());
    app.Start (appStart);
    clientApps.Add (app);
    NS_LOG_INFO ("OnOff application installed in UE nodeId "
                << itPath->second.srcNodeId 
                << " srcIp " << itPath->second.srcIp 
                << " towards UE nodeId " << itPath->second.tgtNodeId 
                << " dstIp "<< itPath->second.tgtIp << ":" <<itPath->first );
    std::cout << "Flow Id: " << itPath->second.id << ", " <<itPath->second.srcIp << " -> " << itPath->second.tgtIp << ":" << itPath->first
              << " start time: " << appStart.GetSeconds ()
              << " s, end time: " <<  stopTrafficTime.GetSeconds () << " s" << std::endl;



    PacketSinkHelper sidelinkSink ("ns3::UdpSocketFactory",
                                  InetSocketAddress (Ipv4Address::GetAny (), itPath->first));
    sidelinkSink.SetAttribute ("EnableSeqTsSizeHeader", BooleanValue (true));
    ApplicationContainer serverApp = sidelinkSink.Install (GetNodeById(ueNodes, itPath->second.tgtNodeId)); // Installed in target UE
    serverApps.Add (serverApp);

    if (bidirectional)
    {
      OnOffHelper sidelinkClient ("ns3::UdpSocketFactory",
                                  InetSocketAddress (itPath->second.srcIp, itPath->first)); // Towards source UE
      sidelinkClient.SetAttribute ("EnableSeqTsSizeHeader", BooleanValue (true));
      sidelinkClient.SetConstantRate (DataRate (dataRateString), udpPacketSize);
      ApplicationContainer app = sidelinkClient.Install (GetNodeById(ueNodes, itPath->second.tgtNodeId)); // Installed in target UE 
      Time appStart = startTrafficTime + Seconds (startTimeRnd->GetValue ());
      app.Start (appStart);
      clientApps.Add (app);
      NS_LOG_INFO ("OnOff application installed in UE nodeId "
              << itPath->second.tgtNodeId 
              << " srcIp " << itPath->second.tgtIp 
              << " towards UE nodeId " << itPath->second.srcNodeId
              << " dstIp "<< itPath->second.srcIp << ":" <<itPath->first );

      std::cout << itPath->second.tgtIp << " -> " << itPath->second.srcIp << ":" << itPath->first
            << " start time: " << appStart.GetSeconds ()
            << " s, end time: " << stopTrafficTime.GetSeconds () << " s" << std::endl;
    
      PacketSinkHelper sidelinkSink ("ns3::UdpSocketFactory",
                                InetSocketAddress (Ipv4Address::GetAny (), itPath->first));
      sidelinkSink.SetAttribute ("EnableSeqTsSizeHeader", BooleanValue (true));
      ApplicationContainer serverApp = sidelinkSink.Install (GetNodeById(ueNodes, itPath->second.srcNodeId)); // Installed in source UE
      serverApps.Add (serverApp);
    
    }
  }
  clientApps.Stop (stopTrafficTime);
  serverApps.Start (Seconds (2.0));
  serverApps.Stop (stopTrafficTime);

  /******************** End application configuration ************************/

  /*
   * Hook the traces, to be used to compute average PIR and to data to be
   * stored in a database
   */

  //Trace receptions; use the following to be robust to node ID changes
  std::ostringstream path;
  path << "/NodeList/*/ApplicationList/*/$ns3::PacketSink/Rx";
  Config::ConnectWithoutContext (path.str (), MakeCallback (&ReceivePacket));
  path.str ("");

  path << "/NodeList/*/ApplicationList/*/$ns3::OnOffApplication/Tx";
  Config::ConnectWithoutContext (path.str (), MakeCallback (&TransmitPacket));
  path.str ("");

  //Database setup
  std::string exampleName = "nr-prose-u2u-multihop";
  SQLiteOutput db (exampleName + ".db");

  UeMacPscchTxOutputStats pscchStats;
  pscchStats.SetDb (&db, "pscchTxUeMac");
  Config::ConnectWithoutContext ("/NodeList/*/DeviceList/*/$ns3::NrUeNetDevice/"
                                 "ComponentCarrierMapUe/*/NrUeMac/SlPscchScheduling",
                                 MakeBoundCallback (&NotifySlPscchScheduling, &pscchStats));

  UeMacPsschTxOutputStats psschStats;
  psschStats.SetDb (&db, "psschTxUeMac");
  Config::ConnectWithoutContext ("/NodeList/*/DeviceList/*/$ns3::NrUeNetDevice/"
                                 "ComponentCarrierMapUe/*/NrUeMac/SlPsschScheduling",
                                 MakeBoundCallback (&NotifySlPsschScheduling, &psschStats));

  UePhyPscchRxOutputStats pscchPhyStats;
  pscchPhyStats.SetDb (&db, "pscchRxUePhy");
  Config::ConnectWithoutContext (
      "/NodeList/*/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/NrUePhy/"
      "NrSpectrumPhyList/*/RxPscchTraceUe",
      MakeBoundCallback (&NotifySlPscchRx, &pscchPhyStats));

  UePhyPsschRxOutputStats psschPhyStats;
  psschPhyStats.SetDb (&db, "psschRxUePhy");
  Config::ConnectWithoutContext (
      "/NodeList/*/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/NrUePhy/"
      "NrSpectrumPhyList/*/RxPsschTraceUe",
      MakeBoundCallback (&NotifySlPsschRx, &psschPhyStats));

  AsciiTraceHelper ascii;
  Ptr<OutputStreamWrapper> RelayNasRxPacketTraceStream = ascii.CreateFileStream ("NrSlRelayNasRxPacketTrace.txt");
  *RelayNasRxPacketTraceStream->GetStream () << "time(s)\tnodeIp\tsrcIp\tdstIp\tdstPort\tsrcLink\tdstLink" << std::endl;
  for (uint32_t i = 0; i < uesNetDev.GetN (); ++i)
  {
    Ptr<EpcUeNas> epcUeNas = uesNetDev.Get (i)->GetObject<NrUeNetDevice> ()->GetNas ();
    epcUeNas->TraceConnect ("NrSlU2uRelayRxPacketTrace", std::to_string(uesNetDev.Get (i)->GetNode ()->GetId ()),
                                    MakeBoundCallback (&TraceSinkRelayNasRxPacketTrace,
                                                        RelayNasRxPacketTraceStream));
  }

  Ptr<OutputStreamWrapper> DiscoveryTraceStream = ascii.CreateFileStream ("DiscoveryRsrpTrace.txt");
  *DiscoveryTraceStream->GetStream () << "time(s)\trxL2Id\tpeerL2Id\trsrp(dBm)" << std::endl;

  for (uint32_t i = 0; i < uesNetDev.GetN (); ++i)
  {
    Ptr<NrSlUeProse> prose = uesNetDev.Get (i)->GetObject<NrUeNetDevice> ()->GetObject <NrSlUeProse> ();
    prose->TraceConnectWithoutContext ("RelayRsrpTrace",
                                        MakeBoundCallback (&TraceSinkDiscoveryRsrpTrace,
                                                          DiscoveryTraceStream));
  }


  std::map<uint32_t, NodePhyStats*> nodePhyStats;
  for (uint32_t i = 0; i < uesNetDev.GetN (); ++i)
  {
        uint32_t nodeId = uesNetDev.Get (i)->GetNode ()->GetId ();
        NodePhyStats* stats = new NodePhyStats(); // Create a new NodePhyStats object
        stats->m_nodeId = nodeId;

        Ptr<LteUeRrc> rrc = uesNetDev.Get (i)->GetObject<NrUeNetDevice> ()->GetRrc ();
        stats->m_l2Id = rrc->GetSourceL2Id (); 
        nodePhyStats[nodeId] = stats;

        std::ostringstream ossTxCtrl;
        ossTxCtrl << "/NodeList/" << nodeId << "/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/NrUePhy/NrSpectrumPhyList/*/TxCtrlTrace";
      Config::ConnectWithoutContext(ossTxCtrl.str(), MakeBoundCallback(&TraceSlPscchTx, stats));

        std::ostringstream ossTxData;
        ossTxData << "/NodeList/" << nodeId << "/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/NrUePhy/NrSpectrumPhyList/*/TxDataTrace";
      Config::ConnectWithoutContext(ossTxData.str(), MakeBoundCallback(&TraceSlPsschTx, stats));

        std::ostringstream ossRxCtrl;
        ossRxCtrl << "/NodeList/" << nodeId << "/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/NrUePhy/NrSpectrumPhyList/*/RxPscchTraceUe";
      Config::ConnectWithoutContext(ossRxCtrl.str(), MakeBoundCallback(&TraceSlPscchRx, stats));

      std::ostringstream ossRxData;
        ossRxData << "/NodeList/" << nodeId << "/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/NrUePhy/NrSpectrumPhyList/*/RxPsschTraceUe";
      Config::ConnectWithoutContext(ossRxData.str(), MakeBoundCallback(&TraceSlPsschRx, stats));

      std::ostringstream ossRxCtrlHd;
      ossRxCtrlHd << "/NodeList/" << nodeId << "/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/NrUePhy/NrSpectrumPhyList/*/SlPscchHalfDuplexDrop";
      Config::ConnectWithoutContext(ossRxCtrlHd.str(), MakeBoundCallback(&TraceSlPscchRxHd, stats));

      std::ostringstream ossRxDataHd;
      ossRxDataHd << "/NodeList/" << nodeId << "/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/NrUePhy/NrSpectrumPhyList/*/SlPsschHalfDuplexDrop";
      Config::ConnectWithoutContext(ossRxDataHd.str(), MakeBoundCallback(&TraceSlPsschRxHd, stats));

      std::ostringstream ossRxDataIgn;
      ossRxDataIgn << "/NodeList/" << nodeId << "/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/NrUePhy/NrSpectrumPhyList/*/SlPsschIgnored";
      Config::ConnectWithoutContext(ossRxDataIgn.str(), MakeBoundCallback(&TraceSlPsschRxIg, stats));
  }

  //Configure FlowMonitor to get traffic flow statistics
  FlowMonitorHelper flowmonHelper;
  NodeContainer endpointNodes;
  endpointNodes.Add (ueNodes);

  Ptr<ns3::FlowMonitor> monitor = flowmonHelper.Install (endpointNodes);
  monitor->SetAttribute ("DelayBinWidth", DoubleValue (0.001));
  monitor->SetAttribute ("JitterBinWidth", DoubleValue (0.001));
  monitor->SetAttribute ("PacketSizeBinWidth", DoubleValue (20));

  Simulator::Stop (stopTrafficTime + Seconds (3.0));
  Simulator::Run ();

  std::cout << "Simulation done!"  << std::endl; 

  //Get per-flow traffic statistics
  monitor->CheckForLostPackets ();
  Ptr<Ipv4FlowClassifier> classifier = DynamicCast<Ipv4FlowClassifier> (flowmonHelper.GetClassifier ());
  FlowMonitor::FlowStatsContainer stats = monitor->GetFlowStats ();

  
  for (auto itPath = flowsInfo.begin (); itPath != flowsInfo.end (); ++itPath)
  {
    for (std::map<FlowId, FlowMonitor::FlowStats>::const_iterator i = stats.begin (); i != stats.end (); ++i)
    {
      Ipv4FlowClassifier::FiveTuple t = classifier->FindFlow (i->first);

      if (itPath->first == t.destinationPort)
      {
        itPath->second.trafficStats.nTxPkts = i->second.txPackets;
        itPath->second.trafficStats.nRxPkts = i->second.rxPackets;
        itPath->second.trafficStats.lossRatio = (double) (i->second.txPackets - i->second.rxPackets) / i->second.txPackets;
        if (i->second.rxPackets > 0)
          {
            itPath->second.trafficStats.meanDelay = 1000 * i->second.delaySum.GetSeconds () / i->second.rxPackets;
            itPath->second.trafficStats.meanJitter =  1000 * i->second.jitterSum.GetSeconds () / i->second.rxPackets;
          }
        else
          {
            itPath->second.trafficStats.meanDelay = 0;
            itPath->second.trafficStats.meanJitter = 0;        
          }
          break;
      }
    }
  }

  //For avg calculations
  double sumLoss = 0;
  double sumDelay = 0;
  double sumJitter = 0;
  double sumNhops = 0;
  double numPathsWithRelays = 0;

  std::ofstream outFile;
  std::string filename = "flowsStats.csv";
  outFile.open (filename.c_str (), std::ofstream::out | std::ofstream::trunc);
  if (!outFile.is_open ())
  {
    std::cerr << "Can't open file " << filename << std::endl;
    return 1;
  }

  outFile.setf (std::ios_base::fixed);
  outFile << "FlowId,srcIp,tgtIp:tgtPort,srcNodeId,tgtNodeId,nHops,nTxPkts,nRxPkts,LossRatio,MeanDelay,MeanJitter\n";
  for (auto itPath = flowsInfo.begin (); itPath != flowsInfo.end (); ++itPath)
  {
    outFile << itPath->second.id << ","
            << itPath->second.srcIp << ","
              << itPath->second.tgtIp << ":" << itPath->first << ","
              << itPath->second.srcNodeId << ","
              << itPath->second.tgtNodeId << ","
              << itPath->second.nHops << ","
              << itPath->second.trafficStats.nTxPkts << ","
              << itPath->second.trafficStats.nRxPkts << ","
              << itPath->second.trafficStats.lossRatio << ","
              << itPath->second.trafficStats.meanDelay << ","
              << itPath->second.trafficStats.meanJitter 
              << std::endl;

    if (itPath->second.nHops > 0)
    {
      sumNhops += itPath->second.nHops;
      sumLoss += itPath->second.trafficStats.lossRatio;
      sumDelay += itPath->second.trafficStats.meanDelay;
      sumJitter += itPath->second.trafficStats.meanJitter;
    }
    if (itPath->second.nHops > 1)
    {
      numPathsWithRelays++;
    }
  }
  outFile.close ();

  std::cout << "Traffic flows statistics: " << std::endl;
  std::ifstream f (filename.c_str ());
  if (f.is_open ())
    {
      std::cout << f.rdbuf ();
    }

  //Simulation stats
  if (pathFindingAlgo == "DistShortest" || pathFindingAlgo == "DistRandom")
  {
    ExportDataForPlotting (ueNodes, selectedRelays, nodePaths, d); 

    std::vector<RelayStats> relayStats = GetRelayStats (selectedRelays, nodePaths);
    std::cout << "Relay Stats:\nRelayNodeID\tPathCount\n";
    for (const auto& relayInfo : relayStats) {
        std::cout << relayInfo.relayNodeId << "\t"
                  << relayInfo.pathCount << std::endl;
    }

    for (const auto& relayInfo : relayStats) {
        g_minPathCount = std::min(g_minPathCount, relayInfo.pathCount);
        g_maxPathCount = std::max(g_maxPathCount, relayInfo.pathCount);
        g_totalPathCount += relayInfo.pathCount;
        if (relayInfo.pathCount > 0) {
              g_minNonZeroPathCount = std::min(g_minNonZeroPathCount, relayInfo.pathCount);
              g_totalNonZeroPathCount += relayInfo.pathCount;
              g_nonZeroCount++;
          }
    }

    if (!relayStats.empty()) 
    {
        g_meanPathCount = static_cast<double>(g_totalPathCount) / relayStats.size();

        std::cout << "Minimum Path Count: " << g_minPathCount << std::endl;
        std::cout << "Maximum Path Count: " << g_maxPathCount << std::endl;
        std::cout << "Mean Path Count: " << g_meanPathCount << std::endl;
        if (g_nonZeroCount > 0) 
        {
          g_meanNonZeroPathCount = static_cast<double>(g_totalNonZeroPathCount) / g_nonZeroCount;

          std::cout << "Minimum Path Count (Non-Zero): " << g_minNonZeroPathCount << std::endl;
          std::cout << "Mean Path Count (Non-Zero): " << g_meanNonZeroPathCount << std::endl;
        } else {
            std::cout << "No non-zero relay path information available." << std::endl;
        }
    } else 
    {
        std::cout << "No relay path information available." << std::endl;
    }
  }
  else
  {
    //This is done in FindPathsFromDiscoveryInfo function
  }
    
  double ratioPathsFound = (double) g_foundNpaths/nPaths; //Ok for distance. discovery?
  std::ofstream outFileSys;

  std::string filenameSys = "simStats.csv";
  outFileSys.open (filenameSys.c_str (), std::ofstream::out | std::ofstream::trunc);
  if (!outFileSys.is_open ())
  {
    std::cerr << "Can't open file " << filename << std::endl;
    return 1;
  }
  outFileSys.setf (std::ios_base::fixed);
  outFileSys << "RngSeed,RngRun,sysNTxPkts,sysNRxPkts,sysLossRatio,";
  outFileSys << "avgNHops,minRelayPathCount,maxRelayPathCount,meanRelayPathCount,minNonZeroRelayPathCount,meanNonZeroRelayPathCount,";
  outFileSys << "avgFlowLossRatio,avgFlowMeanDelay,avgFlowMeanJitter,"
             << "ratioPathsFound\n";
  outFileSys << RngSeedManager::GetSeed ()<< ","
             << RngSeedManager::GetRun ()<< ","
             << txPktCounter << ","
             << rxPktCounter << ","
             << ((double) (txPktCounter - rxPktCounter) / txPktCounter) << ",";

  if (g_foundNpaths > 0)
  {
    outFileSys << sumNhops / g_foundNpaths << ","; //Only paths that were found
  }
  else
  {
    outFileSys << "NaN" << ",";
  }
  if (numPathsWithRelays > 0)
  {
    outFileSys  << g_minPathCount << ","
                << g_maxPathCount << ","
                << g_meanPathCount << ","
                << g_minNonZeroPathCount << ","
                << g_meanNonZeroPathCount<< ",";
  }
  else 
  {
    //No path with relays found in the simulation
     outFileSys << "NaN" << "," << "NaN" << "," << "NaN" << "," << "NaN" << "," << "NaN" << ",";
  }

  if (g_foundNpaths > 0)
  {
    outFileSys << sumLoss / g_foundNpaths << "," //Only paths that were found
              << sumDelay / g_foundNpaths << "," //Only paths that were found
              << sumJitter / g_foundNpaths<< ","; //Only paths that were found
  }
  else
  {
    outFileSys << "NaN" << "," << "NaN" << "," << "NaN" << ",";
  }
  outFileSys << ratioPathsFound 
             << std::endl;
  outFile.close ();

  std::cout << "System statistics: " << std::endl;
  std::ifstream fs (filenameSys.c_str ());
  if (fs.is_open ())
  {
    std::cout << fs.rdbuf ();
  }

  std::cout << "Number of packets relayed by relays:"  << std::endl;
  std::cout << "relayNodeId" << "         " 
            << "relayIp"<< "         " 
            << "srcIp->dstIp"<< "         " 
            << "srcLink->dstLink"<< "\t\t" 
            << "nPackets"  << std::endl;
  for (auto it = g_relayNasPacketCounter.begin (); it != g_relayNasPacketCounter.end (); ++it)
  {
    std::cout << it->first << "\t\t" << it->second << std::endl;
  }

  uint64_t totalnTxCtrl = 0;
  uint64_t totalnTxData = 0;
  uint64_t totalnRxCtrl = 0;
  uint64_t totalnRxData = 0;
  uint64_t totalnCorruptRxCtrl = 0;
  uint64_t totalnCorruptRxData = 0;
  uint64_t totalnHdRxCtrl = 0;
  uint64_t totalnHdRxData = 0;
  uint64_t totalnIgnoredData = 0;

  std::cout << "\n Phy statistics:"  << std::endl;
  std::cout << std::left << std::setw(8) << "Node ID" << "| "
            << std::setw(8) << "L2 ID" << "| "
            << std::setw(8) << "nTxCtrl" << "| "
            << std::setw(8) << "nTxData" << "| "
            << std::setw(8) << "nRxCtrl" << "| "
            << std::setw(8) << "nRxData" << "| "
            << std::setw(16) << "nCorruptRxCtrl" << "| "
            << std::setw(16) << "nCorruptRxData" << "| "
            << std::setw(12) << "nHdRxCtrl" << "| "
            << std::setw(12) << "nHdRxData" << "| "
            << std::setw(12) << "nIgnoredData" << std::endl;
  std::cout << "----------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  for (const auto& [nodeId, stats] : nodePhyStats) 
  {
    std::cout << std::left << std::setw(8) << stats->m_nodeId << "| "
              << std::setw(8) << stats->m_l2Id << "| "
              << std::setw(8) << stats->m_stats.nTxCtrl << "| "
              << std::setw(8) << stats->m_stats.nTxData << "| "
              << std::setw(8) << stats->m_stats.nRxCtrl << "| "
              << std::setw(8) << stats->m_stats.nRxData << "| "
              << std::setw(16) << stats->m_stats.nCorruptRxCtrl << "| "
              << std::setw(16) << stats->m_stats.nCorruptRxData << "| "
              << std::setw(12) << stats->m_stats.nHdRxCtrl << "| "
              << std::setw(12) << stats->m_stats.nHdRxData << "| "
              << std::setw(12) << stats->m_stats.nIgnoredData << std::endl;

    totalnTxCtrl += stats->m_stats.nTxCtrl; 
    totalnTxData += stats->m_stats.nTxData; 
    totalnRxCtrl += stats->m_stats.nRxCtrl;
    totalnRxData += stats->m_stats.nRxData;
    totalnCorruptRxCtrl += stats->m_stats.nCorruptRxCtrl;
    totalnCorruptRxData += stats->m_stats.nCorruptRxData;
    totalnHdRxCtrl += stats->m_stats.nHdRxCtrl;
    totalnHdRxData += stats->m_stats.nHdRxData ;
    totalnIgnoredData += stats->m_stats.nIgnoredData;
  }
  std::cout << "----------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << std::left << std::setw(18) << "Total" << "| "
            << std::setw(8) << totalnTxCtrl << "| "
            << std::setw(8) << totalnTxData << "| "
            << std::setw(8) << totalnRxCtrl << "| "
            << std::setw(8) << totalnRxData << "| "
            << std::setw(16) << totalnCorruptRxCtrl << "| "
            << std::setw(16) << totalnCorruptRxData << "| "
            << std::setw(12) << totalnHdRxCtrl << "| "
            << std::setw(12) << totalnHdRxData << "| "
            << std::setw(12) << totalnIgnoredData << std::endl;
  std::cout << "----------------------------------------------------------------------------------------------------------------------------------------" << std::endl;

  std::cout << "Balance Control: " << totalnTxCtrl - totalnRxCtrl - totalnCorruptRxCtrl - totalnHdRxCtrl <<std::endl;
  std::cout << "Balance Data: " <<totalnTxData - totalnRxData - totalnCorruptRxData - totalnHdRxData - totalnIgnoredData <<std::endl;
  
  std::ofstream outFileSysPhy;
  std::string filenamePhyStats = "simPhyStats.csv";
  outFileSysPhy.open (filenamePhyStats.c_str (), std::ofstream::out | std::ofstream::trunc);
  if (!outFileSysPhy.is_open ())
  {
    std::cerr << "Can't open file " << filename << std::endl;
    return 1;
  }
  outFileSysPhy.setf (std::ios_base::fixed);
  outFileSysPhy << "RngSeed,RngRun,";
  outFileSysPhy << "totalnTxCtrl,totalnTxData,totalnRxCtrl,totalnRxData,totalnCorruptRxCtrl,totalnCorruptRxData,totalnHdRxCtrl,totalnHdRxData,totalnIgnoredData"<< std::endl;;
  outFileSysPhy << RngSeedManager::GetSeed ()<< ","
                << RngSeedManager::GetRun ()<< ",";
  outFileSysPhy << totalnTxCtrl << ","
                << totalnTxData << ","
                << totalnRxCtrl << ","
                << totalnRxData << ","
                << totalnCorruptRxCtrl << ","
                << totalnCorruptRxData << ","
                << totalnHdRxCtrl << ","
                << totalnHdRxData << ","
                << totalnIgnoredData << std::endl;
  outFileSysPhy.close ();

  /*
   * IMPORTANT: Do not forget to empty the database cache, which would
   * dump the data store towards the end of the simulation in to a database.
   */
  pscchStats.EmptyCache ();
  psschStats.EmptyCache ();
  pscchPhyStats.EmptyCache ();
  psschPhyStats.EmptyCache ();

  Simulator::Destroy ();

  return 0;
}
