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
 * \file nr-prose-u2u-multihop-routing.cc
 * \brief A proof of concept of sidelink multihop UE-to-UE relay
 *
 * TODO!
 *
 * \code{.unparsed}
$ ./ns3 run "nr-prose-u2u-multihop-routing --help"
    \endcode
 *
 */

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/applications-module.h"
#include "ns3/mobility-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/nr-module.h"
#include "ns3/lte-module.h"
#include "ns3/stats-module.h"
#include "ns3/config-store-module.h"
#include "ns3/log.h"
#include "ns3/antenna-module.h"
#include <iomanip>
#include <sqlite3.h>
#include "ns3/flow-monitor-module.h"
#include "ns3/olsr-module.h"
#include "ns3/batmand-module.h"
#include "ns3/aodv-module.h"
#include "ns3/dsdv-module.h"
#include "ns3/nr-prose-module.h"
#include "ns3/spectrum-module.h"

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
 * TODO
 */
void
ReceivePacket (Ptr<const Packet> packet, const Address &from , const Address & to, const SeqTsSizeHeader &header)
{
  rxByteCounter += header.GetSize ();
  rxPktCounter++;
//  std::cout <<Simulator::Now().GetSeconds() << " APP -> Rx " <<  header.GetSeq ()<<" " << header.GetSize () <<std::endl;

}

/**
 * \brief Method to listen to the transmitting application trace Tx.
 * \param packet The packet 
 * TODO
 */
void
TransmitPacket (Ptr<const Packet> packet, const Address &from, const Address & to, const SeqTsSizeHeader &header)
{
  txByteCounter += header.GetSize ();
  txPktCounter++;
//  std::cout <<Simulator::Now().GetSeconds() << " APP -> Tx " <<  header.GetSeq () <<" " << header.GetSize () <<std::endl;
}


/* 
 * Global variables to calculate simulation statistics. 
 * Global because depending on the path finding algorithm the calculation is 
 * done either in the main or in a function scheduled for simulation time 
 */
double g_d = 0;

struct TrafficStats
{
  uint32_t nRxPkts=0;
  uint32_t nTxPkts=0;
  double lossRatio=std::numeric_limits<double>::quiet_NaN();
  double meanDelay=std::numeric_limits<double>::quiet_NaN();
  double meanJitter=std::numeric_limits<double>::quiet_NaN();
  double throughput=std::numeric_limits<double>::quiet_NaN();
};

/*
* TODO
*/
struct FlowInfo
{
  uint16_t id;
  Ipv4Address srcIp;
  Ipv4Address tgtIp;
  uint16_t tgtPort;
  uint16_t srcNodeId;
  uint16_t tgtNodeId;
  std::list<Ptr<NetDevice>> ueDevPath;
  uint16_t minNHops = std::numeric_limits<uint16_t>::max();
  uint16_t maxNHops = 0;
  uint16_t nRouteChanges = 0;
  uint16_t nDisruptions = 0;
  SidelinkInfo slInfo;
  TrafficStats trafficStats; //To be filled at the end of the simulation
};


void ExportDataForPlotting(const ns3::NodeContainer& nodes, 
                           const std::list<FlowInfo>& flows,
                           double d) {
  std::ofstream file("nodes_and_paths.csv");

  // Write node information
  for (uint32_t i = 0; i < nodes.GetN(); ++i) {
    Ptr<ns3::Node> node = nodes.Get(i);
    Ptr<ns3::MobilityModel> mobility = node->GetObject<ns3::MobilityModel>();
    ns3::Vector pos = mobility->GetPosition();
    bool isRelay = true;
    file << node->GetId() << "," << pos.x << "," << pos.y << "," << (isRelay ? "relay" : "non-relay") << "\n";
  }

  // Write path information
  for (const auto& flowEntry : flows) {
      file << "Path," << flowEntry.id;
      for (const auto& dev : flowEntry.ueDevPath) {
          file << "," << dev->GetNode()->GetId();
      }
      file << "\n";
  }
  file << "D," << d << "\n";;

  file.close();
}


// Helper function to get node ID from an IP address
uint32_t GetNodeIdFromIp(Ipv4Address ip, NodeContainer nodes) {
    for (uint32_t i = 0; i < nodes.GetN(); ++i) {
        Ptr<Node> node = nodes.Get(i);
        Ptr<Ipv4> ipv4 = node->GetObject<Ipv4>();
        for (uint32_t j = 0; j < ipv4->GetNInterfaces(); ++j) {
            Ipv4InterfaceAddress iface = ipv4->GetAddress(j, 0);
            if (iface.GetLocal() == ip) {
                return node->GetId();
            }
        }
    }
    return 255;  //Return 
}

Ptr<ns3::Node> GetNodeById(const ns3::NodeContainer& nodes, uint32_t nodeId) {
  for (uint32_t i = 0; i < nodes.GetN(); ++i) {
      if (nodes.Get(i)->GetId() == nodeId) {
          return nodes.Get(i);
      }
  }
  return nullptr; // Return nullptr if node is not found
}

Ipv4Address GetIpFromNodeId(uint32_t nodeId, NodeContainer nodes) {
    for (NodeContainer::Iterator i = nodes.Begin(); i != nodes.End(); ++i) {
        Ptr<Node> node = *i;
        if (node->GetId() == nodeId) {
            Ptr<Ipv4> ipv4 = node->GetObject<Ipv4>(); // Get the Ipv4 instance of the node
            if (ipv4 && ipv4->GetNInterfaces() > 1) { // Check if the node has at least one interface
                Ipv4InterfaceAddress iface = ipv4->GetAddress(1, 0); // Get the first interface address, 0 is loopback
                return iface.GetLocal(); // Return the IP address of the first interface
            } else {
                std::cerr << "No valid interfaces found on node " << nodeId << std::endl;
            }
        }
    }
    return Ipv4Address("0.0.0.0"); // Return a default "invalid" address if no node is found or no interface is set up
}

void PrintRoutingTables (NodeContainer ueNodes, uint32_t nodeId=255)
{
    Ptr<OutputStreamWrapper> routingTableStream = Create<OutputStreamWrapper>(&std::cout);

    *routingTableStream->GetStream () << "******************************" << "Simulation time (s): " << Simulator::Now().GetSeconds ()<<  "******************************" << std::endl;
    for ( uint32_t k = 0; k < ueNodes.GetN(); k ++ )
    {
        Ptr <ns3::Ipv4> ipv4 =  ueNodes.Get(k)->GetObject <ns3::Ipv4> ();
        if (!ipv4)
        {
          *routingTableStream->GetStream ()  << "Node " << ueNodes.Get(k)->GetId () << " Does not have an Ipv4 object";
          return;
        }
        if (ueNodes.Get(k)->GetId () == nodeId || nodeId == 255)
        {
          ipv4->GetRoutingProtocol ()->PrintRoutingTable (routingTableStream);
        }
    }
    *routingTableStream->GetStream () <<  "**********************************************************************************************" << std::endl;
}


//Global variable used to count routing info only after traffic
Time g_startTrafficTime;

void
CheckRoute (Ptr<NrSlProseHelper> nrSlProseHelper, FlowInfo& flowInfo, NodeContainer ueNodes, Ptr<OutputStreamWrapper> trace, Time interval, std::string prevPath) {
   
    std::vector<uint32_t> path;
    std::set<uint32_t> visitedNodes; // To keep track of visited nodes for loop detection
    bool loop = false;

    uint32_t currentNodeId = flowInfo.srcNodeId;
    Ipv4Address srcIp = GetIpFromNodeId (flowInfo.srcNodeId, ueNodes);
    Ipv4Address tgtIp = GetIpFromNodeId (flowInfo.tgtNodeId, ueNodes);
//    *trace->GetStream()<< std::fixed << std::setprecision(3)  << Simulator::Now().GetSeconds() << "," << srcIp << "," << flowInfo.srcNodeId << "," << tgtIp << "," << flowInfo.tgtNodeId << ",";

    while (currentNodeId != flowInfo.tgtNodeId) {
        if (visitedNodes.find(currentNodeId) != visitedNodes.end())
        {
            //Loop
            loop = true;
            path.clear();
            break;
        }
        path.push_back(currentNodeId);
        visitedNodes.insert(currentNodeId);

        Ptr<Node> currentNode = ueNodes.Get(currentNodeId);
        Ptr<Ipv4> ipv4 = currentNode->GetObject<Ipv4>();
        Ipv4Header ipv4Header;
        ipv4Header.SetDestination(tgtIp);

        Socket::SocketErrno se;
        Ptr<Ipv4Route> route = ipv4->GetRoutingProtocol()->RouteOutput(Ptr<Packet>(), ipv4Header, nullptr, se);

        if (route == nullptr) {
            path.clear(); // Clear the path as no valid route was found
            break;
        }

        Ipv4Address nextHopIp = route->GetGateway();
        int32_t nextNodeId = GetNodeIdFromIp(nextHopIp, ueNodes);

        if (nextNodeId == 255) {
            path.clear(); // Clear the path as there was an error finding the next node
            break;
        }

        currentNodeId = nextNodeId;
    }

    std::ostringstream currPath;
    if (!path.empty() && currentNodeId == flowInfo.tgtNodeId) 
    {
        path.push_back(flowInfo.tgtNodeId);  // Add the target node at the end only if a valid path is found

        //Create path string for comparison and trace
        for (size_t i = 0; i < path.size(); ++i) {
            if (i > 0) {
                currPath  << "->";
            }
            currPath  << path[i];
        }
    } else {
        if (loop)
        {
           currPath << "loop";
        }
        else
        {
          currPath << "nan";
        }
    }
    
    //Log the latest info on the number of hops  just before traffic starts
    if ((Simulator::Now() > (g_startTrafficTime - Seconds(1)) && Simulator::Now() < g_startTrafficTime) )
    {
      if (!path.empty())
      {
        flowInfo.minNHops = path.size () -1;
        flowInfo.maxNHops = path.size () -1;
      }
    }

    //Check if the path changed
    if (currPath.str () != prevPath)
    {
      std::cout << "Flow " << flowInfo.id << " changed route, previous: " << prevPath << ", new: " << currPath.str () <<std::endl;
      *trace->GetStream()<< std::fixed << std::setprecision(3)  << Simulator::Now().GetSeconds() << "," << srcIp << "," << flowInfo.srcNodeId << "," << tgtIp << "," << flowInfo.tgtNodeId << ","<< currPath.str () << std::endl;

      PrintRoutingTables (ueNodes);
      if (!path.empty())
      {
        //Update flow info data
        if (Simulator::Now() > g_startTrafficTime)
        {
          flowInfo.nRouteChanges ++;
          if (path.size () -1 <= flowInfo.minNHops )
            {
              flowInfo.minNHops = path.size () -1;
            }
            if (path.size () -1 >= flowInfo.maxNHops)
            {
              flowInfo.maxNHops = path.size () -1;
            }
        }
        //Configure new route
        //Clear old dev path
        flowInfo.ueDevPath.clear ();

        //Create dev path
        for (size_t i = 0; i < path.size(); ++i) {
          Ptr<Node> currentNode = ueNodes.Get(path[i]);
          Ptr<NrUeNetDevice> netDevPtr = currentNode->GetDevice (0)->GetObject<NrUeNetDevice> ();
          flowInfo.ueDevPath.push_back (netDevPtr);
        }

        //Call the helper
        nrSlProseHelper->ConfigureU2uRelayPath (flowInfo.srcIp , flowInfo.tgtIp, flowInfo.tgtPort, flowInfo.slInfo, flowInfo.ueDevPath);
      }
    }

    Simulator::Schedule(interval, &CheckRoute, nrSlProseHelper, std::ref(flowInfo), ueNodes, trace, interval, currPath.str ()); //Call it with prevPath = currPath

}

/*
* TODO
*/
uint32_t GetRandomNodeId (ns3::NodeContainer& nodes, Ptr<ns3::UniformRandomVariable> randomVariable) {

  uint32_t randomIndex = randomVariable->GetInteger(0, nodes.GetN() - 1);

  return nodes.Get(randomIndex)->GetId();
}


/*
 * Global variables to count TX/RX routing packets .
 */
uint32_t rxRoutingByteCounter = 0;
uint32_t txRoutingByteCounter = 0;
uint32_t rxRoutingPktCounter = 0;
uint32_t txRoutingPktCounter = 0;

void TxRoutingPacketOlsr(const ns3::olsr::PacketHeader &header, const ns3::olsr::MessageList &messages)
{
  if (Simulator::Now() > g_startTrafficTime)
  {
    txRoutingPktCounter ++;
    txRoutingByteCounter += header.GetPacketLength();
  }
/*  std::cout <<Simulator::Now().GetSeconds ()  << " OLSR TX - Header:  Seq " << header.GetPacketSequenceNumber() << " Len " << header.GetPacketLength() << " serialized size " << header.GetSerializedSize () << std::endl;
  for (const auto &message : messages) {
      std::cout << "Message - Type: " << message.GetMessageType () << " Seq " << message.GetMessageSequenceNumber ()
                << " Originator " <<  message.GetOriginatorAddress() << " Hop count " << +message.GetHopCount () << " serialized size " << message.GetSerializedSize () <<  std::endl;
  }*/
}

void RxRoutingPacketOlsr(const ns3::olsr::PacketHeader &header, const ns3::olsr::MessageList &messages)
{
  if (Simulator::Now() > g_startTrafficTime)
  {
    rxRoutingPktCounter ++;
    rxRoutingByteCounter += header.GetPacketLength();
  }
/*    std::cout <<Simulator::Now().GetSeconds ()  << " OLSR RX - Header:  Seq " << header.GetPacketSequenceNumber() << " Len " << header.GetPacketLength() << std::endl;
  for (const auto &message : messages) {
      std::cout << "Message - Type: " << message.GetMessageType () << " Seq " << message.GetMessageSequenceNumber ()
                << " Originator " <<  message.GetOriginatorAddress() << " Hop count " << +message.GetHopCount () <<   std::endl;
  } */
}


void TxRoutingPacketBatman(std::string nodeId, const ns3::batmand::OGMHeader &header, const ns3::batmand::MessageList &messages)
{
  if (Simulator::Now() > g_startTrafficTime)
  {
    txRoutingPktCounter ++;
  }
  //txRoutingByteCounter += header.GetPacketLength(); //The header is not having the correct lenght....

/* std::cout <<Simulator::Now().GetSeconds ()  << " " << nodeId << " BATMAN TX - Header -"
          << " Len: " << header.GetPacketLength ()  << " Seq: " << +header.GetPacketSequenceNumber ()
          <<  " Originator: " << header.GetOriginatorAddress() << " Forwarder: " << header.GetForwarderAddress ()
          << " Dir link: " << +header.HasDirectLinkFlag ()
          << " Uni: " << +header.HasUnidirectionalLinkFlag ()
          << " TTL: " << +header.GetTtl ()<< " TQ: " << +header.GetTQvalue()
          << " [N messages: " << messages.size ()   <<"]" <<   std::endl; */
  for (const auto &message : messages) {
      /* std::cout << "\t\t\t\t\t\t\t\t\t\t Message - Len: " << message.GetPacketLength ()  << " Seq: " << +message.GetPacketSequenceNumber ()
                << " Originator: " << message.GetOriginatorAddress()
                << " Forwarder: " << message.GetForwarderAddress ()
                << " Dir link: " << +message.HasDirectLinkFlag ()
                << " Uni: " << +message.HasUnidirectionalLinkFlag ()
                << " TTL: " << +message.GetTtl ()
                << " TQ: " << +message.GetTQvalue()
                << std::endl; */
    if (Simulator::Now() > g_startTrafficTime)
    {
      txRoutingByteCounter += message.GetPacketLength ();
    }
  }
}

void RxRoutingPacketBatman(std::string nodeId, const ns3::batmand::OGMHeader &header, const ns3::batmand::MessageList &messages)
{
  if (Simulator::Now() > g_startTrafficTime)
  {
    rxRoutingPktCounter ++;
  }
  //rxRoutingByteCounter += header.GetPacketLength(); //The header is not having the correct lenght....

  /* std::cout <<Simulator::Now().GetSeconds ()  << " " << nodeId  << " BATMAN RX - Header -"
          << " Len: " << header.GetPacketLength ()  << " Seq: " << +header.GetPacketSequenceNumber ()
           <<  " Originator: " << header.GetOriginatorAddress() << " Forwarder: " << header.GetForwarderAddress ()
           << " Dir link: " << +header.HasDirectLinkFlag ()
          << " Uni: " << +header.HasUnidirectionalLinkFlag ()
          << " TTL: " << +header.GetTtl ()<< " TQ: " << +header.GetTQvalue()
          << " [N messages: " << messages.size ()   <<"]" <<   std::endl; */
  for (const auto &message : messages) {
      /* std::cout << "\t\t\t\t\t\t\t\t\t\t Message - Len: " << message.GetPacketLength ()  << " Seq: " << +message.GetPacketSequenceNumber ()
                << " Originator: " << message.GetOriginatorAddress()
                << " Forwarder: " << message.GetForwarderAddress ()
                << " Dir link: " << +message.HasDirectLinkFlag ()
                << " Uni: " << +message.HasUnidirectionalLinkFlag ()
                << " TTL: " << +message.GetTtl ()
                << " TQ: " << +message.GetTQvalue()
                << std::endl; */

    if (Simulator::Now() > g_startTrafficTime)
    {
      rxRoutingByteCounter += message.GetPacketLength ();
    }
  }
}

void SetNodeTxPower(uint32_t nodeId, NodeContainer& nodes, double newTxPower) {
    Ptr<Node> node = nodes.Get(nodeId);

    // Attempt to retrieve the NR UE NetDevice from the node
    Ptr<NetDevice> netDevice = node->GetDevice(0); // Assuming the first device is the NR device
    Ptr<NrUeNetDevice> nrUeNetDevice = DynamicCast<NrUeNetDevice>(netDevice);
    if (nrUeNetDevice == nullptr) {
        std::cerr << "Node with ID " << nodeId << " does not have an NR UE NetDevice" << std::endl;
        return;
    }

    // Get the PHY layer and set the new TX power
    Ptr<NrUePhy> uePhy = nrUeNetDevice->GetPhy(0);
    if (uePhy != nullptr) {
        uePhy->SetTxPower(newTxPower);
        std::cout << "TX Power of Node " << nodeId << " set to " << newTxPower << " dBm" << std::endl;
    } else {
        std::cerr << "No PHY layer found on Node with ID " << nodeId << std::endl;
    }
}

void
CreateDisruption (FlowInfo& flowInfo, NodeContainer ueNodes)
{
  if (!flowInfo.ueDevPath.empty () && flowInfo.ueDevPath.size () > 2)
  {

    //Get path string
    std::ostringstream currPath;
    uint16_t index = 0;
    for (Ptr<NetDevice> netDevice : flowInfo.ueDevPath)
    {
      if (index > 0)
      {
          currPath  << "->";
      }
      currPath  << netDevice->GetNode ()->GetId ();
      index++;
    }
    //Get middle node ID
    uint32_t devIdx =  std::ceil (flowInfo.ueDevPath.size () / 2.0);
    std::list<Ptr<NetDevice>>::iterator it = flowInfo.ueDevPath.begin();
    std::advance(it, (devIdx-1));
    uint32_t nodeId =  (*it)->GetNode()->GetId();

    std::cout << Simulator::Now ().GetSeconds () << "Creating disruption for flow " << flowInfo.id 
              << ". Current path: " << currPath.str () << ", Node ID to disturb: " << nodeId << std::endl;
    flowInfo.nDisruptions ++;
    //Set Tx power to zero
    SetNodeTxPower(nodeId, ueNodes, 0);
  }
  else
  {
    std::cout << Simulator::Now ().GetSeconds () << "Unable to create disruption for flow " << flowInfo.id 
              << ". Current path size: " << flowInfo.ueDevPath.size () << std::endl;
  }
}

uint32_t nGrantDrops = 0;
void TraceGrantDrops (std::string nodeId, uint64_t imsi, uint16_t rnti, uint32_t dstL2Id, uint8_t lcid, uint32_t bytes)
{
  //std::cout << Now ().GetSeconds () << " Grant drop nodeId "<< nodeId << " IMSI " << imsi << " dstL2Id " << dstL2Id << " lcid " << +lcid << std::endl;
  nGrantDrops ++;
}


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
// void
// TraceSlPscchTx (NodePhyStats *stats, const Time duration)
// {
//   if (Simulator::Now() > g_startTrafficTime)
//   {
//     stats->m_stats.nTxCtrl ++;
//   }
// }
void
TraceSlPscchTx (NodePhyStats *stats, const uint64_t dstL2Id)
{
  if (Simulator::Now() > g_startTrafficTime && dstL2Id != 255) //Do not trace broadcast transmissions
  {
    stats->m_stats.nTxCtrl ++;
  }
}

/*
* TODO
*/
// void
// TraceSlPsschTx (NodePhyStats *stats, const Time duration)
// {
//   if (Simulator::Now() > g_startTrafficTime)
//   {
//     stats->m_stats.nTxData ++;
//   }
// }
void
TraceSlPsschTx (NodePhyStats *stats, const uint64_t dstL2Id)
{
  if (Simulator::Now() > g_startTrafficTime && dstL2Id != 255) //Do not trace broadcast transmissions
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
  if (Simulator::Now() > g_startTrafficTime)
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
  if (Simulator::Now() > g_startTrafficTime)
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
  if (Simulator::Now() > g_startTrafficTime)
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
  if (Simulator::Now() > g_startTrafficTime)
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
  if (Simulator::Now() > g_startTrafficTime)
  {
    stats->m_stats.nIgnoredData ++;
  }
}


int
main (int argc, char *argv[])
{
  
  bool useDb = false;
  std::string routingAlgo = "BATMAN";

  
  //Topology parameters
  uint16_t nUes = 3;
  uint16_t d = 60; //meters
  std::string deployment = "SqrRand";
  bool prioToSps = false;

  uint16_t nPaths = 1;

  std::string propagationType = "3GPP";

  // Simulation timeline parameters
  Time startTrafficTime = Seconds (60.0); //Time to start the traffic in the application layer.
  Time trafficDuration = Seconds (5*60.0);

  // NR parameters
  uint16_t numerologyBwpSl = 0; //The numerology to be used in sidelink bandwidth part
  double centralFrequencyBandSl = 5.89e9; // band n47  TDD //Here band is analogous to channel
  double bandwidthBandSl = 40e6; //40 MHz

  double txPower = 23; //dBm
  uint32_t mcs = 5; 

  // SL parameters
  bool sensing = true;
  uint32_t maxNumTx = 1;
  std::string harqType = "No";  

  // Traffic parameters
  uint32_t udpPacketSize = 60; //bytes
  double dataRate = 24; //kbps

  //Scenario parameter
  uint16_t nDisrupt = 0; //Number of route disruptions
  double interval = 1.0; //s

  CommandLine cmd;
  cmd.AddValue ("nUes", "The number of UEs in the topology", nUes);
  cmd.AddValue ("deployment", "The type of deployment. Options: Linear/SqrRand/SqrGrid", deployment);
  cmd.AddValue ("d", "The side of the square when using SqrGrid or SqrRand topologies, or the distance between first and last UE when using linear topology", d);
  cmd.AddValue ("routingAlgo", "The routing algorithm. Options: OLSR/BATMAN ", routingAlgo);
  cmd.AddValue ("mu", "The numerology to be used in sidelink bandwidth part", numerologyBwpSl);
  cmd.AddValue ("sensing", "If true, it enables the sensing based resource selection for SL, otherwise, no sensing is applied", sensing);
  cmd.AddValue ("maxNumTx", "The maximum numeber of transmissions in the PSSCH", maxNumTx);
  cmd.AddValue ("packetSizeBe", "packet size in bytes to be used by the traffic", udpPacketSize);
  cmd.AddValue ("dataRate", "The data rate in kilobits per second", dataRate);
  cmd.AddValue ("nPaths", "The number of paths to create", nPaths);
  cmd.AddValue ("propagationType", "The propagation type to use", propagationType);
  cmd.AddValue ("nDisrupt", "The number of route disruptions during traffic", nDisrupt);
  cmd.AddValue ("interval", "The intervalu used by the routing protocol", interval);
  cmd.AddValue ("harqType", "The type of HARQ used in the SL. Options No/Blind/Feedback1/Feedback2/Feedback4", harqType);

  // Parse the command line
  cmd.Parse (argc, argv);



  g_d = d;
  g_startTrafficTime = startTrafficTime;
  Time stopTrafficTime = startTrafficTime + trafficDuration;

  //Check if the frequency is in the allowed range.
  NS_ABORT_IF (centralFrequencyBandSl > 6e9);

  //Setup large enough buffer size to avoid overflow
  Config::SetDefault ("ns3::LteRlcUm::MaxTxBufferSize", UintegerValue (999999999));
  //Setup -t-Reordering to 0ms to avoid RLC reordering delay
  Config::SetDefault ("ns3::LteRlcUm::ReorderingTimer", TimeValue (MilliSeconds (0)));

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
    double interUeDistance = d / (nUes-1);
    for (uint16_t i = 0; i < nUes; i++)
      {
        listPositionAllocator->Add (Vector (interUeDistance * i, d/2, 1.5));
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
  std::vector<std::unique_ptr<BandwidthPartInfo> > bwpVector;
  CcBwpCreator ccBwpCreator; // This object needs to stay in scope
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
  //bandMask |= NrHelper::INIT_FADING;

  if (propagationType == "3GPP")
  {
      nrHelper->InitializeOperationBand (&bandSl, bandMask);
      allBwps = CcBwpCreator::GetAllBwps ({bandSl});
  }
  else if (propagationType == "AWGN")
  {
      std::unique_ptr<BandwidthPartInfo> bwp1 (new BandwidthPartInfo (1, centralFrequencyBandSl, bandwidthBandSl));
      auto spectrumChannel = CreateObject<MultiModelSpectrumChannel> ();
      Config::SetDefault("ns3::RangePropagationLossModel::MaxRange", DoubleValue(500));
      auto propagationLoss = CreateObject<RangePropagationLossModel> ();
      spectrumChannel->AddPropagationLossModel (propagationLoss);
      bwp1->m_channel = spectrumChannel;
      bwpVector.emplace_back(std::move(bwp1)); // Move ownership to the vector
      auto ref1 = std::ref(bwpVector[0]);
      allBwps.push_back(ref1);
  }
  else
  {
      std::cout << "Unknown propagation type " << propagationType << std::endl;
      exit(1);
  }

  //Packet::EnableChecking (); //BATAMAN DOESN'T WORK IF THIS IS ENABLED!
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
      //t2 = 9; // with T1 = 2, this gives a window of 8 slots with mu = 0 => 8 ms

      break;
    case 1:
      t2 = 33; // with T1 = 2, this gives a window of 32 slots with mu = 1 => 16 ms

      break;
    case 2:
      //t2 = 33; // with T1 = 2, this gives a window of 32 slots with mu = 2 => 8 ms
      t2 = 65; // with T1 = 2, this gives a window of 64 slots with mu = 2 => 16 ms
      //t2 = 17; // with T1 = 2, this gives a window of 16 slots with mu = 2 => 4 ms

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
   * In this example we use NrSlUeMacSchedulerFixedMcs scheduler, which uses
   * fix MCS value and schedules logical channels by priority order first and
   * then by creation order
   */
  nrSlHelper->SetNrSlSchedulerTypeId (NrSlUeMacSchedulerFixedMcs::GetTypeId ());
  nrSlHelper->SetUeSlSchedulerAttribute("Mcs", UintegerValue(mcs));
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
  std::vector<std::bitset<1>> slBitmap = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
  ptrFactory->SetSlTimeResources (slBitmap);
  ptrFactory->SetSlSensingWindow (100); // T0 in ms
  ptrFactory->SetSlSelectionWindow (5);
  ptrFactory->SetSlFreqResourcePscch (10); // PSCCH RBs
  ptrFactory->SetSlSubchannelSize (10);
  ptrFactory->SetSlMaxNumPerReserve (3);
  std::list<uint16_t> resourceReservePeriodList = {0, 20}; // in ms
  ptrFactory->SetSlResourceReservePeriodList (resourceReservePeriodList);
  ptrFactory->SetSlMinTimeGapPsfch (3);
  if (harqType == "Blind")
  {
    ptrFactory->SetSlPsfchPeriod (0); //0 (blind)
  } else if (harqType == "Feedback1")
  {
    ptrFactory->SetSlPsfchPeriod (1); //1 (every slot)
  }else if (harqType == "Feedback1")
  {
    ptrFactory->SetSlPsfchPeriod (2); //2 (every 2 slots)
  }else if (harqType == "Feedback2")
  {
    ptrFactory->SetSlPsfchPeriod (4); //4 (every 4 slots)
  } else
  {
    ptrFactory->SetSlPsfchPeriod (0); //0 (blind)
  }

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
  tddUlDlConfigCommon.tddPattern = "UL|UL|UL|UL|UL|UL|UL|UL|UL|UL|";

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
  BatmandHelper batman;
  Config::SetDefault ("ns3::batmand::RoutingProtocol::OGMInterval", TimeValue (Seconds (interval)));
  Config::SetDefault ("ns3::olsr::RoutingProtocol::HelloInterval", TimeValue (Seconds (interval)));

  OlsrHelper olsr;        
  
  Ipv4ListRoutingHelper list;
  uint16_t routingPort = 0;
  
  if (routingAlgo == "OLSR")
  {
    list.Add (olsr, 10);
    routingPort = 698;
  }
  else if (routingAlgo == "BATMAN")
  {
    list.Add (batman, 10);
    routingPort = 4305;
  }
  else {
    //Error
  }
  InternetStackHelper internet;
  internet.SetRoutingHelper (list); 
  internet.Install (ueNodes);
  internet.AssignStreams( ueNodes, 90000 );

  if (routingAlgo == "OLSR")
  {
    olsr.AssignStreams (ueNodes, 80000);
  }
  else if (routingAlgo == "BATMAN")
  {
    batman.AssignStreams (ueNodes, 80000);
  }

  // Install Ipv4 addresses
  Ipv4AddressHelper address;
  address.SetBase ("10.1.1.0", "255.255.255.0");
  Ipv4InterfaceContainer ueIfaces;
  ueIfaces = address.Assign (uesNetDev);


  //Boradcast bearer for routing messages
  uint32_t dstL2Id = 255;
  Ipv4Address broadcastAddress4 ("10.1.1.255");     //broadcast

  Ptr<LteSlTft> tft;
  SidelinkInfo slInfoRouting;
  slInfoRouting.m_castType = SidelinkInfo::CastType::Broadcast;
  slInfoRouting.m_dynamic = true;
  slInfoRouting.m_dstL2Id = dstL2Id;
  tft = Create<LteSlTft> (LteSlTft::Direction::BIDIRECTIONAL, broadcastAddress4, routingPort, slInfoRouting);
  nrSlHelper->ActivateNrSlBearer (Seconds (2.0), uesNetDev, tft);


  //Create ProSe helper
  Ptr<NrSlProseHelper> nrSlProseHelper = CreateObject<NrSlProseHelper> ();
  // Install ProSe layer and corresponding SAPs in the UEs
  nrSlProseHelper->PrepareUesForProse (uesNetDev);


  //Random variable for getting end nodes
  Ptr<ns3::UniformRandomVariable> randomVariableNodes = CreateObject<ns3::UniformRandomVariable>();
  randomVariableNodes->SetStream (70001); //Fix stream to have same results for each Run indepently of previous random variables creations


  ApplicationContainer clientApps;
  ApplicationContainer serverApps;
  std::list<FlowInfo> flows;

  //For each flow/path configure: 
  //-Routing path trace, 
  //-Configuration of the static route bearers
  //-Traffic flow

  uint16_t tgtPort = 8000;
  for (uint16_t i = 0; i < nPaths; i++)
  {
    //uint32_t srcNodeId = 0;
    //uint32_t tgtNodeId = 8;

    uint32_t srcNodeId;
    uint32_t tgtNodeId;
    do
    {
      srcNodeId = GetRandomNodeId (ueNodes, randomVariableNodes);
      tgtNodeId = GetRandomNodeId (ueNodes, randomVariableNodes);
    } while (srcNodeId == tgtNodeId);


    FlowInfo fInfo;
    fInfo.id = i;
    fInfo.srcNodeId = srcNodeId;
    fInfo.tgtNodeId = tgtNodeId;
    fInfo.srcIp = GetIpFromNodeId(srcNodeId, ueNodes);
    fInfo.tgtIp = GetIpFromNodeId(tgtNodeId, ueNodes);
    fInfo.tgtPort =  tgtPort;

    SidelinkInfo slInfoTraffic;
    slInfoTraffic.m_castType = SidelinkInfo::CastType::Unicast;
    slInfoTraffic.m_dynamic = false;
    slInfoTraffic.m_rri = MilliSeconds (20); 
    if(harqType != "No")
    {
      slInfoTraffic.m_harqEnabled = true;
    }
    fInfo.slInfo = slInfoTraffic;

    flows.push_back (fInfo);


    std::ostringstream routingTracefilename;
    routingTracefilename << "routing-trace_path_" << fInfo.id << ".csv";
    Ptr<OutputStreamWrapper> routingPathStream = Create<OutputStreamWrapper>(routingTracefilename.str(), std::ios::out);
    *routingPathStream->GetStream()  << "Time(s)" << "," << "srcIp" << "," << "srcNodeId" << "," << "tgtIp" << "," << "tgtNodeId" << "," << "PathNodeIds" << std::endl;;
    CheckRoute(nrSlProseHelper, std::ref(flows.back()), ueNodes, routingPathStream, Seconds (0.001),"nan");

    //Pertrubation/Disruption
    //Simulator::Schedule(Seconds (45.0), &SetNodeTxPower, 4, ueNodes, 0.0);

    if (nDisrupt > 0)
    {
      std::cout << "Scheduling disruption(s) for flow " << flows.back().id <<":" << std::endl;
      for (uint16_t n =1 ; n <= nDisrupt; n++)
      {
        Time disruptionTime = startTrafficTime + Seconds (n * (trafficDuration.GetSeconds () / (nDisrupt+1)));
        std::cout <<  " at " << disruptionTime.GetSeconds () << std::endl;
        Simulator::Schedule(disruptionTime, &CreateDisruption, std::ref(flows.back()), ueNodes);

      }
    }
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

    OnOffHelper sidelinkClient ("ns3::UdpSocketFactory",
                                InetSocketAddress (GetIpFromNodeId(tgtNodeId, ueNodes), tgtPort)); //Towards target UE and port
    sidelinkClient.SetAttribute ("EnableSeqTsSizeHeader", BooleanValue (true));
    sidelinkClient.SetConstantRate (DataRate (dataRateString), udpPacketSize);
    ApplicationContainer app = sidelinkClient.Install (GetNodeById(ueNodes, srcNodeId)); // Installed in source UE 
    Time appStartTime = startTrafficTime + Seconds (startTimeRnd->GetValue ());
    app.Start (appStartTime);
    app.Stop (stopTrafficTime);
    clientApps.Add (app);
    std::cout << "Flow " << fInfo.id << ", " << GetIpFromNodeId(srcNodeId, ueNodes) << "(Node ID " <<srcNodeId << ") -> " << GetIpFromNodeId(tgtNodeId, ueNodes)<< "(Node ID " <<tgtNodeId  << "):" << tgtPort
              << " start time: " << appStartTime.GetSeconds ()
              << " s, end time: " <<  stopTrafficTime.GetSeconds () << " s" << std::endl;

    PacketSinkHelper sidelinkSink ("ns3::UdpSocketFactory",
                                  InetSocketAddress (Ipv4Address::GetAny (), tgtPort));
    sidelinkSink.SetAttribute ("EnableSeqTsSizeHeader", BooleanValue (true));
    ApplicationContainer serverApp = sidelinkSink.Install (GetNodeById(ueNodes, tgtNodeId)); // Installed in target UE
    serverApps.Add (serverApp);

    tgtPort ++ ;
  }


  /*
   * Hook the traces
   */

  //Database setup
  std::string exampleName = "nr-prose-u2u-multihop";
  SQLiteOutput db (exampleName + ".db");
  UeMacPscchTxOutputStats pscchStats;
  UeMacPsschTxOutputStats psschStats;
  UePhyPscchRxOutputStats pscchPhyStats;
  UePhyPsschRxOutputStats psschPhyStats;

  if (useDb)
  {
    pscchStats.SetDb (&db, "pscchTxUeMac");
    Config::ConnectWithoutContext ("/NodeList/*/DeviceList/*/$ns3::NrUeNetDevice/"
                                  "ComponentCarrierMapUe/*/NrUeMac/SlPscchScheduling",
                                  MakeBoundCallback (&NotifySlPscchScheduling, &pscchStats));

    psschStats.SetDb (&db, "psschTxUeMac");
    Config::ConnectWithoutContext ("/NodeList/*/DeviceList/*/$ns3::NrUeNetDevice/"
                                  "ComponentCarrierMapUe/*/NrUeMac/SlPsschScheduling",
                                  MakeBoundCallback (&NotifySlPsschScheduling, &psschStats));

    pscchPhyStats.SetDb (&db, "pscchRxUePhy");
    Config::ConnectWithoutContext (
        "/NodeList/*/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/$ns3::BandwidthPartUe/NrUePhy/$ns3::NrUePhy/NrSpectrumPhy/$ns3::NrSpectrumPhy/RxPscchTraceUe",
        MakeBoundCallback (&NotifySlPscchRx, &pscchPhyStats));

    psschPhyStats.SetDb (&db, "psschRxUePhy");
    Config::ConnectWithoutContext (
        "/NodeList/*/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/$ns3::BandwidthPartUe/NrUePhy/$ns3::NrUePhy/NrSpectrumPhy/$ns3::NrSpectrumPhy/RxPsschTraceUe",
        MakeBoundCallback (&NotifySlPsschRx, &psschPhyStats));
  }
 
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


  //Trace traffic
  std::ostringstream path;
  path << "/NodeList/*/ApplicationList/*/$ns3::PacketSink/RxWithSeqTsSize";
  Config::ConnectWithoutContext (path.str (), MakeCallback (&ReceivePacket));
  path.str ("");

  path << "/NodeList/*/ApplicationList/*/$ns3::OnOffApplication/TxWithSeqTsSize";
  Config::ConnectWithoutContext (path.str (), MakeCallback (&TransmitPacket));
  path.str ("");

  //Trace routing messages
  if (routingAlgo == "OLSR")
  {
    path << "/NodeList/*/$ns3::Node/$ns3::olsr::RoutingProtocol/Tx";
    ns3::Config::ConnectWithoutContext(path.str (),ns3::MakeCallback(&TxRoutingPacketOlsr));
    path.str ("");
    path << "/NodeList/*/$ns3::Node/$ns3::olsr::RoutingProtocol/Rx";
    ns3::Config::ConnectWithoutContext(path.str (),ns3::MakeCallback(&RxRoutingPacketOlsr));
    path.str ("");
  }
  else if (routingAlgo == "BATMAN")
  {
    for (uint32_t i = 0; i < nUes; i++)
    {
      std::string nodeId = std::to_string(i);
      path << "/NodeList/"<<nodeId<<"/$ns3::Node/$ns3::batmand::RoutingProtocol/Tx";
      Config::ConnectWithoutContext(path.str (),MakeBoundCallback(&TxRoutingPacketBatman,nodeId));
      path.str ("");
      path << "/NodeList/"<<nodeId<<"/$ns3::Node/$ns3::batmand::RoutingProtocol/Rx";
      Config::ConnectWithoutContext(path.str (),MakeBoundCallback(&RxRoutingPacketBatman, nodeId));
      path.str ("");
    }
  }

  //Trace grant drops
  for (uint32_t i = 0; i < nUes; i++)
  {
    std::string nodeId = std::to_string(i);
    path << "/NodeList/"<<nodeId<<"/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/NrUeMac/NrSlGrantDrop";
    Config::ConnectWithoutContext (path.str (), MakeBoundCallback (&TraceGrantDrops, nodeId));
    path.str ("");

    path << "/NodeList/"<<nodeId<<"/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/$ns3::BandwidthPartUe/NrUeMac/$ns3::NrSlUeMac/NrSlUeMacScheduler/$ns3::NrSlUeMacSchedulerFixedMcs/NrSlGrantDrop";
    Config::ConnectWithoutContext (path.str (), MakeBoundCallback (&TraceGrantDrops, nodeId));
    path.str ("");
  }
  //Trace PHY losses
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
      ossTxCtrl << "/NodeList/" << nodeId << "/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/$ns3::BandwidthPartUe/NrUePhy/$ns3::NrUePhy/NrSpectrumPhy/$ns3::NrSpectrumPhy/TxPscchTrace";
      Config::ConnectWithoutContext(ossTxCtrl.str(), MakeBoundCallback(&TraceSlPscchTx, stats));

      std::ostringstream ossTxData;
      ossTxData << "/NodeList/" << nodeId << "/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/$ns3::BandwidthPartUe/NrUePhy/$ns3::NrUePhy/NrSpectrumPhy/$ns3::NrSpectrumPhy/TxPsschTrace";
      Config::ConnectWithoutContext(ossTxData.str(), MakeBoundCallback(&TraceSlPsschTx, stats));

      std::ostringstream ossRxCtrl;
      ossRxCtrl << "/NodeList/" << nodeId << "/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/$ns3::BandwidthPartUe/NrUePhy/$ns3::NrUePhy/NrSpectrumPhy/$ns3::NrSpectrumPhy/RxPscchTraceUe";
      Config::ConnectWithoutContext(ossRxCtrl.str(), MakeBoundCallback(&TraceSlPscchRx, stats));

      std::ostringstream ossRxData;
      ossRxData << "/NodeList/" << nodeId << "/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/$ns3::BandwidthPartUe/NrUePhy/$ns3::NrUePhy/NrSpectrumPhy/$ns3::NrSpectrumPhy/RxPsschTraceUe";
      Config::ConnectWithoutContext(ossRxData.str(), MakeBoundCallback(&TraceSlPsschRx, stats));

      std::ostringstream ossRxCtrlHd;
      ossRxCtrlHd << "/NodeList/" << nodeId << "/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/$ns3::BandwidthPartUe/NrUePhy/$ns3::NrUePhy/NrSpectrumPhy/$ns3::NrSpectrumPhy/SlPscchHalfDuplexDrop";
      Config::ConnectWithoutContext(ossRxCtrlHd.str(), MakeBoundCallback(&TraceSlPscchRxHd, stats));

      std::ostringstream ossRxDataHd;
      ossRxDataHd << "/NodeList/" << nodeId << "/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/$ns3::BandwidthPartUe/NrUePhy/$ns3::NrUePhy/NrSpectrumPhy/$ns3::NrSpectrumPhy/SlPsschHalfDuplexDrop";
      Config::ConnectWithoutContext(ossRxDataHd.str(), MakeBoundCallback(&TraceSlPsschRxHd, stats));

      std::ostringstream ossRxDataIgn;
      ossRxDataIgn << "/NodeList/" << nodeId << "/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/$ns3::BandwidthPartUe/NrUePhy/$ns3::NrUePhy/NrSpectrumPhy/$ns3::NrSpectrumPhy/SlPsschIgnored";
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

  Config::SetDefault ("ns3::ConfigStore::Filename", StringValue ("ZZZZZZZZZZZZZ.txt"));
  Config::SetDefault ("ns3::ConfigStore::Mode", StringValue ("Save"));
  ConfigStore outputConfig;
  outputConfig.ConfigureDefaults ();
  outputConfig.ConfigureAttributes ();

  Simulator::Run ();

  std::cout << "Simulation done!"  << std::endl; 


  PrintRoutingTables (ueNodes); 
  ExportDataForPlotting (ueNodes, flows, g_d*1.1); 


  //For avg calculations
  double sumLoss = 0;
  double sumDelay = 0;
  double sumJitter = 0;
  double sumTput = 0;
  double sumMinNhops = 0;
  double sumMaxNhops = 0;
  double sumNRouteChanges = 0;
  double numPathsFound = 0;
  double sumNDisruptions = 0;
  //Get per-flow traffic statistics
  monitor->CheckForLostPackets ();
  Ptr<Ipv4FlowClassifier> classifier = DynamicCast<Ipv4FlowClassifier> (flowmonHelper.GetClassifier ());
  FlowMonitor::FlowStatsContainer stats = monitor->GetFlowStats ();
  
  uint32_t ipHeaderSize = 28;
  
  for (auto& flowEntry : flows) 
  {
    for (std::map<FlowId, FlowMonitor::FlowStats>::const_iterator i = stats.begin (); i != stats.end (); ++i)
    {
      Ipv4FlowClassifier::FiveTuple t = classifier->FindFlow (i->first);

      if (flowEntry.tgtPort == t.destinationPort)
      {
        flowEntry.trafficStats.nTxPkts = i->second.txPackets;
        flowEntry.trafficStats.nRxPkts = i->second.rxPackets;
        flowEntry.trafficStats.lossRatio = (double) (i->second.txPackets - i->second.rxPackets) / i->second.txPackets;
        uint32_t rxBytesApp = i->second.rxBytes - (i->second.rxPackets * ipHeaderSize);
        flowEntry.trafficStats.throughput = (double) (rxBytesApp / trafficDuration.GetSeconds ()) * 8 / 1000.0; //kbps
        if (i->second.rxPackets > 0)
          {
            flowEntry.trafficStats.meanDelay = 1000 * i->second.delaySum.GetSeconds () / i->second.rxPackets;
            flowEntry.trafficStats.meanJitter =  1000 * i->second.jitterSum.GetSeconds () / i->second.rxPackets;
          }
        else
          {
            flowEntry.trafficStats.meanDelay = 0;
            flowEntry.trafficStats.meanJitter = 0;        
          }
          break;
      }
    }
  }

  std::ofstream outFileFlows;
  std::string filenameFlows = "flowsStats.csv";
  outFileFlows.open (filenameFlows.c_str (), std::ofstream::out | std::ofstream::trunc);
  if (!outFileFlows.is_open ())
  {
    std::cerr << "Can't open file " << filenameFlows << std::endl;
    return 1;
  }

  outFileFlows.setf (std::ios_base::fixed);
  outFileFlows << "FlowId,srcIp,tgtIp:tgtPort,srcNodeId,tgtNodeId,minNHops,maxNHops,nRouteChanges,nTxPkts,nRxPkts,LossRatio,MeanDelay,MeanJitter,avgThroughput,nDisruptions\n";
  for (auto& flowEntry : flows) 
  {
    outFileFlows << flowEntry.id << ","
                 << flowEntry.srcIp << ","
                 << flowEntry.tgtIp << ":" << flowEntry.tgtPort << ","
                 << flowEntry.srcNodeId << ","
                 << flowEntry.tgtNodeId << ","
                 << flowEntry.minNHops << ","
                 << flowEntry.maxNHops << ","
                 << flowEntry.nRouteChanges << ","
                 << flowEntry.trafficStats.nTxPkts << ","
                 << flowEntry.trafficStats.nRxPkts << ","
                 << flowEntry.trafficStats.lossRatio << ","
                 << flowEntry.trafficStats.meanDelay << ","
                 << flowEntry.trafficStats.meanJitter << ","
                 << flowEntry.trafficStats.throughput << ","
                 << flowEntry.nDisruptions
                 << std::endl;


    if (flowEntry.maxNHops > 0)
    {  
      sumMinNhops += flowEntry.minNHops;
      sumMaxNhops += flowEntry.maxNHops;
      sumNRouteChanges += flowEntry.nRouteChanges;
      sumLoss += flowEntry.trafficStats.lossRatio;
      sumDelay += flowEntry.trafficStats.meanDelay;
      sumJitter += flowEntry.trafficStats.meanJitter;
      sumTput += flowEntry.trafficStats.throughput;
      sumNDisruptions += flowEntry.nDisruptions;
      numPathsFound ++;
    }
  }
  outFileFlows.close ();
  //Print out on standard output:
  std::cout << "Traffic flows statistics: " << std::endl;
  std::ifstream f (filenameFlows.c_str ());
  if (f.is_open ())
    {
      std::cout << f.rdbuf ();
    }

  //System level statistics
  std::ofstream outFileSys;
  std::string filenameSys = "simStats.csv";
  outFileSys.open (filenameSys.c_str (), std::ofstream::out | std::ofstream::trunc);
  if (!outFileSys.is_open ())
  {
    std::cerr << "Can't open file " << filenameSys << std::endl;
    return 1;
  }
  outFileSys.setf (std::ios_base::fixed);
  outFileSys << "RngSeed,RngRun,"; 
  outFileSys << "sysNTxPkts,sysNRxPkts,";
  outFileSys << "sysLossRatio,";
  outFileSys << "sysThroughput,";
  outFileSys << "avgMinNHops,avgMaxNHops,avgFlowLossRatio,avgFlowMeanDelay,avgFlowMeanJitter,avgFlowThroughput,"
             << "ratioPathsFound,"
             << "sysNRouteChanges,sysNTxRoutingPkts,"
             << "sysRoutingOhPkt,sysRoutingOhBytes,"
             << "avgNDisruptions,"
             << "nGrantDrops\n";

  outFileSys << RngSeedManager::GetSeed ()<< ","
             << RngSeedManager::GetRun ()<< ","
             << txPktCounter << ","
             << rxPktCounter << ",";

  if (txPktCounter > 0)
  {
    outFileSys << ((double) (txPktCounter - rxPktCounter) / txPktCounter) << ",";
  }
  else
  {
    outFileSys << "NaN" << ",";
  }

  outFileSys << ((double) (rxByteCounter) / trafficDuration.GetSeconds ()) * 8 / 1000.0 << ","; //kbps 

  if (numPathsFound > 0)
  {
    outFileSys << sumMinNhops / numPathsFound << "," //Only paths that were found
               << sumMaxNhops / numPathsFound << "," //Only paths that were found
               << sumLoss / numPathsFound << "," //Only paths that were found
               << sumDelay / numPathsFound << "," //Only paths that were found
               << sumJitter / numPathsFound<< "," //Only paths that were found
               << sumTput / numPathsFound<< ","; //Only paths that were found
  }
  else
  {
    outFileSys << "NaN"<< "," << "NaN"<< "," "NaN" << "," << "NaN" << "," << "NaN"<< "," << "NaN" << ",";
  }
  outFileSys << numPathsFound / nPaths <<  ","; 
  if (numPathsFound > 0)
  {
    outFileSys << sumNRouteChanges<< "," << txRoutingPktCounter  << ",";
    if (txPktCounter > 0)
    {
      outFileSys << (double) txRoutingPktCounter / txPktCounter << ","
              << (double) txRoutingByteCounter / txByteCounter<< ",";
    }
    else
    {
      outFileSys << "NaN"<< "," << "NaN"<< ",";
    }
    outFileSys << sumNDisruptions/ numPathsFound<< ",";
  }
  else
  {
    outFileSys << "NaN"<< "," << "NaN"<< "," "NaN" << "," << "NaN"<< "," << "NaN"<< ",";;
  }
  outFileSys <<nGrantDrops;
  outFileSys << std::endl;
  outFileSys.close ();

  //Print out on standard output:
  std::cout << "System statistics: " << std::endl;
  std::ifstream fs (filenameSys.c_str ());
  if (fs.is_open ())
  {
    std::cout << fs.rdbuf ();
  }


  //Print NAS stats
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

  //Print PHY stats
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

  //Ouput per node PHY stats to csv
  std::ofstream outFileSysPhyPerNode;
  std::string filenamePhyStatPerNodes = "simPhyStats_perNode.csv";
  outFileSysPhyPerNode.open (filenamePhyStatPerNodes.c_str (), std::ofstream::out | std::ofstream::trunc);
  if (!outFileSysPhyPerNode.is_open ())
  {
    std::cerr << "Can't open file " << filenamePhyStatPerNodes << std::endl;
    return 1;
  }
  outFileSysPhyPerNode.setf (std::ios_base::fixed);
  outFileSysPhyPerNode << "RngSeed,RngRun,";
  outFileSysPhyPerNode << "NodeID" << ","
                       << "L2ID" << ","
                       << "nTxCtrl" << ","
                       <<  "nTxData" << ","
                       << "nRxCtrl" << ","
                       << "nRxData" << ","
                       << "nCorruptRxCtrl" << ","
                       << "nCorruptRxData" << ","
                       << "nHdRxCtrl" << ","
                       << "nHdRxData" << ","
                       << "nIgnoredData\n";
  for (const auto& [nodeId, stats] : nodePhyStats)
  {
    outFileSysPhyPerNode << RngSeedManager::GetSeed ()<< ","
                         << RngSeedManager::GetRun ()<< ",";
    outFileSysPhyPerNode << stats->m_nodeId << ","
                         << stats->m_l2Id << ","
                         << stats->m_stats.nTxCtrl << ","
                         << stats->m_stats.nTxData << ","
                         << stats->m_stats.nRxCtrl << ","
                         << stats->m_stats.nRxData << ","
                         << stats->m_stats.nCorruptRxCtrl << ","
                         << stats->m_stats.nCorruptRxData << ","
                         << stats->m_stats.nHdRxCtrl << ","
                         << stats->m_stats.nHdRxData << ","
                         << stats->m_stats.nIgnoredData << "\n";
  }
  outFileSysPhyPerNode.close ();


  //Output total PHY stats to csv
  std::ofstream outFileSysPhy;
  std::string filenamePhyStats = "simPhyStats.csv";
  outFileSysPhy.open (filenamePhyStats.c_str (), std::ofstream::out | std::ofstream::trunc);
  if (!outFileSysPhy.is_open ())
  {
    std::cerr << "Can't open file " << filenamePhyStats << std::endl;
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



  if (useDb)
  {
    /*
    * IMPORTANT: Do not forget to empty the database cache, which would
    * dump the data store towards the end of the simulation in to a database.
    */
    pscchStats.EmptyCache ();
    psschStats.EmptyCache ();
    pscchPhyStats.EmptyCache ();
    psschPhyStats.EmptyCache ();
  }
  Simulator::Destroy ();

  return 0;
}
