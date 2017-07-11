#include "ReconstructionExample.h"

// podio specific includes
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

#include "datamodel/CaloClusterCollection.h"
#include "datamodel/CaloHitCollection.h"

#include "TVector3.h"

// STL
#include <vector>
#include <iostream>
#include <cmath>

ReconstructionExample::ReconstructionExample(const std::string& aCluserCollName, const std::string& aPosHitCollName, Decoder aDecoder, int aEventToAnalyse, double aEnergy, double aEtaMax, int aNoEta, int aNoPhi, double aDEta, double aDPhi):
  m_clusterCollName(aCluserCollName), m_posHitCollName(aPosHitCollName), m_decoder(aDecoder), m_eventToAnalyse(aEventToAnalyse), m_energy(aEnergy), m_etaMax(aEtaMax), m_noEta(aNoEta), m_noPhi(aNoPhi), m_dEta(aDEta), m_dPhi(aDPhi) {
  Initialize_histos();
}

ReconstructionExample::~ReconstructionExample(){}


void ReconstructionExample::Initialize_histos() {
  hAllCellEnergy = new TH2F("all cells","calo towers (50 GeV e^{-})",m_noEta,-m_etaMax,m_etaMax,m_noPhi,-M_PI-0.5*m_dPhi,M_PI+0.5*m_dPhi);
  hClusterEnergy = new TH2F("cluster","cluster seeds (50 GeV e^{-})",m_noEta,-m_etaMax,m_etaMax,m_noPhi,-M_PI-0.5*m_dPhi,M_PI+0.5*m_dPhi);
  hClusterCellEnergy = new TH2F("cells associated to clusters","cells in reconstructed cluster (50 GeV e^{-})",m_noEta,-m_etaMax,m_etaMax,m_noPhi,-M_PI-0.5*m_dPhi,M_PI+0.5*m_dPhi);
  m_histograms.push_back(hAllCellEnergy);
  m_histograms.push_back(hClusterCellEnergy);
  m_histograms.push_back(hClusterEnergy);
}

void ReconstructionExample::processEvent(podio::EventStore& aStore, int aEventId, bool aVerbose) {
  if (aEventId != m_eventToAnalyse) {
    return;
  }
  // Get the collections
  const fcc::CaloClusterCollection* clusters(nullptr);
  const fcc::CaloHitCollection* cells(nullptr);

  bool testClusters = aStore.get(m_clusterCollName, clusters);
  bool testCells = aStore.get(m_posHitCollName, cells);

  // cells associated to the reconstructed cluster
  std::vector<uint64_t> cellIds;

  // Get clusters reconstructed in an event
  if (testClusters) {
    if (aVerbose) {
      std::cout << "Number of clusters: " << clusters->size() << std::endl;
    }
    //Loop through the collection
    for (const auto iclu = clusters->begin(); iclu != clusters->end(); ++iclu) {
      if (aVerbose) {
        std::cout << "Cluster reconstructed at " << iclu->core().position.x
                  << " , " <<  iclu->core().position.y
                  << " , " <<  iclu->core().position.z
                  << "  with energy " <<  iclu->core().energy << " GeV" << std::endl;
      }
      TVector3 pos (iclu->core().position.x, iclu->core().position.y, iclu->core().position.z);
      float phi = pos.Phi();
      float eta = pos.Eta();
      hClusterEnergy->Fill(eta, phi, iclu->core().energy);
      // get cells hat are associated to that cluster
      int cellsNo = iclu->hits_size();
      for (int icell = 0; icell < cellsNo; icell++) {
        cellIds.push_back(iclu->hits(icell).core().cellId);
      }
    }
  } else {
    std::cout << "No Cluster Collection in the event." << std::endl;
  }

  // Get cells
  if (testCells) {
    if (aVerbose) {
      std::cout << "Number of cells: " << cells->size() << std::endl;
    }
    //Loop through the collection
    for (const auto icell = cells->begin(); icell != cells->end(); ++icell) {
      m_decoder.setValue(icell->core().cellId);
      float etaCell = Eta_offset + m_decoder.value("eta", icell->core().cellId) * m_dEta;
      float phiCell = Phi_offset + m_decoder.value("phi", icell->core().cellId) * m_dPhi;
      if (aVerbose) {
        std::cout << "Cell id: " << icell->core().cellId << " at eta: " << etaCell << " and phi: " << phiCell
                  << "  with energy " <<  icell->core().energy << " GeV" << std::endl;
      }
      hAllCellEnergy->Fill(etaCell, phiCell, icell->core().energy);
      if(std::find(cellIds.begin(), cellIds.end(), icell->core().cellId) != cellIds.end()) {
        hClusterCellEnergy->Fill(etaCell, phiCell, icell->core().energy);
      }
    }
  } else {
    std::cout << "No Cell Collection in the event." << std::endl;
  }
}

void ReconstructionExample::finishLoop(int aNumEvents, bool aVerbose) {}
