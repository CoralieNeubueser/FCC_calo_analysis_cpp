import calo_init
calo_init.add_defaults()
## add arguments relevant only for that script
calo_init.parser.add_argument("--clusterColl", help="Name of the clusters collection (fcc::CaloClusterCollection)", type = str)
calo_init.parser.add_argument("--positionColl", help="Name of the positioned cells collection (fcc::PositionedCaloHitCollection)", type = str)
calo_init.parser.add_argument("--windowSeed", help="Size of the window used for seeding [eta,phi]", type = int, nargs=2)
calo_init.parser.add_argument("--windowPos", help="Size of the window used for berycentre coalculation [eta,phi]", type = int, nargs=2)
calo_init.parser.add_argument("--windowDupl", help="Size of the window used for duplicate removal [eta,phi]", type = int, nargs=2)
calo_init.parser.add_argument("--dEta", help="Size of the tower in eta", type = float, nargs=2)
calo_init.parser.add_argument("--maxEta", help="Maximum eta", type = float)
group = calo_init.parser.add_mutually_exclusive_group()
group.add_argument("--dPhi", help="Size of the tower in phi", type = float)
group.add_argument("--numPhi", help="Number of the towers in phi", type = int)
calo_init.parser.add_argument("--zoom", help="How many bins around centre should be visible", type = int, nargs=2)
calo_init.parser.add_argument("--event", help="Number of an event to draw", type = int)
calo_init.parse_args()
calo_init.print_config()

from math import pi, floor
# set of default parameters
eventToDraw = 0
maxEta = 1.8
nPhi = 512 # artificially increase by 1 (odd number) - to make plots look OK
dEta = 0.001
nEta = int(2*maxEta/dEta + 1)
dPhi = 2*pi/nPhi
nameClusterCollection = "EcalClusters"
nameCellCollection = "caloCells"
zoomEta = 21
zoomPhi = 21
etaWindowSeed = 9
phiWindowSeed = 9
etaWindowPos = 3
phiWindowPos = 3
etaWindowDupl = 5
phiWindowDupl = 5
# get parameters if passed from command line
if calo_init.args.clusterColl:
    nameClusterCollection = calo_init.args.clusterColl
if calo_init.args.positionColl:
    namePositionedHitsCollection = calo_init.args.positionColl
if calo_init.args.maxEta:
    maxEta = calo_init.args.maxEta
if calo_init.args.numPhi:
    nPhi = calo_init.args.numPhi
    dPhi = 2*pi/nPhi
if calo_init.args.dPhi:
    dPhi = calo_init.args.dPhi
    nPhi = int(dPhi*2*pi)
if calo_init.args.windowSeed:
    etaWindowSeed = calo_init.args.windowSeed[0]
    phiWindowSeed = calo_init.args.windowSeed[1]
if calo_init.args.windowPos:
    etaWindowPos = calo_init.args.windowPos[0]
    phiWindowPos = calo_init.args.windowPos[1]
if calo_init.args.windowDupl:
    etaWindowDupl = calo_init.args.windowDupl[0]
    phiWindowDupl = calo_init.args.windowDupl[1]
if calo_init.args.zoom:
    zoomEta = calo_init.args.zoom[0]
    zoomPhi = calo_init.args.zoom[1]
if calo_init.args.event:
    eventToDraw = calo_init.args.event

from ROOT import gSystem
gSystem.Load("libCaloAnalysis")
from ROOT import Decoder, ReconstructionExample, TCanvas, TFile, gStyle, gPad, kGreen, kRed, kBlue, TColor
from draw_functions import *

# use this script for multiple files
# gStyle.SetPalette(56) # kInvertedDarkBodyRadiator
gStyle.SetPalette(73) # kCMYK

for ifile, filename in enumerate(calo_init.filenamesIn):
    energy = calo_init.energy(ifile)
    dec = Decoder("system:4,cryo:1,module:11,type:3,subtype:3,layer:8,eta:15,phi:10")
    
    print "Energy of the initial particle: " + str(energy)
    print "File with reconstruction results: " + filename
    analysis = ReconstructionExample(nameClusterCollection,
                                        nameCellCollection,
                                        dec, # decoder for cell ids
                                        eventToDraw,
                                        energy,
                                        maxEta, # max eta
                                        nEta, # number of bins in eta
                                        nPhi, # number of bins in phi
                                        dEta, # tower size in eta
                                        dPhi) # tower size in phi
    analysis.loop(filename, calo_init.verbose)
    # retrieve histograms to draw them
    hist = []
    hist.append(analysis.hAllCellEnergy)
    hist.append(analysis.hClusterEnergy)
    hist.append(analysis.hClusterCellEnergy)
    print analysis.hClusterEnergy.GetEntries()
    print analysis.hClusterCellEnergy.GetEntries()
    ## Calculate some parameters
    enTotal = hist[0].Integral()
    meanEta = hist[1].ProjectionX().GetMean()
    meanPhi = hist[1].ProjectionY().GetMean()
    hist[2].SetMinimum(0)

    ## Set drawing options
    for h in hist:
        h.GetYaxis().SetRangeUser(meanPhi - (floor(zoomPhi/2) - 0.5 ) * dPhi,meanPhi + (floor(zoomPhi/2) + 0.5 ) * dPhi)
        h.GetXaxis().SetRangeUser(meanEta - (floor(zoomEta/2) - 0.5 ) * dEta,meanEta + (floor(zoomEta/2) + 0.5 ) * dEta)
        h.GetXaxis().SetNdivisions(zoomEta)
        h.GetYaxis().SetNdivisions(zoomPhi)
        h.GetXaxis().SetLabelSize(0.02)
        h.GetYaxis().SetLabelSize(0.02)

    canv = TCanvas('ECal_map_e'+str(energy)+'GeV', 'ECal', 1600, 1400 )
    canv.Divide(2,2)
    pad = canv.cd(1)
    pad.SetGrid()
    draw_hist2d(hist[0])
    pad = canv.cd(2)
    pad.SetLogz()
    pad.SetGrid()
    draw_hist2d(hist[0])
    draw_text(["energy: "+str(round(enTotal,2))+" GeV"])

    pad2 = canv.cd(3)
    draw_hist2d(hist[1])
    pad2.SetGrid()
    draw_text(["energy: "+str(round(hist[1].Integral(),2))+" GeV",
               "           "+str(round(hist[1].Integral()/enTotal*100,1))+" %"])

    pad = canv.cd(4)
    pad.SetGrid()
    draw_hist2d(hist[2])
    draw_text(["energy: "+str(round(hist[2].Integral(),2))+" GeV",
               "           "+str(round(hist[2].Integral()/enTotal*100,1))+" %"])
    for ipad in range(2,5):
        canv.cd(ipad)
        draw_rectangle([meanEta-etaWindowSeed/2.*dEta, meanPhi-phiWindowSeed/2.*dPhi],
                       [meanEta+etaWindowSeed/2.*dEta, meanPhi+phiWindowSeed/2.*dPhi], kRed, 4)
        draw_rectangle([meanEta-etaWindowDupl/2.*dEta, meanPhi-phiWindowDupl/2.*dPhi],
                       [meanEta+etaWindowDupl/2.*dEta, meanPhi+phiWindowDupl/2.*dPhi], kBlue, 3)
        draw_rectangle([meanEta-etaWindowPos/2.*dEta, meanPhi-phiWindowPos/2.*dPhi],
                       [meanEta+etaWindowPos/2.*dEta, meanPhi+phiWindowPos/2.*dPhi], kGreen, 2)
    if calo_init.output(ifile):
        canv.Print(calo_init.output(ifile)+'.root')
        canv.Print(calo_init.output(ifile)+'.png')
    else:
        canv.Print('ECal_map_e'+str(energy)+'GeV.png')
