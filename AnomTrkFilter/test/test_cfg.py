import sys, string
from ROOT import *
#from cmstools import *
import FWCore.ParameterSet.Config as cms
gSystem.Load("libFWCoreFWLite.so")
AutoLibraryLoader.enable()


data = TFile('output/TrkFiltOutEventsFromMETPD.root')
events = EventTree(data.Get("Events"))

print 'events = ', events

trackBranch = events.branch("ctfWithMaterialTracks")
for event in events:
	tracks = trackBranch() # read tracks
	count = 0              # init counter
	for trk in tracks:     # loop over tracks
		if trk.pt() > 5.0:  # cut on pT
			count += 1          

		print "Found ",count," tracks" # print


