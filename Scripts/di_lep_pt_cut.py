import ROOT
import numpy as np

filename = "ntuple.root"
rootfile = ROOT.TFile.Open(filename)

key1 = 'lepton_Pt'
electron_mass = 0.000511
muon_mass = 0.10566

hist_pixel_x = 800
hist_pixel_y = 600

num_hist_x = 1
num_hist_y = 1
canvas = ROOT.TCanvas("canvas", "canvas", hist_pixel_x * num_hist_x, hist_pixel_y * num_hist_y)
canvas.Divide(num_hist_x, num_hist_y)

iev_prev = 0
num_of_electrons = 0
electron_pt = []
electron_eta = []
electron_phi = []
lepton_Pts = []

for entries in rootfile.electron:
	iev = entries.iev
	# move to next particle within an event
	if iev == iev_prev:
		num_of_electrons += 1
		electron_pt.append(entries.electron_pt)
		electron_eta.append(entries.electron_eta)
		electron_phi.append(entries.electron_phi)
	else:
		# reached end of an event, calculate inv mass
		if num_of_electrons >= 2:
			# we select the highest pt leptons
			## find index of max pt electron
			i_max = np.argmax(electron_pt)
			## change max pt to 0 to find index of 2nd max pt  
			electron_pt_max_removed = electron_pt[:]
			electron_pt_max_removed[i_max] = 0
			i_second_max = np.argmax(electron_pt_max_removed)
			lep1 = ROOT.TLorentzVector(0,0,0,0)
			lep2 = ROOT.TLorentzVector(0,0,0,0)
			lep1.SetPtEtaPhiM(electron_pt[i_max], electron_eta[i_max],
							electron_phi[i_max], electron_mass)
			lep2.SetPtEtaPhiM(electron_pt[i_second_max], electron_eta[i_second_max], 
							electron_phi[i_second_max], electron_mass)
			# calculate transverse momentum from these 2 leptons
			x = (np.sqrt(np.square(electron_pt[i_max] * np.cos(electron_phi[i_max]) + electron_pt[i_second_max] * np.cos(electron_phi[i_second_max])) + np.square(electron_pt[i_max] * np.sin(electron_phi[i_max]) + electron_pt[i_second_max]*np.sin(electron_phi[i_second_max]))))
			if x > -50.0:
				lepton_Pts.append(x)
		iev_prev = iev
		num_of_electrons = 1
		electron_pt = [entries.electron_pt]
		electron_eta = [entries.electron_eta]
		electron_phi = [entries.electron_phi]

iev_prev = 0
num_of_muons = 0
muon_pt = []
muon_eta = []
muon_phi = []
lepton_Pts = []

for entries in rootfile.muon:
	iev = entries.iev
	# move to next particle within an event
	if iev == iev_prev:
		num_of_muons += 1
		muon_pt.append(entries.muon_pt)
		muon_eta.append(entries.muon_eta)
		muon_phi.append(entries.muon_phi)
	else:
		# reached end of an event, calculate inv mass
		if num_of_muons >= 2:
			# we select the highest pt leptons
			## find index of max pt muon
			i_max = np.argmax(muon_pt)
			## change max pt to 0 to find index of 2nd max pt  
			muon_pt_max_removed = muon_pt[:]
			muon_pt_max_removed[i_max] = 0
			i_second_max = np.argmax(muon_pt_max_removed)
			lep1 = ROOT.TLorentzVector(0,0,0,0)
			lep2 = ROOT.TLorentzVector(0,0,0,0)
			lep1.SetPtEtaPhiM(muon_pt[i_max], muon_eta[i_max],
							muon_phi[i_max], muon_mass)
			lep2.SetPtEtaPhiM(muon_pt[i_second_max], muon_eta[i_second_max], 
							muon_phi[i_second_max], muon_mass)
			# calculate transverse momentum from these 2 leptons
			x = (np.sqrt(np.square(muon_pt[i_max] * np.cos(muon_phi[i_max]) + muon_pt[i_second_max] * np.cos(muon_phi[i_second_max])) + np.square(muon_pt[i_max] * np.sin(muon_phi[i_max]) + muon_pt[i_second_max] * np.sin(muon_phi[i_second_max]))))
			if x > -50.0:
				lepton_Pts.append(x)
		iev_prev = iev
		num_of_muons = 1
		muon_pt = [entries.muon_pt]
		muon_eta = [entries.muon_eta]
		muon_phi = [entries.muon_phi]

lepton_Pt = ROOT.TH1D(key1, "Di-lepton Transverse Momentum;Transverse Momentum (GeV/c);Number of Events/ GeV", 100, -45, 45)
for x in lepton_Pts:
    lepton_Pt.Fill(x)

canvas.cd(1)
lepton_Pt.Draw()
canvas.SaveAs("pt_dilep_cut.pdf")
