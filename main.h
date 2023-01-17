// This contains:
// 1. The set of LHE histograms I want to plot, saved to a root file
#ifndef LOOPER_H
#define LOOPER_H

#include "Nano.h" // Contains the definition of the "nt" object that reads the NanoAOD
#include "MCTools.h" // Contains the definition of the dumpGenParticleInfos();
#include "ElectronSelections.h" // Contains the definitions of ttH::electronID
#include "MuonSelections.h" // Contains the definitions of ttH::muonID
#include "Config.h" // Contains the definitions of gconf
#include "rooutil.h" // Contains the definitions of various help functions that start with "RooUtil"
#include "cxxopts.h" // Contains how to parse argc and argv 
#include "LHEanalysis.h" // Contains the LHE level data selection class

void readinput(int inputno, char** arguments, std::string& input_path, std::string& output_path, int& n_events);
int main(int argc, char** argv);
