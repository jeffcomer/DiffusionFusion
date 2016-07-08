//////////////////////////////////////////////////////////////////////
// Copyright 2014-2016 Jeffrey Comer
//
// This file is part of DiffusionFusion.
//
// DiffusionFusion is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// DiffusionFusion is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with DiffusionFusion. If not, see http://www.gnu.org/licenses/.
///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef DIFFUSIONFUSION_H
#define DIFFUSIONFUSION_H

#include <ctime>
#include <omp.h>
#include <limits>

#include "useful.H"
#include "CommandLineReader.H"
#include "TrajFilter.H"
#include "TrajCostDesc.H"
#include "Event.H"
#include "Field.H"
#include "Bin.H"
#include "PiecewiseLinear.H"
#include "PiecewiseCubic.H"
#include "PiecewiseBicubic.H"
#include "FieldDesc.H"
#include "HistoryReader.H"
#include "Prior.H"
#include "RandomGsl.H"
#include "MetroMonteCarlo.H"
#include "TrajCostComputer.H"

// Various types of trajectory cost calculators, some experimental.
#include "TrajComer.H"
#include "TrajSmolCrank.H"
#include "TrajFracSmolCrank.H"
#include "TrajReflect.H"
#include "TrajComer2d.H"
#include "TrajReflect2d.H"

// #include "TrajReflect2d.H"
// #include "TrajComerPmf.H"
// #include "TrajSmolCrank.H"
// #include "TrajSmolCrankBias.H"
// #include "TrajFracSmolCrank.H"
// #include "TrajFracSmolCrankDual.H"
// #include "TrajSmolCrankDual.H"

// Mostly experimental
// #include "TrajSmolCrank1.H"
// #include "TrajFoxFunction.H"
// #include "TrajComerPmf2d.H"
// #include "TrajLangevin.H"
// #include "TrajLangevinExp.H"
// #include "TrajLangevinExp2d.H"
// #include "TrajLangevinGamma.H"
// #include "TrajComerLog.H"
// #include "TrajTurkcan.H"
// #include "TrajComer2d2d.H"
// #include "TrajComerPmf2d2d.H"
// #include "TrajComer2d2d2d.H"
// #include "TrajFdSmoluchowski.H"
// #include "TrajSimpleSmol.H"
// #include "TrajSimpleSmolReflect.H"
// #include "TrajSmolReflectBias.H"
// #include "TrajNoise.H"

struct TrajFile {
public:
  String fileName;
  int stride;
  int skip;
  int position;
  int group;
  TrajFilter* filter;
  String biasFile;
  String biasCountFile;
  int biasFullSamples;
  double biasSmoothWidth;
  int eventFirst;
  int eventLast;
  int serialLast;

  static const int posInit = 0;
  static const int posFinal = 1;
  static const int posMid = 2;
  static int getPosIndex(const String& s) {
    if (s == "init") {
      return posInit;
    } else if (s == "final") {
      return posFinal;
    } else if (s == "mid") {
      return posMid;
    } else {
      return -1;
    }
  }

  String toString() const {
    char s[STRLEN];
    snprintf(s, STRLEN, "file %s stride %d skip %d position %d %s", fileName.cs(), stride, skip, position, filter->toString().cs());
    return String(s);
  }
};


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// The class that handles reading the command file and constructing the
// other objects.
class DiffusionFusion {
private:
  int cmdNum;
  CommandLineReader** cmdList;
  
  // Prefix of the output files.
  String outputPrefix;
  String outputTraj;
  int outputPeriod, previewPeriod, updatePeriod;
  // Period for possible trial moves on global variables.
  int globalPeriod;
  // The largest relative error between locally computed costs and the global cost that we permit.
  double relErrMax;

  // Descriptions of the trajectory variables.
  int trajVarNum;
  int* trajVarCol;
  String* trajVarName;
  IndexList colList; // 0-based list of columns corresponding to trajectory vars
  int dispFileMax; // Number of files to display in the log.

  // For bias history files.
  HistoryReader** biasHistoryList;
  int totalTrajFileNum;
  int biasFieldNum;
  const Piecewise1d** biasFieldList;

  // List of events.
  int eventMax; // capacity
  int eventNum;
  Event* event;

  // Fields.
  int fieldNum;
  Field** fieldList;
  FieldDesc* fieldDesc;
  Field** saveList; // Saving the lowest cost field values.
  
  // Priors.
  int priorNum;
  Prior** priorList;

  // The pseudorandom number generator
  Random* rando;

  // Monte Carlo object
  MetroMonteCarlo* monte;
  IndexList mcFieldSel;
  int mcCycles;

  // The objects where the magic happens.
  int trajCostNum;
  TrajCostComputer** trajCostList;
  bool** tcFieldOn;

  // Dump control
  static const int dumpTypeNum = 7;
  static const int dumpEvents = 0;
  static const int dumpFields = 1;
  static const int dumpDistro = 2;
  static const int dumpBest = 3;
  static const int dumpEventCost = 4;
  static const int dumpRefField = 5;
  static const int dumpBias = 6;
  String dumpFile[dumpTypeNum];

public:
  DiffusionFusion(const String& cmdFile, const String& outPreCmdLine) : cmdNum(0), cmdList(NULL), trajVarNum(0), dispFileMax(30), biasHistoryList(NULL), totalTrajFileNum(0), biasFieldNum(0), biasFieldList(NULL), eventMax(1), eventNum(0), event(NULL), fieldNum(0), fieldList(NULL), priorNum(0), priorList(NULL), rando(NULL), monte(NULL), mcCycles(0), trajCostList(NULL) {
    // Read the configuration file.
    cmdNum = countLines(cmdFile.cs(), IndexList(0));
    cmdList = new CommandLineReader*[cmdNum];
    for (int c = 0; c < cmdNum; c++) cmdList[c] = NULL;
    readCommands(cmdFile.cs(), cmdNum);
    printf("Read %d commands from `%s'.\n", cmdNum, cmdFile.cs());

    printf("Organizing commands.\n");
    int mcNum = 0, trajNum = 0, outputNum = 0, dumpNum = 0;
    int cmdMc = 0, cmdTraj = 0;
    IndexList cmdField, cmdPrior, cmdLoad, cmdLoadMany, cmdLoadEvents, cmdTrajCost;
    for (int c = 0; c < cmdNum; c++) {
      String cmd(cmdList[c]->getCommand());
      //printf("COMMAND %s parameters %d options %d\n", cmd.cs(), cmdList[c]->getParamNum(), cmdList[c]->getOptionNum());

      if (cmd == "output") {
	if (cmdList[c]->getParamNum() != 1) {
	  fprintf(stderr, "ERROR DiffusionFusion: Usage of output: `output outputPrefix'.\n");
	  exit(-1);
	}
	outputNum++;
	outputPrefix = cmdList[c]->getParam(0);
      } else if (cmd == "trajectory") {
	trajNum++;
	cmdTraj = c;
      } else if (cmd == "field") {
	cmdField.add(c);
      } else if (cmd == "load") {
	cmdLoad.add(c);
      } else if (cmd == "loadMany") {
	cmdLoadMany.add(c);
      } else if (cmd == "loadEvents") {
	cmdLoadEvents.add(c);
      } else if (cmd == "prior") {
	cmdPrior.add(c);
      } else if (cmd == "trajCost") {
	cmdTrajCost.add(c);
      } else if (cmd == "mc") {
	mcNum++;
	cmdMc = c;
      } else if (cmd == "dump") {
	// Handle dumping.
	if (cmdList[c]->getParamNum() != 2) {
	  fprintf(stderr, "ERROR DiffusionFusion: Usage of dump: `dump events|fields|distro|best dumpFile'.\n");
	  exit(-1);
	}
	String par0 = cmdList[c]->getParam(0);
	if (par0 == "events") {
	  dumpFile[dumpEvents] = cmdList[c]->getParam(1);
	} else if (par0 == "fields") {
	  dumpFile[dumpFields] = cmdList[c]->getParam(1);
	} else if (par0 == "distro") {
	  dumpFile[dumpDistro] = cmdList[c]->getParam(1);
	} else if (par0 == "best") {
	  dumpFile[dumpBest] = cmdList[c]->getParam(1);
	} else if (par0 == "eventCost") {
	  dumpFile[dumpEventCost] = cmdList[c]->getParam(1);
	} else if (par0 == "refField") {
	  dumpFile[dumpRefField] = cmdList[c]->getParam(1);
	} else if (par0 == "bias") {
	  dumpFile[dumpBias] = cmdList[c]->getParam(1);
	} else {
	  fprintf(stderr, "ERROR DiffusionFusion: Unrecognized dump command `%s'.\n", par0.cs());
	  exit(-1);
	}
	dumpNum++;
      } else {
	fprintf(stderr, "ERROR DiffusionFusion: Unrecognized command `%s'.\n", cmd.cs());
	exit(-1);
      }
    }

    // Count the commands that can be used multiple times.
    fieldNum = cmdField.length();
    priorNum = cmdPrior.length();
    trajCostNum = cmdTrajCost.length();

    if ((outputNum==0 && outPreCmdLine.length()==0) || outputNum > 1) {
      fprintf(stderr, "ERROR DiffusionFusion: `%s' must have exactly one `output' command or -o must be set on the command line.\n", cmdFile.cs());
      exit(-1);
    }
    if (trajNum != 1) {
      fprintf(stderr, "ERROR DiffusionFusion: `%s' does not have exactly one `trajectory' command.\n", cmdFile.cs());
      exit(-1);
    }
    if (trajCostNum < 1) {
      fprintf(stderr, "ERROR DiffusionFusion: `%s' must have at least one `trajCost' command.\n", cmdFile.cs());
      exit(-1);
    }
    if (mcNum != 1) {
      fprintf(stderr, "ERROR DiffusionFusion: `%s' does not have exactly one `mc' command.\n", cmdFile.cs());
      exit(-1);
    }

    // Command line 'outputPrefix'
    if ( outPreCmdLine.length() > 0 ) {
      printf("Using output prefix specified on the command line.\n");
      outputPrefix = outPreCmdLine;
    }
    printf("outputPrefix %s\n", outputPrefix.cs());

    // Substitute 'outputPrefix' for '%' in dump files.
    for (int d = 0; d < dumpTypeNum; d++)
      dumpFile[d] = insertOutputPrefix(dumpFile[d]);

    printf("\nCOMMAND_COUNT\n");
    printf("  load %d\n", cmdLoad.length());
    printf("  loadMany %d\n", cmdLoadMany.length());
    printf("  loadEvents %d\n", cmdLoadEvents.length());
    printf("  field %d\n", fieldNum);
    printf("  prior %d\n", priorNum);
    printf("  mc %d\n", mcNum);
    printf("  trajCost %d\n", trajCostNum);
    printf("  dump %d\n", dumpNum);

    // Try to flush the output before loading files.
    fflush(stdout);
    fflush(stderr);    

    printf("\nLoading the trajectory files.\n");
    prepareTrajectory(*cmdList[cmdTraj]);
    prepareData(cmdLoad, cmdLoadMany, cmdLoadEvents);

    // Output trajectory information.
    printEventStats();
    if (dumpFile[dumpEvents].length() > 0) writeEvents(dumpFile[dumpEvents]);
    // We assume that the relevant coordinate is 0.
    if (dumpFile[dumpBias].length() > 0) writeBiases(dumpFile[dumpBias], 0);
 
    printf("\nPreparing fields.\n");
    fieldList = new Field*[fieldNum];
    fieldDesc = new FieldDesc[fieldNum];
    saveList = new Field*[fieldNum];
    // Initialize the fields.
    for (int f = 0; f < fieldNum; f++) {
      fieldList[f] = NULL;
      saveList[f] = NULL;
      prepareField(*cmdList[cmdField.get(f)], f);
    }
    printf("\nFIELDS\n");
    for (int f = 0; f < fieldNum; f++) {
      printf("  %s\n",fieldDesc[f].toString().cs());
      if (dumpFile[dumpFields].length() > 0) {
	char cs[STRLEN];
	snprintf(cs,STRLEN,"%s.%s.dat",dumpFile[dumpFields].cs(),fieldDesc[f].name.cs());
	fieldList[f]->write(String(cs));
	snprintf(cs,STRLEN,"%s.%s.err",dumpFile[dumpFields].cs(),fieldDesc[f].name.cs());
	fieldList[f]->writeErr(String(cs));
      }
    }

    // Prepare the priors.
    priorList = new Prior*[priorNum];
    printf("\nPRIORS\n");
    double* pc = new double[priorNum];
    for (int p = 0; p < priorNum; p++)
      pc[p] = preparePrior(*cmdList[cmdPrior.get(p)], p);
    for (int p = 0; p < priorNum; p++)
      printf("INITIAL_PRIOR_COST %d %.15g\n", p, pc[p]);
    printf("\n");
    delete[] pc;

    // Initialize the random number generator.
    rando = new Random();

    // Prepare what's left.
    prepareMonteCarlo(*cmdList[cmdMc]);

    // Prepare the trajectory cost computers.
    trajCostList = new TrajCostComputer*[trajCostNum];
    tcFieldOn = new bool*[trajCostNum];
    for (int tc = 0; tc < trajCostNum; tc++) {
      tcFieldOn[tc] = new bool[fieldNum];
      prepareTrajCost(*cmdList[cmdTrajCost.get(tc)], tc);
    }

    // Erase all data in the output trajectory file.
    initOutputTraj();
  }

  static void printUsage() {
    fprintf(stdout, "Command usage:\n\n");
    fprintf(stdout, "trajectory [-time varName -col colInd] [-coor varName -col colInd] [-force varName -col colInd]...\n");
    fprintf(stdout, "\t*Note: Replacing the column index with a period `.' makes a variable of all zeros.\n");
    fprintf(stdout, "\t*Note: An important restriction is that the coordinates must be given in the order of the axes of the fields.");
    fprintf(stdout, "\t e.g., the first -coor is assumed to the primary direction in all fields.\n");
    fprintf(stdout, "\nload trajFile [-stride dataStride] [-skip eventSkip] [-minVar colName -minVal val] [-maxVar colName -maxVal val] [-scaleVar colName -scaleVal val] [-periodicVar col -periodicMin qMin -periodicMax qMax] [-pos init|mid|final] [-bias biasHistoryFile] [-biasSmooth width] [-biasCount biasCountHistoryFile] [-biasFullSamples fullSamples] [-group group]\n");
    //fprintf(stdout, "\t*Note: -minVal and -maxVal filters apply to the initial of position of an event\n");
    fprintf(stdout, "\t*Note: -stride changes the time between frames, while -skip\n\t\tthrows out events after the interval has been calculated.\n");
    fprintf(stdout, "\t*Note: Events in the same group (specified by -group) must be loaded by contiguous load commands\n");
    fprintf(stdout, "\nloadMany trajFileFile [same options as load]\n");
    fprintf(stdout, "\t*Note: Using -bias with loadMany probably doesn't make sense -- all files will use the same bias history. No warning will be issued.\n");
    fprintf(stdout, "\nloadEvent [-skip eventSkip] [-offset remainder] eventFile\n");
    fprintf(stdout, "\t*Note: Assumed to already be correctly filtered.\n");
    fprintf(stdout, "\t*Note: loadEvent does not yet support event groups.\n");
    fprintf(stdout, "\nfield fieldName fieldType -step mcStep [-f initialFileName] [-periodic true|false] [-n numEntries] [-fixed fixedNodeFile] [-init initialValue] [-minVal val] [-maxVal val] [-var colName] [-err errorFile] [-outPmf force|prob] [-global true|false]\n");
    fprintf(stdout, "\t*Note: -outPmf force|prob allows you to write the negative integral\n\t\tor -kT log(prob), respectively, in addition writing the\n\t\tfield in the normal way.\n");
    fprintf(stdout, "\t*Note: -global permits parameters that affect all nodes.\n");
    fprintf(stdout, "\nprior scale|known|smooth|couple field [-ref refField] [-err uniformErr] [-grad gradientStd] [-dim gradientDimension] [-start startingFieldIndex] [-end endingFieldIndex] [-couple coupleField] [-std coupleStd] \n");
    fprintf(stdout, "\nmc field0 field1... -n numSteps -output outputPeriod -preview previewPeriod -update updatePeriod -tol relErrMax -seed randomSeed -cycle stepsPerCycle -global globalPeriod\n");
    //fprintf(stdout, "\ntrajCost ccg|ccgPmf|ccg2d|reflect|reflect2d|smolCrank|smolCrankBias|fracSmolCrank|smolCrankDual|fracSmolCrankDual|langevinExp|smoluchowski|simpleSmol|reflectSmol field0 field1... -`varName' costVarName [-kt thermalEnergy] [-leastLocal leastLocalField] [-timestep timestep] [-hop maxHop] [-dim primaryDimension] [-group group]\n");

    fprintf(stdout, "\ntrajCost ccg|reflect|ccg2d|reflect2d|smolCrank|fracSmolCrank field0 field1... -`varName' costVarName [-kt thermalEnergy] [-leastLocal leastLocalField] [-timestep timestep] [-hop maxHop] [-dim primaryDimension] [-group group]\n");
    fprintf(stdout, "\ndump events|fields|distro|best|eventCost|bias dumpFile\n");
    fprintf(stdout, "\t*Note: A '%%' in dumpFile is substituted with 'output'.\n");
  }

  ~DiffusionFusion() {
    if (cmdList != NULL) {
      for (int i = 0; i < cmdNum; i++) if (cmdList[i] != NULL) delete cmdList[i];
      delete[] cmdList;
    }

    if ( trajVarNum > 0) {
      delete[] trajVarCol;
      delete[] trajVarName;
    }

    if (biasHistoryList != NULL) {
      for (int i = 0; i < totalTrajFileNum; i++) if (biasHistoryList[i] != NULL) delete biasHistoryList[i];
      delete[] biasHistoryList;
    }

    if (biasFieldList != NULL) {
      // biasFieldList only points to stuff in biasHistoryList
      // DO NOT deallocate the pointers in biasFieldList, since they
      // were just deallocated above.
      delete[] biasFieldList;
    }

    if (event != NULL) delete event;

    if (fieldList != NULL) {
      for (int i = 0; i < fieldNum; i++)
	if (fieldList[i] != NULL) {
	  delete fieldList[i];
	  delete saveList[i];
	}
      delete[] fieldList;
      delete[] fieldDesc;
      delete[] saveList;
    }
   
    if (priorList != NULL) {
      for (int p = 0; p < priorNum; p++) delete priorList[p];
      delete[] priorList;
    }

    if (rando != NULL) delete rando;
    if (monte != NULL) delete monte;
    if (trajCostList != NULL) {

      for (int tc = 0; tc < trajCostNum; tc++) {
	delete trajCostList[tc];
	delete[] tcFieldOn[tc];
      }
      delete[] trajCostList;
      delete[] tcFieldOn;
    }
  }

  ///////////////////////////////////////////////////////////////////////
  // Run the ugly thing.
  void run() {
    int stepsPerCycle = monte->getStepsPerCycle();
    double costSum = 0.0;
    double costSumSq = 0.0;
    int count = 0;
    int preview = 0;
    TrialMove trialMove;
    char outFile[STRLEN];
    
    // Initial cost.
    //double lastCost = trajCost->calcCost();
    double lastCost = 0.0;
    for (int tc = 0; tc < trajCostNum; tc++)
      lastCost += trajCostList[tc]->updateLocal();
    for (int p = 0; p < priorNum; p++)
      lastCost += priorList[p]->calcCost();
    printf("INITIAL_COST %.15g\n", lastCost);
    double costMin = lastCost;
    double cycleTime = omp_get_wtime();

    // Field move counts and acceptance ratios.
    int* fieldAccept = new int[fieldNum];
    int* fieldCount = new int[fieldNum];
    for (int f = 0; f < fieldNum; f++) {
      fieldCount[f] = 0;
      fieldAccept[f] = 0;
    }

    // Try to write the log before we begin.
    fflush(stdout);
    fflush(stderr);

    // A cycle consists of mcCycles MC steps, where mcCycles is the
    // number field nodes modified by the MC procedure.
    for (int cycle = 1; cycle <= mcCycles; cycle++) {
      // The inner MC step loop.
      for (int step = 0; step < stepsPerCycle; step++) {
	// Make the move.
	trialMove = monte->trialMove(step % globalPeriod);
	fieldCount[trialMove.fieldId]++;

	if (!fieldDesc[trialMove.fieldId].global) {
	  // Trajectory cost.
	  // Calculate the change only locally.
	  double deltaCost = 0.0;

	  // We need not update trajCostComputers that don't depend the field.
	  for (int tc = 0; tc < trajCostNum; tc++) {
	    if (tcFieldOn[tc][trialMove.fieldId])
	      deltaCost += trajCostList[tc]->deltaCost(trialMove);
	  }

	  // Prior cost.
#pragma omp parallel for schedule(dynamic) reduction(+:deltaCost)
	  for (int p = 0; p < priorNum; p++)
	    deltaCost += priorList[p]->deltaCost(trialMove);

	  // Metropolis accept or reject.
	  if (!(deltaCost != deltaCost) && deltaCost <= std::numeric_limits<double>::max() && monte->metropolis(deltaCost)) {
	    // Accept
	    lastCost += deltaCost;
	    fieldAccept[trialMove.fieldId]++;
	  } else {
	    // We reverse the move.
	    monte->reject();
	    // We also need to revert node cost in trajCost.
	    for (int tc = 0; tc < trajCostNum; tc++) {
	      if (tcFieldOn[tc][trialMove.fieldId])
		trajCostList[tc]->revert(trialMove);
	    }
	  }

	} else {
	  // We do things differently for global fields.
	  // Everything is calculated explicitly.
	  double currCost = 0.0;

	  for (int tc = 0; tc < trajCostNum; tc++) {
	    if (tcFieldOn[tc][trialMove.fieldId])
	      currCost += trajCostList[tc]->updateLocal();
	  }

#pragma omp parallel for schedule(dynamic) reduction(+:currCost)
	  for (int p = 0; p < priorNum; p++)
	    currCost += priorList[p]->calcCost();
	  double deltaCost = currCost - lastCost;

	  // Metropolis accept or reject.
	  if (monte->metropolis(deltaCost)) {
	    // Accept
	    lastCost = currCost;
	    fieldAccept[trialMove.fieldId]++;
	  } else {
	    // We reverse the move.
	    monte->reject();
	    // We also need to revert node cost in trajCost.

	    for (int tc = 0; tc < trajCostNum; tc++) {
	      if (tcFieldOn[tc][trialMove.fieldId])
		trajCostList[tc]->updateLocal();
	    }
	  }
	}

	count++;
	costSum += lastCost;
	costSumSq += lastCost*lastCost;
	if (lastCost < costMin) {
	  costMin = lastCost;
	  // Save the best version of each field.
	  for (int f = 0; f < fieldNum; f++) *saveList[f] = *fieldList[f];
	}
      } // end of cycle

      // Output the state.
      if (cycle % outputPeriod == 0) appendOutputTraj(outputTraj, cycle, lastCost);

      // Output information on the progress.
      if (cycle % updatePeriod == 0) {
	fflush(stderr);
	fflush(stdout);

	double costMean = costSum/count;
	double costStd = sqrt( (costSumSq - costSum*costSum/count)/(count-1) );
	printf("\ncycle %d\n", cycle);
	printf("  steps %d\n", count);

	// Timing.
	double tim = omp_get_wtime() - cycleTime;
	cycleTime =  omp_get_wtime();
	if (stepsPerCycle == 1)  printf("  trialMove fieldId %d node %d trialVal %g lastVal %g\n", trialMove.fieldId, trialMove.node, trialMove.trialVal, trialMove.lastVal);
	printf("  updateTime %.4g\n", tim);
	printf("  stepTime %.4g\n", tim/count);

	// Calculate the acceptance ratios.
	int accept = 0;
	for (int f = 0; f < fieldNum; f++) {
	  // If there are any moves for this field, calculate the acceptance ratio.
	  if (fieldCount[f] > 0) {
	    String nam = fieldDesc[f].name;
	    double ratio = double(fieldAccept[f])/fieldCount[f];
	    printf("  acceptRatio %s %.3g\n", nam.cs(), ratio);
	    accept += fieldAccept[f];
	  }
	  fieldCount[f] = 0;
	  fieldAccept[f] = 0;
	}
	printf("  acceptRatio TOTAL %.3g\n", double(accept)/count);

	printf("  costMin %.15g\n", costMin);
	printf("  costMean %.15g\n", costMean);
	printf("  costStd %.15g\n", costStd);

	// Trajectory cost.
	double costUpdate = 0.0;
	for (int tc = 0; tc < trajCostNum; tc++) {
	  double costTraj = trajCostList[tc]->updateLocal();
	  costUpdate += costTraj;
	  printf("  costTraj %d %.15g\n", tc, costTraj);
	}

	// Calculate the prior cost.
	for (int p = 0; p < priorNum; p++) {
	  double costPrior = priorList[p]->calcCost();
	  costUpdate += costPrior;
	  printf("  costPrior %d %.15g\n", p, costPrior);
	}

	// Check that the cost has not drifted too much.
	double costRelErr = fabs((lastCost - costUpdate)/costStd);

	printf("  cost %.15g\n", lastCost);
	printf("  costUpdate %.15g\n", costUpdate);
	printf("  costErr %.15g\n", lastCost - costUpdate);
	printf("  costRelErr %.15g\n", costRelErr);
	if ( costRelErr > relErrMax ) {
	  fprintf(stdout, "ERROR Large difference between locally accumulated cost and global cost.\n");
	  fprintf(stdout, "|costAcc - costGlobal|/costStd = %.15g\n", costRelErr);
	  fprintf(stdout, "Are global fields declared with `-global 1'?\n");
	  fprintf(stdout, "Try setting `-leastLocal' to the grid with the widest grid spacing,\n");
	  fprintf(stdout, "or making the grid spacing of your MC fields identical.\n");
	  fprintf(stdout, "Other options include decreasing `-update' or `-cycle'.\n");
	  fprintf(stdout, "For Smoluchowski solvers, you may need to increase `-hop'.\n");
	  fprintf(stdout, "Could also indicate a bug in the chosen trajCost method.\n");
	  if (stepsPerCycle == 1)  {
	    fprintf(stdout, "fieldId %d node %d trialVal %g lastVal %g\n", trialMove.fieldId, trialMove.node, trialMove.trialVal, trialMove.lastVal);
	    IndexList neigh = fieldList[trialMove.fieldId]->neighbors(trialMove.node);
	    fprintf(stdout,"NEIGH");
	    for (int n = 0; n < neigh.length(); n++) fprintf(stdout," %d", neigh.get(n));
	    fprintf(stdout,"\n");
	  }
	  if (cycle > 3*updatePeriod) {
	    fprintf(stderr, "ERROR Large difference between locally accumulated cost and global cost.\n");
	    exit(-1);
	  }
	}

	lastCost = costUpdate;
	count = 0;
	costSum = 0.0;
	costSumSq = 0.0;
      }

      // Write the current states of the fields.
      if (cycle % previewPeriod == 0) {
	// Write the current states of the fields.
	for (int i = 0; i < mcFieldSel.length(); i++) {
	  int f = mcFieldSel.get(i);
	  snprintf(outFile, STRLEN,"%s.%d.%s", outputPrefix.cs(), preview, fieldDesc[f].name.cs());
	  fieldList[f]->write(String(outFile));

	  // Convert force to pmf.
	  if (fieldDesc[f].outInt) {
	    snprintf(outFile, STRLEN,"%s.%d.pmf", outputPrefix.cs(), preview);
	    fieldList[f]->writeIntegral(String(outFile), -1.0);
	  }
	  // Convert probability to pmf.
	  if (fieldDesc[f].outLog) {
	    snprintf(outFile, STRLEN,"%s.%d.pmf", outputPrefix.cs(), preview);
	    fieldList[f]->writeLog(String(outFile), -trajCostList[0]->getKt());
	  }
	}
	
	preview++;
	// end of preview
      }

    } // Done with all cycles.

    // Write the best field found.
    if (dumpFile[dumpBest].length() > 0) {
      for (int f = 0; f < fieldNum; f++) {
	snprintf(outFile, STRLEN,"%s.%s", dumpFile[dumpBest].cs(), fieldDesc[f].name.cs());
	saveList[f]->write(outFile);

	// Convert force to pmf.
	if (fieldDesc[f].outInt) {
	  snprintf(outFile, STRLEN,"%s.%s.pmf", dumpFile[dumpBest].cs(), fieldDesc[f].name.cs());
	  fieldList[f]->writeIntegral(String(outFile), -1.0);
	}
	// Convert probability to pmf.
	if (fieldDesc[f].outLog) {
	  snprintf(outFile, STRLEN,"%s.%s.pmf", dumpFile[dumpBest].cs(), fieldDesc[f].name.cs());
	  fieldList[f]->writeLog(String(outFile), -trajCostList[0]->getKt());
	}
      }
    }

    // Write the gt distribution for the best file.
    if (dumpFile[dumpDistro].length() > 0) {
      // Copy best field to the current field.
      for (int f = 0; f < fieldNum; f++) *fieldList[f] = *saveList[f];
      // Write the gt distribution.
      writeGtDistro(dumpFile[dumpDistro]);
    }

    // Write the event cost.
    if (dumpFile[dumpEventCost].length() > 0) {
      // Copy best field to the current field.
      for (int f = 0; f < fieldNum; f++) *fieldList[f] = *saveList[f];
      // Copy best field to the current field.
      writeEventCost(dumpFile[dumpEventCost]);
    }

    delete[] fieldAccept;
    delete[] fieldCount;
    printf("Ran %d steps.\n", mcCycles*stepsPerCycle);
  }

  void appendOutputTraj(const String& fileName, int cycle, double cost) const {
    FILE* out = fopen(fileName.cs(), "a");
    if (out == NULL) {
      fprintf(stderr,"ERROR appendOutputTraj: Could not append trajectory file `%s'\n", fileName.cs());
      exit(-1);
    }

    // We first write the cycle number and the cost.
    fprintf(out, "CYCLE %d %.15g\n", cycle, cost);

    // Dump the fields to which the Monte Carlo is applied.
    for (int i = 0; i < mcFieldSel.length(); i++) {
      int f = mcFieldSel.get(i);
      fprintf(out, "FIELD %s %d\n", fieldDesc[f].name.cs(), fieldList[f]->length());
      fieldList[f]->dump(out);
    }
    fclose(out);
  }

  ///////////////////////////////////////////////////////////////////////
  // Initialization
private:
  // Read the trajectory variables.
  void prepareTrajectory(const CommandLineReader& cmd) {
    trajVarNum = cmd.getOptionNum()/2;
    printf("There are %d trajectory variables.\n", trajVarNum);
    trajVarName = new String[trajVarNum];
    trajVarCol = new int[trajVarNum];

    String* coorName = new String[trajVarNum];
    IndexList coorCol;
    String* timeName = new String[trajVarNum];
    IndexList timeCol;
    String* forceName = new String[trajVarNum];
    IndexList forceCol;

    for (int o = 0; o < cmd.getOptionNum(); o++) {
      String opt(cmd.getOption(o));
      String val(cmd.getOptionValue(o));
      String colString;

      if (opt == "coor" || opt == "time" || opt == "force" ) {
	if (o+1 >= cmd.getOptionNum() || cmd.getOption(o+1) != "col") {
	  fprintf(stderr, "ERROR trajectory: -coor/time/force must be followed by -col\n");
	  exit(-1);
	}
	colString = cmd.getOptionValue(o+1);
	int col;

	if ( !String::isNumeral(colString[0]) ) {
	  // "." means that the column taken to be all zeros.
	  // We represent this with a column number of -1.
	  if (colString == ".") col = -1;
	  else {
	    fprintf(stderr, "ERROR trajectory: Column index %s is not a nonnegative number or `.'\n", colString.cs());
	    printUsage();
	    exit(-1);
	  }
	} else col = atoi(colString.cs());
	
	if (opt == "coor") {
	  coorName[coorCol.length()] = val;
	  coorCol.add(col);
	} else if (opt == "time") {
	  timeName[timeCol.length()] = val;
	  timeCol.add(col);
	} else {
	  forceName[forceCol.length()] = val;
	  forceCol.add(col);
	}
	o++; // We already did the next option.
      } else {
	fprintf(stderr, "ERROR trajectory: Unrecognized option `%s'\n", opt.cs());
	exit(-1);
      }
    }
      
    // Load the vars into trajVarName and trajVarCol with the position variables first.
    // The position variables most be first to ease the use of Field::nearestNode()
    int v = 0;
    for (int p = 0; p < coorCol.length(); p++) {
      trajVarName[v] = coorName[p];
      trajVarCol[v] = coorCol.get(p);
      v++;
    }
    for (int p = 0; p < timeCol.length(); p++) {
      trajVarName[v] = timeName[p];
      trajVarCol[v] = timeCol.get(p);
      v++;
    }
    for (int p = 0; p < forceCol.length(); p++) {
      trajVarName[v] = forceName[p];
      trajVarCol[v] = forceCol.get(p);
      v++;
    }

    // Make the column list in the correct order!
    for (int var = 0; var < trajVarNum; var++) colList.add(trajVarCol[var]);

    printf("trajectory");
    for (int var = 0; var < trajVarNum; var++)
      printf(" (%s %d)", trajVarName[var].cs(), trajVarCol[var]);
    printf("\n");

  }

  // Make a big list of trajectory files.
  // Put these files in the TrajFile structure along with their filters with initTrajFile().
  // Read each file with loadAll(), which calls readTrajFile().
  void prepareData(const IndexList& cmdLoad, const IndexList& cmdLoadMany, const IndexList& cmdLoadEvents) {
    int fileNum = 0;
    TrajFile* fileList = NULL;

    // Count the trajectory files.
    IndexList loadManyNum;
    for (int i = 0; i < cmdLoadMany.length(); i++) {
      CommandLineReader* pCmd = cmdList[cmdLoadMany.get(i)];
      if (pCmd->getParamNum() != 1) {
	fprintf(stderr, "ERROR loadMany must have exactly one parameter, the file name.\n");
	printUsage();
	exit(-1);
      }

      int n = countLines(pCmd->getParam(0).cs(), IndexList(0));
      loadManyNum.add(n);
      printf("trajFileFile `%s' contains %d file names.\n", pCmd->getParam(0).cs(), n);
      fileNum += n;
    }
    fileNum += cmdLoad.length();
    printf("We will load a total of %d trajectory files.\n", fileNum);
    printf("TRAJFILENUM %d\n", fileNum);
    fileList = new TrajFile[fileNum];
    
    // Populate the trajectory file list.
    int f = 0;
    // First the "load" commands.
    for (int i = 0; i < cmdLoad.length(); i++) {
      CommandLineReader* pCmd = cmdList[cmdLoad.get(i)];
      if (pCmd->getParamNum() != 1) {
	fprintf(stderr, "ERROR load must have exactly one parameter, the file name.\n");
	printUsage();
	exit(-1);
      }

      String fileName = pCmd->getParam(0);      
      initTrajFile(fileName, *pCmd, fileList[f]);
      f++;
    }
    // Next the "loadMany" commands.
    for (int i = 0; i < cmdLoadMany.length(); i++) {
      CommandLineReader* pCmd = cmdList[cmdLoadMany.get(i)];
      int n = loadManyNum.get(i);
      String* fList = new String[n];
      
      // Read the files from the trajectory file file.
      readFileFile(pCmd->getParam(0), fList, n);

      // Add the files to the fileList.
      for (int j = 0; j < n; j++) {
	initTrajFile(fList[j], *pCmd, fileList[f]);
	f++;
      }

      delete[] fList;
    }

    // Collect the event file names.
    const int eventFileNum = cmdLoadEvents.length();
    printf("We will load a total of %d event files.\n", eventFileNum);
    printf("EVENTFILENUM %d\n", eventFileNum);
    String* eventFileList = new String[eventFileNum];
    int* eventFileSkip = new int[eventFileNum];
    int* eventFileOffset = new int[eventFileNum];
    for (int f = 0; f < eventFileNum; f++) {
      CommandLineReader* pCmd = cmdList[cmdLoadEvents.get(f)];
      if (pCmd->getParamNum() != 1) {
	fprintf(stderr, "ERROR loadEvents must have exactly one parameter, the file name.\n");
	printUsage();
	exit(-1);
      }
      
      // Set the event file and default options.
      eventFileList[f] = pCmd->getParam(0);
      eventFileSkip[f] = 1;
      eventFileOffset[f] = 0;

      // Check out the loadEvents options.
      for (int o = 0; o < pCmd->getOptionNum(); o++) {
	String opt(pCmd->getOption(o));
	String val(pCmd->getOptionValue(o));

	if (opt == "skip") {
	  eventFileSkip[f] = atoi(val);
	  if (eventFileSkip[f] < 1) {
	    fprintf(stderr, "ERROR loadEvents %s: -skip must be greater than or equal to 1\n", eventFileList[f].cs());
	    exit(-1);
	  }
	} else if (opt == "offset") {
	  eventFileOffset[f] = atoi(val);
	  if (eventFileOffset[f] < 0 || eventFileOffset[f] >= eventFileSkip[f] ) {
	    fprintf(stderr, "ERROR loadEvents %s: -offset %d must be 0 <= offset < skip=%d. A -skip command should precede -offset.\n", eventFileList[f].cs(), eventFileOffset[f], eventFileSkip[f]);
	    exit(-1);
	  }
	} else {
	  fprintf(stderr, "ERROR loadEvents %s: Unrecognized option `-%s'\n", eventFileList[f].cs(), opt.cs());
	  exit(-1);
	}
      } // Done with the options.
    } // Done with loadEvents
    
    // Actually load the files, applying filters.
    loadAll(fileList, fileNum, eventFileList, eventFileSkip, eventFileOffset, eventFileNum);

    // This code is to support bias history files.
    // Load the bias history for each file if necessary.
    // We do not yet support bias histories associated with event files.
    int biasHistoryCount = 0;
    for (int f = 0; f < fileNum; f++) if (fileList[f].biasFile.length()>1) biasHistoryCount++;
    printf("%d trajectory files have an associated bias history.\n", biasHistoryCount);

    if (biasHistoryCount > 0) {
      printf("Reading bias history files:\n");
      // While all files might not have a bias history, the array covers all,
      // with nulls where necessary.
      biasHistoryList = new HistoryReader*[fileNum];
      totalTrajFileNum = fileNum;
      // Not all these may be used.
      for (int f = 0; f < fileNum; f++) biasHistoryList[f] = NULL;

      biasFieldNum = 0;
      for (int f = 0; f < fileNum; f++) {
	if (fileList[f].biasFile.length() != 0) {
	  // We have a bias file, load it.
	  biasHistoryList[f] = new HistoryReader(fileList[f].biasFile.cs());

	  // Load the count file if we set fullSamples.
	  // We need to mimic how colvars handles fullSamples.
	  if (fileList[f].biasFullSamples > 0) {
	    if (fileList[f].biasCountFile.length() == 0) {
	      fprintf(stderr, "ERROR: File %d specifies -biasFullSamples > 0 without specifying a count history file with -biasCount\n", f);
	      exit(-1);
	    }
	    
	    int changed = biasHistoryList[f]->enforceFullSamples(fileList[f].biasCountFile.cs(), fileList[f].biasFullSamples);
	    if (f < dispFileMax) printf("Rescaled %d bias forces to imitate the colvars module with fullSamples=%d\n", changed,  fileList[f].biasFullSamples);
	  }

	  // Smooth the data.
	  if (fileList[f].biasSmoothWidth > 0.0) {
	    if (f < dispFileMax) printf("Applying Gaussian smoothing with width %g\n", fileList[f].biasSmoothWidth);
	    biasHistoryList[f]->smooth(fileList[f].biasSmoothWidth);
	    //if (f==0) biasHistoryList[f]->getBiasFrame(0)->write("bias.dat");
	  }

	  // Keep an index of bias fields.
	  biasFieldNum += biasHistoryList[f]->length();

	  if (f < dispFileMax)
	    printf("History file %s associated with trajectory file %s has %d bias frames.\n", fileList[f].biasFile.cs(), fileList[f].fileName.cs(), biasHistoryList[f]->length());
	}
      }

      // Now we must associate each event with a bias.
      biasFieldList = new const Piecewise1d*[biasFieldNum];
      int bfn = 0;
      for (int f = 0; f < fileNum; f++) {
	if (biasHistoryList[f] == NULL) continue;
	// We assume that the history file and the trajectory file
	// represent approximately the same interval.
	int histFrameN = biasHistoryList[f]->length();
	int eventN = fileList[f].eventLast - fileList[f].eventFirst + 1;
	double histEventRatio = histFrameN/double(eventN);
	printf("biasHistory: file %d historyFrames %d events %d histEventRatio %g\n", f, histFrameN, eventN, histEventRatio);

	int lastFrame = -1;
	for (int e = fileList[f].eventFirst; e <= fileList[f].eventLast; e++) {
	  int eventDelta = e - fileList[f].eventFirst;
	  //int histFrame = int(eventDelta*histEventRatio + 0.5);
	  //if (histFrame >= histFrameN) histFrame = histFrameN-1;
	  int histFrame = (histFrameN*eventDelta)/eventN;

	  if (histFrame != lastFrame) {
	    // Add a new bias field if necessary.
	    biasFieldList[bfn] = biasHistoryList[f]->getBiasFrame(histFrame);
	    bfn++;
	    lastFrame = histFrame;
	  }

	  // Link the event to the bias.
	  event[e].bias = bfn-1;
	}
      }
    }

    for (int f = 0; f < fileNum; f++) delete fileList[f].filter;
    delete[] fileList;
    delete[] eventFileList;
  }

  // This function extracts the options and parameters given in the `load' command
  // and applies it to the TrajFile structure.
  void initTrajFile(const String& fileName, const CommandLineReader& cmd, TrajFile& traj) {
    if (cmd.getParamNum() == 0) {
      fprintf(stderr,"ERROR %s has no parameters.\n", cmd.getCommand().cs());
      printUsage();
      exit(-1);
    }

    // Initialize the TrajFile object.
    traj.fileName = fileName.trim();
    traj.stride = 1;
    traj.skip = 1;
    traj.position = TrajFile::posInit;
    traj.group = 0;
    traj.filter = new TrajFilter(trajVarNum);
    traj.biasFile = "";
    traj.biasCountFile = "";
    traj.biasFullSamples = -1;
    traj.biasSmoothWidth = -1;

    // Add the options.
    for (int o = 0; o < cmd.getOptionNum(); o++) {
      String opt(cmd.getOption(o));
      String val(cmd.getOptionValue(o));

      if (opt == "stride") {
	traj.stride = atoi(val);
	if (traj.stride < 1) {
	  fprintf(stderr, "ERROR load %s: -stride must be greater than or equal to 1\n", fileName.cs());
	  exit(-1);
	}
      } else if (opt == "skip") {
	traj.skip = atoi(val);
	if (traj.skip < 1) {
	  fprintf(stderr, "ERROR load %s: -skip must be greater than or equal to 1\n", fileName.cs());
	  exit(-1);
	}
      } else if (opt == "pos") {
	traj.position = TrajFile::getPosIndex(val);
	if (traj.position < 1) {
	  fprintf(stderr, "ERROR load %s: -pos %s unrecognized.\n", fileName.cs(), val.cs());
	  exit(-1);
	}
      } else if (opt == "group") {
	traj.group = atoi(val);
	if (traj.group < 0) {
	  fprintf(stderr, "ERROR load %s: -group %d should be nonnegative.\n", fileName.cs(), traj.group);
	  exit(-1);
	}
      } else if (opt == "bias") {
	traj.biasFile = val;
      } else if (opt == "biasCount") {
	traj.biasCountFile = val;
      } else if (opt == "biasFullSamples") {
	traj.biasFullSamples = atoi(val);
      } else if (opt == "biasSmooth") {
	traj.biasSmoothWidth = strtod(val, NULL);
      } else if (opt == "scaleVar") {
	if (o+1 >= cmd.getOptionNum() || cmd.getOption(o+1) != "scaleVal") {
	  fprintf(stderr, "ERROR load %s: -scaleVar must be followed by -scaleVal\n", fileName.cs());
	  exit(-1);
	}
	int scaleVar = varIndex(val);
	double scaleVal = strtod(cmd.getOptionValue(o+1), NULL);
	traj.filter->addScale(scaleVar, scaleVal);
	o++; // We already did the next option.
      } else if (opt == "minVar") {
	if (o+1 >= cmd.getOptionNum() || cmd.getOption(o+1) != "minVal") {
	  fprintf(stderr, "ERROR load %s: -minVar must be followed by -minVal\n", fileName.cs());
	  exit(-1);
	}
	int minVar = varIndex(val);
	double minVal = strtod(cmd.getOptionValue(o+1), NULL);
	traj.filter->addMin(minVar, minVal);
	o++; // We already did the next option.
      } else if (opt == "maxVar") {
	if (o+1 >= cmd.getOptionNum() || cmd.getOption(o+1) != "maxVal") {
	  fprintf(stderr, "ERROR load %s: -maxVar must be followed by -maxVal\n", fileName.cs());
	  exit(-1);
	}
	int maxVar = varIndex(val);
	double maxVal = strtod(cmd.getOptionValue(o+1), NULL);
	traj.filter->addMax(maxVar, maxVal);
	o++; // We already did the next option.
      } else if (opt == "periodicVar") {
	if (o+2 >= cmd.getOptionNum() || cmd.getOption(o+1) != "periodicMin" || cmd.getOption(o+2) != "periodicMax" ) {
	  fprintf(stderr, "ERROR load %s: -periodicVar must be followed by -periodicMin and -periodicMax\n", fileName.cs());
	  exit(-1);
	}
	int pVar = varIndex(val);
	double pMin = strtod(cmd.getOptionValue(o+1), NULL);
	double pMax = strtod(cmd.getOptionValue(o+2), NULL);
	traj.filter->addPeriodic(pVar, pMin, pMax);
	o+=2; // We already did the next 2 options.
      } else {
	fprintf(stderr, "ERROR load %s: Unrecognized option `-%s'.\n", fileName.cs(), opt.cs());
	printUsage();
	exit(-1);
      }
    }
  }

  // Count the data entries in the files.
  // Initialize the event.
  // Load all files from the list of trajectory files using readTrajFile().
  // readTrajFile() also fills the event buffer and applies the filters.
  // filesList isn't const because we set eventFirst and eventLast.
  void loadAll(TrajFile* fileList, int fileNum, const String* eventFileList, const int* eventFileSkip, const int* eventFileOffset, int eventFileNum) {
    eventMax = 1;
    printf("\nCounting file lines:\n");
    for (int f = 0; f < fileNum; f++) {
      //int n =  countLines(fileList[f].fileName.cs(), colList);
      int n = countLines(fileList[f].fileName.cs());
      eventMax += n/(fileList[f].stride*fileList[f].skip) + 1;
      
      if (f < dispFileMax) printf("file %s lines %d\n", fileList[f].fileName.cs(), n);
      else if (f==dispFileMax) printf("Reading %d more files...\n", fileNum-dispFileMax);
    }
    // Add the event files.
    IndexList eventColList;
    for (int i = 0; i < trajVarNum; i++) eventColList.add(i);
    for (int f = 0; f < eventFileNum; f++) {
      int n = countLines(eventFileList[f].cs());
      eventMax += n/eventFileSkip[f] + 1;
      if (f < dispFileMax) printf("file %s lines %d\n", eventFileList[f].cs(), n);
      else if (f==dispFileMax) printf("Reading %d more files...\n", eventFileNum-dispFileMax);
    }

    // Allocate the event buffer.
    eventNum = 0;
    event = new Event[eventMax];
    printf("Maximum number of events %d\n\n", eventMax);

    // Fill the event list.
    printf("\nReading trajectory files:\n");
    for (int f = 0; f < fileNum; f++) {
      //if (f < dispFileMax) printf("Reading trajectory file `%s'\n", fileList[f].fileName.cs());
      // We keep a record of how each file was mapped into the event buffer (eventFirst, eventLast).
      fileList[f].eventFirst = eventNum;
      fileList[f].serialLast = readTrajFile(fileList[f], colList);
      fileList[f].eventLast = eventNum-1;
      if (f < dispFileMax) printf("trajFile %s events %d first %d last %d\n", fileList[f].fileName.cs(),  eventNum-fileList[f].eventFirst, fileList[f].eventFirst, fileList[f].eventLast);
      else if (f==dispFileMax) printf("Reading %d more files...\n", fileNum-dispFileMax);
    }

    // Read the event files directly.
    printf("\nReading event files:\n");
    for (int f = 0; f < eventFileNum; f++) {

      int n = readEventFile(eventFileList[f].cs(), eventFileSkip[f], eventFileOffset[f]);
      if (f < dispFileMax) printf("eventFile %s events %d\n", eventFileList[f].cs(), n);
      else if (f==dispFileMax) printf("Reading %d more event files...\n", eventFileNum-dispFileMax);
    }

   
  }

  void prepareField(const CommandLineReader& cmdLine, int fieldInd) {
    const int f = fieldInd;
    if (cmdLine.getParamNum() != 2) {
      fprintf(stderr, "ERROR prepareField: Invalid field command for field index %d.\n", f);
      printUsage();
      exit(-1);
    }
    
    // Initialize.
    fieldDesc[f].id = f;
    fieldDesc[f].name = cmdLine.getParam(0);
    fieldDesc[f].type = cmdLine.getParam(1);
    fieldDesc[f].hasMin = false;
    fieldDesc[f].hasMax = false;
    fieldDesc[f].step = -1.0;
    //fieldDesc[f].file; // empty string
    fieldDesc[f].n = 20; // default value
    fieldDesc[f].initVal = 0.0; // default value
    fieldDesc[f].periodic = false; // default value
    fieldDesc[f].outLog = false; // default value
    fieldDesc[f].outInt = false; // default value
    //fieldDesc[f].errFile; // empty string
    fieldDesc[f].global = false; // default value

    // Add the options.
    for (int o = 0; o < cmdLine.getOptionNum(); o++) {
      String opt(cmdLine.getOption(o));
      String val(cmdLine.getOptionValue(o));

      if (opt == "step") {
	fieldDesc[f].step = strtod(val, NULL);
      } else if (opt=="f") {
	fieldDesc[f].file = val.cs();
      } else if (opt=="n") {
	fieldDesc[f].n = atoi(val);
      } else if (opt=="init") {
	fieldDesc[f].initVal = strtod(val, NULL);
      } else if (opt=="periodic") {
	fieldDesc[f].periodic = readBoolean(val);
      } else if (opt=="minVal") {
	fieldDesc[f].hasMin = true;
	fieldDesc[f].minVal = strtod(val, NULL);
      } else if (opt=="maxVal") {
	fieldDesc[f].hasMax = true;
	fieldDesc[f].maxVal = strtod(val, NULL);
      } else if (opt=="fixed") {
	fieldDesc[f].fixedFile = val;
      } else if (opt=="err") {
	fieldDesc[f].errFile = val;
      } else if (opt=="var") {
	fieldDesc[f].varList.add(varIndex(val));
      } else if (opt=="global") {
	fieldDesc[f].global = readBoolean(val);
      } else if (opt=="outPmf") {
	if (val == "prob") fieldDesc[f].outLog = true;
	else if (val == "force") fieldDesc[f].outInt = true;
	else {
	  fprintf(stderr, "ERROR field %s: Unrecognized value for -outPmf, `%s' \n", fieldDesc[f].name.cs(), val.cs());
	  printUsage();
	  exit(-1);
	}

      } else {
	fprintf(stderr, "ERROR field %s: Unrecognized option `-%s' \n", fieldDesc[f].name.cs(), opt.cs());
	printUsage();
	exit(-1);
      }
    }

    if (fieldDesc[f].step < 0.0) {
      fprintf(stderr, "ERROR field %s: You must provide a positive value for `-step mcStep'.\n", fieldDesc[f].name.cs());
      printUsage();
      exit(-1);
    }

    if (fieldDesc[f].type == "cubic" || fieldDesc[f].type == "linear") {
      if (fieldDesc[f].file.length() > 0) {
	printf("Loading field %s from input file `%s'\n", fieldDesc[f].name.cs(), fieldDesc[f].file.cs());
	if (fieldDesc[f].type == "cubic") {
	  fieldList[f] = new PiecewiseCubic(fieldDesc[f].file.cs(), fieldDesc[f].periodic);
	  saveList[f] = new PiecewiseCubic( *(static_cast<const PiecewiseCubic*>(fieldList[f])) );
	} else {
	  fieldList[f] = new PiecewiseLinear(fieldDesc[f].file.cs(), fieldDesc[f].periodic);
	  saveList[f] = new PiecewiseLinear( *(static_cast<const PiecewiseLinear*>(fieldList[f])) );
	}

      } else {
	if (fieldDesc[f].varList.length() == 0) {
	  fprintf(stderr, "ERROR field %s: If you don't specify -f inputFile, you must specify (-var) a variable whose range serves as the range of the field.\n", fieldDesc[f].file.cs());
	  printUsage();
	  exit(-1);
	}
	
	// Get data bounds.
	double boundMin, boundMax;
	eventBounds(fieldDesc[f].varList.get(0), boundMin, boundMax);

	printf("Generating field %s using\n data bounds (%g %g), n=%d, and initVal=%g.\n", fieldDesc[f].name.cs(), boundMin, boundMax, fieldDesc[f].n, fieldDesc[f].initVal);
	double dx = 0.1*(boundMax-boundMin)/fieldDesc[f].n;
	if (fieldDesc[f].type == "cubic") {
	  fieldList[f] = new PiecewiseCubic(boundMin-dx, boundMax+dx, fieldDesc[f].initVal, fieldDesc[f].n, fieldDesc[f].periodic);
	  saveList[f] = new PiecewiseCubic( *(static_cast<const PiecewiseCubic*>(fieldList[f])) );
	} else {
	  fieldList[f] = new PiecewiseLinear(boundMin-dx, boundMax+dx, fieldDesc[f].initVal, fieldDesc[f].n, fieldDesc[f].periodic);
	  saveList[f] = new PiecewiseLinear( *(static_cast<const PiecewiseLinear*>(fieldList[f])) );
	}


      }

      // Load the error field if it's specified.
      if (fieldDesc[f].errFile.length() > 0) {
	PiecewiseCubic errField(fieldDesc[f].errFile, fieldList[f]->getPeriodic());
	const int n = errField.length();
	printf("Loading %d error values from `%s' into field `%s'.\n", n, fieldDesc[f].errFile.cs(), fieldDesc[f].name.cs());

	if (n != fieldList[f]->length()) {
	  fprintf(stderr, "ERROR field -err: ErrField `%s' does not have the same number of elements as `%s'.\n", fieldDesc[f].errFile.cs(), fieldDesc[f].name.cs());
	  exit(-1);
	}
	for (int i = 0; i < n; i++) fieldList[f]->setErr(i, errField.get(i));
      }
      
    } else if(fieldDesc[f].type == "bicubic") {
      if (fieldDesc[f].file.length() > 0) {
	printf("Loading field %s from input file `%s'\n", fieldDesc[f].name.cs(), fieldDesc[f].file.cs());
	fieldList[f] = new PiecewiseBicubic(fieldDesc[f].file.cs(), fieldDesc[f].periodic);
	saveList[f] = new PiecewiseBicubic( *(static_cast<const PiecewiseBicubic*>(fieldList[f])) );
      }
    } else {
      fprintf(stderr, "ERROR prepareField: Invalid field command for field index %d.\n", fieldInd);
      printUsage();
      exit(-1);
    }

    // Load the error field if it's specified.
    if (fieldDesc[f].errFile.length() > 0) {
      PiecewiseBicubic errField(fieldDesc[f].errFile, fieldList[f]->getPeriodic());
      const int n = errField.length();
      printf("Loading %d error values from `%s' into field %s.\n", n, fieldDesc[f].errFile.cs(), fieldDesc[f].name.cs());

      if (n != fieldList[f]->length()) {
	fprintf(stderr, "ERROR field -err: ErrField `%s' does not have the same number of elements as `%s'.\n", fieldDesc[f].errFile.cs(), fieldDesc[f].name.cs());
	exit(-1);
      }
      for (int i = 0; i < n; i++) fieldList[f]->setErr(i, errField.get(i));
    }

    // Determine which nodes are fixed.
    if (fieldDesc[f].fixedFile.length() != 0) {
      printf("Loading file `%s' marking fixed nodes for field %s.\n", fieldDesc[f].fixedFile.cs(), fieldDesc[f].name.cs());
      PiecewiseLinear fixed(fieldDesc[f].fixedFile.cs(),fieldDesc[f].periodic);
      if (fixed.length() != fieldList[f]->length()) {
	fprintf(stderr, "ERROR field -fixed: Fixed node file `%s' does not have the same number of elements as `%s'.\n", fieldDesc[f].fixedFile.cs(), fieldDesc[f].name.cs());
	exit(-1);
      }

      // Build the list of free nodes.
      for (int i = 0; i < fixed.length(); i++)
	if (fixed.get(i) <= 0.0) fieldDesc[f].freeNodeList.add(i);
      printf("Number of fixed nodes: %d\n", fieldList[f]->length()-fieldDesc[f].freeNodeList.length());
    } else {
      // Just put all the nodes in freeNodeList if there is no fixedFile.
      for (int i = 0; i < fieldList[f]->length(); i++) fieldDesc[f].freeNodeList.add(i);
    }

  }

  void prepareTrajCost(const CommandLineReader& cmdLine, int tc) {
    if (cmdLine.getParamNum() < 2) {
      fprintf(stderr, "ERROR prepareTrajCost: Invalid trajCost command.\n");
      printUsage();
      exit(-1);
    }
    String type = cmdLine.getParam(0);

    TrajCostDesc tcd(fieldList,event,eventNum);
    int tcFieldNum = cmdLine.getParamNum() - 1;
    
    // Make a list of the relevant fields for this tc.
    for (int i = 0; i < tcFieldNum; i++) {
      String fieldName = cmdLine.getParam(i+1);

      // Find the field with the matching name.
      int f;
      for (f = 0; f < fieldNum; f++) 
	if (fieldDesc[f].name == fieldName) break;
      if (f == fieldNum) {	
	fprintf(stderr, "ERROR prepareTrajCost: trajCost field `%s' has not been defined by a field command.\n", fieldName.cs());
	printUsage();
	exit(-1);
      }
      
      tcd.fieldSel.add(f);
    }

    // Set some defaults.
    String* optCostVar = new String[cmdLine.getOptionNum()];
    IndexList optTrajVar;

    // Add the options.
    for (int o = 0; o < cmdLine.getOptionNum(); o++) {
      String opt(cmdLine.getOption(o));
      String val(cmdLine.getOptionValue(o));

      if (opt == "kt") tcd.kbt = strtod(val, NULL);
      else if (opt == "timestep") tcd.timestep = strtod(val, NULL);
      else if (opt == "hop") tcd.maxHop = strtod(val, NULL);
      else if (opt == "dim") tcd.dimension = atoi(val);
      else if (opt == "group") tcd.group = atoi(val);
      else if (opt == "leastLocal") {
	// Find the field that correpsponds to the given name.
	int f;
	for (f = 0; f < fieldNum; f++) if (fieldDesc[f].name == val) break;
	if (f == fieldNum) {	
	  fprintf(stderr, "ERROR trajCost -leastLocal `%s' has not been defined by a field command.\n", val.cs());
	  //printUsage();
	  exit(-1);
	}
	tcd.leastLocal = f;
      } else {
	// See if the option value is a trajectory variable name.
	int v;
	for (v = 0; v < trajVarNum; v++) if (trajVarName[v] == val) break;

	if (v == trajVarNum) {
	  fprintf(stderr, "ERROR trajCost %s: Unrecognized option `-%s %s' or trajectory variable `%s' has not been defined.\n", type.cs(), opt.cs(), val.cs(), val.cs());
	  //printUsage();
	  exit(-1);
	}

	// We expect that this is an option of type `-costVar trajVar', but we'll have to check later.
	optCostVar[optTrajVar.length()] = opt;
	optTrajVar.add(v);
      }
    }

    // Default to the first field.
    if (tcd.leastLocal < 0) tcd.leastLocal = tcd.fieldSel.get(0);

    // Initialize the trajectory cost computer.
    // Just the basics.
    if (type == "ccg") {
      trajCostList[tc] = new TrajComer(tcd);
    } else if (type == "reflect") {
      trajCostList[tc] = new TrajReflect(tcd);
    } else if (type == "ccg2d") {
      if (tcd.dimension < 0) {
    	fprintf(stderr,"ERROR trajCost %s: You must define -dim\n", type.cs());
     	exit(-1);
      }
      trajCostList[tc] = new TrajComer2d(tcd);
    } else if (type == "reflect2d") {
      if (tcd.dimension < 0) {
    	fprintf(stderr,"ERROR trajCost %s: You must define -dim\n", type.cs());
     	exit(-1);
      }
      trajCostList[tc] = new TrajReflect2d(tcd);
    } else if (type == "smolCrank") {
      if (tcd.timestep < 0.0 || tcd.maxHop < 0.0) {
    	fprintf(stderr,"ERROR: trajCost smolCrank0 requires a timestep and maxHop. Add '-timestep' and '-hop' to the trajCost command.\n");
    	exit(-1);
      }
      trajCostList[tc] = new TrajSmolCrank(tcd);
    } else if (type == "fracSmolCrank") {
      if (tcd.timestep < 0.0 || tcd.maxHop < 0.0) {
    	fprintf(stderr,"ERROR: trajCost fracSmolCrank requires a timestep and maxHop. Add '-timestep' and '-hop' to the trajCost command.\n");
    	exit(-1);
      }
      trajCostList[tc] = new TrajFracSmolCrank(tcd);
    } else {
      fprintf(stderr, "ERROR trajCost: Unrecognized type `%s'\n", type.cs());
      printUsage();
      exit(-1);
    }


    // // Initialize the trajectory cost computer.
    // if (type == "ccgPmf"  || type == "pmf") {
    //   trajCostList[tc] = new TrajComerPmf(tcd);
    // } else if (type == "ccg") {
    //    trajCostList[tc] = new TrajComer(tcd);
    // } else if (type == "reflect") {
    //    trajCostList[tc] = new TrajReflect(tcd);
    // } else if (type == "noise") {
    //    trajCostList[tc] = new TrajNoise(tcd);
    // } else if (type == "foxFunc") {
    //    trajCostList[tc] = new TrajFoxFunction(tcd);
    // } else if (type == "ccgLog") {
    //   trajCostList[tc] = new TrajComerLog(tcd);
    // } else if (type == "langevin" || type == "bbkPmf") {
    //   trajCostList[tc] = new TrajLangevin(tcd);
    // } else if (type == "langevinGamma") {
    //   trajCostList[tc] = new TrajLangevinGamma(tcd);
    //  } else if (type == "turkcan") {
    //   trajCostList[tc] = new TrajTurkcan(tcd);
    //  } else if (type == "langevinExp") {
    //   trajCostList[tc] = new TrajLangevinExp(tcd);
    //  } else if (type == "langevinExp2d") {
    //   trajCostList[tc] = new TrajLangevinExp2d(tcd);
    // } else if (type == "ccg2d") {
    //   if (tcd.dimension < 0) {
    // 	fprintf(stderr,"ERROR trajCost %s: You must define -dim\n", type.cs());
    // 	exit(-1);
    //   }
    //   trajCostList[tc] = new TrajComer2d(tcd);
    // } else if (type == "reflect2d") {
    //   if (tcd.dimension < 0) {
    // 	fprintf(stderr,"ERROR trajCost %s: You must define -dim\n", type.cs());
    // 	exit(-1);
    //   }
    //   trajCostList[tc] = new TrajReflect2d(tcd);
    // } else if (type == "ccgPmf2d") {
    //    trajCostList[tc] = new TrajComerPmf2d(tcd);
    //    printf("Note that the PMF given by ccgPmf2d is defined only up to a constant function of the second dimension.\n");
    // } else if (type == "ccg2d2d") {
    //    trajCostList[tc] = new TrajComer2d2d(tcd);
    // } else if (type == "ccgPmf2d2d") {
    //    trajCostList[tc] = new TrajComerPmf2d2d(tcd);
    // } else if (type == "ccg2d2d2d") {
    //    trajCostList[tc] = new TrajComer2d2d2d(tcd);
    // } else if (type == "smoluchowski") {
    //   if (tcd.timestep < 0.0) {
    // 	fprintf(stderr,"ERROR: trajCost smoluchowski requires a timestep. Add '-timestep' to the trajCost command.\n");
    // 	exit(-1);
    //   }
    //   trajCostList[tc] = new TrajFdSmoluchowski(tcd, biasFieldList);
    // } else if (type == "simpleSmol") {
    //   if (tcd.timestep < 0.0 || tcd.maxHop < 0.0) {
    // 	fprintf(stderr,"ERROR: trajCost simpleSmol requires a timestep and maxHop. Add '-timestep' and '-hop' to the trajCost command.\n");
    // 	exit(-1);
    //   }
    //   trajCostList[tc] = new TrajSimpleSmol(tcd, biasFieldList);
    // } else if (type == "reflectSmol") {
    //   if (tcd.timestep < 0.0 || tcd.maxHop < 0.0) {
    // 	fprintf(stderr,"ERROR: trajCost reflectSmol requires a timestep and maxHop. Add '-timestep' and '-hop' to the trajCost command.\n");
    // 	exit(-1);
    //   }
    //   trajCostList[tc] = new TrajSimpleSmolReflect(tcd, biasFieldList);
    // } else if (type == "biasSmol") {
    //   if (tcd.timestep < 0.0 || tcd.maxHop < 0.0) {
    // 	fprintf(stderr,"ERROR: trajCost biasSmol requires a timestep and maxHop. Add '-timestep' and '-hop' to the trajCost command.\n");
    // 	exit(-1);
    //   }
    //   trajCostList[tc] = new TrajSmolReflectBias(tcd, biasFieldList);
    // } else if (type == "smolCrankBias") {
    //   if (tcd.timestep < 0.0 || tcd.maxHop < 0.0) {
    // 	fprintf(stderr,"ERROR: trajCost smolCrank requires a timestep and maxHop. Add '-timestep' and '-hop' to the trajCost command.\n");
    // 	exit(-1);
    //   }
    //   trajCostList[tc] = new TrajSmolCrankBias(tcd, biasFieldList);
    // } else if (type == "smolCrank0" || type == "smolCrank") {
    //   if (tcd.timestep < 0.0 || tcd.maxHop < 0.0) {
    // 	fprintf(stderr,"ERROR: trajCost smolCrank0 requires a timestep and maxHop. Add '-timestep' and '-hop' to the trajCost command.\n");
    // 	exit(-1);
    //   }
    //   trajCostList[tc] = new TrajSmolCrank(tcd);
    // } else if (type == "smolCrankDual") {
    //   if (tcd.timestep < 0.0 || tcd.maxHop < 0.0) {
    // 	fprintf(stderr,"ERROR: trajCost smolCrankDual requires a timestep and maxHop. Add '-timestep' and '-hop' to the trajCost command.\n");
    // 	exit(-1);
    //   }
    //   trajCostList[tc] = new TrajSmolCrankDual(tcd);
    // } else if (type == "smolCrank1") {
    //   if (tcd.timestep < 0.0 || tcd.maxHop < 0.0) {
    // 	fprintf(stderr,"ERROR: trajCost smolCrank0 requires a timestep and maxHop. Add '-timestep' and '-hop' to the trajCost command.\n");
    // 	exit(-1);
    //   }
    //   trajCostList[tc] = new TrajSmolCrank1(tcd);
    // } else if (type == "fracSmolCrank") {
    //   if (tcd.timestep < 0.0 || tcd.maxHop < 0.0) {
    // 	fprintf(stderr,"ERROR: trajCost fracSmolCrank requires a timestep and maxHop. Add '-timestep' and '-hop' to the trajCost command.\n");
    // 	exit(-1);
    //   }
    //   trajCostList[tc] = new TrajFracSmolCrank(tcd);
    // } else if (type == "fracSmolCrankDual") {
    //   if (tcd.timestep < 0.0 || tcd.maxHop < 0.0) {
    // 	fprintf(stderr,"ERROR: trajCost fracSmolCrankDual requires a timestep and maxHop. Add '-timestep' and '-hop' to the trajCost command.\n");
    // 	exit(-1);
    //   }
    //   trajCostList[tc] = new TrajFracSmolCrankDual(tcd);
    // } else {
    //   fprintf(stderr, "ERROR trajCost: Unrecognized type `%s'\n", type.cs());
    //   printUsage();
    //   exit(-1);
    // }

    // Initialize the trajectory variables.
    trajCostList[tc]->initCostVars(optCostVar, optTrajVar);

    printf("\nTRAJCOST %s\n", type.cs());
    for (int f = 0; f < tcFieldNum; f++) {
      int ifi = tcd.fieldSel.get(f);
      printf("  field %d %s \"%s\"\n", ifi, trajCostList[tc]->fieldName(f).cs(), fieldDesc[ifi].name.cs());
    }
    for (int v = 0; v < trajCostList[tc]->getTrajVarMin(); v++) {
      int ev = trajCostList[tc]->eventVarIndex(v);
      printf("  var %d %s \"%s\"\n", v, trajCostList[tc]->eventVarName(v).cs(), trajVarName[ev].cs());
    }
    printf("  events %d\n", trajCostList[tc]->getEventNum());

    printf("  kT %.15g\n", trajCostList[tc]->getKt());
    if (tcd.dimension >= 0) printf("  dimension %d\n", tcd.dimension);
    if (tcd.timestep >= 0.0) printf("  timestep %.15g\n", tcd.timestep);
    if (tcd.group >= 0) printf("  group %d\n", tcd.group);

    printf("  leastLocalField %d \"%s\" \n", tcd.leastLocal, fieldDesc[tcd.leastLocal].name.cs());

    // Fill the array to know which fields are relevant.
    for (int f = 0; f < fieldNum; f++) {
      if (tcd.fieldSel.find(f) >= 0) tcFieldOn[tc][f] = true;
      else tcFieldOn[tc][f] = false;
      //printf("trajCost %d field %d %s\n", tc, f, tcFieldOn[tc][f]?"on":"off");
    }

    printf("INITIAL_TRAJCOST %d %.15g\n", tc, trajCostList[tc]->calcCost());

    delete[] optCostVar;
  }

  void prepareMonteCarlo(const CommandLineReader& cmdLine) {
    // Initialize the Monte Carlo object.
    if (cmdLine.getParamNum() == 0) {
      fprintf(stderr, "ERROR mc: The Monte Carlo must act on at least one field.\n");
      printUsage();
      exit(-1);
    }
    
    outputPeriod = 10; // default value
    previewPeriod = 100; // default value
    updatePeriod = 1; // default value
    relErrMax = 1e-3; // default value
    globalPeriod = 10; // default value
    long int seed = 0;
    int stepsPerCycle = -1; // stepsPerCycle defaults to the total nodes in all MC fields

    for (int o = 0; o < cmdLine.getOptionNum(); o++) {
      String opt(cmdLine.getOption(o));
      String val(cmdLine.getOptionValue(o));

      if (opt == "n") {
	mcCycles = atoi(val);
      } else if (opt=="output") {
	outputPeriod = atoi(val);
      } else if (opt=="preview") {
	previewPeriod = atoi(val);
      } else if (opt=="update") {
	updatePeriod = atoi(val);
      } else if (opt=="tol") {
	relErrMax = strtod(val, NULL);
      } else if (opt=="seed") {
	seed = atol(val.cs());
      } else if (opt=="global") {
	globalPeriod = atoi(val);
      } else if (opt=="cycle") {
	stepsPerCycle = atoi(val);
      } else {
	fprintf(stderr, "ERROR mc: Unrecognized option `%s'.\n", opt.cs());
	printUsage();
	exit(-1);

      }
    }

    printf("\nMONTE CARLO");
    for (int p = 0; p < cmdLine.getParamNum(); p++) {
      // Find the field with the name given to the mc command.
      String mcFieldName = cmdLine.getParam(p);
      int f;
      for (f = 0; f < fieldNum; f++) if (fieldDesc[f].name == mcFieldName) break;
      if (f == fieldNum) {
	fprintf(stderr, "ERROR mc: Could not find field with name `%s'.\n", mcFieldName.cs());
	exit(-1);
      }
      
      mcFieldSel.add(f);
      printf(" field %d %s", f, mcFieldName.cs());
    }
    printf("\n");

    // We set the random seed in MC,
    // although the Random object is now owned by the DiffusionFusion object.
    if (seed == 0) seed = (unsigned int)time((time_t *)NULL);
    rando->init(seed);
    printf("Random seed %ld\n", seed);

    monte = new MetroMonteCarlo(fieldList, fieldDesc, mcFieldSel, rando);
    printf("MC fields: %d\n", monte->getFieldNum());
    printf("MC field nodes: %d\n", monte->getStepsPerCycle());
    printf("MC local fields: %d\n", monte->getLocalFieldNum());
    if (stepsPerCycle > 0) monte->setStepsPerCycle(stepsPerCycle);
    printf("MC steps per cycle: %d\n", monte->getStepsPerCycle());
    printf("Period for attempting moves with global fields: %d\n", globalPeriod);
    printf("Period for writing output trajectory: %d\n", outputPeriod);
    printf("Period for writing preview files: %d\n", previewPeriod);
    printf("Period for calculating statistics: %d\n", updatePeriod);    
    printf("Maximum permitted error between the locally accumulated cost and the global cost: %g\n", relErrMax); 
  }

  double preparePrior(const CommandLineReader& cmdLine, int pi) {
    // Initialize the Monte Carlo object.
    if (cmdLine.getParamNum() != 2) {
      fprintf(stderr, "ERROR prior: Prior must be 'prior priorType field [OPTIONS]'\n");
      printUsage();
      exit(-1);
    }

    // Find the field of the prior.
    String priorFieldName = cmdLine.getParam(1);
    int f;
    for (f = 0; f < fieldNum; f++) if (fieldDesc[f].name == priorFieldName) break;
    if (f == fieldNum) {
      fprintf(stderr, "ERROR prior: Could not find field with name `%s'.\n", priorFieldName.cs());
      exit(-1);
    }
    // Construct the object.
    priorList[pi] = new Prior(cmdLine.getParam(0), fieldList[f], f);

    double uniformErr = -1.0;
    String refFieldName;
    String coupleFieldName;

    switch(priorList[pi]->getType()) {

      // Scale invariance prior.
    case Prior::scaleType:
      if (cmdLine.getOptionNum() > 0) {
	fprintf(stderr, "ERROR 'prior scale' takes no options.\n");
	exit(-1);
      }
      break;

      // Known field prior.
    case Prior::knownType:
      for (int o = 0; o < cmdLine.getOptionNum(); o++) {
	String opt(cmdLine.getOption(o));
	String val(cmdLine.getOptionValue(o));
	
	if (opt == "ref") {
	  int f;
	  for (f = 0; f < fieldNum; f++) if (fieldDesc[f].name == val) break;
	  if (f == fieldNum) {
	    fprintf(stderr, "ERROR 'prior known -ref' Could not find field with name `%s'.\n", val.cs());
	    printUsage();
	    exit(-1);
	  }
	  refFieldName = val;
	  priorList[pi]->setRefField(fieldList[f]);
	} else if (opt == "err") {
	  uniformErr = strtod(val, NULL);
	} else if (opt == "start") {
	  if ( !priorList[pi]->setStartIndex(atoi(val)) ) 
	    fprintf(stderr,"ERROR `prior known' (%d): Index %s outside of field bounds\n", pi, val.cs());
	} else if (opt == "end") {
	  if ( !priorList[pi]->setEndIndex(atoi(val)) ) 
	    fprintf(stderr,"ERROR `prior known' (%d): Index %s outside of field bounds\n", pi, val.cs());
	} else {
	  fprintf(stderr, "ERROR 'prior known' Unrecognized option `%s'.\n", opt.cs());
	  exit(-1);
	}
      }

      if (!priorList[pi]->hasRefField()) {
	fprintf(stderr, "ERROR 'prior known' must specify a reference field '-ref'.\n");
	exit(-1);
      }
      
      if (uniformErr >= 0.0) {
	printf("Setting error of the reference field to %g.\n", uniformErr);
	priorList[pi]->setRefError(uniformErr);
      }
      break;

      // Smoothness prior.
    case Prior::smoothType:
      for (int o = 0; o < cmdLine.getOptionNum(); o++) {
	String opt(cmdLine.getOption(o));
	String val(cmdLine.getOptionValue(o));
	
	if (opt == "grad") {
	  priorList[pi]->setGradStd(strtod(val,NULL));
	} else if (opt == "dim") {
	  priorList[pi]->setGradDim(atoi(val));
	  priorList[pi]->setGradDim(atoi(val));
	} else if (opt == "start") {
	  if ( !priorList[pi]->setStartIndex(atoi(val)) ) 
	    fprintf(stderr,"ERROR `prior smooth' (%d): Index %s outside of field bounds\n", pi, val.cs());
	} else if (opt == "end") {
	  if ( !priorList[pi]->setEndIndex(atoi(val)) ) 
	    fprintf(stderr,"ERROR `prior smooth' (%d): Index %s outside of field bounds\n", pi, val.cs());
	} else {
	  fprintf(stderr,"ERROR `prior smooth' (%d): Unrecognized option `%s'\n", pi, opt.cs());
	  printUsage();
	  exit(-1);
	}
      }
      
      if (priorList[pi]->getGradStd() < 0.0) {
	fprintf(stderr,"ERROR `prior smooth' (%d) must include `-grad'\n", pi);
	exit(-1);
      }

      break;

 // Known field prior.
    case Prior::coupleType:
      for (int o = 0; o < cmdLine.getOptionNum(); o++) {
	String opt(cmdLine.getOption(o));
	String val(cmdLine.getOptionValue(o));
	
	if (opt == "couple") {
	  int f;
	  for (f = 0; f < fieldNum; f++) if (fieldDesc[f].name == val) break;
	  if (f == fieldNum) {
	    fprintf(stderr, "ERROR 'prior couple -couple' Could not find field with name `%s'.\n", val.cs());
	    printUsage();
	    exit(-1);
	  }
	  coupleFieldName = val;
	  priorList[pi]->setCoupleField(fieldList[f], f);
	} else if (opt == "std") {
	  priorList[pi]->setCoupleStd(strtod(val,NULL));
	} else if (opt == "start") {
	  if ( !priorList[pi]->setStartIndex(atoi(val)) ) 
	    fprintf(stderr,"ERROR `prior couple' (%d): Index %s outside of field bounds\n", pi, val.cs());
	} else if (opt == "end") {
	  if ( !priorList[pi]->setEndIndex(atoi(val)) ) 
	    fprintf(stderr,"ERROR `prior couple' (%d): Index %s outside of field bounds\n", pi, val.cs());
	} else {
	  fprintf(stderr, "ERROR 'prior couple' Unrecognized option `%s'.\n", opt.cs());
	  exit(-1);
	}
      }

      if (priorList[pi]->getCoupleStd() < 0.0) {
	fprintf(stderr,"ERROR `prior couple' (%d) must include `-std'\n", pi);
	exit(-1);
      }


      if (!priorList[pi]->hasCoupleField()) {
	fprintf(stderr, "ERROR 'prior couple' must specify a field to couple to '-couple'.\n");
	exit(-1);
      }
      break;
    }

    printf("  prior %s %s", priorList[pi]->getTypeName().cs(), fieldDesc[priorList[pi]->getFieldId()].name.cs());
    String text;
    char s[STRLEN];
    if (priorList[pi]->getGradStd() >= 0.0) {
      snprintf(s,STRLEN," gradientStandardDev %.15g",priorList[pi]->getGradStd());
      text.add(s);
      snprintf(s,STRLEN," gradientDimension %d",priorList[pi]->getGradDim());
      text.add(s);
    }
    if (priorList[pi]->hasRefField()) {
      text.add(" referenceField ");
      text.add(refFieldName);
    }
    if (priorList[pi]->hasCoupleField()) {
      text.add(" coupleField ");
      text.add(coupleFieldName);
    }
    printf("%s\n",text.cs());

    if (dumpFile[dumpRefField].length() > 0 && priorList[pi]->hasRefField()) {
      String fileName(dumpFile[dumpRefField]);
      fileName.add(".");
      fileName.add(refFieldName);
      priorList[pi]->writeRefField(fileName);
    }

    return priorList[pi]->calcCost();
  }

  void initOutputTraj() {
    outputTraj = outputPrefix;
    outputTraj.add(".traj");

    // Initialize the file.
    FILE* out = fopen(outputTraj, "w");
    if (out == NULL) {
      fprintf(stderr,"ERROR initOutputTraj Could not write output trajectory file `%s'\n", outputTraj.cs());
      exit(-1);
    }
    fclose(out);
  }

  ///////////////////////////////////////////////////////////////////////
  // Public utilities
  // Gets the column of the trajectory files associated with the column (variable) name.
  int columnIndex(const String& colName) const {
    for (int i = 0; i < trajVarNum; i++) {
      if (colName == trajVarName[i]) return trajVarCol[i];
    }
    fprintf(stderr, "ERROR `%s' is not a trajectory variable\n", colName.cs());
    exit(-1);
  }
  // Gets the variable associated with the column (variable) name.
  int varIndex(const String& varName) const {
    for (int i = 0; i < trajVarNum; i++) {
      if (varName == trajVarName[i]) return i;
    }
    fprintf(stderr, "ERROR `%s' is not a trajectory variable\n", varName.cs());
    exit(-1);
  }

  void eventBounds(int var, double& boundMin, double& boundMax) const {
    if (eventNum == 0) {
      fprintf(stderr, "ERROR eventBounds called with no events.\n");
      exit(-1);
    }

    boundMin = event[0].var[var];
    boundMax = event[0].var[var];
    for (int e = 1; e < eventNum; e++) {
      if (event[e].var[var] < boundMin) boundMin = event[e].var[var];
      if (event[e].var[var] > boundMax) boundMax = event[e].var[var];
    }
  }

  // Inserts the output name at '%' in a fileName.
  String insertOutputPrefix(const String& s) const {
    const int n = s.length();
    String ret;
    for (int i = 0; i < n; i++) {
      if (s[i] == '%') ret.add(outputPrefix);
      else ret.add(s[i]);
    }
    return ret;
  }

  // Returns the value of the boolean or exits if it can't be interpreted.  
  static bool readBoolean(String s) {
    s.lower();
    if (s == String("true") || s == String("yes") || s == String("on") || s==String("1")) return 1;
    if (s == String("false") || s == String("no") || s == String("off") || s==String("0")) return 0;
    
    fprintf(stderr,"ERROR! Invalid boolean value `%s'.\n", s.val());
    exit(-1);
  }

  static bool checkLine(const String& s) {
    String s1 = s.trim();
    // Ignore blank lines.
    if (s1.length() == 0) return false;
    // Ignore comments.
    if (s1[0] == '#') return false;
    return true;
  }

  static int countLines(const char* fileName, IndexList columnList) {
    // Open the file for reading.
    FILE* inp = fopen(fileName,"r");
    if (inp == NULL) {
      fprintf(stderr,"Couldn't open file `%s' for reading.\n", fileName);
      exit(-1);
    }

    const int columns = columnList.length();
    if (columns < 1) {
      fprintf(stderr,"countLines columnList must have at least one item\n");
      exit(-1);
    }
    // Find the biggest column.
    int maxColumn = columnList.get(0);
    for (int i = 1; i < columns; i++)
      if (columnList.get(i) > maxColumn) maxColumn = columnList.get(i);

    char line[STRLEN];
    int l = 0;
    int ln = 0;
    while (fgets(line, STRLEN, inp) != NULL) {
      ln++;
      String s(line);
      if (!checkLine(s)) continue;

      // Check for the right number of columns.
      int tok = s.tokenCount();
      if (tok <= maxColumn) {
	fprintf(stderr, "Warning: line %d of file `%s' has too few columns (%d of %d).\n", ln, fileName, tok, maxColumn);
	continue;
      }

      // This line is valid.
      l++;
    }
    fclose(inp);

    return l;
  }

  // A fast version of countLines with no checks.
  // We are guaranteed to get a value greater than or equal to the number of lines.
  static int countLines(const char* fileName) {
    // Open the file for reading.
    FILE* inp = fopen(fileName,"r");
    if (inp == NULL) {
      fprintf(stderr,"Couldn't open file `%s' for reading.\n", fileName);
      exit(-1);
    }

    char line[STRLEN];
    int l = 0;
    while (fgets(line, STRLEN, inp) != NULL) l++;
    fclose(inp);

    return l;
  }


  static void calcStats(const double* data, int n, double& minVal, double& maxVal, double& mean, double& stdev) {
    minVal = data[0];
    maxVal = data[0];
    double sum = data[0];
    double sumSq = data[0]*data[0];

    for (int i = 1; i < n; i++) {
      if (data[i] < minVal) minVal = data[i];
      if (data[i] > maxVal) maxVal = data[i];
      sum += data[i];
      sumSq += data[i]*data[i];
    }
    mean = sum/n;
    stdev = sqrt( (sumSq - sum*sum/n)/(n-1) );
  }

static bool writeHistogram(const char* fileName, const double* val, int n, const Bin& hist) {
    int total = 0;
    double* count = new double[hist.n];

    // Go through each array item.
    for (int i = 0; i < n; i++) {
      int b = hist.getCell(val[i]);
      if (b >= 0 && b < hist.n) count[b]++;

      total++;
    }

    // Scale to obtain a probability density.
    double scale = 1.0/(hist.dx*total);
    for (int b = 0; b < hist.n; b++) count[b]*=scale;

    // Open the file for writing.
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      fprintf(stderr,"writeHistogram Couldn't open file `%s' for writing.\n", fileName);
      delete[] count;
      return false;
    }

    // Write the distribution.
    for (int b = 0; b < hist.n; b++) {
      double x = hist.getCellCenter(b);
      fprintf(out, "%.14g %.14g\n", x, count[b]);
    }

    fclose(out);
    delete[] count;
    return true;
  }

private:
  // Private utilities
  int readCommands(const char* fileName, int commandNum) {
    // Open the file for reading.
    FILE* inp = fopen(fileName,"r");
    if (inp == NULL) {
      fprintf(stderr,"ERROR Couldn't open file `%s' for reading.\n", fileName);
      exit(-1);
    }

    char line[STRLEN];
    int l = 0;
    int ln = 0;
    while (fgets(line, STRLEN, inp) != NULL) {
      ln++;
      String s(line);
      if (!checkLine(s)) continue;

      // Sanity check to not read more than in the buffer.
      if (l >= commandNum) {
	printf("Warning DiffusionFusion:readCommands Unexpected extra command `%s'.\n", line);
	break;
      }

      // This line is valid.
      cmdList[l] = new CommandLineReader(s, ln);
      l++;
    }
    fclose(inp);

    return l;
  }

  // Read a trajectory file.
  // Return the serial number of the last line read (that wasn't empty or a comment).
  int readTrajFile(const TrajFile& traj, const IndexList& colList) {
    // Open the file for reading.
    FILE* inp = fopen(traj.fileName,"r");
    if (inp == NULL) {
      fprintf(stderr,"ERROR Couldn't open file `%s' for reading.\n", traj.fileName.cs());
      exit(-1);
    }

    const int columns = colList.length();
    if (columns < 1) {
      fprintf(stderr,"ERROR colList must have at least one item\n");
      exit(-1);
    }
    // Find the biggest column.
    int maxColumn = colList.get(0);
    for (int i = 1; i < columns; i++)
      if (colList.get(i) > maxColumn) maxColumn = colList.get(i);

    //printf("%s group %d\n", traj.fileName.cs(), traj.group);

    char line[STRLEN];
    double varLast[Event::varMax];
    double varLast0[Event::varMax];
    int read = -1;
    int ln = 0;
    int eventCounter = 0;
    int serial = 0;

    /////
    // Read all the entries.
    while (fgets(line, STRLEN, inp) != NULL) {
      ln++;
      String s(line);
      if (!checkLine(s)) continue;
      read++;

      // Skip by stride.
      if (read % traj.stride != 0) continue;

      // Check for the right number of columns.
      int tokN = s.tokenCount();
      if (tokN <= maxColumn) {
	fprintf(stderr, "Warning: line %d of file `%s' has too few columns (%d of %d).\n", ln, traj.fileName.cs(), tokN, maxColumn);
	continue;
      }

      // Events can be placed in groups.
      event[eventNum].group = traj.group;

      // Serial tells you how far through the file we are--which
      // can be used to determine the closest bias field in the history file.
      event[eventNum].serial = serial;
      serial++;

      // Get the event data.
      String* tokList = new String[tokN];
      s.tokenize(tokList);
      
      // Read the current variable values and set varLast.
      for (int c = 0; c < columns; c++) {
	int col = colList.get(c);
	// Negative columns are set to zero.
	if (col < 0) {
	  event[eventNum].var[c] = 0.0;
	  event[eventNum].del[c] = 0.0;
	  event[eventNum].del0[c] = 0.0;
	} else {
	  double var = strtod(tokList[col], NULL);

	  // The first, the last, or the middle.
	  switch(traj.position) {
	  case TrajFile::posInit:
	    event[eventNum].var[c] = varLast[c];
	    break;
	  case TrajFile::posFinal:
	    event[eventNum].var[c] = var;
	    break;
	  case TrajFile::posMid:
	    event[eventNum].var[c] = 0.5*(var + varLast[c]);
	    break;
	  }

	  event[eventNum].del[c] = var - varLast[c];
	  event[eventNum].del0[c] = varLast[c] - varLast0[c];
	  //printf("del %g del0 %g\n",event[eventNum].del[c], event[eventNum].del0[c]);
	  varLast0[c] = varLast[c];
	  varLast[c] = var;
	}
      }

      // We throw out the first two events of each file since del and del0
      // weren't set.
      if (read < 2*traj.stride) continue;

      // Now that we've set varLast, we can check if we skip this event.
      // Skip differs from stride in that it doesn't change the interval of the data, just throws out events.
      eventCounter++;
      if (eventCounter % traj.skip != 0) continue;

      // Apply the filters.
      // Scale.
      for (int i = 0; i < traj.filter->getScaleNum(); i++) {
	int v = traj.filter->getScaleVar(i);
	event[eventNum].var[v] *= traj.filter->getScaleVal(i);
	event[eventNum].del[v] *= traj.filter->getScaleVal(i);
	event[eventNum].del0[v] *= traj.filter->getScaleVal(i);
      }
      // Wrap.
      for (int i = 0; i < traj.filter->getPeriodicNum(); i++) {
	int v = traj.filter->getPeriodicVar(i);
	double pMin = traj.filter->getPeriodicMin(i);
	double pMax = traj.filter->getPeriodicMax(i);
	event[eventNum].var[v] = Field::wrapReal(event[eventNum].var[v], pMin, pMax);
	event[eventNum].del[v] = Field::wrapRealDiff(event[eventNum].del[v], pMin, pMax);
	event[eventNum].del0[v] = Field::wrapRealDiff(event[eventNum].del0[v], pMin, pMax);
      }

      bool valid = true;
      // Min filter.
      for (int i = 0; i < traj.filter->getMinNum() && valid; i++) {
	int v = traj.filter->getMinVar(i);
	// Skip events with the initial pos < minVal.
	if (event[eventNum].var[v] < traj.filter->getMinVal(i)) valid = false;
      }
      // Max filter.
      for (int i = 0; i < traj.filter->getMaxNum() && valid; i++) {
	int v = traj.filter->getMaxVar(i);
	// Skip events with the initial pos > maxVal.
	if (event[eventNum].var[v] > traj.filter->getMaxVal(i)) valid = false;
      }

      // Store this event.
      if (valid) {
	eventNum++;

	// Sanity check.
	if (eventNum > eventMax) {
	  fprintf(stderr, "ERROR Overfilled event buffer reading `%s'\n", traj.fileName.cs());
	  fprintf(stderr, "eventNum %d eventMax %d\n", eventNum, eventMax);
	  exit(-1);
	}
      }
    }
    fclose(inp);

    return serial;
  }

  static int readFileFile(const char* fileFile, String* fileList, int readMax) {
    // Open the file for reading.
    FILE* inp = fopen(fileFile,"r");
    if (inp == NULL) {
      fprintf(stderr, "ERROR Couldn't open file `%s' for reading.\n", fileFile);
      exit(-1);
    }

    char line[STRLEN];
    int l = 0;
    while (fgets(line, STRLEN, inp) != NULL) {
      String s(line);
      if (!checkLine(s)) continue;
      //fileList[l] = s.range(0,-2);
      fileList[l] = s;

      l++;
    }
    fclose(inp);

    // Sanity check.
    if (l > readMax) {
      fprintf(stderr, "ERROR Found too many entries in fileFile `%s'\n", fileFile);
      exit(-1);
    }

    return l;
  }

  int readEventFile(const char* eventFile, int skip, int offset) {
    int eventNum0 = eventNum;
    // Open the file for reading.
    FILE* inp = fopen(eventFile,"r");
    if (inp == NULL) {
      fprintf(stderr, "ERROR Couldn't open file `%s' for reading.\n", eventFile);
      exit(-1);
    }

    char line[STRLEN];
    //const int trajVarNumTwo = trajVarNum*2;
    const int trajVarNumThree = trajVarNum*3;
    String tokenList[trajVarNumThree];
    int e = 0;
    while (fgets(line, STRLEN, inp) != NULL) {
      String s(line);
      if (!checkLine(s)) continue;

      // Sanity check.
      if ( eventNum >= eventMax ) {
	fprintf(stderr, "ERROR Found too many entries in eventFile `%s'\n", eventFile);
	exit(-1);
      }
 
      int numTokens = s.tokenCount();
      // +1 is for the serial number
      if (numTokens != trajVarNumThree && numTokens != trajVarNumThree+1) {
	fprintf(stderr,"ERROR Line of eventFile `%s' doesn't have %d entries:\n%s\n", eventFile, trajVarNumThree, line);
	exit(-1);
      }


      s.tokenize(tokenList);
      // Perhaps skip some entries.
      if (e % skip == offset) {
	// Put the event data directly into the event buffer.
	for (int v = 0; v < trajVarNum; v++) {
	  event[eventNum].var[v] = strtod(tokenList[3*v].cs(),NULL);
	  event[eventNum].del[v] = strtod(tokenList[3*v+1].cs(),NULL);
	  event[eventNum].del0[v] = strtod(tokenList[3*v+2].cs(),NULL);
	}
	eventNum++;
      }
      e++;

      // Next line.
    }
    fclose(inp);
    return eventNum - eventNum0;
  }


  // Diagnostic functions
public:
  bool writeEvents(String fileName) {
    // Open the file for writing.
    FILE* out = fopen(fileName.cs(),"w");
    if (out == NULL) {
      fprintf(stderr,"Warning: Couldn't open file `%s' for writing.\n", fileName.cs());
      return false;
    }
    
    fprintf(out, "#");
    for (int v = 0; v < trajVarNum; v++)
      fprintf(out, " %s d%s d%s0", trajVarName[v].cs(), trajVarName[v].cs(), trajVarName[v].cs());
    fprintf(out, " biasIndex");
    fprintf(out, "\n");

    for (int e = 0; e < eventNum; e++) {
      for (int v = 0; v < trajVarNum; v++) {
	fprintf(out, "%.14g %.14g %.14g ", event[e].var[v], event[e].del[v], event[e].del0[v]);
      }
      fprintf(out, "%d %d\n", event[e].group, event[e].bias);
    }

    fclose(out);
    return true;
  }

  bool writeBiases(String fileName, int var) {
    // Open the file for writing.
    FILE* out = fopen(fileName.cs(),"w");
    if (out == NULL) {
      fprintf(stderr,"Warning: Couldn't open file `%s' for writing.\n", fileName.cs());
      return false;
    }
    
    fprintf(out, "# BIAS(X)\n");
    for (int e = 0; e < eventNum; e++) {
      if (event[e].bias >= 0) {
	double bias = biasFieldList[event[e].bias]->computeZeroOrder(event[e].var[var]);
	fprintf(out, "%.14g %.14g %d\n", event[e].var[var], bias, event[e].bias);
      } else fprintf(out, "%.14g %.14g %d\n", event[e].var[var], 0.0, event[e].bias);
    }

    fclose(out);
    return true;
  }

  void printEventStats() {
    const int n = eventNum;

    printf("\nEVENT STATISTICS\n");
    printf("EVENTNUM %d\n", n);
    if (n == 0) return;

    for (int v = 0; v < trajVarNum; v++) {
      double varMin = event[0].var[v];
      double varMax = event[0].var[v];
      double varSum = 0.0;
      double varSumSq = 0.0;
      double delMin = event[0].del[v];
      double delMax = event[0].del[v];
      double delSum = 0.0;
      double delSumSq = 0.0;
      double mean, std;

      // Collect the statistics.
      for (int e = 0; e < n; e++) {
	if (event[e].var[v] < varMin) varMin = event[e].var[v];
	if (event[e].var[v] > varMax) varMax = event[e].var[v];
	varSum += event[e].var[v];
	varSumSq += event[e].var[v]*event[e].var[v];

	if (event[e].del[v] < delMin) delMin = event[e].del[v];
	if (event[e].del[v] > delMax) delMax = event[e].del[v];
	delSum += event[e].del[v];
	delSumSq += event[e].del[v]*event[e].del[v];
      }
      
      printf("  VARIABLE %s", trajVarName[v].cs());
      mean = varSum/n;
      std = sqrt( (varSumSq - varSum*varSum/n)/(n-1) );
      printf(" min %.6g", varMin);
      printf(" max %.6g", varMax);
      printf(" mean %.6g", mean);
      printf(" std %.6g\n", std);

      printf("  DELTA %s", trajVarName[v].cs());
      mean = delSum/n;
      std = sqrt( (delSumSq - delSum*delSum/n)/(n-1) );
      printf(" min %.6g", delMin);
      printf(" max %.6g", delMax);
      printf(" mean %.6g", mean);
      printf(" std %.6g\n", std);
    }
  }

  void writeGtDistro(const String& fileName) {
    for (int tc = 0; tc < trajCostNum; tc++) {
      // Get the gt value for each event.
      const int n =  trajCostList[tc]->getEventNum();
      double* gt = new double[n];
      for (int e = trajCostList[tc]->getEventStart(); e <= trajCostList[tc]->getEventEnd(); e++) {
	trajCostList[tc]->eventCost(e);
	gt[e] = trajCostList[tc]->gtNumer/sqrt(trajCostList[tc]->gtVar);
      }

      // Create bins. Five sigma!
      Bin normalBin(100, -5.0, 5.0);

      if (trajCostNum == 1) {
	writeHistogram(fileName.cs(), gt, n, normalBin);
      } else {
	// Specify which trajCostComputer with a number.
	char exten[128];
	snprintf(exten,128,".%d",tc);
      
	String fileName1;
	int dot = fileName.last('.');
	if (dot < 0) {
	  fileName1 = fileName;
	  fileName1.add(exten);
	} else {
	  fileName1 = fileName.range(0,dot-1);
	  fileName1.add(exten);
	  fileName1.add(fileName.range(dot,-1));
	}
      }

      delete[] gt;
    }
  }

  bool writeEventCost(String fileName) {
    // Open the file for writing.
    FILE* out = fopen(fileName.cs(),"w");
    if (out == NULL) {
      fprintf(stderr,"Warning: Couldn't open file `%s' for writing.\n", fileName.cs());
      return false;
    }
    for (int tc = 0; tc < trajCostNum; tc++) {
      for (int e = trajCostList[tc]->getEventStart(); e <= trajCostList[tc]->getEventEnd(); e++)
	fprintf(out, "%.14g\n", trajCostList[tc]->eventCost(e));
    }

    fclose(out);
    return true;
  }


private:
  // Don't permit.
  DiffusionFusion();
  DiffusionFusion(const DiffusionFusion&);
  void operator=(const DiffusionFusion&);
};

#endif
