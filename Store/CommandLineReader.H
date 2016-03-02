///////////////////////////////////////////////////////////////////////
// Configuration file reader
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef COMMANDLINEREADER_H
#define COMMANDLINEREADER_H

#include "useful.H"

class CommandLineReader {
private:
  String command;
  int paramNum;
  String* param;
  int optionNum;
  String* option;
  String* optionValue;
  int lineNo;

public:
  // Initialize from the command line.
  CommandLineReader(int argc, char** argv) : command(argv[0]), paramNum(0), optionNum(0), lineNo(0) {
    // Count the number of parameters and options.
    for (int i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
	optionNum++;
	if (i == argc-1) {
	  fprintf(stderr, "ERROR CommmandLineReader: Option %s has no value.\n", argv[i]);
	  exit(-1);
	}
	i++;
      }
      else paramNum++;
    }

    // Allocate the parameters and options.
    param = new String[paramNum];
    option = new String[optionNum];
    optionValue = new String[optionNum];

    int on = 0;
    int pn = 0;

    // Store the options and parameters.
    for (int i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
	option[on] = &(argv[i][1]);
	optionValue[on] = argv[i+1];
	i++;
	on++;
      }
      else {
	param[pn] = argv[i];
	pn++;
      }
    }
  }

  CommandLineReader(const char* line, int lineNo0) {
    CommandLineReader(String(line), lineNo0);
  }

  // Initialize with a string.
  CommandLineReader(String line, int lineNo0) : paramNum(0), optionNum(0), lineNo(lineNo0) {
    // Extract the tokens.
    int tokNum = line.tokenCount();
    String* tok = new String[tokNum];
    line.tokenize(tok);

    // Set the command.
    command = tok[0];

    // Count the number of parameters and options.
    for (int i = 1; i < tokNum; i++) {
      if (tok[i][0] == '-') {
	optionNum++;
	if (i == tokNum-1) {
	  fprintf(stderr,"ERROR CommmandLineReader: Option %s has no value.\n", tok[i].cs());
	  exit(-1);
	}
	i++;
      }
      else paramNum++;
    }
    //printf("command %s p %d o %d\n", command.cs(), paramNum, optionNum);

    // Allocate the parameters and options.
    param = new String[paramNum];
    option = new String[optionNum];
    optionValue = new String[optionNum];

    int on = 0;
    int pn = 0;

    // Store the options and parameters.
    for (int i = 1; i < tokNum; i++) {
      if (tok[i][0] == '-') {
	option[on] = tok[i].range(1,-1);
	optionValue[on] = tok[i+1];
	i++;
	on++;
      }
      else {
	param[pn] = tok[i];
	pn++;
      }
    }

    delete[] tok;
  }

  ~CommandLineReader() {
    delete[] param;
    delete[] option;
    delete[] optionValue;
  }
  
  String getCommand() const {
    return command;
  }

  int getParamNum() const {
    return paramNum;
  }

  int getOptionNum() const {
    return optionNum;
  }

  String getParam(int ind) const {
    if (ind >= 0 && ind < paramNum) return param[ind];
    return String();
  }

  String getOption(int ind) const {
    if (ind >= 0 && ind < optionNum) return option[ind];
    return String();
  }

  String getOptionValue(int ind) const {
    if (ind >= 0 && ind < optionNum) return optionValue[ind];
    return String();
  }

  int getLineNumber() const {
    return lineNo;
  }

private:
  // No default copying.
  CommandLineReader();
  CommandLineReader(const CommandLineReader&);
  CommandLineReader& operator=(const CommandLineReader&);
};

#endif
