//
//  Copyright 2015  Taha Masood, Johannes Tausch and Jerome Butler
//
//  Permission to use, copy, and distribute this software and its
//  documentation for any purpose with or without fee is hereby granted,
//  provided that the above copyright notice appear in all copies and
//  that both the copyright notice and this permission notice appear
//  in supporting documentation.
//
//  This software is provided "as is" without express or implied warranty
//  to the extent permitted by applicable law.
//
#include <iostream>
#include <fstream>
#include <sstream>

#include "structure.h"
#include "grating.h"
#include "ssystem.h"

using namespace std;

int readlayer(string layerbuf, structure *epiptr);
void trimspaces(string& str);

int readfile(char* filename, structure *epiptr, grating *gratptr,
	     ssystem *sysptr)
{
  ifstream infile;
  string buffer;
  string laybuf;
  string str;
  string::size_type position;
  string stype;
  string gs;

  double wvl;
  int gl;
  int fl;
  double prd;
  double dc, pl;
  int i;
  int startpos;
  int sssh;
  int cvpts;
  double zgi, zgr;
  double incr;
  int nsteps;
  double amp;
  double ph; // phase of counter propagating wave
  i = 0;
  infile.open(filename);
  if (!infile)
    {
      cout << "The input file does not exist. "
	   << "Program terminating .... " << endl;
      return 1;
    }
  
  while (getline(infile, buffer) != NULL)
    {
      //delete leading blank spaces and tabs
      trimspaces(buffer);

      if (buffer[0] == '#')
	{
	  //Its a comment
	}
      else if (buffer.find("WVL") != string::npos)
	{
	  if ((position = buffer.find('=')) != string::npos)
	    {
	      startpos = static_cast<unsigned int> (position);
	      str = buffer.substr((startpos + 1));
	      istringstream ins(str);
	      ins >> wvl;
	      epiptr->setwavelength(wvl);
	      cout << "Wavelength = " << wvl << " um" << endl;
	    }
	}
      else if (buffer.find("GS") != string::npos)
	{
	  if ((position = buffer.find('=')) != string::npos)
	    {
	      startpos = static_cast<unsigned int> (position);
	      str = buffer.substr((startpos + 1));
	      istringstream ins(str);
	      ins >> gs;
	      gratptr->setgratingshape(gs);
	      cout << "Grating shape = " << gs << endl;
	    }
	}
      else if (buffer.find("GL") != string::npos)
	{
	  if ((position = buffer.find('=')) != string::npos)
	    {
	      startpos = static_cast<unsigned int> (position);
	      str = buffer.substr((startpos + 1));
	      istringstream ins(str);
	      ins >> gl;
	      gratptr->setgllayer(gl);
	      cout << "Grating Layer = " << gl << endl;
	    }
	}
      else if (buffer.find("FL") != string::npos)
	{
	  if ((position = buffer.find('=')) != string::npos)
	    {
	      startpos = static_cast<unsigned int> (position);
	      str = buffer.substr((startpos + 1));
	      istringstream ins(str);
	      ins >> fl;
	      gratptr->setfllayer(fl);
	    }
	}
      else if (buffer.find("PRD") != string::npos)
	{
	  if ((position = buffer.find('=')) != string::npos)
	    {
	      startpos = static_cast<unsigned int> (position);
	      str = buffer.substr((startpos + 1));
	      istringstream ins(str);
	      ins >> prd;
	      gratptr->setperiod(prd);
	      cout << "Grating period = " << prd << " um" << endl;
	    }
	}
      else if (buffer.find("DC") != string::npos)
	{
	  if ((position = buffer.find('=')) != string::npos)
	    {
	      startpos = static_cast<unsigned int> (position);
	      str = buffer.substr((startpos + 1));
	      istringstream ins(str);
	      ins >> dc;
	      gratptr->setdc(dc);
	      cout << "Duty cycle = " << dc << endl;
	    }
	}
      else if (buffer.find("DX") != string::npos)
	{
	  if ((position = buffer.find('=')) != string::npos)
	    {
	      startpos = static_cast<unsigned int> (position);
	      str = buffer.substr((startpos + 1));
	      istringstream ins(str);
	      ins >> pl;
	      gratptr->setdx(pl);
	    }
	}
      else if (buffer.find("SH") != string::npos)
	{
	  if ((position = buffer.find('=')) != string::npos)
	    {
	      startpos = static_cast<unsigned int> (position);
	      str = buffer.substr((startpos + 1));
	      istringstream ins(str);
	      ins >> sssh;
	      gratptr->setssspchrmncs(sssh);
	    }
	}
      else if (buffer.find("NCP") != string::npos)
	{
	  if ((position = buffer.find('=')) != string::npos)
	    {
	      startpos = static_cast<unsigned int> (position);
	      str = buffer.substr((startpos + 1));
	      istringstream ins(str);
	      ins >> cvpts;
	      gratptr->setcvertpnts(cvpts);
	    }
	}
      else if (buffer.find("ZG") != string::npos)
	{
	  if ((position = buffer.find('=')) != string::npos)
	    {
	      startpos = static_cast<unsigned int> (position);
	      str = buffer.substr((startpos + 1));
	      istringstream ins(str);
	      ins >> zgr >> zgi;
	      epiptr->setguess(zgi, zgr);
	      cout << "Initial n_eff guess = (" << zgr << "," <<
		zgi << ")" << endl;
	    }
	}
      else if (buffer.find("STYPE") != string::npos)
	{
	  if ((position = buffer.find('=')) != string::npos)
	    {
	      startpos = static_cast<unsigned int> (position);
	      str = buffer.substr((startpos + 1));
	      istringstream ins(str);
	      ins >> stype;
	      epiptr->setstructtype(stype);
	      cout << "Structure type = " << stype << endl;
	    }
	}
      else if (buffer.find("INCR") != string::npos)
	{
	  if ((position = buffer.find('=')) != string::npos)
	    {
	      startpos = static_cast<unsigned int> (position);
	      str = buffer.substr((startpos + 1));
	      istringstream ins(str);
	      ins >> incr;
	      sysptr->setloopincr(incr);
	    }
	}
      else if (buffer.find("NSTEPS") != string::npos)
	{
	  if ((position = buffer.find('=')) != string::npos)
	    {
	      startpos = static_cast<unsigned int> (position);
	      str = buffer.substr((startpos + 1));
	      istringstream ins(str);
	      ins >> nsteps;
	      sysptr->setnumsteps(nsteps);
	    }
	}
      // counter propagating wave amplitude
      else if (buffer.find("AMP") != string::npos)
	{
	  if ((position = buffer.find('=')) != string::npos)
	    {
	      startpos = static_cast<unsigned int> (position);
	      str = buffer.substr((startpos + 1));
	      istringstream ins(str);
	      ins >> amp;
	      sysptr->setcpwamplitude(amp);
	    }
	}
      // counter propagating wave phase
      else if (buffer.find("PH") != string::npos)
	{
	  if ((position = buffer.find('=')) != string::npos)
	    {
	      startpos = static_cast<unsigned int> (position);
	      str = buffer.substr((startpos + 1));
	      istringstream ins(str);
	      ins >> ph;
	      sysptr->setcpwphase(ph);
	    }
	}
      else if (buffer.find("LAYER_START") != string::npos)
	{
	  if (getline(infile, laybuf) != NULL)
	    {
	      while (laybuf.find("LAYER_END") == string::npos)
		{
		  trimspaces(laybuf);
		  if (laybuf[0] == '#')
		    {
		      //Its a comment
		    }
		  else
		    {
		      readlayer(laybuf, epiptr);
		    }
		  getline(infile, laybuf);
		}
	    }
	}
      i++;
    }
  infile.close();
  	  
  return 0;
}


void trimspaces(string& str)
{
  size_t startpos, endpos;
  
  startpos = str.find_first_not_of(" \t");
  endpos = str.find_last_not_of(" \t");

  if ((string::npos == startpos) || (string::npos == endpos))
    {
      str = "";
    }
  else
    {
      str = str.substr(startpos, endpos-startpos+1);
    }
}
