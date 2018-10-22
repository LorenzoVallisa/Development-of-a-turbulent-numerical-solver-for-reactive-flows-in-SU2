#include "../include/reacting_model_library.hpp"
#include "../include/not_setup_exception.hpp"

#include <experimental/filesystem>
#include <locale>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace Framework {

su2double ReactingModelLibrary::Density(const su2double& temp,const su2double& pressure) {
  su2double mol_mass = 0.0; // Moles over mass of mixture
  for(auto i=0;i<nSpecies;++i)
    mol_mass += Ys[i]/mMasses[i];

  return pressure/(Rgas*mol_mass*temp);
}

void ReactingModelLibrary::SetMolarFractions(const RealVec& xs) {
  Xs = xs;

  for(auto i=0;i<nSpecies;++i) {
    if(Xs[i] < 0.0)
      Xs[i] = 0.0;

    assert(Xs[i]<=1.0);
  }
}

void ReactingModelLibrary::SetRgas(const RealVec& Ys,const RealVec& Ri) {
  Rgas = 0.0;
  assert(Ys.size()==nSpecies);
  assert(Ri.size()==nSpecies);
  for(auto i=0;i<nSpecies;++i)
    Rgas += Ys[i]*Ri[i];
}

void ReactingModelLibrary::SetMoleculesIDs(std::vector<unsigned>& Ids) {
  Ids.reserve(nSpecies);

  for(auto i=0;i<nSpecies;++i) {
    if(Atoms[i]>1)
      Ids.push_back(i);
  }
}

void ReactingModelLibrary::SetMassFractions(const RealVec& ys) {
  Ys = ys;

  for(auto i=0;i<nSpecies;++i) {
    if(Ys[i] < 0.0)
      Ys[i] = 0.0;

    assert(Ys[i]<=1.0);
  }
}

/* This file set the order of the species */

void ReactingModelLibrary::ReadDataMixture(const std::string& f_name) {

  std::string line;
  unsigned n_line = 0;

  std::ifstream mixfile(f_name);
  if(mixfile.is_open()) {
    while(mixfile.good() and !mixfile.eof()) {
      std::getline(mixfile,line);
      // We avoid clearly reading comments in the file
      if(line[0]!='/') {
        if(n_line==0 and std::isdigit(line[0])) {
          std::istringstream curr_line(line);
          curr_line>>nSpecies;
          n_line++;
          // Resize vector that wil be often used
          mMasses.resize(nSpecies);
          Ri.resize(nSpecies);
          Ys.resize(nSpecies);
          Xs.resize(nSpecies);
          Viscosities.resize(nSpecies);
          Internal_Energies.resize(nSpecies);
          Enthalpies.resize(nSpecies);
          Heat_Capacities.resize(nSpecies);
          CPs.resize(nSpecies);
          CVs.resize(nSpecies);
          Thermal_Conductivities.resize(nSpecies);
          Sound_Speeds.resize(nSpecies);
          Formation_Enthalpies.resize(nSpecies);
        }
      }
        else if(n_line!=0 and std::isalpha(line[0])) {
          Species_Names[n_line-1] = line;
          n_line++;
      }
    }
    mixfile.close();
  }
  else {
    std::cerr<<"Unable to open the mixture file: "<<f_name<<std::endl;
    std::exit(1);
  }

  assert(Species_Names.size() == nSpecies);
}

void ReactingModelLibrary::ReadDataChem(const std::string& f_name) {

  std::string line;
  unsigned n_line = 0;

  std::ifstream chemfile(f_name);
  if(chemfile.is_open()) {
    while(chemfile.good() and !chemfile.eof()) {
      std::getline(chemfile,line);
      // We avoid clearly reading comments in the file
      if(line[0]!='/') {
        if(n_line==0 and std::isdigit(line[0])) {
          n_line++;
          std::istringstream curr_line(line);
          curr_line>>nReactions;
        }
        else if(n_line!=0 and std::isalpha(line[0])) {
          ReadReactSpecies(Stoich_Coeffs,line);
          readChemCoefs(Chems_Coeffs,line);
        }
    }
    chemfile.close();
  }
  else {
    std::cerr<<"Unable to open the chemical file: "<<f_name<<std::endl;
    std::exit(1);
  }
}

void ReactingModelLibrary::ReadDataTransp(const std::string& f_name,const unsigned short& iSpecies) {

  std::string line;
  unsigned n_line = 0;

  std::ifstream speciesfile(f_name);

  if(speciesfile.is_open()) {
    while(speciesfile.good() and !speciesfile.eof()) {
      std::getline(speciesfile,line);

      if(line[0]!='/' and std::isdigit(line[0])) {
        std::istringstream curr_line(line);
        n_line++;
        if((n_line==1 && (curr_line>>Atoms[iSpecies]).fail())
         || (n_line==2 && (curr_line>>mMasses[iSpecies]).fail())
         || (n_line==3 && (curr_line>>Viscosities[iSpecies]).fail())
         || (n_line==4 && (curr_line>>Thermal_Conductivities[iSpecies]).fail()))
         {
           std::cerr<<"Wrong format of species file: "<<f_name<<std::endl;
           std::exit(1);
         }
      }
    }
    speciesfile.close();
  }
  else {
    std::cerr<<"Unable to open the species file: "<<f_name<<std::endl;
    std::exit(1);
  }
}

void ReactingModelLibrary::ReadDataThermo(const std::string& f_name,const unsigned short& iSpecies) {
  std::string line;
  unsigned n_line = 0;

  std::ifstream thermofile(f_name);

  if(thermofile.is_open()) {
    while(thermofile.good() and !thermofile.eof()) {
      std::getline(thermofile,line);

      if(line[0]!='/' and std::isdigit(line[0])) {
        std::istringstream curr_line(line);
        n_line++;
        if((n_line==1 && (curr_line>>Formation_Enthalpies[iSpecies]).fail())
        || (n_line==2 && (curr_line>>CPs[iSpecies]).fail())
        || (n_line==3 && (curr_line>>CVs[iSpecies]).fail())
        || (n_line==4 && (curr_line>>Enthalpies[iSpecies]).fail())
        || (n_line==5 && (curr_line>>Internal_Energies[iSpecies]).fail())
        || (n_line==6 && (curr_line>>Heat_Capacities[iSpecies]).fail())
        || (n_line==7 && (curr_line>>Sound_Speeds[iSpecies]).fail()))
        {
          std::cerr<<"Wrong format of thermo file: "<<f_name<<std::endl;
          std::exit(1);
        }
      }
    }
    thermofile.close();
  }
  else {
    std::cerr<<"Unable to open the thermo file: "<<f_name<<std::endl;
    std::exit(1);
  }
}

void ReactingModelLibrary::Setup(void) {
  if(!Lib_Setup) {
    Lib_Setup = true;

    AF = 0.192;
    T_ref = 298.16;
    Le = 1.0;

    // if nobody has configured the library path, we try to do it here with a default value
    if(Lib_Path=="") {
      std::cout<<"Library path set to default"<<std::endl;
      auto base_dir = std::experimental::filesystem::current_path().string();
      Lib_Path = base_dir + "../../Common/include";
    }

    std::vector<std::string> list_file;
    std::ifstream curr_file(File_Names);
    if(curr_file.is_open()) {
      while(curr_file.good() and !curr_file.eof()) {
        std::string curr_line;
        std::getline(curr_file,curr_line);
        if(curr_line[0]!='/')
          list_file.push_back(curr_line);
      }
    }
    else {
      std::cerr<<"Unable to open the specified file with all file names."<<std::endl;
      std::exit(1);
    }

    assert(list_file.size()==nSpecies+2);

    std::string file_mix = list_file[0];
    ReadDataMixture(file_mix);
    std::cout<<"Mixture Data read"<<std::endl;

    std::string file_chem = list_file[1];
    ReadDataChem(file_chem);
    std::cout<<"Chemical Reactions read"<<std::endl;

    for(auto i=0;i<nSpecies;++i) {
      std::string file_species = list_file[i*2+2];
      ReadDataTransp(file_species,i);
      std::cout<<"Transport Data Species "<<i<<" read"<<std::endl;
      std::string file_thermo = list_file[i*2+3];
      ReadDataThermo(file_thermo,i);
      std::cout<<"Thermo Data Species "<<i<<" read"<<std::endl;
    }

    SetRiGas(Ri);
    //SetMassFractions(Ys);
    SetRgas(Ys,Ri);
    std::cout<<"Library set."<<std::endl;
  }

  else
    throw Common::NotSetup("Trying to setup again without calling unsetup first.");
}

void ReactingModelLibrary::Unsetup(void) {
  if(Lib_Setup) {
    Species_Names.clear();
    Atoms.clear();
    Ri.clear();
    mMasses.clear();
    Ys.clear();
    Xs.clear();
    Viscosities.clear();
    Internal_Energies.clear();
    Enthalpies.clear();
    Heat_Capacities.clear();
    CPs.clear();
    CVs.clear();
    Thermal_Conductivities.clear();
    Sound_Speeds.clear();
    Formation_Enthalpies.clear();
    Stoich_Coeffs.clear();
    Chem_Coeffs.clear();

    Lib_Setup = false;
  }

  else
    throw Common::NotSetup("Trying to unsetup without calling setup first.");

}

}
