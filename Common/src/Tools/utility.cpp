#include "../include/Tools/utility.hpp"

#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

namespace Utility {
  //
  //
  /*--- This function reads a reaction from line ---*/
  void Parse_Terms(std::string& line, unsigned n_reac, bool is_rev, bool is_reac, const MyMap& map_names, RealMatrix& stoich_coeff,
                   RealMatrix& stoich_coeff_exp_reac, RealMatrix& stoich_coeff_exp_prod) {
     auto size = line.size();
     std::size_t idx = 0;
     std::string coeff, symbol, exp_coeff;

     /*--- Find first stoichiometric coefficient or first letter of species name ---*/
     while(!std::isdigit(line.at(idx)) && !std::isalpha(line.at(idx)))
       idx++;

     /*--- Extract leading coefficient ---*/
     while(std::isdigit(line.at(idx)) || std::ispunct(line.at(idx))) {
       coeff += line.at(idx);
       idx++;
     }
     do {
       /*--- Extract chemical symbol ---*/
       while((std::isalpha(line[idx]) || std::isdigit(line[idx])) && idx < size) {
         symbol += line.at(idx);
         idx++;
       }
     } while(idx < size && !std::ispunct(line.at(idx)) && !std::isspace(line.at(idx)) && line.at(idx) != '+');

     auto it = map_names.find(symbol);
     SU2_Assert(it != map_names.end(),std::string("The symbol " + symbol + " is not in the mixture list"));

     /*--- Retrieve the index of the species detected ---*/
     unsigned short idx_species = it->second;

     /*--- Saving stoichiometric coefficient ---*/
     double coefficient;
     if(coeff.empty())
       coefficient = 1.0;
     else {
       std::istringstream str_to_coeff(coeff);
       str_to_coeff>>coefficient;
     }
     stoich_coeff(idx_species,n_reac - 1) += coefficient;

     /*--- Saving possibly coefficient at the exponent of coencentration in reaction rate ---*/
     if(idx < size) {
       if(std::ispunct(line.at(idx))) {
         idx++;
         while(idx < size && (std::isdigit(line.at(idx)) || std::ispunct(line.at(idx)))) {
           exp_coeff += line.at(idx);
           idx++;
         }
       }
     }

     double exp_coefficient = 0.0;
     if(!exp_coeff.empty()) {
       std::istringstream str_to_expcoeff(exp_coeff);
       str_to_expcoeff>>exp_coefficient;
       if(is_reac)
        stoich_coeff_exp_reac(n_reac - 1,idx_species) += exp_coefficient;
       else
        if(is_rev)
          stoich_coeff_exp_prod(n_reac - 1,idx_species) += exp_coefficient;
     }
     else {
       if(is_reac)
         stoich_coeff_exp_reac(n_reac - 1,idx_species) += stoich_coeff(idx_species,n_reac - 1);
     }

     /*--- Any more terms to extract? ---*/
     if(idx == size)
       return;

     std::string sub = line.substr((idx + 1), (size - idx));
     if(sub.empty())
       return;
     else
       Parse_Terms(sub,n_reac,is_rev,is_reac,map_names,stoich_coeff,stoich_coeff_exp_reac,stoich_coeff_exp_prod);
  }

  //
  //
  /*--- This function reads extra terms to compute correctly forward and backward rates ---*/
  void Parse_ExtraTerms(std::string& line, unsigned idx_reac, const MyMap& map_names, RealMatrix& stoich_coeff_exp) {
    std::size_t size = line.size();
    std::size_t idx = 0;
    std::string symbol, exp_coeff;
    double exp_coefficient;

    while(!std::isalpha(line.at(idx)))
      idx++;

    /*--- Extract the symbol of species ---*/
    while(!std::ispunct(line.at(idx))) {
      symbol += line.at(idx);
      idx++;
    }

    /*--- Check if the species in the mixture ---*/
    auto it = map_names.find(symbol);
    SU2_Assert(it != map_names.end(),std::string("The symbol " + symbol + " is not in the mixture list"));

    /*--- Retrieve the index of the species detected ---*/
    unsigned short idx_species = it->second;

    /*--- Jump the division symbol between species name and exponent value ---*/
    idx++;
    while(idx < size) {
      exp_coeff += line.at(idx);
      idx++;
    }
    SU2_Assert(!exp_coeff.empty(), "No exponent detected after the species in the extra term for rates");
    std::istringstream str_to_expcoeff(exp_coeff);
    str_to_expcoeff>>exp_coefficient;

    /*--- Save the detected exponent ---*/
    if(stoich_coeff_exp(idx_reac, idx_species) != 0.0)
      std::cerr<<"Warning: You are going to overwrite a term you already read"<<std::endl;
    stoich_coeff_exp(idx_reac, idx_species) = exp_coefficient;

    /*--- Any more terms to extract? ---*/
    if(idx == size)
      return;

    std::string sub = line.substr((idx + 1), (size - idx));
    if(sub.empty())
      return;
    else
      Parse_ExtraTerms(sub, idx_reac, map_names, stoich_coeff_exp);
   }

} /*--- End of namespace Utility ---*/
