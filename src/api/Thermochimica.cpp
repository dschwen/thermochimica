#include "thermochimica/Thermochimica.hpp"
#include "thermochimica/util/ErrorCodes.hpp"
#include <algorithm>
#include <cctype>
#include <string>

namespace Thermochimica {

// Element symbols for lookup
static const char* kElementSymbols[] = {
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
};

int getAtomicNumber(const std::string& symbol) {
    std::string upper = symbol;
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);

    // Handle special cases
    if (upper == "E-" || upper == "E(-)" || upper.substr(0, 2) == "E(") {
        return 0;  // Electron
    }

    // Search element symbols
    for (int i = 1; i <= 118; ++i) {
        std::string elemUpper = kElementSymbols[i];
        std::transform(elemUpper.begin(), elemUpper.end(), elemUpper.begin(), ::toupper);
        if (upper == elemUpper) {
            return i;
        }
    }

    return -1;  // Not found
}

std::string getElementSymbol(int atomicNumber) {
    if (atomicNumber >= 1 && atomicNumber <= 118) {
        return kElementSymbols[atomicNumber];
    }
    return "";
}

const char* getErrorMessage(int errorCode) {
    return ErrorCode::getMessage(errorCode);
}

} // namespace Thermochimica
