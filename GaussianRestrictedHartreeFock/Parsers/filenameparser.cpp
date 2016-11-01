#include "filenameparser.h"
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>

using std::cout;
using std::endl;
using std::string;
using boost::regex;
using boost::regex_replace;

string FileNameParser::findBasisFile(string atom, string basis) {
    string path = m_path + atom + "/";

    // Replace [ * ] --> s
    basis = regex_replace(basis, regex("\\*"), "s");
    // Replace [ + ] --> p
    basis = regex_replace(basis, regex("\\+"), "p");
    // Replace [ ( ] --> _
    basis = regex_replace(basis, regex("\\("), "_");
    // Replace [ ) ] --> _
    basis = regex_replace(basis, regex("\\)"), "_");
    // Replace [ whitespace ] --> _
    basis = regex_replace(basis, regex(" "),   "_");

    return path + basis;
}
