//
// Created by Xin Sui on 2/28/18.
//

#ifndef PROTEIN_REGISTRATION_CONFIGURATION_H
#define PROTEIN_REGISTRATION_CONFIGURATION_H

#include <boost/filesystem.hpp>

struct ConfigData{
    boost::filesystem::path data_directory;
    boost::filesystem::path out_file_name;
    unsigned long N;
    unsigned long dSSD_cutoff;
    unsigned long SSD_cutoff;
};

ConfigData configuration (const boost::filesystem::path& );

#endif //PROTEIN_REGISTRATION_CONFIGURATION_H
