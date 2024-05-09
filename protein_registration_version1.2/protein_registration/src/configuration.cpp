//
// Created by Xin Sui on 2/28/18.
//

#include "configuration.h"
#include <boost/algorithm/string/trim.hpp>
#include <iostream>
#include <string>

ConfigData configuration (const boost::filesystem::path& cfg_file){

    if ( !boost::filesystem::exists(cfg_file) ) {
        std::cerr << cfg_file << " does not exist. " << std::endl;
    }

    ConfigData ret;

    boost::filesystem::ifstream cfg;
    cfg.open(cfg_file, std::ios::in);

    std::string line;
//    while(getline(cfg, line))
//    {
        // line with #
//        if ( line.find_first_of('#') == std::string::npos )
//            break;
//
//    }
    getline(cfg, line);
    ret.data_directory = line;
    if ( !boost::filesystem::is_directory(ret.data_directory) ) {
        std::cerr << ret.data_directory << " does not exist. " << std::endl;
    }

    getline(cfg, line);
    ret.out_file_name = line;

    cfg >> ret.N >> ret.dSSD_cutoff >> ret.SSD_cutoff;

    return ret;



}