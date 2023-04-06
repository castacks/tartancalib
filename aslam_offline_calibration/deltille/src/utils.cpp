/**
* Copyright (C) 2017-present, Facebook, Inc.
*
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#if defined(WITH_BOOST)
#include <algorithm>

#include <boost/filesystem.hpp>

#include <checkerboard_detector/utils.h>

std::vector<std::string> getFilesWithExtension(std::string dname_,
                                               std::string ext_) {
  namespace fs = boost::filesystem;
  if (!fs::is_directory(dname_)) {
    throw std::runtime_error(
        "getFilesWithExtension: expecting a directory name");
  }

  std::vector<std::string> ret;
  for (auto &entry :
       boost::make_iterator_range(fs::directory_iterator(dname_))) {
    if (ext_ == fs::extension(entry))
      ret.push_back(entry.path().string());
  }

  std::sort(ret.begin(), ret.end());
  return ret;
}

ProgramOptions::ProgramOptions(std::string name) : _desc(name) {
  _desc.add_options()("help,h", "Print this help message");
}

void ProgramOptions::parse(int argc, char **argv) {
  namespace po = boost::program_options;
  try {
    po::store(po::parse_command_line(argc, argv, _desc), _vm);
  } catch (const std::exception &ex) {
    std::cerr << "Error parsing command line: " << ex.what() << std::endl;
    throw ex;
  }

  po::notify(_vm);

  if (hasOption("help"))
    printHelpAndExit();
}

bool ProgramOptions::hasOption(std::string option) const {
  return _vm.count(option) > 0;
}

void ProgramOptions::printHelp() const { std::cout << _desc << std::endl; }

void ProgramOptions::printHelpAndExit(int exit_code) const {
  std::cout << _desc << std::endl;
  exit(exit_code);
}

ProgramOptions &ProgramOptions::addOption(std::string name, std::string msg) {
  _desc.add_options()(name.c_str(), msg.c_str());
  return *this;
}

void Timer::start() { _start_time = std::chrono::high_resolution_clock::now(); }

auto Timer::stop() -> Milliseconds {
  auto t_now = std::chrono::high_resolution_clock::now();
  auto ret = std::chrono::duration_cast<Milliseconds>(t_now - _start_time);
  _start_time = std::chrono::high_resolution_clock::now();
  return ret;
}

auto Timer::elapsed() const -> Milliseconds {
  auto t_now = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<Milliseconds>(t_now - _start_time);
}

#endif
