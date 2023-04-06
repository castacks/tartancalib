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

#ifndef CHECKERBOARD_DETECTOR_UTILS_H
#define CHECKERBOARD_DETECTOR_UTILS_H

#include <chrono>
#include <iostream>
#include <string>
#include <vector>

template <class... T> inline void UNUSED(const T &...) {}

template <class T> inline void doNotOptimizeAway(T &&d) {
  asm volatile("" : "+r"(d));
}

#if defined(WITH_BOOST)

#include <boost/program_options.hpp>
#include <boost/program_options/value_semantic.hpp>

std::vector<std::string> getFilesWithExtension(std::string directory_name,
                                               std::string extension);

class ProgramOptions {
public:
  ProgramOptions(std::string name = "ProgramOptions");

  /**
   */
  ProgramOptions &addOption(std::string name, std::string msg);

  template <typename T>
  inline ProgramOptions &addOption(std::string name, T *value,
                                   std::string msg = "") {
    return this->operator()(name, value, msg);
  }

  /**
   */
  template <typename T>
  inline ProgramOptions &addOptionMultiToken(std::string name,
                                             std::vector<T> *v,
                                             std::string msg = "") {
    _desc.add_options()(
        name.c_str(),
        boost::program_options::value<std::vector<T>>(v)->multitoken(),
        msg.c_str());
    return *this;
  }

  void parse(int argc, char **argv);
  void printHelp() const;
  void printHelpAndExit(int exit_code = -1) const;
  bool hasOption(std::string) const;

  inline ProgramOptions &operator()(std::string name, std::string msg) {
    return this->addOption(name, msg);
  }

  template <class T>
  inline ProgramOptions &operator()(std::string name, T v, std::string msg) {
    _desc.add_options()(name.c_str(),
                        boost::program_options::value<T>()->default_value(v),
                        msg.c_str());
    return *this;
  }

  template <class T>
  inline ProgramOptions &operator()(std::string name, T *v, std::string msg) {
    _desc.add_options()(name.c_str(), boost::program_options::value<T>(v),
                        msg.c_str());
    return *this;
  }

  inline ProgramOptions &operator()(std::string name, const char *v,
                                    std::string msg) {
    return this->operator()(name, std::string(v), msg);
  }

  template <typename T> inline T get(std::string name) const {
    try {
      return _vm[name].template as<T>();
    } catch (const std::exception &ex) {
      std::cerr << "Error: " << ex.what() << std::endl;
      throw ex;
    }
  }

private:
  boost::program_options::options_description _desc;
  boost::program_options::variables_map _vm;
}; // ProgramOptions

/**
 * A class for simple wall clock timing
 */
class Timer {
public:
  typedef std::chrono::milliseconds Milliseconds;

public:
  inline Timer() { start(); }

  /**
   * starts the timer
   */
  void start();

  /**
   * \return elapsed time since the last call to start() and resets the
   * start point
   */
  Milliseconds stop();

  /**
   * \return elapsed time since the last call to start()
   */
  Milliseconds elapsed() const;

private:
  std::chrono::high_resolution_clock::time_point _start_time;
}; // Timer

/**
 * Runs 'Func' for 'N_rep" times and reports the average elapsed time in
 * milliseconds
 */
template <class Func, class... Args>
inline double TimeCode(int N_rep, Func &&f, Args... args) {
  Timer timer;
  for (int i = 0; i < N_rep; ++i)
    f(args...);
  auto t = timer.stop().count() / static_cast<double>(N_rep);
  doNotOptimizeAway(t);
  return t;
}

#endif

#endif
