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

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include <deltille/target_detector.h>

#include <opencv2/highgui.hpp>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>


namespace fs = boost::filesystem;
namespace po = boost::program_options;

using namespace std;

/**
 */
bool writeCornersToFile(std::ostream &os,
                        const std::vector<CalibrationCorner> &corners,
                        string filename, const cv::Size &image_size,
                        bool write_ordered_only = true) {

  auto num_corners = corners.size();
  if(write_ordered_only) {
    num_corners = 0;
    for(auto& c : corners) {
      num_corners += c.isValid() && c.isOrdered;
    }
  }

  cout << "writing " << num_corners << "\n";

  os << "filename: " << filename << endl;
  os << "width: " << image_size.width << endl;
  os << "height: " << image_size.height << endl;
  os << "num_corners: " << num_corners << endl;
  os << "encoding: ascii" << endl;

  auto p = os.precision();
  os.precision(numeric_limits<double>::max_digits10);
  for (auto &c : corners) {
    if (!c.isValid() || (write_ordered_only && !c.isOrdered)) {
      continue;
    }

    os << c << endl;
  }

  os.precision(p);
  return os.good();
}

/**
 */
class DataSource {
public:
  virtual ~DataSource() {}

  bool getImage(cv::Mat &image, int index = -1) {
    if (this->get_image(image, index)) {
      convert_to_grayscale(image);
      convert_type(image);
      return true;
    } else {
      return false;
    }
  }

  const std::string &getLastFilename() const { return _last_file_name; }

private:
  virtual bool get_image(cv::Mat &, int) = 0;

  void convert_to_grayscale(cv::Mat &image) const {
    if (image.channels() == 3) {
      cv::cvtColor(image, image, cv::COLOR_BGR2GRAY);
    } else if (image.channels() == 4) {
      cv::cvtColor(image, image, cv::COLOR_BGRA2GRAY);
    }
  }

  void convert_type(cv::Mat &image) const {
    if (image.depth() == cv::DataDepth<uint16_t>::value) {
      double max_val = 0.0;
      cv::minMaxLoc(image, nullptr, &max_val);
      image.convertTo(image, CV_MAKETYPE(cv::DataDepth<float>::value, 1),
                      255.0 * (1.0 / max_val));
    }
  }

protected:
  string _last_file_name;
};

/**
 */
class ImageListDataSource : public DataSource {
public:
  ImageListDataSource(vector<string> &&filenames)
      : _filenames(move(filenames)) {}

private:
  bool get_image(cv::Mat &I, int f_i) override {
    if (f_i < 0)
      f_i = _counter++;

    if (std::size_t(f_i) < _filenames.size()) {
      this->_last_file_name = _filenames[f_i];
      return !(I = cv::imread(_filenames[f_i],
                              cv::IMREAD_ANYDEPTH | cv::IMREAD_GRAYSCALE))
                  .empty();
    } else {
      return false;
    }
  }

private:
  int _counter{0};
  vector<string> _filenames;
};

void RunDetector(DataSource *data_source, string target_dsc_fn,
                 const fs::path &output_dir, bool debug) {
  TargetDetector target_detector(target_dsc_fn);

  string debug_path;
  if (debug) {
    debug_path = (output_dir / fs::path("deltille_out")).string();

    if (!fs::exists(debug_path))
      if (!fs::create_directory(debug_path)) {
        throw runtime_error("failed to create debug output directory");
      }

    if (!fs::is_directory(debug_path)) {
      throw invalid_argument("invalid debug output path");
    }
  }

  if (debug) {
    cv::namedWindow("detection", cv::WINDOW_NORMAL);
  }

  auto do_write_corners = !output_dir.string().empty();

  // NOTE: if debug == false, this could be run in parallel
  cv::Mat I, debug_image;
  for (int i = 0; data_source->getImage(I, i); ++i) {
    if (!I.empty()) {
      vector<CalibrationCorner> corners;
      target_detector.run(I, corners, debug ? &debug_image : nullptr);

      if (do_write_corners) {
        auto filename = data_source->getLastFilename();
        auto basename = fs::path(filename).stem();
        auto output_filename =
            output_dir / fs::change_extension(basename, ".orpc");
 
        std::ofstream file(output_filename.string());
        if (file.is_open()) {
          writeCornersToFile(file, corners, filename, I.size(), true);
        } else {
          cerr << "Failed to open: " << output_filename << " for writing"
               << endl;
        }
      }

      if (debug) {
        cv::imshow("detection", debug_image);
        int k = cv::waitKey(100);
        if (k == ' ') // pause on space
          k = cv::waitKey(0);
        if (k == 'q' || k == 27) // quite on 'q' or ESC
          break;

        auto fn_stem = fs::path(data_source->getLastFilename()).stem();
        auto fn = fs::change_extension(fs::path(debug_path) / fn_stem, ".png");
        cout << "Writing detection result to : " << fn << endl;
        cv::imwrite(fn.string(), debug_image);
      }
    }
  }
}

int main(int argc, char **argv) {
  string target_dsc_fn;

  vector<string> files;
  fs::path output_dir;

  po::options_description desc(argv[0]);
  desc.add_options()("help,h", "Produce this help message")(
      "target,t", po::value<string>(&target_dsc_fn)->required(),
      "Target *.dsc file")("files,f",
                           po::value<vector<string>>(&files)->multitoken(),
                           "List of image files")(
      "output,o", po::value<string>(), "Output directory")(
      "debug,d", "draw the detected corners and store the result to disk");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (const exception &ex) {
    cerr << ex.what() << endl << endl;
    cerr << desc << endl;
    return 1;
  }

  if (vm.count("help")) {
    cout << desc << endl;
    return 1;
  }

  if (!fs::exists(target_dsc_fn) || !fs::is_regular_file(target_dsc_fn)) {
    throw invalid_argument("invalid target *.dsc file '" + target_dsc_fn + "'");
  }

  if (vm.count("output")) {
    output_dir = fs::path(vm["output"].as<string>());
    if (!fs::exists(output_dir)) {
      if (!fs::create_directory(output_dir)) {
        throw invalid_argument("invalid output directory " +
                               output_dir.string());
      }
    } else {
      if (!fs::is_directory(output_dir)) {
        throw invalid_argument("argument to --output is NOT a directory");
      }
    }
  }

  if (!files.empty()) {
    ImageListDataSource data_source(move(files));
    RunDetector(&data_source, target_dsc_fn, output_dir, vm.count("debug"));
  }

  return 0;
}
